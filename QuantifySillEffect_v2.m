% QuantifySillEffect_v1

% Mike Wolovick, 6/16/2017

% This script quantifies the effect of the artificial sill, relative to the
% warming run. It computes the prevented sea level rise, prevented
% grounding line retreat, sea level delay, and grounding line delay.

% Prevented sea level rise and prevented grounding line retreat both have
% units of meter-years.  However, for sea level rise the 'meters' refers to
% eustatic sea level, while for the grounding line 'meters' refers to
% grounding line position.

% Prevented sea level rise and prevented grounding line retreat are both
% computed as the integral of the difference between the sill run time
% series and the warming run time series.  

% For narrow flowband boundaries, a correction factor is applied to sea 
% level rise to account for the unmodeled catchment.  This correction 
% factor is just the ratio of the wide flowband area to the narrow flowband
% area.  

% If a run crashes, both VAF and GL are assumed to drop to zero immediately 
% after the crash.  

% Ice melt has the density of fresh water for the purpose of computing sea 
% level rise.

% The script assumes that all of the time series have the same dt.

% v2:  Now the script also calculates:
%    Maximum prevented sea level rise (at any instant in time)
%    Maximum prevented GL retreat (at any instant in time)
%    Retreat delay
%    GL recovery (max advance after sill built)
%    Sea level recovery (max VAF increase after sill built)
% The above are calculated to produce values with understandable units
% instead of "meter-years".  Delays are only meaningful when both the sill
% run and the warming run experience retreat and mass loss.  When one or
% both runs advances, then the delay is set to NaN.
% The script outputs figures visualizing how it measured all of the above
% so that I can check that the measurements are reasonable.
% In addition, the script outputs a table of all model runs in csv format
% so I can easily make a table for the supplement.  
% Text files output:
%    Basic experiment (minus Jakobshavn)
%    Thwaites v3a (with column for sill type)
%    Helheim/Kanger v3a
% All numbers in the text file output are rounded to two significant
% figures.

clear variables
tic

%% Parameters:

% Files and folders:
timeseriesfile='/net/mjw/FjordSillPaper/ModelOutput/ParamSweepRuns_v3_v3a.mat';
inputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/';
outputfile='/net/mjw/FjordSillPaper/ModelOutput/QuantifiedSillEffect_v2.mat';
figfolder='/net/mjw/FjordSillPaper/Figures/SillQuantification_v2/';
% timeseriesfile='/home/wolovick/Dropbox/FjordSillPaper/ModelOutput/ParamSweepRuns_v3_v3a.mat';
% inputfolder='/home/wolovick/Dropbox/FjordSillPaper/ModelInput/';
% outputfile='/home/wolovick/Dropbox/FjordSillPaper/ModelOutput/QuantifiedSillEffect_v2.mat';
% outputtextfile='/home/wolovick/Dropbox/FjordSillPaper/ModelOutput/QuantifiedSillEffect_v2.csv';
% figfolder='/home/wolovick/Dropbox/FjordSillPaper/Figures/SillQuantification_v2/';
figsuffix='sillquantification_v2';

% Surface area of the ocean:
oceanarea=3.6e14;               % m^2

% Delay measurement settings:
icbuffer=10;                    % integer 
minretreat=1e3;                 % m
minretreatratio=.1;             % unitless

% Margins:
verttextbuffer=.075;            % unitless
titlebuffer=.06;                % unitless
horztextbuffer=.11;             % unitless
horznotextbuffer=.03;           % unitless

% Font size: 
fontsize=12;                    % points

% Pagesize and resolution:
pagesize=[10,7];                % [1x2] inches
resolution=300;                 % dpi

%% Work:

% Hardcoded plot arrangement:
numvertplots=2;
numhorzplots=1;

% Compute subplot positions: 
availablehorzspace=1-horztextbuffer-numhorzplots*horznotextbuffer;
availablevertspace=1-verttextbuffer-numvertplots*titlebuffer;
Boxes=cell(numvertplots,numhorzplots);
% Loop through boxes: 
for ii=1:numvertplots
    for jj=1:numhorzplots
        Boxes{ii,jj}=[horztextbuffer+(jj-1)*(horznotextbuffer+availablehorzspace/numhorzplots),... 
            verttextbuffer+(numvertplots-ii)*(titlebuffer+availablevertspace/numvertplots),...
            availablehorzspace/numhorzplots,availablevertspace/numvertplots];
    end
end

% Load the time series file:
load(timeseriesfile,'TimeSeries','Parameters')

% Compute number of model runs:
numruns=length(TimeSeries);

% Compute number of sill runs:
numsillruns=0;
IsSillRun=false(numruns,1);
for ii=1:numruns
    % Locate underscores in filename:
    underscorenum=strfind(TimeSeries(ii).filename,'_');
    % Determine if this is a sill run:
    if strcmp(TimeSeries(ii).filename(underscorenum(2)+1),'s')
        IsSillRun(ii)=1;
        numsillruns=numsillruns+1;
    end
end

% Compute dt:
dt=TimeSeries(1).Time(2)-TimeSeries(1).Time(1);

% Pre-allocate variables:
SillRunFileName=cell(numsillruns,1);
WarmingRunFileName=cell(numsillruns,1);
SillRunParameters=repmat(Parameters(1),[numsillruns,1]);
Flowband=cell(numsillruns,1);
SillKey=zeros(numsillruns,1);
M=zeros(numsillruns,1);
Calving=cell(numsillruns,1);
Isv3a=false(numsillruns,1);
IntegratedPreventedSeaLevelRise=zeros(numsillruns,1);
IntegratedPreventedGLRetreat=zeros(numsillruns,1);
MaxPreventedSeaLevelRise=zeros(numsillruns,1);
MaxPreventedGLRetreat=zeros(numsillruns,1);
RetreatDelay=zeros(numsillruns,1);
WarmingRunCollapse=false(numsillruns,1);
SillRunCollapse=false(numsillruns,1);
GLRecovery=zeros(numsillruns,1);
SeaLevelRecovery=zeros(numsillruns,1);
SillBuiltUnderGroundedIce=false(numsillruns,1);

% Make notes:
NOTE_units='Integrated prevented sea level rise is in meter-years (meters of eustatic sea level), integrated prevented GL retreat is also in meter-years, but meters of GL position.  All grounding line distances are in meters.';
NOTE_flowbandwidth='For narrow flowbands ("NameB"), a correction factor has been applied to sea level rise to account for the difference in flowband area between the wide and narrow models.';
NOTE_sillkey='0=normal sill run, 1,2,3=s1,s2,s3 (the thwaites runs with sill erosion)';

% Loop through all runs:
sillruncounter=1;
for ii=1:numruns                    
    
    % Check if this is a sill run:
    if IsSillRun(ii)==0
        continue
    end
    
    % Locate underscores in filename:
    underscorenum=strfind(TimeSeries(ii).filename,'_');
    
    % Parse filename:
    flowband=TimeSeries(ii).filename(1:underscorenum(1)-1);
    sliding=TimeSeries(ii).filename(underscorenum(3)+1:underscorenum(4)-1);
    calving=TimeSeries(ii).filename(underscorenum(4)+1:underscorenum(5)-1);
    versionkey=TimeSeries(ii).filename(underscorenum(5)+1:end-4);
    
    % Find matching warming run:
    thiswarmingrun=NaN;
    for jj=1:numruns
        % Locate underscores in filename:
        underscorenum=strfind(TimeSeries(jj).filename,'_');
        % Check that this is a warming run:
        if strcmp(TimeSeries(jj).filename(underscorenum(2)+1),'w')==0
            continue
        end
        % Check that flowband matches:
        if strcmp(TimeSeries(jj).filename(1:underscorenum(1)-1),flowband)==0
            continue
        end
        % Check that sliding matches:
        if strcmp(TimeSeries(jj).filename(underscorenum(3)+1:underscorenum(4)-1),sliding)==0
            continue
        end
        % Check that calving matches:
        if strcmp(TimeSeries(jj).filename(underscorenum(4)+1:underscorenum(5)-1),calving)==0
            continue
        end
        % Check that version key matches:
        if strcmp(TimeSeries(jj).filename(underscorenum(5)+1:end-4),versionkey)==0
            % Check if this is the v3a sill runs that use v3 warming runs
            % as their reference:
            if strcmp(versionkey,'v3a')==0 || strcmp(flowband(1:min(length(flowband),10)),'Jakobshavn')
                continue
            end
        end
        % Everything matches, this is the warming run:
        thiswarmingrun=jj;
        break
    end
    
    % Check if the warming run was located:
    if isnan(thiswarmingrun)
        error('Unable to locate warming run.')
    end
    
    % Re-Locate underscores in filename:
    underscorenum=strfind(TimeSeries(ii).filename,'_');
    
    % Check if we need a correction factor for the narrow flowbands:
    if strcmp(TimeSeries(ii).filename(underscorenum(1)-1),'B')
        % Load narrow width:
        load([inputfolder,TimeSeries(ii).filename(1:underscorenum(2)-1),'.mat'],'Width_input','lasticeind','steplength')
        % Compute narrow area:
        narrowarea=sum(Width_input(1:lasticeind+1)*steplength);
        % Compute wide filename:
        widefilename=[TimeSeries(ii).filename(1:underscorenum(1)-2),'A',TimeSeries(ii).filename(underscorenum(1):underscorenum(2)-1)];
        % Load wide width:
        load([inputfolder,widefilename,'.mat'],'Width_input','lasticeind','steplength')
        % Compute wide area:
        widearea=sum(Width_input(1:lasticeind+1)*steplength);
        % Compute correction factor:
        correctionfactor=widearea/narrowarea;
    else
        % Correction factor is 1:
        correctionfactor=1;
    end
    
    % Assign file names and parameters:
    SillRunFileName{sillruncounter}=TimeSeries(ii).filename;
    WarmingRunFileName{sillruncounter}=TimeSeries(thiswarmingrun).filename;
    SillRunParameters(sillruncounter)=Parameters(ii);
    
    % Assign experiment variables:
    Flowband{sillruncounter}=flowband;
    M(sillruncounter)=str2double(sliding(2:end));
    Calving{sillruncounter}=calving;
    
    % Assign sill key:
    if strcmp(TimeSeries(ii).filename(underscorenum(2)+1:underscorenum(3)-1),'s')
        SillKey(sillruncounter)=0;
    else
        SillKey(sillruncounter)=str2double(TimeSeries(ii).filename(underscorenum(3)-1));
    end
    
    % Assign version key:
    if strcmp(versionkey,'v3a')
        Isv3a(sillruncounter)=1;
    end
    
%     if strcmp(flowband,'ThwaitesC')==0 || strcmp(sliding,'m10')==0 || strcmp(calving,'chm')==0 || SillKey(sillruncounter)~=1
%         sillruncounter=sillruncounter+1;
%         continue
%     end
    
    % Compute numrecords:
    numrecords=length(TimeSeries(ii).Time);
    
    % Compute conversion factor between vaf and sea level:
    vafslfactor=correctionfactor*(Parameters(ii).ThermalParameters.rho_i/Parameters(ii).ThermalParameters.rho_fw)/oceanarea;
    
    % Check whether either run collapsed:
    lastind_sill=find(TimeSeries(ii).VAF~=0,1,'last');
    lastind_warming=find(TimeSeries(thiswarmingrun).VAF~=0,1,'last');
    if lastind_sill<numrecords
        SillRunCollapse(sillruncounter)=1;
    end
    if lastind_warming<numrecords
        WarmingRunCollapse(sillruncounter)=1;
    end
    
    % Make sure VAF vector doesn't overrun time vector: (kludge shit)
    if TimeSeries(ii).Time(lastind_sill)==0
        TimeSeries(ii).VAF(lastind_sill)=0;
        lastind_sill=lastind_sill-1;
    end
    if TimeSeries(thiswarmingrun).Time(lastind_warming)==0
        TimeSeries(thiswarmingrun).VAF(lastind_warming)=0;
        lastind_warming=lastind_warming-1;
    end
    
    % Compute integrated prevented sea level rise:
    IntegratedPreventedSeaLevelRise(sillruncounter)=vafslfactor*dt*sum(TimeSeries(ii).VAF-TimeSeries(thiswarmingrun).VAF);
    
    % Compute maximum prevented sea level rise:
    MaxPreventedSeaLevelRise(sillruncounter)=vafslfactor*max(TimeSeries(ii).VAF-TimeSeries(thiswarmingrun).VAF);
    
    % Compute integrated prevented grounding line retreat:
    IntegratedPreventedGLRetreat(sillruncounter)=dt*sum(TimeSeries(ii).X_gl-TimeSeries(thiswarmingrun).X_gl);
    
    % Compute maximum prevented grounding line retreat:
    MaxPreventedGLRetreat(sillruncounter)=max(TimeSeries(ii).X_gl-TimeSeries(thiswarmingrun).X_gl);
    
    % Compute grounding line recovery due to sill:
    sillstartind=find(TimeSeries(ii).Time>=Parameters(ii).SillParameters.sillstarttime_yr,1,'first');
    sillstopind=find(TimeSeries(ii).Time>=Parameters(ii).SillParameters.sillstarttime_yr+Parameters(ii).SillParameters.sillconstructiontime_yr,1,'first');
    glmaxind=sillstartind-1+find(TimeSeries(ii).X_gl(sillstartind:end)==max(TimeSeries(ii).X_gl(sillstartind:end)),1,'first');
    glmin=min(TimeSeries(ii).X_gl(sillstartind:glmaxind));
    GLRecovery(sillruncounter)=max(TimeSeries(ii).X_gl(sillstartind:end)-glmin);
    
    % Compute sea level recovery due to sill:
    slmaxind=sillstartind-1+find(TimeSeries(ii).VAF(sillstartind:end)==max(TimeSeries(ii).VAF(sillstartind:end)),1,'first');
    slmin=min(TimeSeries(ii).VAF(sillstartind:slmaxind))*vafslfactor;
    SeaLevelRecovery(sillruncounter)=max(TimeSeries(ii).VAF(sillstartind:end)*vafslfactor-slmin);
    %     slmin=min(TimeSeries(ii).VAF(sillstartind:end))*vafslfactor;
    %     if slmin~=0
    %         slminind=sillstartind-1+find(abs(TimeSeries(ii).VAF(sillstartind:end)-slmin/vafslfactor)/(slmin/vafslfactor)<1e-14,1,'first');
    %         slmax=max(TimeSeries(ii).VAF(slminind:end)*vafslfactor);
    %         SeaLevelRecovery(sillruncounter)=slmax-slmin;
    %     end
    
    % Check whether the sill was built under grounded ice:
    if max(TimeSeries(ii).X_gl(sillstartind:sillstopind))>Parameters(ii).SillParameters.sillcentx-Parameters(ii).SillParameters.sillwidth
        SillBuiltUnderGroundedIce(sillruncounter)=1;
    end
    
    
    % COMPUTE DELAY:
    
    % Check whether any runs collapsed:
    if SillRunCollapse(sillruncounter) && WarmingRunCollapse(sillruncounter) % both runs collapsed
        
        % Uce point-based measurement at collapse:
        % Check which grounding line was more advanced at collapse:
        if TimeSeries(ii).X_gl(lastind_sill)>TimeSeries(thiswarmingrun).X_gl(lastind_warming) % sill run more advanced
            % Locate last crossing of sill-run collapse point in warming run time series:
            % determine level and time at sill-run collapse:
            gllevel=TimeSeries(ii).X_gl(lastind_sill);
            silltime_gl=TimeSeries(ii).Time(lastind_sill);
            % locate index:
            thisind=find((TimeSeries(thiswarmingrun).X_gl(1:end-1)<=gllevel&TimeSeries(thiswarmingrun).X_gl(2:end)>gllevel)|...
                (TimeSeries(thiswarmingrun).X_gl(1:end-1)>gllevel&TimeSeries(thiswarmingrun).X_gl(2:end)<=gllevel),1,'last');
            % interpolate exact time:
            warmingtime_gl=TimeSeries(thiswarmingrun).Time(thisind)+dt*(gllevel-TimeSeries(thiswarmingrun).X_gl(thisind))/...
                (TimeSeries(thiswarmingrun).X_gl(thisind+1)-TimeSeries(thiswarmingrun).X_gl(thisind));
            % compute delay:
            RetreatDelay(sillruncounter)=silltime_gl-warmingtime_gl;
        else % warming run more advanced
            % Locate last crossing of warming-run collapse point in sill run time series:
            % determine level and time at warming-run collapse:
            gllevel=TimeSeries(thiswarmingrun).X_gl(lastind_warming);
            warmingtime_gl=TimeSeries(thiswarmingrun).Time(lastind_warming);
            % locate index:
            thisind=find((TimeSeries(ii).X_gl(1:end-1)<=gllevel&TimeSeries(ii).X_gl(2:end)>gllevel)|...
                (TimeSeries(ii).X_gl(1:end-1)>gllevel&TimeSeries(ii).X_gl(2:end)<=gllevel),1,'last');
            % interpolate exact time:
            silltime_gl=TimeSeries(ii).Time(thisind)+dt*(gllevel-TimeSeries(ii).X_gl(thisind))/...
                (TimeSeries(ii).X_gl(thisind+1)-TimeSeries(ii).X_gl(thisind));
            % compute delay:
            RetreatDelay(sillruncounter)=silltime_gl-warmingtime_gl;
        end
        
    elseif SillRunCollapse(sillruncounter) && WarmingRunCollapse(sillruncounter)==0 % only the sill run collapsed
        
        % Use collapse time to define (negative) retreat delay:
        RetreatDelay(sillruncounter)=TimeSeries(ii).Time(lastind_sill)-Parameters(ii).TimingParameters.runtime_yr;
        
        % Define display variables:
        gllevel=TimeSeries(ii).X_gl(lastind_sill);
        silltime_gl=TimeSeries(ii).Time(lastind_sill);
        warmingtime_gl=Parameters(ii).TimingParameters.runtime_yr;
        
    elseif SillRunCollapse(sillruncounter)==0 && WarmingRunCollapse(sillruncounter) % only the warming run collapsed
        
        % Use collapse time to define minimum retreat delay:
        RetreatDelay(sillruncounter)=Parameters(ii).TimingParameters.runtime_yr-TimeSeries(thiswarmingrun).Time(lastind_warming);
        
        % Define display variables:
        gllevel=TimeSeries(thiswarmingrun).X_gl(lastind_warming);
        silltime_gl=Parameters(ii).TimingParameters.runtime_yr;
        warmingtime_gl=TimeSeries(thiswarmingrun).Time(lastind_warming);
        
    else % both runs did not collapse
        
        % Check whether delays are meaningful in the absence of a collapse:
        if TimeSeries(thiswarmingrun).VAF(end)>TimeSeries(thiswarmingrun).VAF(1) || ... % did the warming run gain mass?
                TimeSeries(thiswarmingrun).X_gl(end)>TimeSeries(thiswarmingrun).X_gl(1)-minretreat || ... % did the warming run not retreat enough?
                TimeSeries(ii).VAF(end)>TimeSeries(ii).VAF(1) || ... % did the sill run gain mass?
                TimeSeries(ii).X_gl(end)>TimeSeries(ii).X_gl(1)-minretreat || ... % did the sill run not retreat enough?
                TimeSeries(ii).X_gl(1)-TimeSeries(ii).X_gl(end)<minretreatratio*(TimeSeries(thiswarmingrun).X_gl(1)-TimeSeries(thiswarmingrun).X_gl(end)) || ... % did the sill run not retreat enough in comparison to the warming run?
                MaxPreventedGLRetreat(sillruncounter)<minretreatratio*(TimeSeries(thiswarmingrun).X_gl(1)-TimeSeries(thiswarmingrun).X_gl(end))  % was the sill run grounding line history not sufficiently different from the warming run history?
            
            % Set delay to NaN:
            RetreatDelay(sillruncounter)=NaN;
            
        else % delays are meaningful
            
            % Uce cross-correlation analysis:
            
            % Compute cross-correlation of gl rate:
            [Xcor,Lags]=xcorr(diff(TimeSeries(thiswarmingrun).X_gl(icbuffer:end))-mean(diff(TimeSeries(thiswarmingrun).X_gl(icbuffer:end))),...
                diff(TimeSeries(ii).X_gl(sillstopind+icbuffer:end))-mean(diff(TimeSeries(ii).X_gl(sillstopind+icbuffer:end))));
            
            % Determine retreat delay from maximum correlation:
            delayind=sillstopind-Lags(find(Xcor==max(Xcor),1,'first'));
            RetreatDelay(sillruncounter)=dt*delayind;
            
            % Quality control:
            % Figure out indices of comparison:
            if delayind>=0
                firstind_sill=sillstopind+delayind;
                lastind_sill=numrecords;
                firstind_warming=sillstopind;
                lastind_warming=numrecords-delayind;
            else
                firstind_sill=sillstopind;
                lastind_sill=numrecords+delayind;
                firstind_warming=sillstopind-delayind;
                lastind_warming=numrecords;
            end
            % Ensure this delay improves the fit between the two raw time series
            if abs(mean(TimeSeries(ii).X_gl(firstind_sill:lastind_sill)-TimeSeries(thiswarmingrun).X_gl(firstind_warming:lastind_warming)))>...
                    abs(mean(TimeSeries(ii).X_gl(sillstopind:end)-TimeSeries(thiswarmingrun).X_gl(sillstopind:end)))
                RetreatDelay(sillruncounter)=NaN;
            end
            
        end
        
    end
    
    if ii<541
        sillruncounter=sillruncounter+1;
        continue
    end
    
    % CREATE FIGURE:
    % Delete extra handle:
    if exist('h3','var')
        clear h3
    end
    % Call figure:
    figure(1)
    % Call first subplot:
    subplot('Position',Boxes{1})
    % Plot zero line:
    hold off
    plot([0,Parameters(ii).TimingParameters.runtime_yr],[0,0],'--k')
    hold on
    % Compute ylims:
    ylims_gl=[min([TimeSeries(ii).X_gl(TimeSeries(ii).X_gl~=0)-TimeSeries(ii).X_gl(1);TimeSeries(thiswarmingrun).X_gl(TimeSeries(thiswarmingrun).X_gl~=0)-TimeSeries(thiswarmingrun).X_gl(1)]),...
        max([TimeSeries(ii).X_gl(TimeSeries(ii).X_gl~=0)-TimeSeries(ii).X_gl(1);TimeSeries(thiswarmingrun).X_gl(TimeSeries(thiswarmingrun).X_gl~=0)-TimeSeries(thiswarmingrun).X_gl(1)])]/1e3;
    % Plot sill extent:
    sillstarttime=Parameters(ii).SillParameters.sillstarttime_yr;
    sillendtime=Parameters(ii).SillParameters.sillstarttime_yr+Parameters(ii).SillParameters.sillconstructiontime_yr;
    patch('Vertices',[[sillstarttime;sillstarttime;Parameters(ii).TimingParameters.runtime_yr;Parameters(ii).TimingParameters.runtime_yr],...
        (Parameters(ii).SillParameters.sillcentx-TimeSeries(ii).X_gl(1)+Parameters(ii).SillParameters.sillwidth*[-1;1;1;-1])/1e3],...
        'Faces',1:4,'FaceVertexCData',[.75,.75,.75],'EdgeColor','none','FaceColor','flat')
    patch('Vertices',[[sillstarttime;sillstarttime;Parameters(ii).TimingParameters.runtime_yr;Parameters(ii).TimingParameters.runtime_yr],...
        (Parameters(ii).SillParameters.sillcentx-TimeSeries(ii).X_gl(1)+.5*Parameters(ii).SillParameters.sillwidth*[-1;1;1;-1])/1e3],...
        'Faces',1:4,'FaceVertexCData',[.5,.5,.5],'EdgeColor','none','FaceColor','flat')
    plot([sillstarttime,Parameters(ii).TimingParameters.runtime_yr],(Parameters(ii).SillParameters.sillcentx-TimeSeries(ii).X_gl(1))*[1,1]/1e3,...
        'Color','k','LineWidth',1)
    % Plot sill construction time:
    patch('Vertices',[[sillstarttime;sillstarttime;sillendtime;sillendtime],[ylims_gl(1);ylims_gl(2);ylims_gl(2);ylims_gl(1)]],...
        'Faces',1:4,'FaceVertexCData',[.5,.5,.5],'EdgeColor','none','FaceColor','flat')
    % Plot grounding line:
    h1=plot(TimeSeries(ii).Time(TimeSeries(ii).X_gl~=0),(TimeSeries(ii).X_gl(TimeSeries(ii).X_gl~=0)-TimeSeries(ii).X_gl(1))/1e3,'Color','g','LineWidth',2);
    h2=plot(TimeSeries(thiswarmingrun).Time(TimeSeries(thiswarmingrun).X_gl~=0),(TimeSeries(thiswarmingrun).X_gl(TimeSeries(thiswarmingrun).X_gl~=0)-TimeSeries(thiswarmingrun).X_gl(1))/1e3,'Color','r','LineWidth',2);
    % Plot model crash point:
    lastind=find(TimeSeries(ii).X_gl~=0,1,'last');
    if lastind~=numrecords
        plot(TimeSeries(ii).Time(lastind),(TimeSeries(ii).X_gl(lastind)-TimeSeries(ii).X_gl(1))/1e3,'Color','k','MarkerFaceColor','g','Marker','p','MarkerSize',5)
    end
    lastind=find(TimeSeries(thiswarmingrun).X_gl~=0,1,'last');
    if lastind~=numrecords
        plot(TimeSeries(thiswarmingrun).Time(lastind),(TimeSeries(thiswarmingrun).X_gl(lastind)-TimeSeries(thiswarmingrun).X_gl(1))/1e3,'Color','k','MarkerFaceColor','r','Marker','p','MarkerSize',5)
    end
    % Plot GL delay:
    if exist('gllevel','var')
        plot([warmingtime_gl,silltime_gl],(gllevel-TimeSeries(ii).X_gl(1))*[1,1]/1e3,...
            'Color','k','MarkerFaceColor','k','Marker','o','MarkerSize',3,'LineWidth',1)
        text(mean([warmingtime_gl,silltime_gl]),(gllevel-TimeSeries(ii).X_gl(1))/1e3,...
            [num2str(round(RetreatDelay(sillruncounter))),' yr'],'Color','k','VerticalAlignment','Bottom','HorizontalAlignment','Center')
        clear warmingtime_gl silltime_gl gllevel
    elseif isnan(RetreatDelay(sillruncounter))==0
        h3=plot(TimeSeries(ii).Time(TimeSeries(ii).X_gl~=0)-RetreatDelay(sillruncounter),(TimeSeries(ii).X_gl(TimeSeries(ii).X_gl~=0)-TimeSeries(ii).X_gl(1))/1e3,'Color','g','LineWidth',1,'LineStyle','--');
    end
    % Plot max prevented GL retreat:
    if MaxPreventedGLRetreat(sillruncounter)>0
        thisind=find(abs(TimeSeries(ii).X_gl-TimeSeries(thiswarmingrun).X_gl-MaxPreventedGLRetreat(sillruncounter))/MaxPreventedGLRetreat(sillruncounter)<1e-14,1,'first');
        plot(TimeSeries(ii).Time(thisind)*[1,1],...
            [TimeSeries(ii).X_gl(thisind)-TimeSeries(ii).X_gl(1),TimeSeries(thiswarmingrun).X_gl(thisind)-TimeSeries(thiswarmingrun).X_gl(1)]/1e3,...
            'Color','k','MarkerFaceColor','k','Marker','o','MarkerSize',3,'LineWidth',1)
        text(TimeSeries(ii).Time(thisind),mean([TimeSeries(ii).X_gl(thisind)-TimeSeries(ii).X_gl(1),max(ylims_gl(1)*1e3,TimeSeries(thiswarmingrun).X_gl(thisind)-TimeSeries(thiswarmingrun).X_gl(1))]/1e3),...
            [num2str(round(MaxPreventedGLRetreat(sillruncounter)/1e3)),' km'],'Color','k','VerticalAlignment','Bottom','HorizontalAlignment','Center','Rotation',-90)
    end
    % Plot immediate GL recovery:
    if GLRecovery(sillruncounter)~=0
        glminind=sillstartind-1+find(abs(TimeSeries(ii).X_gl(sillstartind:glmaxind)-glmin)/glmin<1e-14,1,'first');
        glmaxtime=TimeSeries(ii).Time(glmaxind);
        glmintime=TimeSeries(ii).Time(glminind);
        plot([glmintime,glmaxtime],(glmin-TimeSeries(ii).X_gl(1))*[1,1]/1e3,'Color',[.5,.5,.5],'LineWidth',.5)
        plot(glmaxtime*[1,1],[TimeSeries(ii).X_gl(glmaxind)-TimeSeries(ii).X_gl(1),glmin-TimeSeries(ii).X_gl(1)]/1e3,...
            'Color','k','MarkerFaceColor','k','Marker','o','MarkerSize',3,'LineWidth',1)
        text(glmaxtime,mean([TimeSeries(ii).X_gl(glmaxind)-TimeSeries(ii).X_gl(1),glmin-TimeSeries(ii).X_gl(1)]/1e3),...
            [num2str(round(GLRecovery(sillruncounter)/1e3)),' km'],'Color','k','VerticalAlignment','Bottom','HorizontalAlignment','Center','Rotation',90)
    end
    % Fix up axes and title:
    xlim([0,Parameters(ii).TimingParameters.runtime_yr])
    set(gca,'XTickLabel',[])
    ylim(ylims_gl)
    ylabel('Change in Position (km)','FontSize',fontsize)
    set(gca,'FontSize',fontsize)
    if strcmp(calving,'ch')
        title({['Flowband: ',flowband,', sliding: m=',num2str(M(sillruncounter)),', calving=c(H), Sill Type: ',num2str(SillKey(sillruncounter))];'Grounding Line Measurements'},'interpreter','none','FontSize',fontsize)
    elseif strcmp(calving,'chu')
        title({['Flowband: ',flowband,', sliding: m=',num2str(M(sillruncounter)),', calving=c(H,u), Sill Type: ',num2str(SillKey(sillruncounter))];'Grounding Line Measurements'},'interpreter','none','FontSize',fontsize)
    elseif strcmp(calving,'chm')
        title({['Flowband: ',flowband,', sliding: m=',num2str(M(sillruncounter)),', calving=c(H,m), Sill Type: ',num2str(SillKey(sillruncounter))];'Grounding Line Measurements'},'interpreter','none','FontSize',fontsize)
    else
        title({['Flowband: ',flowband,', sliding: m=',num2str(M(sillruncounter)),', Sill Type: ',num2str(SillKey(sillruncounter))];'Grounding Line Measurements'},'interpreter','none','FontSize',fontsize)
    end
    % Create legend:
    if exist('h3','var')
        legend([h1;h2;h3],'Sill Run','Warming Run',['Sill Run - ',num2str(round(RetreatDelay(sillruncounter))),' yr'],'location','Best')
    else
        legend([h1;h2],'Sill Run','Warming Run','location','Best')
    end
    % Call second subplot:
    subplot('Position',Boxes{2})
    % Plot zero line:
    hold off
    plot([0,Parameters(ii).TimingParameters.runtime_yr],[0,0],'--k')
    hold on
    % Compute ylims:
    ylims_sl=[min([TimeSeries(ii).VAF(TimeSeries(ii).VAF~=0)-TimeSeries(ii).VAF(1);TimeSeries(thiswarmingrun).VAF(TimeSeries(thiswarmingrun).VAF~=0)-TimeSeries(thiswarmingrun).VAF(1)]),...
        max([TimeSeries(ii).VAF(TimeSeries(ii).VAF~=0)-TimeSeries(ii).VAF(1);TimeSeries(thiswarmingrun).VAF(TimeSeries(thiswarmingrun).VAF~=0)-TimeSeries(thiswarmingrun).VAF(1)])]*100*vafslfactor;
    % Plot sill construction time:
    patch('Vertices',[[sillstarttime;sillstarttime;sillendtime;sillendtime],[ylims_sl(1);ylims_sl(2);ylims_sl(2);ylims_sl(1)]],...
        'Faces',1:4,'FaceVertexCData',[.5,.5,.5],'EdgeColor','none','FaceColor','flat')
    % Plot volume above flotation:
    plot(TimeSeries(ii).Time(TimeSeries(ii).VAF~=0),(TimeSeries(ii).VAF(TimeSeries(ii).VAF~=0)-TimeSeries(ii).VAF(1))*100*vafslfactor,'Color','g','LineWidth',2)
    plot(TimeSeries(thiswarmingrun).Time(TimeSeries(thiswarmingrun).VAF~=0),(TimeSeries(thiswarmingrun).VAF(TimeSeries(thiswarmingrun).VAF~=0)-TimeSeries(thiswarmingrun).VAF(1))*100*vafslfactor,'Color','r','LineWidth',2)
    % Plot model crash point:
    lastind=find(TimeSeries(ii).VAF~=0,1,'last');
    if lastind~=numrecords
        plot(TimeSeries(ii).Time(lastind),(TimeSeries(ii).VAF(lastind)-TimeSeries(ii).VAF(1))*100*vafslfactor,'Color','k','MarkerFaceColor','g','Marker','p','MarkerSize',5)
    end
    lastind=find(TimeSeries(thiswarmingrun).VAF~=0,1,'last');
    if lastind~=numrecords
        plot(TimeSeries(thiswarmingrun).Time(lastind),(TimeSeries(thiswarmingrun).VAF(lastind)-TimeSeries(thiswarmingrun).VAF(1))*100*vafslfactor,'Color','k','MarkerFaceColor','r','Marker','p','MarkerSize',5)
    end
    % Plot max prevented sea level rise:
    if MaxPreventedSeaLevelRise(sillruncounter)>0
        thisind=find(abs((TimeSeries(ii).VAF-TimeSeries(thiswarmingrun).VAF)*vafslfactor-MaxPreventedSeaLevelRise(sillruncounter))/MaxPreventedSeaLevelRise(sillruncounter)<1e-14,1,'first');
        plot(TimeSeries(ii).Time(thisind)*[1,1],...
            [TimeSeries(ii).VAF(thisind)-TimeSeries(ii).VAF(1),TimeSeries(thiswarmingrun).VAF(thisind)-TimeSeries(thiswarmingrun).VAF(1)]*100*vafslfactor,...
            'Color','k','MarkerFaceColor','k','Marker','o','MarkerSize',3,'LineWidth',1)
        lastpow10=10^floor(log10(MaxPreventedSeaLevelRise(sillruncounter)));
        roundedslrise=lastpow10*.1*round(10*MaxPreventedSeaLevelRise(sillruncounter)/lastpow10); % round to 2 sig figs
        text(TimeSeries(ii).Time(thisind),mean([TimeSeries(ii).VAF(thisind)-TimeSeries(ii).VAF(1),max(ylims_sl(1)/(100*vafslfactor),TimeSeries(thiswarmingrun).VAF(thisind)-TimeSeries(thiswarmingrun).VAF(1))]*100*vafslfactor),...
            [num2str(100*roundedslrise),' cm'],'Color','k','VerticalAlignment','Bottom','HorizontalAlignment','Center','Rotation',-90)
    end
    % Plot sea level recovery:
    if SeaLevelRecovery(sillruncounter)~=0
        slminind=sillstartind-1+find(abs(TimeSeries(ii).VAF(sillstartind:slmaxind)*vafslfactor-slmin)/slmin<1e-14,1,'first');
        slmaxtime=TimeSeries(ii).Time(slmaxind);
        slmintime=TimeSeries(ii).Time(slminind);
        plot([slmintime,slmaxtime],(slmin/vafslfactor-TimeSeries(ii).VAF(1))*[1,1]*100*vafslfactor,'Color',[.5,.5,.5],'LineWidth',.5)
        plot(slmaxtime*[1,1],[TimeSeries(ii).VAF(slmaxind)-TimeSeries(ii).VAF(1),slmin/vafslfactor-TimeSeries(ii).VAF(1)]*100*vafslfactor,...
            'Color','k','MarkerFaceColor','k','Marker','o','MarkerSize',3,'LineWidth',1)
        lastpow10=10^floor(log10(SeaLevelRecovery(sillruncounter)));
        roundedslrecovery=lastpow10*.1*round(10*SeaLevelRecovery(sillruncounter)/lastpow10); % round to 2 sig figs
        text(slmaxtime,mean([TimeSeries(ii).VAF(slmaxind)-TimeSeries(ii).VAF(1),slmin/vafslfactor-TimeSeries(ii).VAF(1)]*100*vafslfactor),...
            [num2str(roundedslrecovery*100),' cm'],'Color','k','VerticalAlignment','Bottom','HorizontalAlignment','Center','Rotation',90)
    end
    % Make a note of correction factor:
    if strcmp(flowband(end),'B')
        text(.5*Parameters(ii).TimingParameters.runtime_yr,ylims_sl(1),'Note: ice volume corrected for difference in area between narrow and wide flowbands',...
            'Color','k','VerticalAlignment','Bottom','HorizontalAlignment','Center')
    end
    % Fix up axes and title:
    xlim([0,Parameters(ii).TimingParameters.runtime_yr])
    xlabel('Time (yr)','FontSize',fontsize)
    ylim(ylims_sl)
    ylabel('Ice Volume Change (cm sleq)','FontSize',fontsize)
    title('Volume Above Flotation Measurements','FontSize',fontsize)
    set(gca,'FontSize',fontsize)
    % Save figure:
    set(gcf,'PaperSize',pagesize)
    set(gcf,'PaperPosition',[0,0,pagesize])
    figname=[figfolder,TimeSeries(ii).filename(1:end-4),'_',figsuffix,'.png'];
    print('-dpng',figname,['-r',num2str(resolution)])
    
    
    % Advance counter:
    sillruncounter=sillruncounter+1;
    
end

%% Save Outputs:

% Save results to matfile:
save(outputfile,'SillRunFileName','WarmingRunFileName','SillRunParameters',...
    'Flowband','M','Calving','SillKey','Isv3a',...
    '*PreventedSeaLevelRise','*PreventedGLRetreat','RetreatDelay','*Recovery','*RunCollapse','SillBuiltUnderGroundedIce',...
    'NOTE_*','oceanarea','icbuffer','minretreat','minretreatratio')


%% Create main experiment text file:
% Generate filename:
thisfilename=[outputfile(1:end-4),'_mainexperiment.csv'];
% Generate column headers:
TableHeaders={'Flowband','Sliding exponent','Calving rule','Max prevented SLR (cm)','Max prevented GL retreat (km)',...
    'Integrated prevented SLR (cm*yr)','Integrated prevented GL retreat (km*yr)','Initial SLR recovery (cm)','Initial GL recovery (km)',...
    'Retreat delay (yr)','Warming run collapse?','Sill run collapse?'};
% Generate table body as cell array:
% Pre-allocate:
TableBody=cell(sum(Isv3a==0),12);
% Generate each row of the table:
rowcounter=1;
for ii=1:numsillruns
    % Check whether to add this run to the table:
    if Isv3a(ii) || strcmp(Flowband{ii}(1:end-1),'Jakobshavn')
        continue
    end
    % Add experiment description:
    TableBody{rowcounter,1}=Flowband{ii};
    TableBody{rowcounter,2}=M(ii);
    if strcmp(Calving{ii},'ch')
        TableBody{rowcounter,3}='c(H)';
    elseif strcmp(Calving{ii},'chm')
        TableBody{rowcounter,3}='c(H,m)';
    elseif strcmp(Calving{ii},'chu')
        TableBody{rowcounter,3}='c(H,u)';
    end
    % Add numerical experiment results: (2 sig figs)
    TableBody{rowcounter,4}=roundsd(MaxPreventedSeaLevelRise(ii)*100,2);
    TableBody{rowcounter,5}=roundsd(MaxPreventedGLRetreat(ii)/1e3,2);
    TableBody{rowcounter,6}=roundsd(IntegratedPreventedSeaLevelRise(ii)*100,2);
    TableBody{rowcounter,7}=roundsd(IntegratedPreventedGLRetreat(ii)/1e3,2);
    TableBody{rowcounter,8}=roundsd(SeaLevelRecovery(ii)*100,2);
    TableBody{rowcounter,9}=roundsd(GLRecovery(ii)/1e3,2);
    % Deal with < and > symbols for delays: 
    if WarmingRunCollapse(ii) && SillRunCollapse(ii)==0
        TableBody{rowcounter,10}=['>',num2str(roundsd(RetreatDelay(ii),2))];
    elseif WarmingRunCollapse(ii)==0 && SillRunCollapse(ii)
        TableBody{rowcounter,10}=['<',num2str(roundsd(RetreatDelay(ii),2))];
    else
        TableBody{rowcounter,10}=roundsd(RetreatDelay(ii),2);
    end
    % Flag collapses:
    if WarmingRunCollapse(ii)
        TableBody{rowcounter,11}='y';
    else
        TableBody{rowcounter,11}='n';
    end
    if SillRunCollapse(ii)
        TableBody{rowcounter,12}='y';
    else
        TableBody{rowcounter,12}='n';
    end
    % Advance rowcounter:
    rowcounter=rowcounter+1;
end
% Combine to form table:
Table=[TableHeaders;TableBody(1:rowcounter-1,:)];
% Create file:
cell2csv(thisfilename,Table,';')


% Create rapid sill experiment text file:
% Generate filename:
thisfilename=[outputfile(1:end-4),'_rapidsillexperiment.csv'];
% Generate column headers:
TableHeaders={'Flowband','Sliding exponent','Calving rule','Max prevented SLR (cm)','Max prevented GL retreat (km)',...
    'Integrated prevented SLR (cm*yr)','Integrated prevented GL retreat (km*yr)','Initial SLR recovery (cm)','Initial GL recovery (km)',...
    'Retreat delay (yr)','Warming run collapse?','Sill run collapse?'};
% Generate table body as cell array:
% Pre-allocate:
TableBody=cell(6,12);
% Generate each row of the table:
rowcounter=1;
for ii=1:numsillruns
    % Check whether to add this run to the table:
    if Isv3a(ii)==0 || (strcmp(Flowband{ii}(1:end-1),'Helheim')==0 && strcmp(Flowband{ii}(1:end-1),'Kanger')==0)
        continue
    end
    % Add experiment description:
    TableBody{rowcounter,1}=Flowband{ii};
    TableBody{rowcounter,2}=M(ii);
    if strcmp(Calving{ii},'ch')
        TableBody{rowcounter,3}='c(H)';
    elseif strcmp(Calving{ii},'chm')
        TableBody{rowcounter,3}='c(H,m)';
    elseif strcmp(Calving{ii},'chu')
        TableBody{rowcounter,3}='c(H,u)';
    end
    % Add numerical experiment results: (2 sig figs)
    TableBody{rowcounter,4}=roundsd(MaxPreventedSeaLevelRise(ii)*100,2);
    TableBody{rowcounter,5}=roundsd(MaxPreventedGLRetreat(ii)/1e3,2);
    TableBody{rowcounter,6}=roundsd(IntegratedPreventedSeaLevelRise(ii)*100,2);
    TableBody{rowcounter,7}=roundsd(IntegratedPreventedGLRetreat(ii)/1e3,2);
    TableBody{rowcounter,8}=roundsd(SeaLevelRecovery(ii)*100,2);
    TableBody{rowcounter,9}=roundsd(GLRecovery(ii)/1e3,2);
    % Deal with < and > symbols for delays: 
    if WarmingRunCollapse(ii) && SillRunCollapse(ii)==0
        TableBody{rowcounter,10}=['>',num2str(roundsd(RetreatDelay(ii),2))];
    elseif WarmingRunCollapse(ii)==0 && SillRunCollapse(ii)
        TableBody{rowcounter,10}=['<',num2str(roundsd(RetreatDelay(ii),2))];
    else
        TableBody{rowcounter,10}=roundsd(RetreatDelay(ii),2);
    end
    % Flag collapses:
    if WarmingRunCollapse(ii)
        TableBody{rowcounter,11}='y';
    else
        TableBody{rowcounter,11}='n';
    end
    if SillRunCollapse(ii)
        TableBody{rowcounter,12}='y';
    else
        TableBody{rowcounter,12}='n';
    end
    % Advance rowcounter:
    rowcounter=rowcounter+1;
end
% Combine to form table:
Table=[TableHeaders;TableBody(1:rowcounter-1,:)];
% Create file:
cell2csv(thisfilename,Table,';')


% Create sill erosion experiment text file:
% Generate filename:
thisfilename=[outputfile(1:end-4),'_sillerosionexperiment.csv'];
% Generate column headers:
TableHeaders={'Flowband','Sliding exponent','Calving rule','Sill strength','Max prevented SLR (cm)','Max prevented GL retreat (km)',...
    'Integrated prevented SLR (cm*yr)','Integrated prevented GL retreat (km*yr)','Initial SLR recovery (cm)','Initial GL recovery (km)',...
    'Retreat delay (yr)','Warming run collapse?','Sill run collapse?'};
% Generate table body as cell array:
% Pre-allocate:
TableBody=cell(sum(SillKey~=0),13);
% Generate each row of the table:
rowcounter=1;
for ii=1:numsillruns
    % Check whether to add this run to the table:
    if SillKey(ii)==0 || SillBuiltUnderGroundedIce(ii)
        continue
    end
    % Add experiment description:
    TableBody{rowcounter,1}=Flowband{ii};
    TableBody{rowcounter,2}=M(ii);
    if strcmp(Calving{ii},'ch')
        TableBody{rowcounter,3}='c(H)';
    elseif strcmp(Calving{ii},'chm')
        TableBody{rowcounter,3}='c(H,m)';
    elseif strcmp(Calving{ii},'chu')
        TableBody{rowcounter,3}='c(H,u)';
    end
    if SillKey(ii)==1
        TableBody{rowcounter,4}='weak';
    elseif SillKey(ii)==2
        TableBody{rowcounter,4}='medium';
    elseif SillKey(ii)==3
        TableBody{rowcounter,4}='strong';
    end
    % Add numerical experiment results: (2 sig figs)
    TableBody{rowcounter,5}=roundsd(MaxPreventedSeaLevelRise(ii)*100,2);
    TableBody{rowcounter,6}=roundsd(MaxPreventedGLRetreat(ii)/1e3,2);
    TableBody{rowcounter,7}=roundsd(IntegratedPreventedSeaLevelRise(ii)*100,2);
    TableBody{rowcounter,8}=roundsd(IntegratedPreventedGLRetreat(ii)/1e3,2);
    TableBody{rowcounter,9}=roundsd(SeaLevelRecovery(ii)*100,2);
    TableBody{rowcounter,10}=roundsd(GLRecovery(ii)/1e3,2);
    % Deal with < and > symbols for delays: 
    if WarmingRunCollapse(ii) && SillRunCollapse(ii)==0
        TableBody{rowcounter,11}=['>',num2str(roundsd(RetreatDelay(ii),2))];
    elseif WarmingRunCollapse(ii)==0 && SillRunCollapse(ii)
        TableBody{rowcounter,11}=['<',num2str(roundsd(RetreatDelay(ii),2))];
    else
        TableBody{rowcounter,11}=roundsd(RetreatDelay(ii),2);
    end
    % Flag collapses:
    if WarmingRunCollapse(ii)
        TableBody{rowcounter,12}='y';
    else
        TableBody{rowcounter,12}='n';
    end
    if SillRunCollapse(ii)
        TableBody{rowcounter,13}='y';
    else
        TableBody{rowcounter,13}='n';
    end
    % Advance rowcounter:
    rowcounter=rowcounter+1;
end
% Combine to form table:
Table=[TableHeaders;TableBody(1:rowcounter-1,:)];
% Create file:
cell2csv(thisfilename,Table,';')


%% Final display:
disp('Done!')
toc


%% Cutting Room Floor:

% Point-wise measurement of delay:

% % Pre-allocate levels:
% VAFlevels=linspace(min(TimeSeries(thiswarmingrun).VAF(TimeSeries(thiswarmingrun).VAF~=0)),...
%     max(TimeSeries(thiswarmingrun).VAF),numsamples_delay)';
% GLlevels=linspace(min(TimeSeries(thiswarmingrun).X_gl(TimeSeries(thiswarmingrun).X_gl~=0)),...
%     max(TimeSeries(thiswarmingrun).X_gl),numsamples_delay)';
% % Check whether delays are meaningful: (was there enough retreat? was it losing mass at the end?)
% if TimeSeries(ii).VAF(end)>TimeSeries(ii).VAF(1) || TimeSeries(thiswarmingrun).VAF(end)>TimeSeries(thiswarmingrun).VAF(1) || TimeSeries(ii).VAF(end)>TimeSeries(ii).VAF(end-1)
%     RetreatDelay(sillruncounter)=NaN;
% end
% if TimeSeries(ii).X_gl(end)>TimeSeries(ii).X_gl(1)-.1*(TimeSeries(thiswarmingrun).X_gl(1)-TimeSeries(thiswarmingrun).X_gl(end)) || TimeSeries(thiswarmingrun).X_gl(end)>TimeSeries(thiswarmingrun).X_gl(1)
%     RetreatDelay(sillruncounter)=NaN;
% end
% % Loop through levels:
% for jj=1:numsamples_delay
%     
%     % Locate last crossing of these levels in sill run:
%     if isnan(RetreatDelay(sillruncounter))==0
%         sillind_vaf=find((TimeSeries(ii).VAF(1:end-1)<=VAFlevels(jj)&TimeSeries(ii).VAF(2:end)>VAFlevels(jj))|...
%             (TimeSeries(ii).VAF(1:end-1)>VAFlevels(jj)&TimeSeries(ii).VAF(2:end)<=VAFlevels(jj)),1,'last');
%         if sillind_vaf<=sillstartind
%             sillind_vaf=[];
%         end
%     else
%         sillind_vaf=[];
%     end
%     if isnan(RetreatDelay(sillruncounter))==0
%         sillind_gl=find((TimeSeries(ii).X_gl(1:end-1)<=GLlevels(jj)&TimeSeries(ii).X_gl(2:end)>GLlevels(jj))|...
%             (TimeSeries(ii).X_gl(1:end-1)>GLlevels(jj)&TimeSeries(ii).X_gl(2:end)<=GLlevels(jj)),1,'last');
%         if sillind_gl<=sillstartind
%             sillind_gl=[];
%         end
%     else
%         sillind_gl=[];
%     end
%     % Interpolate times of last crossing in sill run:
%     if isempty(sillind_vaf)==0
%         silltime_vaf=TimeSeries(ii).Time(sillind_vaf)+dt*(VAFlevels(jj)-TimeSeries(ii).VAF(sillind_vaf))/(TimeSeries(ii).VAF(sillind_vaf+1)-TimeSeries(ii).VAF(sillind_vaf));
%     end
%     if isempty(sillind_gl)==0
%         silltime_gl=TimeSeries(ii).Time(sillind_gl)+dt*(GLlevels(jj)-TimeSeries(ii).X_gl(sillind_gl))/(TimeSeries(ii).X_gl(sillind_gl+1)-TimeSeries(ii).X_gl(sillind_gl));
%     end
%     
%     % Locate last crossing of these levels in warming run:
%     if isnan(RetreatDelay(sillruncounter))==0
%         warmingind_vaf=find((TimeSeries(thiswarmingrun).VAF(1:end-1)<=VAFlevels(jj)&TimeSeries(thiswarmingrun).VAF(2:end)>VAFlevels(jj))|...
%             (TimeSeries(thiswarmingrun).VAF(1:end-1)>VAFlevels(jj)&TimeSeries(thiswarmingrun).VAF(2:end)<=VAFlevels(jj)),1,'last');
%         if warmingind_vaf<=sillstartind
%             warmingind_vaf=[];
%         end
%     else
%         warmingind_vaf=[];
%     end
%     if isnan(RetreatDelay(sillruncounter))==0
%         warmingind_gl=find((TimeSeries(thiswarmingrun).X_gl(1:end-1)<=GLlevels(jj)&TimeSeries(thiswarmingrun).X_gl(2:end)>GLlevels(jj))|...
%             (TimeSeries(thiswarmingrun).X_gl(1:end-1)>GLlevels(jj)&TimeSeries(thiswarmingrun).X_gl(2:end)<=GLlevels(jj)),1,'last');
%         if warmingind_gl<=sillstartind
%             warmingind_gl=[];
%         end
%     else
%         warmingind_gl=[];
%     end
%     % Interpolate times of last crossing in warming run:
%     if isempty(warmingind_vaf)==0
%         warmingtime_vaf=TimeSeries(thiswarmingrun).Time(warmingind_vaf)+dt*(VAFlevels(jj)-TimeSeries(thiswarmingrun).VAF(warmingind_vaf))/(TimeSeries(thiswarmingrun).VAF(warmingind_vaf+1)-TimeSeries(thiswarmingrun).VAF(warmingind_vaf));
%     end
%     if isempty(warmingind_gl)==0
%         warmingtime_gl=TimeSeries(thiswarmingrun).Time(warmingind_gl)+dt*(GLlevels(jj)-TimeSeries(thiswarmingrun).X_gl(warmingind_gl))/(TimeSeries(thiswarmingrun).X_gl(warmingind_gl+1)-TimeSeries(thiswarmingrun).X_gl(warmingind_gl));
%     end
%     
%     % Compute delays:
%     if isempty(sillind_vaf)==0 && isempty(warmingind_vaf)==0
%         RetreatDelay(sillruncounter)=max(RetreatDelay(sillruncounter),silltime_vaf-warmingtime_vaf);
%         if RetreatDelay(sillruncounter)==silltime_vaf-warmingtime_vaf
%             finalsilltime_vaf=silltime_vaf;
%             finalwarmingtime_vaf=warmingtime_vaf;
%             finalvaflevel=VAFlevels(jj);
%         end
%     end
%     if isempty(sillind_gl)==0 && isempty(warmingind_gl)==0
%         RetreatDelay(sillruncounter)=max(RetreatDelay(sillruncounter),silltime_gl-warmingtime_gl);
%         if RetreatDelay(sillruncounter)==silltime_gl-warmingtime_gl
%             finalsilltime_gl=silltime_gl;
%             finalwarmingtime_gl=warmingtime_gl;
%             finalgllevel=GLlevels(jj);
%         end
%     end
%     
% end