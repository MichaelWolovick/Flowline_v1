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

% Ice is assumed to melt into fresh water.

% The script assumes that all of the time series have the same dt.

clear all
tic

%% Parameters:

% Files and folders:
timeseriesfile='/net/mjw/FjordSillPaper/ModelOutput/ParamSweepRuns_v3_v3a.mat';
inputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/';
outputfile='/net/mjw/FjordSillPaper/ModelOutput/PreventedSeaLevelRise_v1.mat';

% Surface area of the ocean:
oceanarea=3.6e14;                   % m^2


%% Work:

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
PreventedSeaLevelRise=zeros(numsillruns,1);
PreventedGLRetreat=zeros(numsillruns,1);

% Make notes:
NOTE_units='Prevented sea level rise is in meter-years (meters of eustatic sea level), prevented GL retreat is also in meter-years, but meters of GL position.';
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
    % Compute prevented sea level rise:
    PreventedSeaLevelRise(sillruncounter)=(Parameters(ii).ThermalParameters.rho_i/Parameters(ii).ThermalParameters.rho_fw)*...
        correctionfactor*dt*sum(TimeSeries(ii).VAF-TimeSeries(thiswarmingrun).VAF)/oceanarea;
    % Compute prevented grounding line retreat:
    PreventedGLRetreat(sillruncounter)=dt*sum(TimeSeries(ii).X_gl-TimeSeries(thiswarmingrun).X_gl);
    % Advance counter:
    sillruncounter=sillruncounter+1;
end

% Save results:
save(outputfile,'SillRunFileName','WarmingRunFileName','SillRunParameters',...
    'Flowband','M','Calving','SillKey','Isv3a',...
    'PreventedSeaLevelRise','PreventedGLRetreat','NOTE_*','oceanarea')

% Final display:
disp('Done!')
toc