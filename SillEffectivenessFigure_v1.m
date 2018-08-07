% SillEffectivenessFigure_v1

% Mike Wolovick, 1/28/2017

% This figure shows the effectiveness of multiple sill scenarios.  Each
% scenario gets a bar plot showing the proportion of model runs that
% regrounded (not including never retreated).  The runs that are
% regrounded are colored according to the sea level rise rate relative to 
% the rate in the warming run.

% Relative sea level rise rate:

% RelRiseRate= lossratesill / lossratewarming

% This figure shows the probability of regrounding the ice shelf as a
% function of sill volume for several different values of the angle of
% repose.

% This figure also shows a delay time for sea level rise.  That delay time
% is calculated in the same way that GL retreat delay was calculated in
% QuantifySillEffect_v2

% Is it ligit to say that, for the experiment that blocks half of the warm
% water, that half the sill length was built?  Technically the 50% blockage
% and 100% blockage experiments have the same sill geoemtry, just different
% assumptions about how good the sill will be at blocking warm water.  I 
% think it's a fair assumption to say that blocking more water will
% always require more material, the problem is that the choice of a factor
% of two for the difference is arbitrary.  

clearvars
close all
tic

%% Parameters:

% Files and folders:
inputfile='/home/wolovick/Dropbox/FjordSillPaper/ModelOutput/ParamSweepRuns_v3_v3a.mat';
figname='/home/wolovick/Dropbox/FjordSillPaper/Figures/SillEffectivenessFigure_v2a.png';

% Experiment description:
flowbandnames={'ThwaitesA_v2';...
    'ThwaitesC_v2'}; 
mkeys={'m01';'m03';'m10'};
ckeys={'ch';'chu';'chm'};
scenariokeys={'s5';'s4';'s6';'s'};  % in increasing order of sill volume

% Scenario names:
% scenarionames={{'0% Water';'Blockage,';'Small Sill'},{'50% Water';'Blockage,';'Small Sill'},...
%     {'100% Water';'Blockage,';'Small Sill'},{'100% Water';'Blockage,';'Large Sill'}};
scenarionames={{'Isolated';'Pinning';'Points'},{'50% Water';'Blockage,';'Small Sill'},...
    {'100% Water';'Blockage,';'Small Sill'},{'100% Water';'Blockage,';'Large Sill'}};

% Analysis settings:   
compare2warming=1;              % logical
mingljump=25e3;                 % m
lossratewindow=30;              % yr

% Ticks and limits:
lossratelims=[-20,20];          % [1x2] percent
lossratetick=10;                % percent
percenttick=25;                 % percent

% Other Aesthetic settings:
neverregroundedcolor=[.5,.5,.5];% [R,G,B]
fontsize=14;

% Margins:
verttextbuffer=.14;             % unitless
titlebuffer=.06;                % unitless
horztextbuffer=.115;            % unitless
horznotextbuffer=.118;          % unitless

% Output size and resolution:
pagesize=[8,5];                 % [1x2] inches
resolution=300;                 % dpi

%% Compute Sill Effectiveness:

% Load input file:
load(inputfile)

% Compute numbers of things:
numflowbands=length(flowbandnames);
numm=length(mkeys);
numc=length(ckeys);
numsamples=numflowbands*numm*numc;
numruns=length(TimeSeries);
numscenarios=length(scenariokeys);
if length(scenarionames)~=numscenarios
    error('"scenarionames" must be same length as "scenariokeys"')
end

% Compute dt:
dt=TimeSeries(1).Time(2)-TimeSeries(1).Time(1);
numrecords=length(TimeSeries(1).Time);

% Compute loss rate window in indices:
lossratewindow_ind=round(lossratewindow/dt);

% Pre-allocate:
Regrounded=false(numsamples,numscenarios);
NeverRetreated=false(numsamples,numscenarios);
RelRiseRate=zeros(numsamples,numscenarios);
CouldHaveRegrounded=false(numsamples,numscenarios);

% Compute regrounding fraction for each sill:
% Set the counters:
fcounter=1;
mcounter=1;
ccounter=1;
% Loop through samples:
for thissample=1:numsamples
    
    % Identify these key strings:
    thisflowband=flowbandnames{fcounter};
    thism=mkeys{mcounter};
    thisc=ckeys{ccounter};
    % Advance counters:
    fcounter=fcounter+1;
    if fcounter==length(flowbandnames)+1
        fcounter=1;
        mcounter=mcounter+1;
        if mcounter==length(mkeys)+1
            mcounter=1;
            ccounter=ccounter+1;
            if ccounter>length(ckeys) && thissample<numsamples
                error('Problem with assigning the sample key strings')
            end
        end
    end
    
    % Loop through scenarios:
    for thisscenario=1:numscenarios
        
        % Locate the model run:
        thiswarmingrun=[];
        % Loop through all model runs:
        for thisrun=1:numruns
            % Parse filename:
            % Find underscores:
            underscorenums=strfind(TimeSeries(thisrun).filename,'_');
            % Check flowband:
            if strcmp(TimeSeries(thisrun).filename(1:underscorenums(2)-1),thisflowband)
                % Check sliding:
                if strcmp(TimeSeries(thisrun).filename(underscorenums(3)+1:underscorenums(4)-1),thism)
                    % Check calving:
                    if strcmp(TimeSeries(thisrun).filename(underscorenums(4)+1:underscorenums(5)-1),thisc)
                        % Check scenario:
                        if strcmp(TimeSeries(thisrun).filename(underscorenums(2)+1:underscorenums(3)-1),scenariokeys{thisscenario})
                            % This is the one!!!
                            break
                        elseif strcmp(TimeSeries(thisrun).filename(underscorenums(2)+1:underscorenums(3)-1),'w')
                            thiswarmingrun=thisrun;
                        end
                    end
                end
            end
            % Check that we found this model run:
            if thisrun==numruns
                error('Unable to locate this model run')
            end
            % Check that we located a warming run:
            if thisrun==numruns && isempty(thiswarmingrun)
                error('Unable to locate warming run')
            end
        end
        
        
        % Define sill starting point:
        sillstartx=Parameters(thisrun).SillParameters.sillcentx-2*Parameters(thisrun).SillParameters.sillwidth;
        sillstarttime=Parameters(thisrun).SillParameters.sillstarttime_yr;
        
        % Identify valid points:
        ValidPoints=TimeSeries(thisrun).Time~=0&TimeSeries(thisrun).X_gl~=0;
        
        % Check that grounding line was behind the sill when
        % construction began:
        if interp1(TimeSeries(thisrun).Time(ValidPoints),TimeSeries(thisrun).X_gl(ValidPoints),sillstarttime)>sillstartx
            % Move on to the next scenario:
            NeverRetreated(thissample,thisscenario)=1;
            RelRiseRate(thissample,thisscenario)=NaN;
            continue
        end
        
        % Check whether the glacier never even retreated:
        if min(TimeSeries(thisrun).X_gl(ValidPoints)-TimeSeries(thisrun).X_gl(1))>-mingljump
            % Move on to the next scenario:
            NeverRetreated(thissample,thisscenario)=1;
            RelRiseRate(thissample,thisscenario)=NaN;
            continue
            % Check whether the glacier retreated and then regrounded:
        elseif max(diff(TimeSeries(thisrun).X_gl(ValidPoints&TimeSeries(thisrun).Time>sillstarttime)))>mingljump  
            % Flag that is regrounded:
            Regrounded(thissample,thisscenario)=1;
        else % the glacier never regrounded
            % Check whether the glacier *could have* regrounded:
            if max(TimeSeries(thisrun).Domainwidth(ValidPoints&TimeSeries(thisrun).Time>sillstarttime))>=Parameters(thisrun).SillParameters.sillcentx
                CouldHaveRegrounded(thissample,thisscenario)=1;
            end
            % Move on to the next scenario:
            RelRiseRate(thissample,thisscenario)=NaN;
            continue
        end
        
        % Find various indices:
        lastind_warming=find(TimeSeries(thiswarmingrun).Time~=0,1,'last');
        lastind_sill=find(TimeSeries(thisrun).Time~=0,1,'last');
        regroundingind=find(diff(TimeSeries(thisrun).X_gl(ValidPoints))>mingljump&...
            TimeSeries(thisrun).Time(ValidPoints&[false(1);true(numrecords-1,1)])>sillstarttime,1,'first');  % index is just before regrounding
        
        % Compute loss rate for comparison:
        if compare2warming
            % Get linear best-fit:
            A=[ones(lossratewindow_ind+1,1),TimeSeries(thiswarmingrun).Time(lastind_warming-lossratewindow:lastind_warming)];
            equation=(A'*A)\A'*TimeSeries(thiswarmingrun).VAF(lastind_warming-lossratewindow:lastind_warming);
        else
            % Get linear best-fit:
            A=[ones(lossratewindow_ind+1,1),TimeSeries(thisrun).Time(regroundingind-lossratewindow:regroundingind)];
            equation=(A'*A)\A'*TimeSeries(thisrun).VAF(regroundingind-lossratewindow:regroundingind);
        end
        lossratebefore=-equation(2);
        
        % Compute loss rate after regrounding:
        B=[ones(lastind_sill-regroundingind,1),TimeSeries(thisrun).Time(regroundingind+1:lastind_sill)];
        equation=(B'*B)\B'*TimeSeries(thisrun).VAF(regroundingind+1:lastind_sill);
        lossratesill=-equation(2);
        
        % Compute relative change in loss rate:
        RelRiseRate(thissample,thisscenario)=lossratesill/lossratebefore;
        
    end
    
end

% Determine number of runs that retreated:
numretreated=sum(NeverRetreated(:,1)==0);
numneverretreated=numsamples-numretreated;

% Create dummy index vector:
DummyInd=linspace(1,numsamples,numsamples)';

% Loop through scenarios:
for ii=1:numscenarios
    % Sort by whether or not the run retreated:
    [NeverRetreated(:,ii),SortedInd]=sort(NeverRetreated(:,ii));
    RelRiseRate(:,ii)=RelRiseRate(SortedInd,ii);
    CouldHaveRegrounded(:,ii)=CouldHaveRegrounded(SortedInd,ii);
    Regrounded(:,ii)=Regrounded(SortedInd,ii);
    % Sort the runs that retreated by their relative rise rate:
    [RelRiseRate(1:numretreated,ii),SortedInd]=sort(RelRiseRate(1:numretreated,ii));
    CouldHaveRegrounded(1:numretreated,ii)=CouldHaveRegrounded(SortedInd,ii);
    Regrounded(1:numretreated,ii)=Regrounded(SortedInd,ii);
    % Sort this scenario:
    %RelRiseRate(:,ii)=[sort(RelRiseRate(isnan(RelRiseRate(:,ii))==0,ii),'ascend');NaN*ones(sum(isnan(RelRiseRate(:,ii))),1)];
    % Check the numretreated variable:
    if sum(NeverRetreated(:,ii)==0)~=numretreated
        error('Different values of numretreated in each scenario')
    end
end

%error('end here')

%% Make a Figure:

% Hardcoded plot arrangement:
numvertplots=1;
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

% Make figure:
figure(1)

% Make faded jet colormap:
FadedJet=colormap('jet');
FadedJet(3*round(size(FadedJet,1)/8)+1:round(size(FadedJet,1)/2),1)=linspace(0,1,round(size(FadedJet,1)/8))';
FadedJet(3*round(size(FadedJet,1)/8)+1:round(size(FadedJet,1)/2),3)=1;
FadedJet(round(size(FadedJet,1)/2)+1:5*round(size(FadedJet,1)/8),1)=1;
FadedJet(round(size(FadedJet,1)/2)+1:5*round(size(FadedJet,1)/8),3)=linspace(1,0,round(size(FadedJet,1)/8))';

% Create RGB matrix:
RGB=zeros(numretreated,numscenarios,3);

% Interpolate RGB values:
for ii=1:3
    % Interpolate:
    RGB(:,:,ii)=interp1(linspace(0,1,size(FadedJet,1))',FadedJet(:,ii),...
        max(0,min(1,(100*RelRiseRate(1:numretreated,:)-lossratelims(1))/(lossratelims(2)-lossratelims(1)))));
end

% Set RGB values of never regrounded runs:
for ii=1:numscenarios
    RGB(isnan(RelRiseRate(1:numretreated,ii)),ii,:)=...
        repmat(permute(neverregroundedcolor,[1,3,2]),[sum(isnan(RelRiseRate(1:numretreated,ii))),1,1]);
end

% Create dummy axes vectors:
DummyX=linspace(.5/numscenarios,1-.5/numscenarios,numscenarios);
DummyY=100*linspace(.5/numretreated,1-.5/numretreated,numretreated)';

% Create subplot:
subplot('Position',Boxes{1})

% Show image:
image(DummyX,DummyY,RGB)
set(gca,'ydir','normal')

% Set up axes, etc:
xlim([0,1])
ylim([0,100])
set(gca,'YTick',[0:percenttick:100])
%xlabel('Scenario')
set(gca,'XTick',DummyX)
set(gca,'XTickLabel',[])
for ii=1:numscenarios
    text(DummyX(ii),-.05,scenarionames{ii},'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fontsize)
end
ylabel('Percent of Model Runs','FontSize',fontsize)
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out')
title('Sill Performance','FontSize',fontsize+2)

% Text box for regrounding labels:
text(mean(DummyX(1:2)),85,{'Never';'Regrounded'},'Color','k','BackgroundColor','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fontsize)
text(mean(DummyX(2:3)),50,{'Successfully';'Regrounded'},'Color','k','BackgroundColor','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fontsize)

% Set up colorbar:
caxis(lossratelims)
colormap(FadedJet)
hc=colorbar;
set(get(hc,'Ylabel'),'String',{'Post-Regrounding';'Sea Level Rise Rate';'(% of no-intervention rise rate)'})
set(get(hc,'Ylabel'),'FontSize',fontsize)
set(hc,'FontSize',fontsize)
set(hc,'YTick',[lossratelims(1):lossratetick:lossratelims(2)])

% Save figure:
set(gcf,'PaperSize',pagesize)
set(gcf,'PaperPosition',[0,0,pagesize])
print('-dpng',figname,'-r300')


% FInal Display:
disp('Done!')
toc

%% Cutting Room Floor:

% % COMPUTE DELAY:
% 
% % Find last good inds of each run:
% lastind_sill=find(TimeSeries(thisrun).Time~=0,1,'last');
% lastind_warming=find(TimeSeries(thiswarmingrun).Time~=0,1,'last');
% % Check whether any runs collapsed:
% if TimeSeries(thisrun).Time(end)==0 && TimeSeries(thiswarmingrun).Time(end)==0 % both runs collapsed
%     
%     % Uce point-based measurement at collapse:
%     % Check which terminal VAF was higher:
%     if TimeSeries(thisrun).VAF(lastind_sill)>TimeSeries(thiswarmingrun).VAF(lastind_warming) % sill run larger
%         % Locate last crossing of sill-run collapse point in warming run time series:
%         % determine level and time at sill-run collapse:
%         vaflevel=TimeSeries(thisrun).VAF(lastind_sill);
%         silltime=TimeSeries(thisrun).Time(lastind_sill);
%         % locate index:
%         thisind=find((TimeSeries(thiswarmingrun).VAF(1:end-1)<=vaflevel&TimeSeries(thiswarmingrun).VAF(2:end)>vaflevel)|...
%             (TimeSeries(thiswarmingrun).VAF(1:end-1)>vaflevel&TimeSeries(thiswarmingrun).VAF(2:end)<=vaflevel),1,'last');
%         % interpolate exact time:
%         warmingtime=TimeSeries(thiswarmingrun).Time(thisind)+dt*(vaflevel-TimeSeries(thiswarmingrun).VAF(thisind))/...
%             (TimeSeries(thiswarmingrun).VAF(thisind+1)-TimeSeries(thiswarmingrun).VAF(thisind));
%     else % warming run larger
%         % Locate last crossing of warming-run collapse point in sill run time series:
%         % determine level and time at warming-run collapse:
%         vaflevel=TimeSeries(thiswarmingrun).VAF(lastind_warming);
%         warmingtime=TimeSeries(thiswarmingrun).Time(lastind_warming);
%         % locate index:
%         thisind=find((TimeSeries(thisrun).VAF(1:end-1)<=vaflevel&TimeSeries(thisrun).VAF(2:end)>vaflevel)|...
%             (TimeSeries(thisrun).VAF(1:end-1)>vaflevel&TimeSeries(thisrun).VAF(2:end)<=vaflevel),1,'last');
%         % interpolate exact time:
%         silltime=TimeSeries(thisrun).Time(thisind)+dt*(vaflevel-TimeSeries(thisrun).VAF(thisind))/...
%             (TimeSeries(thisrun).VAF(thisind+1)-TimeSeries(thisrun).VAF(thisind));
%     end
%     
% elseif TimeSeries(thisrun).Time(end)==0 && TimeSeries(thiswarmingrun).Time(end)~=0 % only the sill run collapsed
%     
%     % Determine level and time at sill-run collapse:
%     vaflevel=TimeSeries(thisrun).VAF(lastind_sill);
%     silltime=TimeSeries(thisrun).Time(lastind_sill);
%     
%     % Check whether the warming run crossed that level:
%     if TimeSeries(thiswarmingrun).VAF(end)<vaflevel
%         % locate index:
%         thisind=find((TimeSeries(thiswarmingrun).VAF(1:end-1)<=vaflevel&TimeSeries(thiswarmingrun).VAF(2:end)>vaflevel)|...
%             (TimeSeries(thiswarmingrun).VAF(1:end-1)>vaflevel&TimeSeries(thiswarmingrun).VAF(2:end)<=vaflevel),1,'last');
%         % interpolate exact time:
%         warmingtime=TimeSeries(thiswarmingrun).Time(thisind)+dt*(vaflevel-TimeSeries(thiswarmingrun).VAF(thisind))/...
%             (TimeSeries(thiswarmingrun).VAF(thisind+1)-TimeSeries(thiswarmingrun).VAF(thisind));
%     else
%         % Warming run time is runtime:
%         warmingtime=Parameters(thiswarmingrun).TimingParameters.runtime_yr;
%     end
%     
% elseif TimeSeries(thisrun).Time(end)~=0 && TimeSeries(thiswarmingrun).Time(end)==0 % only the warming run collapsed
%     
%     % determine level and time at warming-run collapse:
%     vaflevel=TimeSeries(thiswarmingrun).VAF(lastind_warming);
%     warmingtime=TimeSeries(thiswarmingrun).Time(lastind_warming);
%     
%     % Check whether the sill run crossed that level:
%     if TimeSeries(thisrun).VAF(end)<vaflevel
%         % locate index:
%         thisind=find((TimeSeries(thisrun).VAF(1:end-1)<=vaflevel&TimeSeries(thisrun).VAF(2:end)>vaflevel)|...
%             (TimeSeries(thisrun).VAF(1:end-1)>vaflevel&TimeSeries(thisrun).VAF(2:end)<=vaflevel),1,'last');
%         % interpolate exact time:
%         silltime=TimeSeries(thisrun).Time(thisind)+dt*(vaflevel-TimeSeries(thisrun).VAF(thisind))/...
%             (TimeSeries(thisrun).VAF(thisind+1)-TimeSeries(thisrun).VAF(thisind));
%     else
%         % Sill run time is runtime:
%         silltime=Parameters(thisrun).TimingParameters.runtime_yr;
%     end
%     
% else % both runs did not collapse
%     
%     % Set times to NaN:
%     silltime=NaN;
%     warmingtime=NaN;
%     
% end
% 
% % Compute delay:
% SeaLevelDelay(thissample,thisscenario)=silltime-warmingtime;

% % Make figure:
% figure(1)
% 
% % First Subplot:
% subplot(2,1,1)
% 
% % Plot regrounding against volume:
% semilogx(SillVolumes,100*repmat(FractionRegrounded+FractionNeverRetreated,[1,numangles]),'k','Marker','s','MarkerFaceColor','k')
% hold on
% 
% % Plot other structures:
% plot(repmat(othervolumes,[2,1]),[0,0;100,100],'--k')
% 
% % Set up axes, etc:
% ylim([0,100])
% ylabel({'Percent Regrounded or';'Never Retreated'})
% xlim([1e-2,1e3])
% %xlabel('Sill Volume (km^3)')
% title('Success Rate')
% 
% % Second subplot:
% subplot(2,1,2)
% 
% % Plot delay time against volume:
% semilogx(SillVolumes,repmat(max(0,MedianDelay),[1,numangles]),'k','Marker','s','MarkerFaceColor','k')
% hold on
% 
% % Plot other structures:
% plot(repmat(othervolumes,[2,1]),[0,0;1000,1000],'--k')
% 
% % Set up axes, etc:
% ylim([0,1000])
% ylabel('Years')
% xlim([1e-2,1e3])
% xlabel('Sill Volume (km^3)')
% title('Median Collapse Delay')