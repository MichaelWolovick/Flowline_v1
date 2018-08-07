% PlotFlowlineTimeSeries_v1

% Mike Wolovick, 1/11/2017

% This script takes the output of GetFlowlineTimeSeries (which compiled the
% time series from many model runs) and plots a time series figure for each
% run.  The figure shows change in grounding line position and change in
% volume above flotation.

% Each plot has three model runs for the three scenarios.

% v1a: modified for climate experiment

clearvars
close all
tic

%% Parameters:

% File and folder names:
inputfile='/net/mjw/FjordSillPaper/ModelOutput/ClimateExperiment1.mat';
figfolder='/net/mjw/FjordSillPaper/Figures/';

% Experiment description: (n=numexperiments, m(n)=numsteps in experiments)
experimentnames={'Accumulation','Ablation','CDW Upwelling','CDW Warming'};  % 1xm cell array of strings
scenariokeys={{'c';'ac5';'ac10';'ac20'},...
    {'c';'ab3';'ab6';'ab9'},...
    {'c';'up1';'up2';'up3';'up4'},...
    {'c';'cdw.5';'cdw1';'cdw2'}}; % 1xm cell array of nx1 cell arrays of strings

% Sampling description:
flowbandnames={'ThwaitesA_v2';...
    'ThwaitesB_v2';...
    'ThwaitesC_v2'};
mkeys={'m03';'m10';'mva'};
ckeys={'ch';'chu';'chm'};

% Line colors:
linecolors={'k';'b';'g';'r';'m'};       % mx1 cell array of strings or RGB (max numsteps)

% Conversion between VAF and sea level:
oceansurfacearea=3.6e14;        % m^2

% Margins:
verttextbuffer=.05;             % unitless
titlebuffer=.05;                % unitless
horztextbuffer=.11;             % unitless
horznotextbuffer=.03;           % unitless

% Font size:
fontsize=11;                    % points

% Pagesize and resolution:
pagesize=[10,8];                % [1x2] inches
resolution=300;                 % dpi

%% Work:

% Hardcoded plot arrangement:
numvertplots=2;
numhorzplots=1;

% Compute subplot positions:
availablehorzspace=1-2*horztextbuffer-(numhorzplots-1)*horznotextbuffer; % 2 y-axes
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

% Load all of the time series:
load(inputfile)

% Compute number of plots:
numplots=length(experimentnames)*length(flowbandnames)*length(mkeys)*length(ckeys);

% Compute number of model runs:
numruns=length(TimeSeries);

% Compute sample size for all experiments:
numsamples=length(flowbandnames)*length(mkeys)*length(ckeys);

% Compute number of experiments:
numexperiments=length(experimentnames);

% Compute number of steps in each experiment:
numsteps=zeros(1,numexperiments);
for ii=1:numexperiments
    numsteps(ii)=length(scenariokeys{ii});
end

% Get runtime, numrecords:
runtime=Parameters(1).TimingParameters.runtime_yr;
numrecords=length(TimeSeries(1).Time);

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
    
    % Loop through experiments:
    for thisexperiment=1:numexperiments
        % Loop through steps:
        for thisstep=1:numsteps(thisexperiment)
            
            % Locate the model run:
            % Loop through all model runs:
            for thisrun=1:numruns
                % Parse filename:
                % Find underscores:
                underscorenums=strfind(TimeSeries(thisrun).filename,'_');
                % Check flowband:
                if strcmp(TimeSeries(thisrun).filename(1:underscorenums(2)-1),thisflowband)
                    % Check scenario:
                    if strcmp(TimeSeries(thisrun).filename(underscorenums(2)+1:underscorenums(3)-1),scenariokeys{thisexperiment}{thisstep})
                        % Check sliding:
                        if strcmp(TimeSeries(thisrun).filename(underscorenums(3)+1:underscorenums(4)-1),thism)
                            % Check calving:
                            if strcmp(TimeSeries(thisrun).filename(underscorenums(4)+1:underscorenums(5)-1),thisc)
                                % This is the one!!!
                                break
                            end
                        end
                    end
                end
                % Check that we found a model run:
                if thisrun==numruns
                    error('Unable to locate model run')
                end
            end
            
            % Call figure:
            figure(1)
                
            % Call first subplot:
            subplot('Position',Boxes{1})
            
            % Plot zero line:
            if thisstep==1
                hold off
                plot([0,runtime],[0,0],'--k')
                hold on
            end
            
            % Identify valid data:
            ValidData=TimeSeries(thisrun).Time~=0;
            
            % Compute ylims:
            if thisstep==1
                ylims_gl=[0,0];
            end
            ylims_gl(1)=min(ylims_gl(1),min(TimeSeries(thisrun).X_gl(ValidData)-TimeSeries(thisrun).X_gl(1))/1e3);
            ylims_gl(2)=max(ylims_gl(2),max(TimeSeries(thisrun).X_gl(ValidData)-TimeSeries(thisrun).X_gl(1))/1e3);
            
            % Plot grounding line:
            plot(TimeSeries(thisrun).Time(ValidData),(TimeSeries(thisrun).X_gl(ValidData)-TimeSeries(thisrun).X_gl(1))/1e3,'Color',linecolors{thisstep},'LineWidth',2)
            
            % Fix up axes and title:
            xlim([0,runtime])
            set(gca,'XTickLabel',[])
            ylim(ylims_gl)
            ylabel('Change in Position (km)','FontSize',fontsize)
            set(gca,'FontSize',fontsize)
            title({['Flowband: ',thisflowband,', sliding: ',thism,', calving: ',thisc,', Climate Experiment : ',experimentnames{thisexperiment}];'Grounding Line Position'},'interpreter','none','FontSize',fontsize)
            
            % Call second subplot:
            subplot('Position',Boxes{2})
            yyaxis left
            
            % Plot zero line:
            if thisstep==1
                hold off
                plot([0,runtime],[0,0],'--k')
                hold on
            end
            
            % Compute ylims:
            if thisstep==1
                ylims_vaf=[0,0];
            end
            ylims_vaf(1)=min(ylims_vaf(1),min(TimeSeries(thisrun).VAF(ValidData)-TimeSeries(thisrun).VAF(1))/1e12);
            ylims_vaf(2)=max(ylims_vaf(2),max(TimeSeries(thisrun).VAF(ValidData)-TimeSeries(thisrun).VAF(1))/1e12);
            
            % Plot volume above flotation:
            if thisstep==1
                h=zeros(numsteps(thisexperiment),1);
            end
            h(thisstep)=plot(TimeSeries(thisrun).Time(ValidData),(TimeSeries(thisrun).VAF(ValidData)-TimeSeries(thisrun).VAF(1))/1e12,...
                'Color',linecolors{thisstep},'LineWidth',2,'LineStyle','-','Marker','none');
            
            % Fix up axes and title:
            xlim([0,runtime])
            xlabel('Time (yr)','FontSize',fontsize)
            ylim(ylims_vaf)
            ylabel('Change in Volume (10^3 km^3)','FontSize',fontsize,'Color','k')
            set(gca,'YColor','k')
            title('Volume Above Flotation','FontSize',fontsize)
            set(gca,'FontSize',fontsize)
            
            % Create legend:
            if thisstep==numsteps(thisexperiment)
                if abs(ylims_vaf(1))>abs(ylims_vaf(2))
                    legend(h,scenariokeys{thisexperiment},'location','SouthWest')
                else
                    legend(h,scenariokeys{thisexperiment},'location','NorthWest')
                end
            end
            
            % Create 2nd y-axis for sea level equivalent:
            if thisstep==numsteps(thisexperiment)
                % Create 2nd y-axis:
                yyaxis right
                % Set ticks and limits:
                sllims=100*ylims_vaf*1e12/oceansurfacearea;
                ylim(sllims)
                slticks=get(gca,'YTick');
                set(gca,'YTickLabel',-slticks)
                set(gca,'YColor','k')
                % Set label:
                ylabel('Sea Level Rise (cm)','FontSize',fontsize,'Color','k')
            end
            
            % Save figure:
            if thisstep==numsteps(thisexperiment)
                set(gcf,'PaperSize',pagesize)
                set(gcf,'PaperPosition',[0,0,pagesize])
                figname=[figfolder,'ClimateExperimentTimeSeries_',experimentnames{thisexperiment},'_',thisflowband,'_',thism,'_',thisc,'.png'];
                print('-dpng',figname,['-r',num2str(resolution)])
            end
            
        end
    end
    
    
    
    
    
end


% Final display:
disp('Done!')
toc

