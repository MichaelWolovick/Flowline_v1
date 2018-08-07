% PlotFlowlineTimeSeries_v1

% Mike Wolovick, 1/11/2017

% This script takes the output of GetFlowlineTimeSeries (which compiled the
% time series from many model runs) and plots a time series figure for each
% run.  The figure shows change in grounding line position and change in
% volume above flotation.

% Each plot has three model runs for the three scenarios.

% v1a: modified for climate experiment

% v1b:  modified for sill partial blocking experiment

clearvars
close all
tic

%% Parameters:

% File and folder names:
inputfile='/net/mjw/FjordSillPaper/ModelOutput/ParamSweepRuns_v3_v3a.mat';
figfolder='/net/mjw/FjordSillPaper/Figures/';

% Experiment description: (n=numscenarios)
scenariokeys={'s6';'s4';'s5';'w';'c'}; % 1xn cell array strings
versionkeys={'v3a';'v3a';'v3a';'v3';'v3'}; % 1xn cell array strings
legendtext={'100% Blocking';'50% Blocking';'0% Blocking';'Warming Climate';'Constant Climate'}; % 1xn cell array of strings

% Sampling description:
flowbandnames={'ThwaitesA_v2';...
    'ThwaitesC_v2'};
mkeys={'m01';'m03';'m10'};
ckeys={'ch';'chu';'chm'};

% Line colors:
linecolors={'b';'g';'m';'r';'k'};       % nx1 cell array of strings or RGB

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
numplots=length(flowbandnames)*length(mkeys)*length(ckeys);

% Compute total number of model runs:
numruns=length(TimeSeries);

% Compute sample size for each scenario:
numsamples=length(flowbandnames)*length(mkeys)*length(ckeys);

% Compute number of scenarios:
numscenarios=length(scenariokeys);

% Get runtime, numrecords:
runtime=Parameters(1).TimingParameters.runtime_yr;
numrecords=length(TimeSeries(1).Time);

% Set the counters:
fcounter=1;
mcounter=1;
ccounter=1;
patchcounter=1;
patchhandles=zeros(2,2);

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
    
    % Loop through samples:
    for thisscenario=1:numscenarios
        
        % Locate the model run:
        % Loop through all model runs:
        for thisrun=1:numruns
            % Parse filename:
            % Find underscores:
            underscorenums=strfind(TimeSeries(thisrun).filename,'_');
            % Check flowband:
            if strcmp(TimeSeries(thisrun).filename(1:underscorenums(2)-1),thisflowband)
                % Check scenario:
                if strcmp(TimeSeries(thisrun).filename(underscorenums(2)+1:underscorenums(3)-1),scenariokeys{thisscenario})
                    % Check sliding:
                    if strcmp(TimeSeries(thisrun).filename(underscorenums(3)+1:underscorenums(4)-1),thism)
                        % Check calving:
                        if strcmp(TimeSeries(thisrun).filename(underscorenums(4)+1:underscorenums(5)-1),thisc)
                            % Check version key:
                            if strcmp(TimeSeries(thisrun).filename(underscorenums(5)+1:end-4),versionkeys{thisscenario})
                                % This is the one!!!
                                break
                            end
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
        if thisscenario==1
            hold off
            plot([0,runtime],[0,0],'--k')
            hold on
        end
        
        % Plot sill position: (both sill positions)
        if strcmp(scenariokeys{thisscenario},'s') || strcmp(scenariokeys{thisscenario},'s4')
            sillstarttime=Parameters(thisrun).SillParameters.sillstarttime_yr;
            sillendtime=Parameters(thisrun).SillParameters.sillstarttime_yr+Parameters(thisrun).SillParameters.sillconstructiontime_yr;
            sillstartx=Parameters(thisrun).SillParameters.sillcentx-2*Parameters(thisrun).SillParameters.sillwidth;
            sillendx=Parameters(thisrun).SillParameters.sillcentx+2*Parameters(thisrun).SillParameters.sillwidth;
            patchhandles(patchcounter,1)=patch('Vertices',[[sillstarttime;sillstarttime;sillendtime;sillendtime],(1e10)*[-1;1;1;-1]],'Faces',1:4,'FaceVertexCData',[.5,.5,.5],'EdgeColor','none','FaceColor','flat');
            patchhandles(patchcounter,2)=patch('Vertices',[[sillstarttime;sillstarttime;runtime;runtime],([sillstartx;sillendx;sillendx;sillstartx]-TimeSeries(thisrun).X_gl(1))/1e3],'Faces',1:4,'FaceVertexCData',[.5,.5,.5],'EdgeColor','none','FaceColor','flat');
            uistack(patchhandles(patchcounter,:),'bottom')
            patchcounter=patchcounter+1;
            if patchcounter==3
                patchcounter=1;
            end
        end
        
        % Identify valid data:
        ValidData=TimeSeries(thisrun).Time~=0;
        
        % Compute ylims:
        if thisscenario==1
            ylims_gl=[0,0];
        end
        ylims_gl(1)=min(ylims_gl(1),min(TimeSeries(thisrun).X_gl(ValidData)-TimeSeries(thisrun).X_gl(1))/1e3);
        ylims_gl(2)=max(ylims_gl(2),max(TimeSeries(thisrun).X_gl(ValidData)-TimeSeries(thisrun).X_gl(1))/1e3);
        
        % Plot grounding line:
        plot(TimeSeries(thisrun).Time(ValidData),(TimeSeries(thisrun).X_gl(ValidData)-TimeSeries(thisrun).X_gl(1))/1e3,'Color',linecolors{thisscenario},'LineWidth',2)
        
        % Fix up axes and title:
        xlim([0,runtime])
        set(gca,'XTickLabel',[])
        ylim(ylims_gl)
        ylabel('Change in Position (km)','FontSize',fontsize)
        set(gca,'FontSize',fontsize)
        title({['Sill Partial Blocking Experiment.  Flowband: ',thisflowband,', sliding: ',thism,', calving: ',thisc];'Grounding Line Position'},'interpreter','none','FontSize',fontsize)
        
        % Call second subplot:
        subplot('Position',Boxes{2})
        yyaxis left
        
        % Plot zero line:
        if thisscenario==1
            hold off
            plot([0,runtime],[0,0],'--k')
            hold on
        end
        
        % Compute ylims:
        if thisscenario==1
            ylims_vaf=[0,0];
        end
        ylims_vaf(1)=min(ylims_vaf(1),min(TimeSeries(thisrun).VAF(ValidData)-TimeSeries(thisrun).VAF(1))/1e12);
        ylims_vaf(2)=max(ylims_vaf(2),max(TimeSeries(thisrun).VAF(ValidData)-TimeSeries(thisrun).VAF(1))/1e12);
        
        % Plot volume above flotation:
        if thisscenario==1
            h=zeros(numscenarios,1);
        end
        h(thisscenario)=plot(TimeSeries(thisrun).Time(ValidData),(TimeSeries(thisrun).VAF(ValidData)-TimeSeries(thisrun).VAF(1))/1e12,...
            'Color',linecolors{thisscenario},'LineWidth',2,'LineStyle','-','Marker','none');
        
        % Fix up axes and title:
        xlim([0,runtime])
        xlabel('Time (yr)','FontSize',fontsize)
        ylim(ylims_vaf)
        ylabel('Change in Volume (10^3 km^3)','FontSize',fontsize,'Color','k')
        set(gca,'YColor','k')
        title('Volume Above Flotation','FontSize',fontsize)
        set(gca,'FontSize',fontsize)
        
        % Create legend:
        if thisscenario==numscenarios
            if abs(ylims_vaf(1))>abs(ylims_vaf(2))
                legend(h,legendtext,'location','SouthEast')
            else
                legend(h,legendtext,'location','NorthEast')
            end
        end
        
        % Create 2nd y-axis for sea level equivalent:
        if thisscenario==numscenarios
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
        if thisscenario==numscenarios
            set(gcf,'PaperSize',pagesize)
            set(gcf,'PaperPosition',[0,0,pagesize])
            figname=[figfolder,'PartialBlockingExperiment_',thisflowband,'_',thism,'_',thisc,'.png'];
            print('-dpng',figname,['-r',num2str(resolution)])
        end
        
    end
    
    
    
    
    
end


% Final display:
disp('Done!')
toc

