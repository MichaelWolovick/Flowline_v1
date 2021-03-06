% FjordSillResultsFigure_v1

% Mike Wolovick, 1/23/2017

% This figure plots model results for the artificial sill paper.  The
% figure shows a series of small snapshots overlayed on time series of
% grounding line position and volume above flotation.  The snapshot times
% are chosen manually to capture important configurations.

% The snapshots are simple displays relative to other snapshot figures I
% have made, because they are going to be shrunk down and overlain on the
% time series.  The only color scale in the snapshots is for ocean
% temperature, with solid colors for ice, rock, sky, and sill.  The
% snapshots are output as separate figures which I will overlay on the time
% series in photoshop.

% The snapshots are set to a target vertical exageration, so that snapshots
% of different glaciers can be readily compared.

% Each snapshot must specify both the time and the model run it is taken
% from.

% If you specify a zoom-in area for the snapshots, the last snapshot shows
% the whole domain with the box overlain (should be time zero).

% v2:  The snapshots are subplots within the same figure as the time
% series.  

% v3:  Only shows a timeseries of VAF, not GL position.


clear all
close all
tic

%% Parameters:

% File names and paths:
inputfolder='/media/wolovick/My Passport/GFDL/Work/FjordSillPaper/ModelOutput/';
inputfiles={'ThwaitesC_v2_s4_m03_ch_v3a.mat';...
    'ThwaitesC_v2_w_m03_ch_v3.mat';...
    'ThwaitesC_v2_c_m03_ch_v3.mat'}; % cell array of strings
% inputfiles={'ThwaitesA_v2_s1_m03_chu_v3a.mat';...
%     'ThwaitesA_v2_s2_m03_chu_v3a.mat';...
%     'ThwaitesA_v2_s3_m03_chu_v3a.mat'}; % cell array of strings
% inputfiles={'ThwaitesC_v2_s_m03_chm_v3.mat';...
%     'ThwaitesC_v2_w_m03_chm_v3.mat';...
%     'ThwaitesC_v2_c_m03_chm_v3.mat'}; % cell array of strings
% inputfiles={'ThwaitesC_v2_s_m01_ch_v3.mat';...
%     'ThwaitesC_v2_w_m01_ch_v3.mat';...
%     'ThwaitesC_v2_c_m01_ch_v3.mat'}; % cell array of strings
outputfolder='/home/wolovick/Dropbox/FjordSillPaper/Figures/';
mainfigname='ResultsFigureThwaites_v10b';

% Title:
titletext={'Initial Conditions';...
    'Collapse Underway';...
    'Sill Built, Initial Regrounding';...
    'Inner Grounding Line Recovery';...
    'Inner Grounding Line Retreat'}; % cell array of strings for snapshots
% titletext={'Marine Ice Sheet Retreat Begins';...
%     'Initial Regrounding';...
%     'Weak Sill Delayed Collapse';...
%     'Medium Sill Collapse Paused';...
%     'Strong Sill Recovery'}; % cell array of strings for snapshots
% titletext={'Initial Conditions';...
%     'Retreat Across Overdeepening';...
%     'Readvance and Shelf Growth';...
%     'Retreat Back Across Overdeepening';...
%     'Final Conditions'}; % cell array of strings for snapshots
% titletext={'Initial Conditions';...
%     'Marine Ice Sheet Retreat';...
%     'Initial Regrounding';...
%     'Continued Regrounding';...
%     'Final Conditions'}; % cell array of strings for snapshots
% titletext={'Initial Conditions';...
%     'Marine Ice Sheet Retreat';...
%     'Initial Regrounding';...
%     'Sparse, Ephemeral Grounding';...
%     'Final Conditions'}; % cell array of strings for snapshots

% Time series legend:
legendtext={'Warming Climate + Sill';'Warming Climate';'Constant Climate'}; % nx1 cell array of strings
%legendtext={'Weak Sill (1 kPa)';'Medium Sill (10 kPa)';'Strong Sill (100 kPa)'}; % nx1 cell array of strings
%legendtext={'100% Blockage';'50% Blockage';'0% Blockage'}; % nx1 cell array of strings
legendpos='NorthEast';          % valid legend location string
legendfontsize=12;                % points

% Target Vertical exaggeration:
%targetexag=25;                       % unitless
%targetexag=60;                       % unitless
targetexag=50;                       % unitless

% Snapshot controls:
% HelheimA_v2_m03_chm_v3:
% snapshottimes=[0;90;130;280;999];    % mx1 years
% snapshotruns=[3;2;1;1;1];            % mx1 indices into inputfiles
% snapshotlabelalignments=[7;1;5;5;3]; % mx1 integers [1-8] (clockwise from upper right, relative to the letter)
% KangerB_v2_m10_chm_v3:
% snapshottimes=[0;180;115;353;999];    % mx1 years
% snapshotruns=[1;2;3;3;3];            % mx1 indices into inputfiles
% snapshotlabelalignments=[7;1;5;5;3]; % mx1 integers [1-8] (clockwise from upper right, relative to the letter)
% HelheimA_v2_m03_chm_v3:
% snapshottimes=[0;90;130;280;999];    % mx1 years
% snapshotruns=[3;2;1;1;1];            % mx1 indices into inputfiles
% snapshotlabelalignments=[7;1;5;5;3]; % mx1 integers [1-8] (clockwise from upper right, relative to the letter)
% ThwaitesC_v2_m01_ch_v3:
% snapshottimes=[0;110;160;500;1000];  % mx1 years
% snapshotruns=[3;2;1;1;1];            % mx1 indices into inputfiles
% snapshotlabelalignments=[7;1;5;5;3]; % mx1 integers [1-8] (clockwise from upper right, relative to the letter)
% ThwaitesC_v2_m03_chm_v3:
% snapshottimes=[0;110;116;150;1000];  % mx1 years
% snapshotruns=[3;2;1;1;1];            % mx1 indices into inputfiles
% snapshotlabelalignments=[5;3;7;5;3]; % mx1 integers [1-8] (clockwise from upper right, relative to the letter)
% ThwaitesA_v2_m10_chu_v3a:
% snapshottimes=[99;172;315;800;800];  % mx1 years
% snapshotruns=[1;1;1;2;3];            % mx1 indices into inputfiles
% snapshotlabelalignments=[1;3;5;1;3]; % mx1 integers [1-8] (clockwise from upper right, relative to the letter)
% snapshotlabelfontsize=14;       % points
% ThwaitesC_v2_m03_ch_v3a:
snapshottimes=[0;99;111;255;650];  % mx1 years
snapshotruns=[3;2;1;1;1];            % mx1 indices into inputfiles
snapshotlabelalignments=[7;1;1;5;5]; % mx1 integers [1-8] (clockwise from upper right, relative to the letter)
snapshotlabelfontsize=16;       % points

% Line settings for time series:
tslinewidth=1;                  % points
tslinecolor={'g';'r';'b'};      % nx1 cell array (color string or RGB)
%tslinecolor={'r';'g';'b'};      % nx1 cell array (color string or RGB)
tslinestyle={'-';'-';'-'};      % nx1 cell array (linestyle strings)
snapshotmarker='o';             % valid marker string
snapshotmarkersize=5;           % points

% Axes ticks and limits:
tlims='auto';                   % 1x2 yr or 'auto'
ttick=200;                      % yr
% vaflims=[-1500,100];            % 1x2 km^3 or 10^3 km^3
% vaftick=500;                    % km^3 or 10^3 km^3
% vaflims=[-10,60];               % 1x2 km^3 or 10^3 km^3
% vaftick=20;                     % km^3 or 10^3 km^3
% vaf103=0;                       % logical (controls whether to use units of km^3 or 10^3 km^3)
% vafrel=1;                       % logical (toggles whether to show VAF of the sill run relative to the warming run. Assumes warming run is 2 and sill run is 1).
% xlims=[150,355];                % 1x2 km or 'auto'
% xtick=25;                       % km
% elevlims=[-1000,2500];          % 1x2 m
% elevtick=500;                   % m
vaflims=[-100,0];               % 1x2 km^3 or 10^3 km^3
vaftick=50;                     % km^3 ) 10^3 km^3
% vaflims=[-75,25];               % 1x2 km^3 or 10^3 km^3
% vaftick=25;                     % km^3 ) 10^3 km^3
% sltick=5;                       % cm
sltick=10;                       % cm
vaf103=1;                       % logical (controls whether to use units of km^3 or 10^3 km^3)
vafrel=0;                       % logical (toggles whether to show VAF of the sill run relative to the warming run. Assumes warming run is 2 and sill run is 1).
xlims='auto';                   % 1x2 km or 'auto'
xtick=100;                      % km
elevlims=[-1500,2000];          % 1x2 m
elevtick=1000;                  % m
elevkm=1;                       % logical (doesn't affect values of elevlims or elevticks, only labels)

% Conversion between VAF and sea level:
oceansurfacearea=3.6e14;        % m^2

% Colormaps and limits:
oceantempcmap='parula';         % valid colormap string or 'fadedjet'
%oceantemplims=[0,6];            % [1x2] deg C (can be 'auto')
oceantemplims=[-1.25,1];        % [1x2] deg C (can be 'auto')

% Plot which parts of the IC profile?
ploticsurface=1;                % logical
ploticbottom=0;                 % logical
ploticfront=0;                  % logical
ploticgl=0;                     % logical

% Line widths (snapshots):
bedlinewidth=.5;                % points
surflinewidth=.5;               % points
gllinewidth=.5;                 % points
frontlinewidth=.5;              % points
elalinewidth=0;                % points (0 deactivates)
inputlinewidth=1;               % points
borderlinewidth=1;              % points
zoomboxlinewidth=.5;            % points

% Block colors in the snapshots:
skycolor=[.75,.75,1];           % [R,G,B]
skyzsize=10;                    % integer
bedcolor=[.63,.32,.18];         % [R,G,B]
icecolor=[.45,.62,.78];         % [R,G,B]
sillcolor=[.5,.5,.5];           % [R,G,B]
shelfcolor=[.875,0,.875];       % [R,G,B]
oceanzsize=50;                  % integer

% Other Settings:
ticklength=.02;                 % unitless
fontsize=18;                    % points
tstickdir='in';                 % 'in' or 'out'
tsminorticksx='on';             % 'on' or 'off'
tsminorticksy='on';             % 'on' or 'off'
snapshottickdir='out';          % 'in' or 'out'
snapshotminorticksx='off';      % 'on' or 'off'
snapshotminorticksy='off';      % 'on' or 'off'

% Margins:
verttextbuffer=.081;            % unitless
titlebuffer=.05;                % unitless
horztextbuffer=.066;            % unitless
horznotextbuffer=.02;           % unitless

% Output size and resolution:
%pagesize=[8,12];                % [1x2] inches
pagewidth=16;                   % inches
resolution=300;                 % dpi

%% Preparation:

% Turn off this warning:
warning('off','MATLAB:hg:patch:RGBColorDataNotSupported')

% Determine number of snapshots and model runs:
numsnapshots=length(snapshottimes);
numruns=length(inputfiles);

% Check lengths of other parameters:
if size(tslinecolor,1)~=numruns
    error('Parameter "tslinecolor" must be the same length as "inputfiles".')
end
if size(tslinestyle,1)~=numruns
    error('Parameter "tslinestyle" must be the same length as "inputfiles".')
end
if size(legendtext,1)~=numruns
    error('Parameter "legendtext" must be the same length as "inputfiles".')
end
if length(snapshotruns)~=numsnapshots
    error('Parameter "snapshotruns" must be the same length as "snapshottimes".')
end
if length(titletext)~=numsnapshots
    error('Parameter "titletext" must have length equal to numsnapshots.')
end
if length(snapshotlabelalignments)~=numsnapshots
    error('Parameter "snapshotlabelalignments" must be the same length as "snapshottimes".')
end

% Compute plot arrangement:
numvertplots=floor(numsnapshots/2)+1;
numhorzplots=2;

% Compute subplot positions (snapshots): (2 y-axes)
availablehorzspace=1-2*horztextbuffer-(numhorzplots-1)*horznotextbuffer;
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

% Modify top left box for 2 y-axes (3rd y-axis on figure overall)
Boxes{1,1}(3)=Boxes{1,1}(3)-(horztextbuffer-horznotextbuffer);

% Make top left box centered horizontally in its column:
%Boxes{1,1}(3)=Boxes{1,1}(3)-(horztextbuffer-horznotextbuffer);
%Boxes{1,1}(1)=Boxes{1,1}(1)+(horztextbuffer-horznotextbuffer);

% Modify top left box for it's own x-axis:
Boxes{1,1}(2)=Boxes{1,1}(2)+verttextbuffer;
Boxes{1,1}(4)=Boxes{1,1}(4)-verttextbuffer;

% Make a figure:
figure(1)

% Make faded jet colormap:
FadedJet=colormap('jet');
FadedJet(3*round(size(FadedJet,1)/8)+1:round(size(FadedJet,1)/2),1)=linspace(0,1,round(size(FadedJet,1)/8))';
FadedJet(3*round(size(FadedJet,1)/8)+1:round(size(FadedJet,1)/2),3)=1;
FadedJet(round(size(FadedJet,1)/2)+1:5*round(size(FadedJet,1)/8),1)=1;
FadedJet(round(size(FadedJet,1)/2)+1:5*round(size(FadedJet,1)/8),3)=linspace(1,0,round(size(FadedJet,1)/8))';

% Assign colormaps:
if strcmp(oceantempcmap,'fadedjet')
    OceanTempcmap=FadedJet;
else
    OceanTempcmap=colormap(oceantempcmap);
end

% Determine glacier name and parameter values:
underscorenum=strfind(inputfiles{1},'_');
glaciername=inputfiles{1}(1:underscorenum(2)-1);
mvalue=str2double(inputfiles{1}(underscorenum(3)+2:underscorenum(3)+3));
ctype=inputfiles{1}(underscorenum(4)+(2:3));

% Generate calving string:
if strcmp(ctype,'h_')
    calvingstring='\dot{c}(H)';
elseif strcmp(ctype,'hu')
    calvingstring='\dot{c}(u,H)';
elseif strcmp(ctype,'hm')
    calvingstring='\dot{c}(\dot{m},H)';
else
    error('Calving Law Not Recognized')
end

% Check that all runs come from the same glacier:
for runnum=1:numruns
    underscorenum=strfind(inputfiles{runnum},'_');
    thisglaciername=inputfiles{runnum}(1:underscorenum(2)-1);
    if strcmp(thisglaciername,glaciername)==0
        error('This script is expecting all of the model runs to be for the same glacier flowband.')
    end
end

% Define VAF divisor:
if vaf103
    vafdivisor=1e12;
else
    vafdivisor=1e9;
end

% Load VAF reference profile from the warming run:
if vafrel
    load([inputfolder,inputfiles{2}],'ModelTimeSeries')
    VAFref=ModelTimeSeries.VAF;
end

%% Make Time Series Figure:

% Pre-allocate line handles:
tshandles=zeros(numruns,1);

% Loop through model runs:
for runnum=1:numruns
    
    % Load time series:
    load([inputfolder,inputfiles{runnum}],'ModelTimeSeries','*Parameters')
    
    % Get runtime:
    runtime=TimingParameters.runtime_yr;
    
    % Get tlims:
    if strcmp(tlims,'auto')
        tlims=[0,runtime];
    end
    
    % Identify valid data:
    if SillParameters.dosill
        ValidData=ModelTimeSeries.X_gl~=0; %&ModelTimeSeries.Time>=SillParameters.sillstarttime_yr;
    else
        ValidData=ModelTimeSeries.X_gl~=0;
    end
    lastvaliddata=find(ValidData,1,'last');
    
    % Call figure:
    figure(1)
    
    % Call first subplot:
    subplot('Position',Boxes{1,1})
    
    % Plot zero line:
    plot([0,tlims(2)],[0,0],'k')
    hold on
    
    %     % Compute running ylims:
    %     if runnum==1
    %         vaflims=[min(ModelTimeSeries.VAF(ValidData)-ModelTimeSeries.VAF(1)),max(ModelTimeSeries.VAF(ValidData)-ModelTimeSeries.VAF(1))];
    %     else
    %         vaflims=[min(vaflims(1),min(ModelTimeSeries.VAF(ValidData)-ModelTimeSeries.VAF(1))),max(vaflims(2),max(ModelTimeSeries.VAF(ValidData)-ModelTimeSeries.VAF(1)))];
    %     end
    
    % Plot volume above flotation:
    if vafrel
        if runnum==1
            % Plot VAF change from warming run:
            tshandles(runnum)=plot(ModelTimeSeries.Time(ValidData),(ModelTimeSeries.VAF(ValidData)-VAFref(ValidData))/vafdivisor,...
                'Color',tslinecolor{runnum},'LineWidth',tslinewidth,'LineStyle',tslinestyle{runnum});
        else
            % Dummy plot:
            tshandles(runnum)=plot(tlims(1)-[100,99],[0,0],'Color',tslinecolor{runnum},'LineWidth',tslinewidth,'LineStyle',tslinestyle{runnum});
        end
    else
        % Plot VAF change since start:
        tshandles(runnum)=plot(ModelTimeSeries.Time(ValidData),(ModelTimeSeries.VAF(ValidData)-ModelTimeSeries.VAF(1))/vafdivisor,...
            'Color',tslinecolor{runnum},'LineWidth',tslinewidth,'LineStyle',tslinestyle{runnum});
    end
    
    % Plot snapshot times:
    if vafrel==0
        for ii=1:numsnapshots
            %         % Check if this is a unique time:
            %         if ii>1
            %             if isempty(find(snapshottimes(1:ii-1)==snapshottimes(ii),1,'last'))
            %                 % Plot vertical line for this snapshot:
            %                 plot(snapshottimes(ii)*[1,1],vaflims,'--k')
            %             end
            %         end
            if snapshotruns(ii)==runnum
                % Interpolate VAF at this snapshot time:
                thisvaf=interp1([0;ModelTimeSeries.Time(ValidData);runtime],...
                    [0;ModelTimeSeries.VAF(ValidData)-ModelTimeSeries.VAF(1);ModelTimeSeries.VAF(lastvaliddata)-ModelTimeSeries.VAF(1)],...
                    snapshottimes(ii),'linear');
                % Plot this snapshot marker:
                plot(snapshottimes(ii),thisvaf/vafdivisor,'LineStyle','none','Marker',snapshotmarker,'MarkerSize',snapshotmarkersize,'Color','k','MarkerFaceColor',tslinecolor{runnum})
                % Label this snapshot marker:
                % Check the text alignment:
                if snapshotlabelalignments(ii)==1
                    text(snapshottimes(ii),thisvaf/vafdivisor,char(65+ii),'HorizontalAlignment','right','VerticalAlignment','top','Color','k','FontSize',snapshotlabelfontsize)
                elseif snapshotlabelalignments(ii)==2
                    text(snapshottimes(ii),thisvaf/vafdivisor,char(65+ii),'HorizontalAlignment','right','VerticalAlignment','middle','Color','k','FontSize',snapshotlabelfontsize)
                elseif snapshotlabelalignments(ii)==3
                    text(snapshottimes(ii),thisvaf/vafdivisor,char(65+ii),'HorizontalAlignment','right','VerticalAlignment','bottom','Color','k','FontSize',snapshotlabelfontsize)
                elseif snapshotlabelalignments(ii)==4
                    text(snapshottimes(ii),thisvaf/vafdivisor,char(65+ii),'HorizontalAlignment','center','VerticalAlignment','bottom','Color','k','FontSize',snapshotlabelfontsize)
                elseif snapshotlabelalignments(ii)==5
                    text(snapshottimes(ii),thisvaf/vafdivisor,char(65+ii),'HorizontalAlignment','left','VerticalAlignment','bottom','Color','k','FontSize',snapshotlabelfontsize)
                elseif snapshotlabelalignments(ii)==6
                    text(snapshottimes(ii),thisvaf/vafdivisor,char(65+ii),'HorizontalAlignment','left','VerticalAlignment','middle','Color','k','FontSize',snapshotlabelfontsize)
                elseif snapshotlabelalignments(ii)==7
                    text(snapshottimes(ii),thisvaf/vafdivisor,char(65+ii),'HorizontalAlignment','left','VerticalAlignment','top','Color','k','FontSize',snapshotlabelfontsize)
                elseif snapshotlabelalignments(ii)==8
                    text(snapshottimes(ii),thisvaf/vafdivisor,char(65+ii),'HorizontalAlignment','center','VerticalAlignment','top','Color','k','FontSize',snapshotlabelfontsize)
                else
                    error('Parameter "snapshotlabelalignments" must be an integer from 1 to 8.')
                end
            end
        end
    end
    
    % Plot sill construction time:
    if runnum==1
        sillstarttime=SillParameters.sillstarttime_yr;
        sillendtime=SillParameters.sillstarttime_yr+SillParameters.sillconstructiontime_yr;
    end
    if runnum==numruns
        silltimehandle2=patch('Vertices',[[sillstarttime;sillstarttime;sillendtime;sillendtime],[vaflims(1);vaflims(2);vaflims(2);vaflims(1)]],...
            'Faces',1:4,'FaceVertexCData',[.5,.5,.5],'EdgeColor','none','FaceColor','flat');
        uistack(silltimehandle2,'bottom')
    end
    
    % Make a text box describing this model experiment:
    if runnum==numruns
        text(.98*tlims(2),vaflims(1)+.13*(vaflims(2)-vaflims(1)),...
            {['Flowband: ',inputfiles{1}(underscorenum(1)-1)];['Sliding: r=',num2str(mvalue)];['Calving: $',calvingstring,'$']},...
            'HorizontalAlignment','Right','VerticalAlignment','Bottom','Color','k','Interpreter','Latex','FontSize',legendfontsize)
    end
    
    % Fix up axes and title:
    xlim(tlims)
    set(gca,'XTick',ttick*(0:1:floor(tlims(2)/ttick)))
    xlabel('Model Time (yr)','FontSize',fontsize)
    ylim(vaflims)
    set(gca,'YTick',vaftick*(ceil(vaflims(1)/vaftick):1:floor(vaflims(2)/vaftick)))
    if vaf103
        %ylabel({'Change in';'Volume (10^3 km^3)'},'FontSize',fontsize)
        ylabel('\Delta VAF (10^3 km^3)','FontSize',fontsize)
    else
        %ylabel({'Change in';'Volume (km^3)'},'FontSize',fontsize)
        ylabel('\Delta VAF (km^3)','FontSize',fontsize)
    end
    if vafrel
        title('a) Volume Above Flotation Relative to Warming Run','FontSize',fontsize)
    else
        title('a) Change in Volume Above Flotation','FontSize',fontsize)
    end
    set(gca,'FontSize',fontsize)
    set(gca,'TickDir',tstickdir)
    set(gca,'XMinorTick',tsminorticksx)
    set(gca,'YMinorTick',tsminorticksy)
    
    % Create 2nd y-axis for sea level equivalent:
    if runnum==numruns
        % Create 2nd y-axis:
        yyaxis right
        % Set ticks and limits:
        sllims=100*vaflims*vafdivisor/oceansurfacearea;
        ylim(sllims)
        set(gca,'YTick',sltick*(ceil(sllims(1)/sltick):1:floor(sllims(2)/sltick)))
        set(gca,'YTickLabel',-sltick*(ceil(sllims(1)/sltick):1:floor(sllims(2)/sltick)))
        set(gca,'YColor','k')
        set(gca,'TickDir',tstickdir)
        set(gca,'YMinorTick',tsminorticksy)
        % Set label:
        ylabel('Sea Level Rise (cm)','FontSize',fontsize,'Color','k')
    end
    
    % Create legend:
    if runnum==numruns
        legend(tshandles,legendtext,'location',legendpos,'FontSize',legendfontsize)
    end
    
end

% Draw:
drawnow

%% Make Snapshots:

% Loop through snapshots:
for snapshot=1:numsnapshots
    
    % Load and unpack model output:
    load([inputfolder,inputfiles{snapshotruns(snapshot)}],'*Parameters','*_input','Zhat_*','DZhat_*','numrecords','numdigits')
    unpack(ModelParameters)
    unpack(GridParameters)
    unpack(SideParameters)
    unpack(PlumeParameters)
    unpack(TimingParameters)
    unpack(TracerParameters)
    unpack(ThermalParameters)
    unpack(BottomParameters)
    unpack(OtherParameters)
    if exist('SillParameters','var')
        unpack(SillParameters)
        if exist('sillblockingfraction','var')==0
            sillblockingfraction=1;
        end
    else
        dosill=0;
    end
    if exist('interpmethod','var')==0
        interpmethod='linear';
    end
    
    % Define ocean temperature limits:
    if strcmp(oceantemplims,'auto')
        if oceantimedependence
            theseoceantemplims=[min(min(Temperature_input)),max(max(Temperature_input))]-tmelt;
        else
            theseoceantemplims=[min(Temperature_input(:,1)),max(Temperature_input(:,1))]-tmelt;
        end
    else
        theseoceantemplims=oceantemplims;
    end
    
    % Define maximum domain width:
    maxdomainwidth=X_input(end);
    
    % Define input xsize:
    inputxsize=length(X_input);
    
    % Set x-limits:
    if strcmp(xlims,'auto') && snapshot==1
        xlims=[0,maxdomainwidth]/1000;
    end
    
    % Compute page size from target exaggeration:
    if snapshot==1
        pageheight=pagewidth*(Boxes{2}(3)/Boxes{2}(4))*targetexag*.001*(elevlims(2)-elevlims(1))/(xlims(2)-xlims(1));
    end
    
    % Get this model record:
    if snapshottimes(snapshot)==0
        % This record is the IC:
        thisrecord=0;
        % Load IC:
        load([inputfolder,inputfiles{snapshotruns(snapshot)}],'InitialConditions')
        ModelRecord=InitialConditions;
        clear InitialConditions
        % Assign things that aren't recorded in the IC:
        % time is zero:
        ModelRecord.time_yr=0;
        % Surface melt is annual melt:
        ModelRecord.MeltRate_u=ModelRecord.AnnualMelt_u;
        % interpolate basal drag:
        if strcmp(slidingstressscale,'file')
            X_lr=linspace(0,ModelRecord.domainwidth,xsize+1);
            ModelRecord.Drag_lrd=interp1(X_input,SlidingStressScale_input,X_lr,interpmethod);
        else
            ModelRecord.Drag_lrd=slidingstressscale*ones(1,xsize+1);
        end
        % interpolate velocity:
        ModelRecord.U_lr=interp1(X_input,VelMag_yr_input/secondsperyear,X_lr,interpmethod);
        % set plume melt rate to zero:
        ModelRecord.MeltRate_c_plume=zeros(1,xsize_plume+zsize_plume);
        % set calving rate to zero:
        ModelRecord.calvingrate_r=0;
    elseif snapshottimes(snapshot)==runtime_yr
        % This record is the final state:
        load([inputfolder,inputfiles{snapshotruns(snapshot)}],'FinalConditions')
        if exist('FinalConditions','var')
            ModelRecord=FinalConditions;
            clear FinalConditions
        else
            error('Unable to locate FinalConditions variable')
        end
    else
        % Determine this record number:
        thisrecord=round(snapshottimes(snapshot)/recordinterval_yr+.5);  % assumes that record time is the middle of the record
        % Load this model record:
        prefix='0'*ones(1,numdigits(end)-floor(log10(thisrecord))-1);
        idnumber=num2str(thisrecord);
        load([inputfolder,inputfiles{snapshotruns(snapshot)}],['ModelRecord_',prefix,idnumber])
        % Check if model record exists:
        if exist(['ModelRecord_',prefix,idnumber],'var')
            % Rename model record:
            eval(['ModelRecord=ModelRecord_',prefix,idnumber,';'])
            clear(['ModelRecord_',prefix,idnumber])
        else
            error(['Unable to locate ModelRecord_',prefix,idnumber])
        end
    end
    
    % Assign sea level:
    if strcmp(rightbctype,'front') && strcmp(rightbcparam,'file')==0
        sealevel=rightbcparam;
    elseif strcmp(rightbctype,'front')
        sealevel=interp1(Time_yr_input,RightBCParam_input,snapshottimes(snapshot),interpmethod,RightBCParam_input(1));
    else
        sealevel=min(BedElev_input)-1;
    end
    
    % Define sill profile:
    if dosill
        if isfield(ModelRecord,'SillThick_input')
            SillThick_input=ModelRecord.SillThick_input;
        else
            % Compute sill thickness profile:
            if sillstarttime_yr==0 && sillconstructiontime_yr==0 && snapshottimes(snapshot)==0
                SillThick_input=sillamp*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2));
            else
                SillThick_input=sillamp*max(0,min(1,(snapshottimes(snapshot)-sillstarttime_yr)/(sillconstructiontime_yr+1)))*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2)); % note stabilizer
            end
        end
    else
        SillThick_input=zeros(size(X_input));
    end
    
    % Define initial domain:
    dx=ModelRecord.domainwidth/xsize;
    X_lr=linspace(0,ModelRecord.domainwidth,xsize+1);
    X_c=linspace(.5*dx,ModelRecord.domainwidth-.5*dx,xsize);
    BedElev_lr=interp1(X_input,BedElev_input+SillThick_input,X_lr,interpmethod);
    BedElev_c=.5*(BedElev_lr(1:end-1)+BedElev_lr(2:end));
    
    % Compute initial profiles:
    HydroHead_input=BedElev_input+(rho_i/rho_sw)*Icethick_input;
    SurfElev_input=zeros(1,inputxsize);
    SurfElev_input(HydroHead_input>=sealevel)=BedElev_input(HydroHead_input>=0)+Icethick_input(HydroHead_input>=0);
    SurfElev_input(HydroHead_input<sealevel)=(1-rho_i/rho_sw)*Icethick_input(HydroHead_input<0);
    IceBottom_input=SurfElev_input-Icethick_input;
    
    % Compute grounding line of the input profiles:
    lastgroundedind_input=find(HydroHead_input>=0,1,'last');
    lasticeind_input=find(isnan(HydroHead_input)==0,1,'last');
    if lastgroundedind_input<lasticeind_input
        x_gl_input=interp1(HydroHead_input(lastgroundedind_input-4:min(lastgroundedind_input+4,lasticeind_input)),...
            X_input(lastgroundedind_input-4:min(lastgroundedind_input+4,lasticeind_input)),sealevel,'pchip');
        bedelev_gl_input=interp1(X_input,BedElev_input,x_gl_input,interpmethod,BedElev_input(end));
        icethick_gl_input=(sealevel-bedelev_gl_input)*(rho_sw/rho_i);
        hasshelf_input=1;
    else
        x_gl_input=X_input(lastgroundedind_input);
        bedelev_gl_input=interp1(X_input,BedElev_input,x_gl_input,interpmethod,BedElev_input(end));
        icethick_gl_input=Icethick_input(lastgroundedind_input);
        hasshelf_input=0;
    end
    surfelev_gl_input=bedelev_gl_input+icethick_gl_input;
    
    % Define vertical grid:
    Zhat_ud=linspace(0,1,zsize+1)'/maxdensify+(1-1/maxdensify)*linspace(0,1,zsize+1)'.^densifypower;
    Zhat_c=.5*(Zhat_ud(1:end-1)+Zhat_ud(2:end));
    DZhat_c=Zhat_ud(2:end)-Zhat_ud(1:end-1);
    DZhat_ud=[Zhat_c(1);Zhat_c(2:end)-Zhat_c(1:end-1);1-Zhat_c(end)];  % first and last are half-cells
    
    % Calculate surface elevation:
    SurfElev_c=ModelRecord.IceBottom_c+ModelRecord.Icethick_c;
    SurfElev_lr=[1.5*SurfElev_c(1)-.5*SurfElev_c(2),.5*(SurfElev_c(1:end-1)+SurfElev_c(2:end)),1.5*SurfElev_c(end)-.5*SurfElev_c(end-1)]; % linear extrapolation
    
    % Compute ELA:
    [MonotonicSurface,MonotonicSurfInd]=sort(SurfElev_c);
    ela=interp1(ModelRecord.Accum_u(MonotonicSurfInd)-ModelRecord.MeltRate_u(MonotonicSurfInd),MonotonicSurface,0,'linear');
    
    % Interpolate ice thickness to grid edges:
    Icethick_lr=[1.5*ModelRecord.Icethick_c(1)-.5*ModelRecord.Icethick_c(2),.5*(ModelRecord.Icethick_c(1:end-1)+ModelRecord.Icethick_c(2:end)),1.5*ModelRecord.Icethick_c(end)-.5*ModelRecord.Icethick_c(end-1)];
    
    % Interpolate ice bottom to grid edges:
    IceBottom_lr=[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(2),.5*(ModelRecord.IceBottom_c(1:end-1)+ModelRecord.IceBottom_c(2:end)),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)];
    
    % Interpolate geometry onto the grounding line: (hydrostatic assumption)
    lastgroundedind=find(X_c<=ModelRecord.x_gl,1,'last');
    bedelev_gl=interp1(X_input,BedElev_input+SillThick_input,ModelRecord.x_gl,interpmethod);
    if lastgroundedind<xsize
        icethick_gl=(sealevel-bedelev_gl)*(rho_sw/rho_i);
        hasshelf=1;
    else
        icethick_gl=Icethick_lr(end);
        hasshelf=0;
    end
    surfelev_gl=bedelev_gl+icethick_gl;
    
    % Compute all grounding lines (in addition to the most seaward one):
    % Compute hydraulic head:
    HydroHead_c=BedElev_c+(rho_i/rho_sw)*ModelRecord.Icethick_c;
    % Locate sea level crossings:
    sealevelcrossinginds=find(diff(sign(HydroHead_c-sealevel))); % index of last grid cell before the crossing
    % Check for multiple grounding lines:
    if length(sealevelcrossinginds)>1
        % Flag that there are multiple grounding lines:
        hasmultiplegls=1;
        % Compute number of grounding lines:
        numgls=length(sealevelcrossinginds);
        % Pre-allocate multiple grounding lines:
        X_gls=zeros(numgls,1);
        SurfElev_gls=zeros(numgls,1);
        BedElev_gls=zeros(numgls,1);
        % Loop through grounding lines:
        for thisgl=1:numgls
            % Locate this grounding line precisely:
            X_gls(thisgl)=interp1(HydroHead_c(max(sealevelcrossinginds(thisgl)-4,1):min(sealevelcrossinginds(thisgl)+4,xsize)),...
                X_c(max(sealevelcrossinginds(thisgl)-4,1):min(sealevelcrossinginds(thisgl)+4,xsize)),sealevel,'pchip',ModelRecord.domainwidth);
            % Interpolate bed topography onto this grounding line:
            BedElev_gls(thisgl)=interp1(X_input,BedElev_input+SillThick_input,X_gls(thisgl),'pchip');
            % Compute surface elevation from hydrostatic assumption:
            SurfElev_gls(thisgl)=sealevel+(1-rho_i/rho_sw)*(sealevel-BedElev_gls(thisgl));
        end
    else
        hasmultiplegls=0;
    end
    
    
    % Create color patches:
    
    % ICE PATCH:;
    IcePatch=struct('Vertices',[[X_lr,X_lr]'/1000,[IceBottom_lr,SurfElev_lr]'],...
        'Faces',[linspace(1,xsize,xsize)',linspace(2,xsize+1,xsize)',linspace(xsize+3,2*xsize+2,xsize)',linspace(xsize+2,2*xsize+1,xsize)'],...
        'FaceVertexCData',zeros(xsize,3));
    % Loop through columns to set ice patch color:
    for ii=1:xsize
        if HydroHead_c(ii)>=sealevel
            IcePatch.FaceVertexCData(ii,:)=icecolor;
        else
            IcePatch.FaceVertexCData(ii,:)=shelfcolor;
        end
    end
    
    % OCEAN PATCH:
    % Create the ocean grid:
    % Vertical grid:
    Zhat_ud_ocean=linspace(0,1,oceanzsize+1)';
    Zhat_c_ocean=.5*(Zhat_ud_ocean(1:end-1)+Zhat_ud_ocean(2:end));
    % Interpolate ocean properties in time:
    if oceantimedependence
        ThisTemperature_input=interp1(Time_yr_input,Temperature_input',snapshottimes(snapshot),'linear');
        if snapshottimes(snapshot)>Time_yr_input(end)
            ThisTemperature_input=Temperature_input(:,end);
        end
    else
        ThisTemperature_input=Temperature_input(:,1);
    end
    % Pre-allocate ocean temperature grid:
    Temp_lrud_ocean=zeros(oceanzsize+1,inputxsize);
    % Loop through grid columns to assign ocean temperature:
    ThisTemperature_input1=ThisTemperature_input;
    for ii=1:inputxsize
        % Find sill depth:
        silldepth=sealevel-max(BedElev_input(min(ii+1,inputxsize):end)+SillThick_input(min(ii+1,inputxsize):end));
        abovesillind=find(Depth_input<=silldepth,1,'last');
        % Interpolate water properties:
        if abovesillind>1
            % Compute partial blocking:
            ThisTemperature_input1(Depth_input>silldepth)=...
                sillblockingfraction*ThisTemperature_input(abovesillind)+(1-sillblockingfraction)*ThisTemperature_input(Depth_input>silldepth);
            % Interpolate:
            Temp_lrud_ocean(:,ii)=interp1(Depth_input,ThisTemperature_input1-tmelt,sealevel-(BedElev_input(ii)+SillThick_input(ii)+Zhat_ud_ocean*(sealevel-(BedElev_input(ii)+SillThick_input(ii)))),'linear',ThisTemperature_input1(end)-tmelt);
        else
            Temp_lrud_ocean(:,ii)=tmelt;
        end
    end
    % Create patch face indexing for ocean patch:
    Faces=zeros((inputxsize-1)*oceanzsize,4);
    Faces(1:oceanzsize,:)=[linspace(1,oceanzsize,oceanzsize)',linspace(2,oceanzsize+1,oceanzsize)',linspace(oceanzsize+3,2*oceanzsize+2,oceanzsize)',linspace(oceanzsize+2,2*oceanzsize+1,oceanzsize)'];
    for ii=2:inputxsize-1
        Faces((ii-1)*oceanzsize+1:ii*oceanzsize,:)=repmat(Faces((ii-1)*oceanzsize,:)+1,[oceanzsize,1])+repmat(linspace(1,oceanzsize,oceanzsize)',[1,4]);
    end
    % Create ocean patch structure:
    OceanPatch=struct('Vertices',[reshape(repmat(X_input/1000,[oceanzsize+1,1]),[],1),reshape(repmat(BedElev_input+SillThick_input,[oceanzsize+1,1])+repmat(Zhat_ud_ocean,[1,inputxsize]).*repmat(sealevel-(BedElev_input+SillThick_input),[oceanzsize+1,1]),[],1),zeros(inputxsize*(oceanzsize+1),1)],...
        'Faces',Faces,'FaceVertexCData',zeros(inputxsize*(oceanzsize+1),3));
    % Interpolate RGB data for ocean patch structure:
    OceanPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,1),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
    OceanPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,2),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
    OceanPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,3),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
    
    % SKY PATCH:
    % Interpolate surface elevation onto input grid:
    ThisSurfElev_input=interp1(X_lr,SurfElev_lr,X_input,'linear',sealevel);
    % Create patch face indexing for sky patch:
    Faces=zeros((inputxsize-1)*skyzsize,4);
    Faces(1:skyzsize,:)=[linspace(1,skyzsize,skyzsize)',linspace(2,skyzsize+1,skyzsize)',linspace(skyzsize+3,2*skyzsize+2,skyzsize)',linspace(skyzsize+2,2*skyzsize+1,skyzsize)'];
    for ii=2:inputxsize-1
        Faces((ii-1)*skyzsize+1:ii*skyzsize,:)=repmat(Faces((ii-1)*skyzsize,:)+1,[skyzsize,1])+repmat(linspace(1,skyzsize,skyzsize)',[1,4]);
    end
    % Create sky patch structure:
    %SkyPatch=struct('Vertices',zeros(skyzsize+3,2),'Faces',Faces,'FaceVertexCData',repmat(skycolor,[(skyzsize+1)*inputxsize,1])+repmat([1,1,1]-skycolor,[(skyzsize+1)*inputxsize,1]).*repmat(reshape(repmat(linspace(0,1,skyzsize+1)',[1,inputxsize]),[],1),[1,3]));
    SkyPatch=struct('Vertices',zeros((skyzsize+1)*inputxsize,2),'Faces',Faces,'FaceVertexCData',zeros((skyzsize+1)*inputxsize,3));
    % Assign values to sky patch:
    SkyPatch.Vertices=[reshape(repmat(X_input/1000,[skyzsize+1,1]),[],1),reshape(repmat(ThisSurfElev_input,[skyzsize+1,1])+repmat(linspace(0,1,skyzsize+1)',[1,inputxsize]).*repmat(elevlims(2)-ThisSurfElev_input,[skyzsize+1,1]),[],1)];
    SkyPatch.FaceVertexCData=repmat(skycolor,[(skyzsize+1)*inputxsize,1])+repmat([1,1,1]-skycolor,[(skyzsize+1)*inputxsize,1]).*repmat((SkyPatch.Vertices(:,2)-sealevel)/(elevlims(2)-sealevel),[1,3]);
    
    % BED PATCH:
    % Create bed patch:
    BedPatch=struct('Vertices',[[X_input,maxdomainwidth,0]'/1000,[BedElev_input,elevlims(1),elevlims(1)]'],...
        'Faces',linspace(1,length(X_input)+2,length(X_input)+2),'FaceVertexCData',bedcolor);
    
    % SILL PATCH:
    if dosill
        SillPatch=struct('Vertices',[[X_input,fliplr(X_input)]'/1000,[BedElev_input,fliplr(BedElev_input+SillThick_input)]'],...
            'Faces',linspace(1,2*length(X_input),2*length(X_input)),'FaceVertexCData',sillcolor);
    end
    
    % Make a subplot:
    subplot('Position',Boxes{snapshot+1})
    
    % Plot ocean patch:
    patch(OceanPatch,'FaceColor','interp','EdgeColor','none')
    hold on
    
    % Plot bed patch:
    patch(BedPatch,'FaceColor','flat','EdgeColor','none')
    
    % Plot sill patch and line:
    if dosill
        patch(SillPatch,'FaceColor','flat','EdgeColor','none')
        plot(X_input/1000,BedElev_input+SillThick_input,'Color','k','LineWidth',bedlinewidth)
    end
    
    % Plot sky:
    patch(SkyPatch,'FaceColor','interp','EdgeColor','none')
    
    % Plot ice patch:
    patch(IcePatch,'FaceColor','flat','EdgeColor','none')
    
    % Plot bed line:
    plot(X_input/1000,BedElev_input,'Color','k','LineWidth',bedlinewidth)
    
    % Plot ice surface:
    plot([0,X_c,ModelRecord.domainwidth]/1000,[SurfElev_lr(1),SurfElev_c,SurfElev_lr(end)],'Color','k','LineWidth',surflinewidth)
    
    % Plot ice bottom:
    plot([0,X_c,ModelRecord.domainwidth]/1000,[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(1),ModelRecord.IceBottom_c,1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],'Color','k','LineWidth',bedlinewidth)
    
    % Plot ice front:
    plot(ModelRecord.domainwidth*[1,1]/1000,[SurfElev_lr(end),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],'Color','k','LineWidth',frontlinewidth)
    
    % Plot grounding line:
    plot(ModelRecord.x_gl*[1,1]/1000,[bedelev_gl,surfelev_gl],'Color','k','LineWidth',gllinewidth)
    
    % Plot multiple grounding lines:
    if hasmultiplegls
        for thisgl=1:numgls
            plot(X_gls(thisgl)*[1,1]/1000,[BedElev_gls(thisgl),SurfElev_gls(thisgl)],'Color','k','LineWidth',gllinewidth)
        end
    end
    
    % Plot input profiles:
    if snapshottimes(snapshot)>0
        if ploticsurface
            plot(X_input(1:lasticeind_input)/1000,SurfElev_input(1:lasticeind_input),'Color','k','LineWidth',inputlinewidth,'LineStyle','--')
        end
        if ploticbottom
            plot(X_input(1:lasticeind_input)/1000,IceBottom_input(1:lasticeind_input),'Color','k','LineWidth',inputlinewidth,'LineStyle','--')
        end
        if ploticgl
            plot(x_gl_input*[1,1]/1000,[bedelev_gl_input,surfelev_gl_input],'Color','k','LineWidth',inputlinewidth,'LineStyle','--')
        end
        if ploticfront
            plot(X_input(lasticeind_input)*[1,1]/1000,[IceBottom_input(lasticeind_input),SurfElev_input(lasticeind_input)],'Color','k','LineWidth',inputlinewidth,'LineStyle','--')
        end
    end
    
    % Plot ELA:
    if elalinewidth>0
        elainterceptx=X_c(find(SurfElev_c>=ela,1,'last'));
        plot([elainterceptx/1000,xlims(2)],ela*[1,1],'Color','k','LineStyle','-','LineWidth',elalinewidth)
        text(.5*(xlims(2)+elainterceptx/1000),ela,'ELA','Color','k','FontSize',fontsize,'VerticalAlignment','bottom','HorizontalAlignment','center')
    end
    
    % Organize axes:
    set(gca,'XLim',xlims)
    set(gca,'XTick',xtick*[ceil(xlims(1)/xtick):1:floor(xlims(2)/xtick)])
    if snapshot==numsnapshots || snapshot==floor(numsnapshots/2)
        xlabel('Distance (km)','FontSize',fontsize)
    else
        set(gca,'XTickLabel',[])
    end
    set(gca,'YLim',elevlims)
    set(gca,'YTick',elevtick*[ceil(elevlims(1)/elevtick):1:floor(elevlims(2)/elevtick)])
    if elevkm
        set(gca,'YTickLabel',elevtick*[ceil(elevlims(1)/elevtick):1:floor(elevlims(2)/elevtick)]/1000)
        ylabel('Elevation (km)','FontSize',fontsize)
    else
        ylabel('Elevation (m)','FontSize',fontsize)
    end
    
    set(gca,'TickDir',snapshottickdir)
    set(gca,'XMinorTick',snapshotminorticksx)
    set(gca,'YMinorTick',snapshotminorticksy)
    %set(gca,'TickLength',[ticklength,1])
    set(gca,'LineWidth',borderlinewidth)
    set(gca,'Box','on')
    set(gca,'FontSize',fontsize)
    if snapshot>floor(numsnapshots/2)
        set(gca,'YAxisLocation','right')
    end
    
    % Make title:
    title([char(97+snapshot),') ',titletext{snapshot},', t=',num2str(round(snapshottimes(snapshot))),' yr'],'FontSize',fontsize)
    %title([char(65+snapshot),') ',titletext{snapshot},', t=',num2str(round(snapshottimes(snapshot))),' yr'],'FontSize',fontsize)
    
    % Draw:
    drawnow
    
    %     % Text label in the corner:
    %     if dozoom && snapshot<numsnapshots
    %         text(zoomxlims(2),zoomelevlims(2),[char(98+snapshot),')'],'HorizontalAlignment','right','VerticalAlignment','top','Color','k','FontSize',snapshotfontsize)
    %     else
    %         text(xlims(2),elevlims(2),[char(98+snapshot),')'],'HorizontalAlignment','right','VerticalAlignment','top','Color','k','FontSize',snapshotfontsize)
    %     end
    
    % Construct figure name:
    % figname=[outputfolder,mainfigname,'_snapshot',num2str(snapshot),'.png'];
    
    % Save figure:
    %     if dozoom && snapshot<numsnapshots
    %         set(gcf,'PaperSize',[zoompagewidth,zoompageheight])
    %         set(gcf,'PaperPosition',[0,0,zoompagewidth,zoompageheight])
    %     else
    %         set(gcf,'PaperSize',[snapshotpagewidth,snapshotpageheight])
    %         set(gcf,'PaperPosition',[0,0,snapshotpagewidth,snapshotpageheight])
    %     end
    %    print('-dpng',figname,['-r',num2str(resolution)])
    
end

% Save figure:
set(gcf,'PaperSize',[pagewidth,pageheight])
set(gcf,'PaperPosition',[0,0,pagewidth,pageheight])
print('-dpng',[outputfolder,mainfigname,'.png'],['-r',num2str(resolution)])

% Final Display:
disp('Done!')
toc