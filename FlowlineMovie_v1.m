%function FTM_2D_display_v16(inputfile)

% FlowlineMovie_v1

% Mike Wolovick, 5/17/2016

% This script makes a movie displaying the results of Flowline_v1.  The
% script is based on FTM_2D_display_v16

% The movie shows several colored ribbons along the boundaries representing
% BC, colored ocean temperature, and vertical lines that move with the ice
% flow.  The vertical lines advect with column-average velocity and their
% advection is computed here (not in the model itself).  New lines are
% spawned whenever the spacing between old lines doubles.

clear variables
close all
tic

%% Parameters:

% File names and folders:
inputfolder='/work/Michael.Wolovick/FjordSillPaper/ModelOutput/';
inputfiles={'ThwaitesC_v2_s4_m03_ch_v3a.mat';...
    'ThwaitesC_v2_s5_m03_ch_v3a.mat';...
    'ThwaitesC_v2_s6_m03_ch_v3a.mat';...
    'ThwaitesC_v2_w_m03_ch_v3.mat';...
    'ThwaitesC_v2_s_m03_ch_v3.mat'};
outputfolder='/net/mjw/FjordSillPaper/ModelOutput_movies/';
outputsuffix='_movie.avi';

% Short model description:
titletext={'50% Water Blocking, Small Sill';...
    '0% Water Blocking, Small Sill';...
    '100% Water Blocking, Small Sill';...
    'No Sill, Warming Climate';...
    '100% Water Blocking, Large Sill'};  % cell array of strings, no trailing comma necessary
titleinterpreter='none';        % 'tex', 'latex', or 'none'
titlefontsize=12;               % points
fontsize=12;                    % points

% Model times:
starttime=0;                    % yr
reftime=0;                      % yr
endtime=1000;                   % yr

% Movie settings:
fps=18;                         % 1/s
moviespeed=4;                   % yr/s
icholdtime=1;                   % s
finalholdtime=.5;               % s

% Ticks and limits:
xtick=50;                       % km
elevlims=[-1500,2000];          % m
elevtick=500;                   % m

% Color ribbon thickness:
ribbonthick=0;                  % m (measured vertically, 0 deactives)

% Vertical lines for flow:
numvertlines=25;                % integer
vertlinesstoragebuffer=2;       % unitless>1

% Colormaps and limits:
smbcmap='fadedjet';             % valid colormap string or 'fadedjet'
oceantempcmap='parula';         % valid colormap string or 'fadedjet'
oceanmeltcmap='jet';            % valid colormap string or 'fadedjet'
basaldragcmap='cool';           % valid colormap string or 'fadedjet'
smblims=[-5,2];                 % [1x2] m/yr
oceantemplims='auto';           % [1x2] deg C or 'auto'
meltcalvinglims=[1,1e4];        % [1x2] m/yr
basaldraglims=[0,200];          % [1x2] kPa 

% Line widths:
bedlinewidth=2;                 % points
surflinewidth=2;                % points
gllinewidth=2;                  % points
frontlinewidth=2;               % points
vertlinewidth=1;                % points
otherlinewidth=1;               % points
inputlinewidth=1;               % points

% Block color settings:
skycolor=[.75,.75,1];           % [R,G,B]
skyzsize=10;                    % integer
bedcolor=[.63,.32,.18];         % [R,G,B]
icecolor=[.45,.62,.78];         % [R,G,B]
sillcolor=[.5,.5,.5];           % [R,G,B]
shelfcolor=[.875,0,.875];       % [R,G,B]
%oceanxsize=200;                 % integer
oceanzsize=50;                  % integer

% Size:
figbox=[108,164,1414,843];      % [left, bottom, width, height], in pixels, corresponds to 'Position', not 'OuterPosition'
plotbox=[.075,.16,.87,.75];     % [left, bottom, width, height], unitless

%% Preparation:

% Check that title text is compatable with input files:
if length(inputfiles)~=length(titletext)
    error('"inputfiles" and "titletext" must be the same length.')
end

% Loop through model runs:
for runnum=1:length(inputfiles)

% Assign this output file:
outputfile=[outputfolder,inputfiles{runnum}(1:end-4),outputsuffix];

% Display this input file:
disp(inputfiles{runnum})

% Determine glacier name and parameter values:
underscorenum=strfind(inputfiles{runnum},'_');
glaciername=inputfiles{runnum}(1:underscorenum(2)-1);
mvalue=str2double(inputfiles{runnum}(underscorenum(3)+2:underscorenum(3)+3));
ctype=inputfiles{runnum}(underscorenum(4)+(2:3));

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

% Make a figure:
figure(1)

% Make faded jet colormap:
FadedJet=colormap('jet');
FadedJet(3*round(size(FadedJet,1)/8)+1:round(size(FadedJet,1)/2),1)=linspace(0,1,round(size(FadedJet,1)/8))';
FadedJet(3*round(size(FadedJet,1)/8)+1:round(size(FadedJet,1)/2),3)=1;
FadedJet(round(size(FadedJet,1)/2)+1:5*round(size(FadedJet,1)/8),1)=1;
FadedJet(round(size(FadedJet,1)/2)+1:5*round(size(FadedJet,1)/8),3)=linspace(1,0,round(size(FadedJet,1)/8))';

% Assign colormaps:
if strcmp(smbcmap,'fadedjet')
    SMBcmap=FadedJet;
else
    SMBcmap=colormap(smbcmap);
end
if strcmp(oceantempcmap,'fadedjet')
    OceanTempcmap=FadedJet;
else
    OceanTempcmap=colormap(oceantempcmap);
end
if strcmp(oceanmeltcmap,'fadedjet')
    OceanMeltcmap=FadedJet;
else
    OceanMeltcmap=colormap(oceanmeltcmap);
end
if strcmp(basaldragcmap,'fadedjet')
    BasalDragcmap=FadedJet;
else
    BasalDragcmap=colormap(basaldragcmap);
end

% Convert color limits to log scale:
if runnum==1
    meltcalvinglims=log(meltcalvinglims);
    % basaldraglims=log(basaldraglims);
end

% Load and unpack model output:
load([inputfolder,inputfiles{runnum}],'*Parameters','*_input','Zhat_*','DZhat_*','numrecords','numdigits')
unpack(GridParameters)
unpack(SideParameters)
unpack(PlumeParameters)
unpack(TimingParameters)
unpack(TracerParameters)
unpack(ThermalParameters)
unpack(BottomParameters)
unpack(OtherParameters)
unpack(ModelParameters)
if exist('SillParameters','var')
    unpack(SillParameters)
    if exist('sillblockingfraction','var')==0
        sillblockingfraction=1;
    end
else
    dosill=0;
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
xlims=[0,maxdomainwidth]/1000;

% Get first model record:
if starttime==0
    % Start record is IC:
    startrecord=0;
    % Load IC:
    load([inputfolder,inputfiles{runnum}],'InitialConditions')
    ModelRecord=InitialConditions;
    clear InitialConditions
    % Assign things that aren't recorded in the IC:
    % time is zero:
    ModelRecord.time_yr=0;
    % zero surface melt:
    ModelRecord.MeltRate_u=zeros(1,xsize);
    % interpolate basal drag:
    if strcmp(slidingstressscale,'file')
        X_lr=linspace(0,ModelRecord.domainwidth,xsize+1);
        ModelRecord.Drag_lrd=interp1(X_input,SlidingStressScale_input,X_lr,'linear');
    else
        ModelRecord.Drag_lrd=slidingstressscale*ones(1,xsize+1);
    end
    % set plume melt rate to bottom of color limits:
    ModelRecord.MeltRate_c_plume=meltcalvinglims(1)*ones(1,xsize_plume+zsize_plume);
    % set calving rate to zero:
    ModelRecord.calvingrate_r=0;
else
    % Determine start record:
    startrecord=round(starttime/recordinterval_yr+.5);  % assumes that record time is the middle of the record
    % Load this model record:
    prefix='0'*ones(1,numdigits(end)-floor(log10(startrecord))-1);
    idnumber=num2str(startrecord);
    load([inputfolder,inputfiles{runnum}],['ModelRecord_',prefix,idnumber])
    % Check if model record exists:
    if exist(['ModelRecord_',prefix,idnumber],'var')
        % Rename model record:
        eval(['ModelRecord=ModelRecord_',prefix,idnumber,';'])
        clear(['ModelRecord_',prefix,idnumber])
    else
        error('Starting record does not exist.')
    end
end

% Determine end record:
endrecord=round(endtime/recordinterval_yr+.5)-1;  % assumes that record time is the middle of the record

% Assign sea level:
if strcmp(rightbctype,'front') && strcmp(rightbcparam,'file')==0
    sealevel=rightbcparam;
elseif strcmp(rightbctype,'front')
    sealevel=interp1(Time_yr_input,RightBCParam_input,starttime,'linear',RightBCParam_input(1));
else
    sealevel=min(BedElev_input)-1;
end

% Define sill profile:
if dosill
    if isfield(ModelRecord,'SillThick_input')
        SillThick_input=ModelRecord.SillThick_input;
    else
        % Compute sill thickness profile:
        if sillstarttime_yr==0 && sillconstructiontime_yr==0 && starttime==0
            SillThick_input=sillamp*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2));
        else
            SillThick_input=sillamp*max(0,min(1,(starttime-sillstarttime_yr)/(sillconstructiontime_yr+1)))*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2)); % note stabilizer
        end
    end
else
    SillThick_input=zeros(size(X_input));
end

% Define initial domain:
dx=ModelRecord.domainwidth/xsize;
X_lr=linspace(0,ModelRecord.domainwidth,xsize+1);
X_c=linspace(.5*dx,ModelRecord.domainwidth-.5*dx,xsize);
BedElev_lr=interp1(X_input,BedElev_input+SillThick_input,X_lr,interpstyle);
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
    bedelev_gl_input=interp1(X_input,BedElev_input,x_gl_input,interpstyle,BedElev_input(end));
    icethick_gl_input=(sealevel-bedelev_gl_input)*(rho_sw/rho_i);
    hasshelf_input=1;
else
    x_gl_input=X_input(lastgroundedind_input);
    bedelev_gl_input=interp1(X_input,BedElev_input,x_gl_input,interpstyle,BedElev_input(end));
    icethick_gl_input=Icethick_input(lastgroundedind_input);
    hasshelf_input=0;
end
surfelev_gl_input=bedelev_gl_input+icethick_gl_input;

% Define vertical grid:
Zhat_ud=linspace(0,1,zsize+1)'/maxdensify+(1-1/maxdensify)*linspace(0,1,zsize+1)'.^densifypower;
Zhat_c=.5*(Zhat_ud(1:end-1)+Zhat_ud(2:end));
DZhat_c=Zhat_ud(2:end)-Zhat_ud(1:end-1);
DZhat_ud=[Zhat_c(1);Zhat_c(2:end)-Zhat_c(1:end-1);1-Zhat_c(end)];  % first and last are half-cells

% Define size of input vectors:
inputxsize=length(X_input);

% Compute record skipping:
recordsperframe=ceil(moviespeed/(fps*recordinterval_yr));

% Compute time per frame and true movie speed:
timeperframe=recordinterval_yr*recordsperframe;  % yr
truemoviespeed=timeperframe*fps;

% Calculate surface elevation:
SurfElev_c=ModelRecord.IceBottom_c+ModelRecord.Icethick_c;
SurfElev_lr=[1.5*SurfElev_c(1)-.5*SurfElev_c(2),.5*(SurfElev_c(1:end-1)+SurfElev_c(2:end)),1.5*SurfElev_c(end)-.5*SurfElev_c(end-1)]; % linear extrapolation

% Interpolate ice thickness to grid edges:
Icethick_lr=[1.5*ModelRecord.Icethick_c(1)-.5*ModelRecord.Icethick_c(2),.5*(ModelRecord.Icethick_c(1:end-1)+ModelRecord.Icethick_c(2:end)),1.5*ModelRecord.Icethick_c(end)-.5*ModelRecord.Icethick_c(end-1)];

% Interpolate ice bottom to grid edges:
IceBottom_lr=[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(2),.5*(ModelRecord.IceBottom_c(1:end-1)+ModelRecord.IceBottom_c(2:end)),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)];

% Interpolate geometry onto the grounding line: (hydrostatic assumption)
lastgroundedind=find(X_c<=ModelRecord.x_gl,1,'last');
bedelev_gl=interp1(X_lr,BedElev_lr,ModelRecord.x_gl,interpstyle,BedElev_lr(end));
if lastgroundedind<xsize
    icethick_gl=(sealevel-bedelev_gl)*(rho_sw/rho_i);
    hasshelf=1;
else
    icethick_gl=Icethick_lr(end);
    hasshelf=0;
end
surfelev_gl=bedelev_gl+icethick_gl;

% Create initial distribution of vertical lines:
totalnumvertlines=round(numvertlines*vertlinesstoragebuffer);
numactivevertlines=numvertlines;
dx_vertlines=ModelRecord.domainwidth/numvertlines;
X_vertlines=[linspace(0,ModelRecord.domainwidth,numvertlines),NaN*zeros(1,totalnumvertlines-numvertlines)];

% Interpolate ice bottom and surface to the vertical lines:
Z_vertlines=[interp1(X_lr,IceBottom_lr,X_vertlines,interpstyle);interp1(X_lr,SurfElev_lr,X_vertlines,interpstyle)];

% Compute vertical exaggeration:
length_h=figbox(3)*plotbox(3);
length_v=figbox(4)*plotbox(4);
density_h=(xlims(2)-xlims(1))/length_h;
density_v=.001*(elevlims(2)-elevlims(1))/length_v;
exag=density_h/density_v;

% Compute surface and bed gradients and angles:
SurfGrad_lr=gradient(SurfElev_lr,dx);
BedGrad_c=[.5*(BedElev_c(1)-BedElev_lr(1)),(BedElev_lr(2:end)-BedElev_lr(1:end-1)),.5*(BedElev_lr(end)-BedElev_c(end))]/dx; % also includes domain edges
SurfAngle_lr=atan(exag*SurfGrad_lr); % angle on exaggerated plot
BedAngle_c=atan(exag*BedGrad_c); % angle on exaggerated plot

% Compute all grounding lines (in addition to the most seaward one):
% Compute hydraulic head:
HydroHead_c=BedElev_c+(rho_i/rho_sw)*ModelRecord.Icethick_c;
% Locate sea level crossings:
sealevelcrossinginds=find(diff(sign(HydroHead_c-sealevel))); % index of last grid cell before the crossing
% Check for multiple grounding lines:
if length(sealevelcrossinginds)>1
    % Flag that there are multiple grounding lines:
    hasmultiplegls=1;
    hasmultglslast=1;
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
    hasmultglslast=0;
end


% ICE PATCH:
% % Create ice patch:
% IcePatch=struct('Vertices',[[X_lr,fliplr(X_lr)]'/1000,[IceBottom_lr,fliplr(SurfElev_lr)]'],...
%     'Faces',linspace(1,2*(xsize+1),2*(xsize+1)),'FaceVertexCData',icecolor);
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


% SMB PATCH:
% Compute surface ribbon coordinates:
SurfRibbonTopZ_lr=SurfElev_lr+ribbonthick*cos(SurfAngle_lr);
SurfRibbonTopX_lr=X_lr-exag*ribbonthick*sin(SurfAngle_lr);
% Create surface ribbon patch structure:
SurfRibbonPatch=struct('Vertices',[[X_lr,SurfRibbonTopX_lr]'/1000,...
    [SurfElev_lr,SurfRibbonTopZ_lr]'],...
    'Faces',[linspace(1,xsize,xsize)',linspace(2,xsize+1,xsize)',linspace(xsize+3,2*xsize+2,xsize)',linspace(xsize+2,2*xsize+1,xsize)'],...
    'FaceVertexCData',zeros(xsize,3));
% Flip SMB sign and scale to be symmetric about 0:
ScaledSMB=(ModelRecord.Accum_u-ModelRecord.MeltRate_u)*secondsperyear;
ScaledSMB(ScaledSMB<0)=max(-1,ScaledSMB(ScaledSMB<0)/(-smblims(1))); % assumes smblims(1) is negative
ScaledSMB(ScaledSMB>0)=min(1,ScaledSMB(ScaledSMB>0)/smblims(2));
ScaledSMB=-ScaledSMB;
% Interpolate surface patch RGB data:
SurfRibbonPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(SMBcmap,1))',SMBcmap(:,1),(ScaledSMB'+1)/2,interpstyle)));
SurfRibbonPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(SMBcmap,1))',SMBcmap(:,2),(ScaledSMB'+1)/2,interpstyle)));
SurfRibbonPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(SMBcmap,1))',SMBcmap(:,3),(ScaledSMB'+1)/2,interpstyle)));


% BASAL DRAG PATCH:
% Compute bed ribbon coordinates:
BedRibbonBotZ_c=[BedElev_lr(1),BedElev_c,BedElev_lr(end)]-ribbonthick*cos(BedAngle_c); % also includes domain edges
BedRibbonBotX_c=[0,X_c,ModelRecord.domainwidth]+exag*ribbonthick*sin(BedAngle_c); % also includes domain edges
% Create bed ribbon patch structure:
if hasshelf
    BedRibbonPatch=struct('Vertices',[[0,X_c(1:lastgroundedind+1),BedRibbonBotX_c(1:lastgroundedind+2)]'/1000,...
        [BedElev_lr(1),BedElev_c(1:lastgroundedind+1),BedRibbonBotZ_c(1:lastgroundedind+2)]'],...
        'Faces',[linspace(1,lastgroundedind+1,lastgroundedind+1)',linspace(2,lastgroundedind+2,lastgroundedind+1)',linspace(lastgroundedind+4,2*lastgroundedind+4,lastgroundedind+1)',linspace(lastgroundedind+3,2*lastgroundedind+3,lastgroundedind+1)'],...
        'FaceVertexCData',zeros(lastgroundedind+1,3));
else
    BedRibbonPatch=struct('Vertices',[[0,X_c(1:lastgroundedind),ModelRecord.domainwidth,BedRibbonBotX_c(1:lastgroundedind+2)]'/1000,...
        [BedElev_lr(1),BedElev_c(1:lastgroundedind),BedElev_lr(end),BedRibbonBotZ_c(1:lastgroundedind+2)]'],...
        'Faces',[linspace(1,lastgroundedind+1,lastgroundedind+1)',linspace(2,lastgroundedind+2,lastgroundedind+1)',linspace(lastgroundedind+4,2*lastgroundedind+4,lastgroundedind+1)',linspace(lastgroundedind+3,2*lastgroundedind+3,lastgroundedind+1)'],...
        'FaceVertexCData',zeros(lastgroundedind+1,3));
end
% Interpolate bed patch RGB data:
BedRibbonPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(BasalDragcmap,1))',BasalDragcmap(:,1),max(0,min(1,(ModelRecord.Drag_lrd(1:lastgroundedind+1)'/1000-basaldraglims(1))/(basaldraglims(2)-basaldraglims(1)))),interpstyle)));
BedRibbonPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(BasalDragcmap,1))',BasalDragcmap(:,2),max(0,min(1,(ModelRecord.Drag_lrd(1:lastgroundedind+1)'/1000-basaldraglims(1))/(basaldraglims(2)-basaldraglims(1)))),interpstyle)));
BedRibbonPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(BasalDragcmap,1))',BasalDragcmap(:,3),max(0,min(1,(ModelRecord.Drag_lrd(1:lastgroundedind+1)'/1000-basaldraglims(1))/(basaldraglims(2)-basaldraglims(1)))),interpstyle)));

% OCEAN PATCH:
% Create the ocean grid:
% Vertical grid:
Zhat_ud_ocean=linspace(0,1,oceanzsize+1)';
Zhat_c_ocean=.5*(Zhat_ud_ocean(1:end-1)+Zhat_ud_ocean(2:end));
% Interpolate ocean properties in time:
if oceantimedependence
    ThisTemperature_input=interp1(Time_yr_input,Temperature_input',starttime,interpstyle);
    if starttime<Time_yr_input(1)
        ThisTemperature_input=Temperature_input(:,1);
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
OceanPatch=struct('Vertices',[reshape(repmat(X_input/1000,[oceanzsize+1,1]),[],1),reshape(repmat(BedElev_input,[oceanzsize+1,1])+repmat(Zhat_ud_ocean,[1,inputxsize]).*repmat(sealevel-BedElev_input,[oceanzsize+1,1]),[],1),zeros(inputxsize*(oceanzsize+1),1)],...
    'Faces',Faces,'FaceVertexCData',zeros(inputxsize*(oceanzsize+1),3));
% Interpolate RGB data for ocean patch structure:
OceanPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,1),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
OceanPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,2),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
OceanPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,3),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
    
% % OCEAN PATCH:
% % Create the ocean grid:
% % Vertical grid:
% Zhat_ud_ocean=linspace(0,1,oceanzsize+1)';
% Zhat_c_ocean=.5*(Zhat_ud_ocean(1:end-1)+Zhat_ud_ocean(2:end));
% % Check whether bed topography at the grounding line is above sea level:
% ind1=find(BedElev_input+SillThick_input<sealevel&X_input>ModelRecord.x_gl,1,'first');
% ind2=find(X_input>ModelRecord.x_gl,1,'first');
% if ind1==ind2
%     % The ocean starts at the grounding line:
%     oceanstartx=ModelRecord.x_gl;
% else
%     % The oceans starts where the bed goes below sea level:
%     oceanstartx=interp1(BedElev_input(ind1-1:ind1)+SillThick_input(ind1-1:ind1),X_input(ind1-1:ind1),sealevel,interpstyle);
% end
% % Create ocean horizontal grid:
% X_lr_ocean=linspace(oceanstartx,maxdomainwidth,oceanxsize+1);
% % Create ocean top and bottom:
% BedElev_lr_ocean=interp1(X_input,BedElev_input+SillThick_input,X_lr_ocean);
% % Interpolate ocean properties in time:
% if oceantimedependence
%     ThisTemperature_input=interp1(Time_yr_input,Temperature_input',starttime,interpstyle);
%     if starttime>Time_yr_input(end)
%         ThisTemperature_input=Temperature_input(:,end);
%     end
% else
%     ThisTemperature_input=Temperature_input(:,1);
% end
% % Pre-allocate:
% Temp_lrud_ocean=zeros(oceanzsize+1,oceanxsize+1);
% % Loop through grid columns to assign ocean temperature:
% for ii=1:oceanxsize+1
%     % Find sill depth:
%     silldepth=sealevel-max(BedElev_lr_ocean(min(ii+1,oceanxsize+1):end));
%     abovesillind=find(Depth_input<=silldepth,1,'last');
%     % Interpolate water properties:
%     Temp_lrud_ocean(:,ii)=interp1(Depth_input(1:abovesillind),ThisTemperature_input(1:abovesillind)-tmelt,sealevel-(BedElev_lr_ocean(ii)+Zhat_ud_ocean*(sealevel-BedElev_lr_ocean(ii))),interpstyle,ThisTemperature_input(abovesillind)-tmelt);
% end
% % Create patch face indexing for ocean patch:
% Faces=zeros(oceanxsize*oceanzsize,4);
% Faces(1:oceanzsize,:)=[linspace(1,oceanzsize,oceanzsize)',linspace(2,oceanzsize+1,oceanzsize)',linspace(oceanzsize+3,2*oceanzsize+2,oceanzsize)',linspace(oceanzsize+2,2*oceanzsize+1,oceanzsize)'];
% for ii=2:oceanxsize
%     Faces((ii-1)*oceanzsize+1:ii*oceanzsize,:)=repmat(Faces((ii-1)*oceanzsize,:)+1,[oceanzsize,1])+repmat(linspace(1,oceanzsize,oceanzsize)',[1,4]);
% end
% % Create ocean patch structure:
% OceanPatch=struct('Vertices',[reshape(repmat(X_lr_ocean/1000,[oceanzsize+1,1]),[],1),reshape(repmat(BedElev_lr_ocean,[oceanzsize+1,1])+repmat(Zhat_ud_ocean,[1,oceanxsize+1]).*repmat(sealevel-BedElev_lr_ocean,[oceanzsize+1,1]),[],1),zeros((oceanxsize+1)*(oceanzsize+1),1)],...
%     'Faces',Faces,'FaceVertexCData',zeros((oceanxsize+1)*(oceanzsize+1),3));
% % Interpolate RGB data for ocean patch structure:
% OceanPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,1),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
% OceanPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,2),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
% OceanPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,3),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));


% SKY PATCH:
% Interpolate surface elevation onto input grid:
SurfElev_input=interp1(X_lr,SurfElev_lr,X_input,interpstyle,sealevel);
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
SkyPatch.Vertices=[reshape(repmat(X_input/1000,[skyzsize+1,1]),[],1),reshape(repmat(SurfElev_input,[skyzsize+1,1])+repmat(linspace(0,1,skyzsize+1)',[1,inputxsize]).*repmat(elevlims(2)-SurfElev_input,[skyzsize+1,1]),[],1)];
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


% Compute monotonic floating ice bottom:
% Set the monotonic bottom to the real ice bottom:
MonotonicIceBottom_c=ModelRecord.IceBottom_c;
% Check if we're running the plume model and have a shelf:
if hasshelf
    % Detect downward sloping ice bottom in the shelf:
    DownwardSlopingShelf_lr=X_lr>ModelRecord.x_gl&([ModelRecord.IceBottom_c,1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)]-[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(2),ModelRecord.IceBottom_c])<0;
    % Check if there are any cells that need to be fixed:
    if sum(DownwardSlopingShelf_lr)>0
        % Identify cells that need to be fixed:
        needfixinginds=find(DownwardSlopingShelf_lr);
        % Pre-allocate a record of whether they have been fixed:
        hasbeenfixed=false(size(needfixinginds));
        % Loop through from the calving front to the grounding line:
        for thisind=length(needfixinginds):-1:1
            % Check if this grid cell has been fixed:
            if hasbeenfixed(thisind)
                continue
            end
            % Identify real index of this grid edge:
            d2=needfixinginds(thisind);
            % Identify downstream ice bottom:
            if d2<xsize+1
                downstreamicebottomz=ModelRecord.IceBottom_c(d2);
            else
                downstreamicebottomz=1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1);
            end
            % Identify next ice bottom below that point:
            nextind=find(ModelRecord.IceBottom_c(1:d2-1)<downstreamicebottomz,1,'last');
            if isempty(nextind)
                continue
            end
            nexticebottom=ModelRecord.IceBottom_c(nextind);
            % Create a linear profile to the next ice bottom:
            if d2<xsize+1
                MonotonicIceBottom_c(nextind:d2-1)=linspace(nexticebottom,downstreamicebottomz-(downstreamicebottomz-nexticebottom)/(d2-nextind+1),d2-nextind);
            else
                MonotonicIceBottom_c(nextind:d2-1)=linspace(nexticebottom,downstreamicebottomz-.5*(downstreamicebottomz-nexticebottom)/(d2-nextind+1),d2-nextind);
            end
            % Flag cells that have been fixed:
            hasbeenfixed(needfixinginds>=nextind+1)=1;
        end
    end
end

% OCEAN MELT PATCH:
% Compute plume model grid size:
dx_plume=(ModelRecord.domainwidth-ModelRecord.x_gl)/xsize_plume;
% Create plume model horizontal grid:
X_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
Z_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
X_lr_plume(1:xsize_plume+1)=linspace(ModelRecord.x_gl,ModelRecord.domainwidth,xsize_plume+1);
X_lr_plume(xsize_plume+2:xsize_plume+zsize_plume+1)=ModelRecord.domainwidth;
X_c_plume=.5*(X_lr_plume(1:end-1)+X_lr_plume(2:end));
% Compute plume model elevation:
if hasshelf
    % Interpolate elevation under the floating shelf:
    Z_lr_plume(1:xsize_plume+1)=interp1([ModelRecord.x_gl,X_c(X_c>ModelRecord.x_gl),ModelRecord.domainwidth],...
        [bedelev_gl,MonotonicIceBottom_c(X_c>ModelRecord.x_gl),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],...
        X_lr_plume(1:xsize_plume+1),interpstyle);
else
    Z_lr_plume(1:xsize_plume+1)=1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1);
end
% Compute elevation along the vertical face:
Z_lr_plume(xsize_plume+1:xsize_plume+zsize_plume+1)=linspace(1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1),sealevel,zsize_plume+1);
% Interpolate to grid centers:
Z_c_plume=.5*(Z_lr_plume(1:end-1)+Z_lr_plume(2:end));
% Compute plume grid gradient and angle:
PlumeGrad_lr_plume=gradient(Z_lr_plume,dx_plume); % also includes domain edges
PlumeAngle_lr_plume=atan(exag*PlumeGrad_lr_plume); % angle on exaggerated plot
% Force vertical face to be vertical:
PlumeAngle_lr_plume(xsize_plume+2:end)=pi/2;
% Fix the angle just before the vertical face:
PlumeAngle_lr_plume(xsize_plume+1)=PlumeAngle_lr_plume(xsize_plume);
% Compute plume ribbon coordinates: (one extra point to accomodate right
% angle turn at ice front)
PlumeRibbonBotZ_lr=[Z_lr_plume-ribbonthick*cos(PlumeAngle_lr_plume),Z_lr_plume(xsize_plume+1)];
PlumeRibbonBotX_lr=[X_lr_plume+exag*ribbonthick*sin(PlumeAngle_lr_plume),X_lr_plume(xsize_plume+1)+exag*ribbonthick];
% Reduce under-shelf ribbon to zero size if there is no shelf:
if hasshelf==0
    PlumeRibbonBotZ_lr(1:xsize_plume+1)=Z_lr_plume(1:xsize_plume+1);
    PlumeRibbonBotX_lr(1:xsize_plume+1)=X_lr_plume(1:xsize_plume+1);
end
% Create plume ribbon patch structure:
PlumePatch=struct('Vertices',[[X_lr_plume,PlumeRibbonBotX_lr]'/1000,...
    [Z_lr_plume,PlumeRibbonBotZ_lr]'],...
    'Faces',[linspace(1,xsize_plume+zsize_plume,xsize_plume+zsize_plume)',...
    linspace(2,xsize_plume+zsize_plume+1,xsize_plume+zsize_plume)',...
    linspace(xsize_plume+zsize_plume+3,2*(xsize_plume+zsize_plume)+2,xsize_plume+zsize_plume)',... 
    [linspace(xsize_plume+zsize_plume+2,2*xsize_plume+zsize_plume+1,xsize_plume),2*(xsize_plume+zsize_plume)+3,linspace(2*xsize_plume+zsize_plume+3,2*(xsize_plume+zsize_plume)+1,zsize_plume-1)]'],...% inserted extra point at right angle
    'FaceVertexCData',zeros(xsize_plume+zsize_plume,3));
% Create combined melt and calving vector:
MeltPlusCalving_c_plume=ModelRecord.MeltRate_c_plume;
MeltPlusCalving_c_plume(xsize_plume+1:end)=MeltPlusCalving_c_plume(xsize_plume+1:end)+ModelRecord.calvingrate_r;
% Interpolate plume patch RGB data:
PlumePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,1),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltcalvinglims(1))))-meltcalvinglims(1))/(meltcalvinglims(2)-meltcalvinglims(1)))),interpstyle)));
PlumePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,2),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltcalvinglims(1))))-meltcalvinglims(1))/(meltcalvinglims(2)-meltcalvinglims(1)))),interpstyle)));
PlumePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,3),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltcalvinglims(1))))-meltcalvinglims(1))/(meltcalvinglims(2)-meltcalvinglims(1)))),interpstyle)));


%% IC Frame:

% Make figure:
figure(1)
set(gcf,'Position',figbox)

% Make subplot:
subplot('Position',plotbox)

% Plot bed patch:
hbedpatch=patch(BedPatch,'FaceColor','flat','EdgeColor','none');
hold on

% Plot ocean:
hocean=patch(OceanPatch,'FaceColor','interp','EdgeColor','none');

% Plot sill patch and line:
if dosill
    hsillpatch=patch(SillPatch,'FaceColor','flat','EdgeColor','none');
    hsillline=plot(X_input/1000,BedElev_input+SillThick_input,'Color','k','LineWidth',bedlinewidth);
end

% Plot sky:
hsky=patch(SkyPatch,'FaceColor','interp','EdgeColor','none');

% Plot ice patch:
hicepatch=patch(IcePatch,'FaceColor','flat','EdgeColor','none');

% Plot color ribbons:
if ribbonthick>0
    % Plot surface mass balance:
    hsurfribbon=patch(SurfRibbonPatch,'FaceColor','flat','EdgeColor','none');
    % Plot basal drag:
    hbedribbon=patch(BedRibbonPatch,'FaceColor','flat','EdgeColor','none');
    % Plot plume:
    hplumeribbon=patch(PlumePatch,'FaceColor','flat','EdgeColor','none');
    hplumeedge(1)=plot(PlumeRibbonBotX_lr(1:xsize_plume+1)/1000,PlumeRibbonBotZ_lr(1:xsize_plume+1),'Color','k','LineWidth',otherlinewidth);
    hplumeedge(2)=plot(PlumeRibbonBotX_lr(xsize_plume+2:end)/1000,PlumeRibbonBotZ_lr(xsize_plume+2:end),'Color','k','LineWidth',otherlinewidth);
end

% Plot input profiles:
hinput=zeros(4,1);
hinput(1)=plot(X_input(1:lasticeind_input)/1000,SurfElev_input(1:lasticeind_input),'Color','k','LineWidth',inputlinewidth,'LineStyle','--');
hinput(2)=plot(X_input(1:lasticeind_input)/1000,IceBottom_input(1:lasticeind_input),'Color','k','LineWidth',inputlinewidth,'LineStyle','--');
hinput(3)=plot(x_gl_input*[1,1]/1000,[bedelev_gl_input,surfelev_gl_input],'Color','k','LineWidth',inputlinewidth,'LineStyle','--');
hinput(4)=plot(X_input(lasticeind_input)*[1,1]/1000,[IceBottom_input(lasticeind_input),SurfElev_input(lasticeind_input)],'Color','k','LineWidth',inputlinewidth,'LineStyle','--');

% Plot bed line:
hbedline=plot(X_input/1000,BedElev_input,'Color','k','LineWidth',bedlinewidth);

% Plot ice surface:
hsurf=plot([0,X_c,ModelRecord.domainwidth]/1000,[SurfElev_lr(1),SurfElev_c,SurfElev_lr(end)],'Color','k','LineWidth',surflinewidth);

% Plot ice bottom:
hbot=plot([0,X_c,ModelRecord.domainwidth]/1000,[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(1),ModelRecord.IceBottom_c,1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],'Color','k','LineWidth',bedlinewidth);

% Plot ice front:
hfront=plot(ModelRecord.domainwidth*[1,1]/1000,[SurfElev_lr(end),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],'Color','k','LineWidth',frontlinewidth);

% Plot grounding line:
hgl=plot(ModelRecord.x_gl*[1,1]/1000,[bedelev_gl,surfelev_gl],'Color','k','LineWidth',gllinewidth);

% Plot multiple grounding lines:
if hasmultiplegls
    hmultgls=zeros(numgls,1);
    for thisgl=1:numgls
        hmultgls(thisgl)=plot(X_gls(thisgl)*[1,1]/1000,[BedElev_gls(thisgl),SurfElev_gls(thisgl)],'Color','k','LineWidth',gllinewidth);
    end
end

% Plot vertical lines:
hvertlines=plot(repmat(X_vertlines/1000,[2,1]),Z_vertlines,'Color','k','LineWidth',vertlinewidth);

% Make a text box describing this model experiment:
htext=text(xlims(1)+.98*(xlims(2)-xlims(1)),elevlims(1)+.98*(elevlims(2)-elevlims(1)),...
    {['Flowband: ',inputfiles{runnum}(underscorenum(1)-1)];['Sliding: m=',num2str(mvalue)];['Calving: $',calvingstring,'$']},...
    'HorizontalAlignment','Right','VerticalAlignment','Top','Color','k','Interpreter','Latex');

% Organize axes:
set(gca,'XLim',xlims)
set(gca,'XTick',xtick*[ceil(xlims(1)/xtick):1:floor(xlims(2)/xtick)])
set(gca,'YLim',elevlims)
set(gca,'YTick',elevtick*[ceil(elevlims(1)/elevtick):1:floor(elevlims(2)/elevtick)])
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')
set(gca,'TickDir','out')
set(gca,'TickLength',[.005,.025])
set(gca,'LineWidth',otherlinewidth)
set(gca,'Box','on')
%set(gca,'Grid','none')
xlabel('Distance (km)','FontSize',fontsize)
ylabel('Elevation (m)','FontSize',fontsize)
set(gca,'FontSize',fontsize)

% Turn off visibility:
set(gcf,'visible','off') 

% Make title:
title([titletext{runnum},', time= ',num2str(floor(ModelRecord.time_yr-reftime)),' years, ',num2str(floor(12*(ModelRecord.time_yr-reftime-floor(ModelRecord.time_yr-reftime)))),' months'],'interpreter',titleinterpreter,'FontSize',titlefontsize)

%% Make IC frame:
drawnow
%frame=getframe(gcf);
addpath('/home/mjw/Documents/MATLAB/export_fig')
frame=export_fig(gcf,'-nocrop', '-a1');
Movie=VideoWriter(outputfile);
Movie.FrameRate=fps;
open(Movie);
writeVideo(Movie,frame);

% Hold IC frame:
for ii=1:icholdtime*fps
    writeVideo(Movie,frame);
end

%% Make Movie:

% Loop through model records:
for recordnum=startrecord+recordsperframe:recordsperframe:endrecord
    
    % Load this model record:
    prefix='0'*ones(1,numdigits(end)-floor(log10(recordnum))-1);
    idnumber=num2str(recordnum);
    load([inputfolder,inputfiles{runnum}],['ModelRecord_',prefix,idnumber])
    
    % Check if model record exists:
    if exist(['ModelRecord_',prefix,idnumber],'var')
        % Rename model record:
        eval(['ModelRecord=ModelRecord_',prefix,idnumber,';'])
        clear(['ModelRecord_',prefix,idnumber])
    else
        % Break from loop:
        break
    end
    
    % Assign sea level:
    if strcmp(rightbctype,'front') && strcmp(rightbcparam,'file')==0
        sealevel=rightbcparam;
    elseif strcmp(rightbctype,'front')
        sealevel=interp1(Time_yr_input,RightBCParam_input,ModelRecord.time_yr,interpstyle,RightBCParam_input(1));
    else
        sealevel=min(BedElev_input)-1;
    end
    
    % Define sill profile:
    if dosill
        if isfield(ModelRecord,'SillThick_input')
            SillThick_input=ModelRecord.SillThick_input;
        else
            % Compute sill thickness profile:
            if sillstarttime_yr==0 && sillconstructiontime_yr==0 && ModelRecord.time_yr==0
                SillThick_input=sillamp*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2));
            else
                SillThick_input=sillamp*max(0,min(1,(ModelRecord.time_yr-sillstarttime_yr)/(sillconstructiontime_yr+1)))*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2)); % note stabilizer
            end
        end
    else
        SillThick_input=zeros(size(X_input));
    end
    
    % Define domain:
    dx=ModelRecord.domainwidth/xsize;
    X_lr=linspace(0,ModelRecord.domainwidth,xsize+1);
    X_c=linspace(.5*dx,ModelRecord.domainwidth-.5*dx,xsize);
    BedElev_lr=interp1(X_input,BedElev_input+SillThick_input,X_lr,interpstyle);
    BedElev_c=.5*(BedElev_lr(1:end-1)+BedElev_lr(2:end));
    
    % Calculate surface elevation:
    SurfElev_c=ModelRecord.IceBottom_c+ModelRecord.Icethick_c;
    SurfElev_lr=[1.5*SurfElev_c(1)-.5*SurfElev_c(2),.5*(SurfElev_c(1:end-1)+SurfElev_c(2:end)),1.5*SurfElev_c(end)-.5*SurfElev_c(end-1)]; % linear extrapolation
    
    % Interpolate ice thickness to grid edges:
    Icethick_lr=[1.5*ModelRecord.Icethick_c(1)-.5*ModelRecord.Icethick_c(2),.5*(ModelRecord.Icethick_c(1:end-1)+ModelRecord.Icethick_c(2:end)),1.5*ModelRecord.Icethick_c(end)-.5*ModelRecord.Icethick_c(end-1)];
    
    % Interpolate ice bottom to grid edges:
    IceBottom_lr=[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(2),.5*(ModelRecord.IceBottom_c(1:end-1)+ModelRecord.IceBottom_c(2:end)),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)];
    
    % Compute surface and bed gradients and angles:
    SurfGrad_lr=gradient(SurfElev_lr,dx);
    BedGrad_c=[.5*(BedElev_c(1)-BedElev_lr(1)),(BedElev_lr(2:end)-BedElev_lr(1:end-1)),.5*(BedElev_lr(end)-BedElev_c(end))]/dx; % also includes domain edges
    SurfAngle_lr=atan(exag*SurfGrad_lr); % angle on exaggerated plot
    BedAngle_c=atan(exag*BedGrad_c); % angle on exaggerated plot
    
    % Interpolate geometry onto the grounding line: (hydrostatic assumption)
    lastgroundedind=find(X_c<=ModelRecord.x_gl,1,'last');
    bedelev_gl=interp1(X_lr,BedElev_lr,ModelRecord.x_gl,interpstyle,BedElev_lr(end));
    if lastgroundedind<xsize
        icethick_gl=(sealevel-bedelev_gl)*(rho_sw/rho_i);
        hasshelf=1;
    else
        icethick_gl=Icethick_lr(end);
        hasshelf=0;
    end
    surfelev_gl=bedelev_gl+icethick_gl;
    
    % Compute column-average velocity:
    Ubar_lr=sum(ModelRecord.U_lr.*repmat(DZhat_c,[1,xsize+1]),1);
    
    % Advect vertical lines:
    X_vertlines=X_vertlines+interp1(X_lr,Ubar_lr*secondsperyear,X_vertlines,interpstyle)*timeperframe;
    
    % Flag vertical lines that have crossed the calving front:
    X_vertlines(X_vertlines>ModelRecord.domainwidth)=NaN;
    
    % Compute spacing between vertical lines:
    DX_vertlines=[X_vertlines(1),X_vertlines(2:end)-X_vertlines(1:end-1)]; % includes spacing between first line and left edge
    
    % Spawn new vertical lines:
    if sum(DX_vertlines>2*dx_vertlines)>0
        % Identify too-large gaps:
        GapTooBig=find(DX_vertlines>2*dx_vertlines);
        % Loop backwards through gaps:
        for jj=length(GapTooBig):-1:1
            % Check to see if this is the first gap:
            if GapTooBig(jj)==1
                % Compute midpoint location:
                newx=.5*X_vertlines(1);
                % Put the new vertical line in order in the list:
                X_vertlines=[newx,X_vertlines(2:end-1)]; % throw out the last one to keep the list the same size
            else
                % Compute midpoint location:
                newx=.5*(X_vertlines(GapTooBig(jj)-1)+X_vertlines(GapTooBig(jj)));
                % Put the new vertical line in order in the list:
                X_vertlines=[X_vertlines(1:GapTooBig(jj)-1),newx,X_vertlines(GapTooBig(jj):end-1)]; % throw out the last one to keep the list the same size
            end
        end
    end
    
    % Interpolate ice bottom and surface to the vertical lines:
    Z_vertlines=[interp1(X_lr,IceBottom_lr,X_vertlines,interpstyle);interp1(X_lr,SurfElev_lr,X_vertlines,interpstyle)];
    
    % Compute all grounding lines (in addition to the most seaward one):
    % Compute hydraulic head:
    HydroHead_c=BedElev_c+(rho_i/rho_sw)*ModelRecord.Icethick_c;
    % Locate sea level crossings:
    sealevelcrossinginds=find(diff(sign(HydroHead_c-sealevel))); % index of last grid cell before the crossing
    % Record whether the last frame had multiple gls:
    hasmultglslast=hasmultiplegls;
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
    
    
    % ICE PATCH:
    % Update ice patch:
    %IcePatch.Vertices=[[X_lr,fliplr(X_lr)]'/1000,[IceBottom_lr,fliplr(SurfElev_lr)]'];
    % ICE PATCH:;
    % Update ice patch:
    IcePatch.Vertices=[[X_lr,X_lr]'/1000,[IceBottom_lr,SurfElev_lr]'];
    % Loop through columns to set ice patch color:
    for ii=1:xsize
        if HydroHead_c(ii)>=sealevel
            IcePatch.FaceVertexCData(ii,:)=icecolor;
        else
            IcePatch.FaceVertexCData(ii,:)=shelfcolor;
        end
    end
    
    
    % SMB PATCH:
    % Compute surface ribbon coordinates:
    SurfRibbonTopZ_lr=SurfElev_lr+ribbonthick*cos(SurfAngle_lr);
    SurfRibbonTopX_lr=X_lr-exag*ribbonthick*sin(SurfAngle_lr);
    % Update surface ribbon patch structure:
    SurfRibbonPatch.Vertices=[[X_lr,SurfRibbonTopX_lr]'/1000,[SurfElev_lr,SurfRibbonTopZ_lr]'];
    % Flip SMB sign and scale to be symmetric about 0:
    ScaledSMB=(ModelRecord.Accum_u-ModelRecord.MeltRate_u)*secondsperyear;
    ScaledSMB(ScaledSMB<0)=max(-1,ScaledSMB(ScaledSMB<0)/(-smblims(1))); % assumes smblims(1) is negative
    ScaledSMB(ScaledSMB>0)=min(1,ScaledSMB(ScaledSMB>0)/smblims(2));
    ScaledSMB=-ScaledSMB;
    % Interpolate surface patch RGB data:
    SurfRibbonPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(SMBcmap,1))',SMBcmap(:,1),(ScaledSMB'+1)/2,interpstyle)));
    SurfRibbonPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(SMBcmap,1))',SMBcmap(:,2),(ScaledSMB'+1)/2,interpstyle)));
    SurfRibbonPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(SMBcmap,1))',SMBcmap(:,3),(ScaledSMB'+1)/2,interpstyle)));
    
    
    % BASAL DRAG PATCH:
    % Compute bed ribbon coordinates:
    BedRibbonBotZ_c=[BedElev_lr(1),BedElev_c,BedElev_lr(end)]-ribbonthick*cos(BedAngle_c); % also includes domain edges
    BedRibbonBotX_c=[0,X_c,ModelRecord.domainwidth]+exag*ribbonthick*sin(BedAngle_c); % also includes domain edges
    % Create new bed ribbon patch structure:
    if hasshelf
        BedRibbonPatch=struct('Vertices',[[0,X_c(1:lastgroundedind+1),BedRibbonBotX_c(1:lastgroundedind+2)]'/1000,...
            [BedElev_lr(1),BedElev_c(1:lastgroundedind+1),BedRibbonBotZ_c(1:lastgroundedind+2)]'],...
            'Faces',[linspace(1,lastgroundedind+1,lastgroundedind+1)',linspace(2,lastgroundedind+2,lastgroundedind+1)',linspace(lastgroundedind+4,2*lastgroundedind+4,lastgroundedind+1)',linspace(lastgroundedind+3,2*lastgroundedind+3,lastgroundedind+1)'],...
            'FaceVertexCData',zeros(lastgroundedind+1,3));
    else
        BedRibbonPatch=struct('Vertices',[[0,X_c(1:lastgroundedind),ModelRecord.domainwidth,BedRibbonBotX_c(1:lastgroundedind+2)]'/1000,...
            [BedElev_lr(1),BedElev_c(1:lastgroundedind),BedElev_lr(end),BedRibbonBotZ_c(1:lastgroundedind+2)]'],...
            'Faces',[linspace(1,lastgroundedind+1,lastgroundedind+1)',linspace(2,lastgroundedind+2,lastgroundedind+1)',linspace(lastgroundedind+4,2*lastgroundedind+4,lastgroundedind+1)',linspace(lastgroundedind+3,2*lastgroundedind+3,lastgroundedind+1)'],...
            'FaceVertexCData',zeros(lastgroundedind+1,3));
    end
    % Interpolate bed patch RGB data:
    BedRibbonPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(BasalDragcmap,1))',BasalDragcmap(:,1),max(0,min(1,(ModelRecord.Drag_lrd(1:lastgroundedind+1)'/1000-basaldraglims(1))/(basaldraglims(2)-basaldraglims(1)))),interpstyle)));
    BedRibbonPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(BasalDragcmap,1))',BasalDragcmap(:,2),max(0,min(1,(ModelRecord.Drag_lrd(1:lastgroundedind+1)'/1000-basaldraglims(1))/(basaldraglims(2)-basaldraglims(1)))),interpstyle)));
    BedRibbonPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(BasalDragcmap,1))',BasalDragcmap(:,3),max(0,min(1,(ModelRecord.Drag_lrd(1:lastgroundedind+1)'/1000-basaldraglims(1))/(basaldraglims(2)-basaldraglims(1)))),interpstyle)));
    
    
%     % OCEAN PATCH:
%     % Check whether bed topography at the grounding line is above sea level:
%     ind1=find(BedElev_input+SillThick_input<sealevel&X_input>ModelRecord.x_gl,1,'first');
%     ind2=find(X_input>ModelRecord.x_gl,1,'first');
%     if ind1==ind2
%         % The ocean starts at the grounding line:
%         oceanstartx=ModelRecord.x_gl;
%     else
%         % The oceans starts where the bed goes below sea level:
%         oceanstartx=interp1(BedElev_input(ind1-1:ind1)+SillThick_input(ind1-1:ind1),X_input(ind1-1:ind1),sealevel,interpstyle);
%     end
%     % Create ocean horizontal grid:
%     X_lr_ocean=linspace(oceanstartx,maxdomainwidth,oceanxsize+1);
%     % Create ocean top and bottom:
%     BedElev_lr_ocean=interp1(X_input,BedElev_input+SillThick_input,X_lr_ocean);
%     % Interpolate ocean properties in time:
%     if oceantimedependence
%         ThisTemperature_input=interp1(Time_yr_input,Temperature_input',ModelRecord.time_yr,interpstyle);
%         if ModelRecord.time_yr>Time_yr_input(end)
%             ThisTemperature_input=Temperature_input(:,end);
%         end
%     else
%         ThisTemperature_input=Temperature_input(:,1);
%     end
%     % Assign ocean grid temperature:
%     for jj=1:oceanxsize+1
%         % Find sill depth:
%         silldepth=sealevel-max(BedElev_lr_ocean(min(jj+1,oceanxsize+1):end));
%         abovesillind=find(Depth_input<=silldepth,1,'last');
%         % Interpolate water properties:
%         Temp_lrud_ocean(:,jj)=interp1(Depth_input(1:abovesillind),ThisTemperature_input(1:abovesillind)-tmelt,sealevel-(BedElev_lr_ocean(jj)+Zhat_ud_ocean*(sealevel-BedElev_lr_ocean(jj))),interpstyle,ThisTemperature_input(abovesillind)-tmelt);
%     end
%     % Update ocean patch structure:
%     OceanPatch.Vertices=[reshape(repmat(X_lr_ocean/1000,[oceanzsize+1,1]),[],1),reshape(repmat(BedElev_lr_ocean,[oceanzsize+1,1])+repmat(Zhat_ud_ocean,[1,oceanxsize+1]).*repmat(sealevel-BedElev_lr_ocean,[oceanzsize+1,1]),[],1),zeros((oceanxsize+1)*(oceanzsize+1),1)];
%     % Interpolate RGB data for ocean patch structure:
%     OceanPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,1),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
%     OceanPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,2),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
%     OceanPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,3),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
    
    
    % OCEAN PATCH:
    % Interpolate ocean properties in time:
    if oceantimedependence
        ThisTemperature_input=interp1(Time_yr_input,Temperature_input',ModelRecord.time_yr,interpstyle);
        if ModelRecord.time_yr>Time_yr_input(end)
            ThisTemperature_input=Temperature_input(:,end);
        end
    else
        ThisTemperature_input=Temperature_input(:,1);
    end
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
    % Update ocean patch vertices:
    OceanPatch.Vertices(:,2)=reshape(repmat(BedElev_input+SillThick_input,[oceanzsize+1,1])+...
        repmat(Zhat_ud_ocean,[1,inputxsize]).*repmat(sealevel-(BedElev_input+SillThick_input),[oceanzsize+1,1]),[],1);
    % Interpolate RGB data for ocean patch structure:
    OceanPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,1),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
    OceanPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,2),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
    OceanPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,3),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),interpstyle)));
    

    
    % SKY PATCH:
    % Interpolate surface elevation onto input grid:
    SurfElev_input=interp1(X_lr,SurfElev_lr,X_input,interpstyle,sealevel);
    % Update sky patch:
    SkyPatch.Vertices=[reshape(repmat(X_input/1000,[skyzsize+1,1]),[],1),reshape(repmat(SurfElev_input,[skyzsize+1,1])+repmat(linspace(0,1,skyzsize+1)',[1,inputxsize]).*repmat(elevlims(2)-SurfElev_input,[skyzsize+1,1]),[],1)];
    SkyPatch.FaceVertexCData=repmat(skycolor,[(skyzsize+1)*inputxsize,1])+repmat([1,1,1]-skycolor,[(skyzsize+1)*inputxsize,1]).*repmat((SkyPatch.Vertices(:,2)-sealevel)/(elevlims(2)-sealevel),[1,3]);
    
    
    % SILL PATCH:
    if dosill
        SillPatch.Vertices(:,2)=[BedElev_input,fliplr(BedElev_input+SillThick_input)]';
    end
    
    
    % Compute monotonic floating ice bottom:
    % Set the monotonic bottom to the real ice bottom:
    MonotonicIceBottom_c=ModelRecord.IceBottom_c;
    % Check if we're running the plume model and have a shelf:
    if hasshelf
        % Detect downward sloping ice bottom in the shelf:
        DownwardSlopingShelf_lr=X_lr>ModelRecord.x_gl&([ModelRecord.IceBottom_c,1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)]-[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(2),ModelRecord.IceBottom_c])<0;
        % Check if there are any cells that need to be fixed:
        if sum(DownwardSlopingShelf_lr)>0
            % Identify cells that need to be fixed:
            needfixinginds=find(DownwardSlopingShelf_lr);
            % Pre-allocate a record of whether they have been fixed:
            hasbeenfixed=false(size(needfixinginds));
            % Loop through from the calving front to the grounding line:
            for thisind=length(needfixinginds):-1:1
                % Check if this grid cell has been fixed:
                if hasbeenfixed(thisind)
                    continue
                end
                % Identify real index of this grid edge:
                d2=needfixinginds(thisind);
                % Identify downstream ice bottom:
                if d2<xsize+1
                    downstreamicebottomx=X_c(d2);
                    downstreamicebottomz=ModelRecord.IceBottom_c(d2);
                else
                    downstreamicebottomx=ModelRecord.domainwidth;
                    downstreamicebottomz=1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1);
                end
                % Identify next ice bottom below that point:
                nextind=find(ModelRecord.IceBottom_c(1:d2-1)<downstreamicebottomz,1,'last');
                if isempty(nextind)
                    nextind=lastgroundedind;
                    nexticebottom=downstreamicebottomz-minbotgrad*(downstreamicebottomx-X_c(nextind));
                else
                    nexticebottom=ModelRecord.IceBottom_c(nextind);
                end
                % Create a linear profile to the next ice bottom:
                if d2<xsize+1
                    MonotonicIceBottom_c(nextind:d2-1)=linspace(nexticebottom,downstreamicebottomz-(downstreamicebottomz-nexticebottom)/(d2-nextind+1),d2-nextind);
                else
                    MonotonicIceBottom_c(nextind:d2-1)=linspace(nexticebottom,downstreamicebottomz-.5*(downstreamicebottomz-nexticebottom)/(d2-nextind+1),d2-nextind);
                end
                % Flag cells that have been fixed:
                hasbeenfixed(needfixinginds>=nextind+1)=1;
            end
        end
    end
    
    % OCEAN MELT PATCH:
    % Compute plume model grid size:
    dx_plume=(ModelRecord.domainwidth-ModelRecord.x_gl)/xsize_plume;
    % Create plume model horizontal grid:
    X_lr_plume(1:xsize_plume+1)=linspace(ModelRecord.x_gl,ModelRecord.domainwidth,xsize_plume+1);
    X_lr_plume(xsize_plume+2:xsize_plume+zsize_plume+1)=ModelRecord.domainwidth;
    X_c_plume=.5*(X_lr_plume(1:end-1)+X_lr_plume(2:end));
    % Compute plume model elevation:
    if hasshelf
        % Interpolate elevation under the floating shelf:
        Z_lr_plume(1:xsize_plume+1)=interp1([ModelRecord.x_gl,X_c(X_c>ModelRecord.x_gl),ModelRecord.domainwidth],...
            [bedelev_gl,MonotonicIceBottom_c(X_c>ModelRecord.x_gl),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],...
            X_lr_plume(1:xsize_plume+1),interpstyle);
    else
        Z_lr_plume(1:xsize_plume+1)=1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1);
    end
    % Compute elevation along the vertical face:
    Z_lr_plume(xsize_plume+1:xsize_plume+zsize_plume+1)=linspace(1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1),sealevel,zsize_plume+1);
    % Interpolate to grid centers:
    Z_c_plume=.5*(Z_lr_plume(1:end-1)+Z_lr_plume(2:end));
    % Compute plume grid gradient and angle:
    PlumeGrad_lr_plume=gradient(Z_lr_plume,dx_plume); % also includes domain edges
    PlumeAngle_lr_plume=atan(exag*PlumeGrad_lr_plume); % angle on exaggerated plot
    % Force vertical face to be vertical:
    PlumeAngle_lr_plume(xsize_plume+2:end)=pi/2;
    % Fix the angle just before the vertical face:
    PlumeAngle_lr_plume(xsize_plume+1)=PlumeAngle_lr_plume(xsize_plume);
    % Compute plume ribbon coordinates: (one extra point to accomodate right
    % angle turn at ice front)
    PlumeRibbonBotZ_lr=[Z_lr_plume-ribbonthick*cos(PlumeAngle_lr_plume),Z_lr_plume(xsize_plume+1)];
    PlumeRibbonBotX_lr=[X_lr_plume+exag*ribbonthick*sin(PlumeAngle_lr_plume),X_lr_plume(xsize_plume+1)+exag*ribbonthick];
    % Reduce under-shelf ribbon to zero size if there is no shelf:
    if hasshelf==0
        PlumeRibbonBotZ_lr(1:xsize_plume+1)=Z_lr_plume(1:xsize_plume+1);
        PlumeRibbonBotX_lr(1:xsize_plume+1)=X_lr_plume(1:xsize_plume+1);
    end
    % Upddate plume patch structure:
    PlumePatch.Vertices=[[X_lr_plume,PlumeRibbonBotX_lr]'/1000,[Z_lr_plume,PlumeRibbonBotZ_lr]'];
    % Create combined melt and calving vector:
    MeltPlusCalving_c_plume=ModelRecord.MeltRate_c_plume;
    MeltPlusCalving_c_plume(xsize_plume+1:end)=MeltPlusCalving_c_plume(xsize_plume+1:end)+ModelRecord.calvingrate_r;
    % Interpolate plume patch RGB data:
    PlumePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,1),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltcalvinglims(1))))-meltcalvinglims(1))/(meltcalvinglims(2)-meltcalvinglims(1)))),interpstyle)));
    PlumePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,2),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltcalvinglims(1))))-meltcalvinglims(1))/(meltcalvinglims(2)-meltcalvinglims(1)))),interpstyle)));
    PlumePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,3),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltcalvinglims(1))))-meltcalvinglims(1))/(meltcalvinglims(2)-meltcalvinglims(1)))),interpstyle)));
    
    
    
    % Call figure:
    %figure(1)
    
    % Delete (some of the) old patches:
    if ribbonthick>0
        delete([hocean;hsky;hicepatch;hsurfribbon;hbedribbon;hplumeribbon;hvertlines])
    else
        delete([hocean;hsky;hicepatch;hvertlines])
    end
    if hasmultglslast
        delete(hmultgls)
    end
    
    % Replot ocean:
    hocean=patch(OceanPatch,'FaceColor','flat','EdgeColor','none');
    
    % Replot sill patch and line:
    if dosill
        delete(hsillpatch)
        hsillpatch=patch(SillPatch,'FaceColor','flat','EdgeColor','none');
        set(hsillline,'YData',BedElev_input+SillThick_input)
    end
    
    % Replot sky:
    hsky=patch(SkyPatch,'FaceColor','interp','EdgeColor','none');
    
    % Replot ice patch:
    hicepatch=patch(IcePatch,'FaceColor','flat','EdgeColor','none');
    
    % Replot color ribbons:
    if ribbonthick>0
        % Replot surface mass balance:
        hsurfribbon=patch(SurfRibbonPatch,'FaceColor','flat','EdgeColor','none');
        % Replot basal drag:
        hbedribbon=patch(BedRibbonPatch,'FaceColor','flat','EdgeColor','none');
        % Replot plume:
        hplumeribbon=patch(PlumePatch,'FaceColor','flat','EdgeColor','none');
        set(hplumeedge(1),'XData',PlumeRibbonBotX_lr(1:xsize_plume+1)/1000)
        set(hplumeedge(1),'YData',PlumeRibbonBotZ_lr(1:xsize_plume+1))
        set(hplumeedge(2),'XData',PlumeRibbonBotX_lr(xsize_plume+2:end)/1000)
        set(hplumeedge(2),'YData',PlumeRibbonBotZ_lr(xsize_plume+2:end))
    end
    
    % Update ice surface:
    set(hsurf,'XData',[0,X_c,ModelRecord.domainwidth]/1000)
    set(hsurf,'YData',[SurfElev_lr(1),SurfElev_c,SurfElev_lr(end)])
    
    % Update ice bottom:
    set(hbot,'XData',[0,X_c,ModelRecord.domainwidth]/1000)
    set(hbot,'YData',[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(1),ModelRecord.IceBottom_c,1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)])
    
    % Update ice front:
    set(hfront,'XData',ModelRecord.domainwidth*[1,1]/1000)
    set(hfront,'YData',[SurfElev_lr(end),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)])
    
    % Update grounding line:
    set(hgl,'XData',ModelRecord.x_gl*[1,1]/1000)
    set(hgl,'YData',[bedelev_gl,surfelev_gl])
    
    % Plot multiple grounding lines:
    if hasmultiplegls
        hmultgls=zeros(numgls,1);
        for thisgl=1:numgls
            hmultgls(thisgl)=plot(X_gls(thisgl)*[1,1]/1000,[BedElev_gls(thisgl),SurfElev_gls(thisgl)],'Color','k','LineWidth',gllinewidth);
        end
    end
    
    % Replot vertical lines:
    hvertlines=plot(repmat(X_vertlines/1000,[2,1]),Z_vertlines,'Color','k','LineWidth',vertlinewidth);
    
    % Rearrange plot order:
    if dosill
        uistack([htext;hsurf;hbot;hfront;hgl;hsillline;hinput],'top')
    else
        uistack([htext;hsurf;hbot;hfront;hgl;hinput],'top')
    end
    
    % Replot title:
    title([titletext{runnum},', time= ',num2str(floor(ModelRecord.time_yr-reftime)),' years, ',num2str(floor(12*(ModelRecord.time_yr-reftime-floor(ModelRecord.time_yr-reftime)))),' months'],'interpreter',titleinterpreter,'FontSize',titlefontsize)
    
    % Make frame:
    drawnow
    %frame=getframe(gcf);
    frame=export_fig(gcf,'-nocrop', '-a1');
    writeVideo(Movie,frame);
    
end

% Hold final frame:
for recordnum=1:finalholdtime*fps
    writeVideo(Movie,frame);
end

% Finish up:
close(Movie)
close all

end

disp('Done!')
toc


%% How to assemble the patch structure necessary to display model variables:

% In all cases, insert the variable you want to plot in "FaceVertexCData"
% and plot using:
% patch(Patch,'EdgeColor','color','FaceColor','flat')

% % Variable at grid centers:
% % Create patch face indexing:
% Faces=zeros(xsize*zsize,4);
% Faces(1:zsize,:)=[linspace(1,zsize,zsize)',linspace(2,zsize+1,zsize)',linspace(zsize+3,2*zsize+2,zsize)',linspace(zsize+2,2*zsize+1,zsize)'];
% for recordnum=2:xsize
%     Faces((recordnum-1)*zsize+1:recordnum*zsize,:)=repmat(Faces((recordnum-1)*zsize,:)+1,[zsize,1])+repmat(linspace(1,zsize,zsize)',[1,4]);
% end
% % Create patch structure:
% Patch=struct('Vertices',[reshape(repmat(X_lr/1000,[zsize+1,1]),[],1),reshape(repmat(BedElev_lr+WaterThick_lr,[zsize+1,1])+repmat(Zhat_ud,[1,xsize+1]).*repmat(SurfElev_lr-(BedElev_lr+WaterThick_lr),[zsize+1,1]),[],1),zeros((xsize+1)*(zsize+1),1)],'Faces',Faces,'FaceVertexCData',reshape(Variable_c,[],1));

% % Variable at grid left/right edges:
% % Create patch face indexing:
% Faces=zeros((xsize+1)*zsize,4);
% Faces(1:zsize,:)=[linspace(1,zsize,zsize)',linspace(2,zsize+1,zsize)',linspace(zsize+3,2*zsize+2,zsize)',linspace(zsize+2,2*zsize+1,zsize)'];
% for recordnum=2:xsize+1
%     Faces((recordnum-1)*zsize+1:recordnum*zsize,:)=repmat(Faces((recordnum-1)*zsize,:)+1,[zsize,1])+repmat(linspace(1,zsize,zsize)',[1,4]);
% end
% % Create patch structure:
% Patch=struct('Vertices',[reshape(repmat([0,X_c,domainwidth]/1000,[zsize+1,1]),[],1),reshape([bedelev_l+Zhat_ud*icethick_l,GridElev_ud,bedelev_r+Zhat_ud*icethick_r],[],1),zeros((xsize+2)*(zsize+1),1)],'Faces',Faces,'FaceVertexCData',reshape(Variable_lr,[],1));

% % Variable at grid up/down edges:
% Create patch face indexing:
% Faces=zeros(xsize*(zsize+1),4);
% Faces(1:zsize+1,:)=[linspace(1,zsize+1,zsize+1)',linspace(2,zsize+2,zsize+1)',linspace(zsize+4,2*(zsize+1)+2,zsize+1)',linspace(zsize+3,2*(zsize+1)+1,zsize+1)'];
% for recordnum=2:xsize
%     Faces((recordnum-1)*(zsize+1)+1:recordnum*(zsize+1),:)=repmat(Faces((recordnum-1)*(zsize+1),:)+1,[zsize+1,1])+repmat(linspace(1,zsize+1,zsize+1)',[1,4]);
% end
% % Create patch structure:
% Patch=struct('Vertices',[reshape(repmat(X_lr/1000,[zsize+2,1]),[],1),reshape([BedElev_lr+WaterThick_lr;repmat(BedElev_lr+WaterThick_lr,[zsize,1])+repmat(Zhat_c,[1,xsize+1]).*repmat(Icethick_lr,[zsize,1]);SurfElev_lr],[],1),zeros((xsize+1)*(zsize+2),1)],'Faces',Faces,'FaceVertexCData',reshape(Variable_ud,[],1));

% % Variable at grid corners:
% % Create patch face indexing:
% Faces=zeros((xsize+1)*(zsize+1),4);
% Faces(1:zsize+1,:)=[linspace(1,zsize+1,zsize+1)',linspace(2,zsize+2,zsize+1)',linspace(zsize+4,2*(zsize+1)+2,zsize+1)',linspace(zsize+3,2*(zsize+1)+1,zsize+1)'];
% for recordnum=2:xsize+1
%     Faces((recordnum-1)*(zsize+1)+1:recordnum*(zsize+1),:)=repmat(Faces((recordnum-1)*(zsize+1),:)+1,[zsize+1,1])+repmat(linspace(1,zsize+1,zsize+1)',[1,4]);
% end
% % Create patch structure:
% Patch=struct('Vertices',[reshape(repmat([0,X_c,domainwidth]/1000,[zsize+2,1]),[],1),reshape([bedelev_l,BedElev_c+WaterThick_c,bedelev_r;bedelev_l+Zhat_c*icethick_l,repmat(BedElev_c+WaterThick_c,[zsize,1])+repmat(Zhat_c,[1,xsize]).*repmat(Icethick_c,[zsize,1]),bedelev_r+Zhat_c*icethick_r;surfelev_l,SurfElev_c,surfelev_r],[],1),zeros((xsize+2)*(zsize+2),1)],'Faces',Faces,'FaceVertexCData',reshape(Variable_lrud,[],1));

