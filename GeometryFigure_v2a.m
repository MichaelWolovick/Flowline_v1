% GeometryFigure_v1

% Mike Wolovick, 7/1/2016

% This script makes a figure showing bed and surface geometry for all six
% glaciers.

% v2:  This script only makes a figure for Thwaites, but it shows a
% map view of bed geometry with flowband boundaries overlain in addition to
% profiles of the geometry.  There may also be an inset showing all of
% Antarctica

% v2a:  no longer shows narrow boundaries, but puts distance contours on
% both maps.

clear all
close all
tic

%% Parameters:

% File names and paths:
% gridfile='/net/mjw/CombinedGrids/AntarcticaCombinedGrids_v1.mat';
% inputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/';
gridfile='/home/wolovick/Dropbox/CombinedGrids/AntarcticaCombinedGrids_v1.mat';
inputfolder='/home/wolovick/Dropbox/FjordSillPaper/ModelInput/';
inputfiles={'ThwaitesA_v2.mat';...
    'ThwaitesB_v2.mat';...
    'ThwaitesC_v2.mat'};
% colormapfile='/home/mjw/Documents/GMA_combined.mat';
% figname='/home/mjw/Documents/FjordSillPaper/Figures/ThwaitesSettingFigure_v4a.png';
% insetfigname='/home/mjw/Documents/FjordSillPaper/Figures/ThwaitesSettingFigure_v4a_inset.png';
colormapfile='/home/wolovick/Dropbox/GMA_combined.mat';
figname='/home/wolovick/Dropbox/FjordSillPaper/Figures/ThwaitesSettingFigure_v4c.png';
insetfigname='/home/wolovick/Dropbox/FjordSillPaper/Figures/ThwaitesSettingFigure_v4c_inset.png';

% Densities:
rho_i=917;                       % kg/m^3
rho_sw=1028;                     % kg/m^3

% Ticks and limits:
xlims=[-1700,-1050];             % [1x2] km
ytick=200;                       % km
ylims=[-750,-100];               % [1x2] km
xtick=200;                       % km
zlims=[-1500,1999.99];           % [1x2] m
ztick=500;                       % m
zkm=1;                           % logical (doesn't affect values of zlims or zticks, only labels)
vlims=[1,1e3];                   % [1x2] m/yr
vcmap='jet';                     % valid colormap string
bedlims=[-2000,2000];            % [1x2] m
bedcmap='file';                  % valid colormap string or 'file'
distcontour=50;                  % km (0 deactivates)

% Surface smoothing wavelength:  (for balance velocity and direction)
surfwavelength=2e4;              % m

% Minimum velocity error: (reported errors do not include baseline errors)
minvelerror=3;                   % m/yr

% Velocity error factor: (for weighting function)
velerrorfactor=5;                % unitless (>1 downweights obs)

% Number of distance smoothing iterations:
numsmoothingiterations=5;        % integer (0 deactivates)

% Vertical Exageration: (for profiles)
exag=150;                        % unitless

% Aesthetic settings:
widelincolor='b';                % valid color string or [R,G,B]
narrowcolor='g';                 % valid color string or [R,G,B]
widefluxcolor='r';               % valid color string or [R,G,B]
linewidth=.75;                   % points
fontsize=14;                     % points
ticklength=.02;                  % unitless

% Lighting settings: (hillshaded bedelev)
lightangle=75;                   % degrees above horizontal
diffusefraction=.33;             % unitless 
totallight=1;                    % unitless

% Margins:
verttextbuffer=.08;              % unitless
titlebuffer=.06;                 % unitless
horztextbuffer=.07;              % unitless
horznotextbuffer=.075;           % unitless

% Output size and resolution:
pagesize=[11,5];                 % [1x2] inches
resolution=300;                  % dpi
insetpagesize=[2,2];             % [1x2] inches

%% Preparation:

% Communicate:
disp('Loading inputs, doing initial preparation...')

% Hardcoded plot arrangement:
numvertplots=1;
numhorzplots=3;

% Compute subplot positions:
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

% Pre-allocate boundary picks structure:
Boundaries=struct('X_boundary',cell(1,2),'Y_boundary',cell(1,2),...
    'X_gate',cell(1,2),'Y_gate',cell(1,2),'pathlength',cell(1,2));

% Load wide boundaries:
load([inputfolder,inputfiles{1}],'X_boundary','Y_boundary')
Boundaries(1).X_boundary=X_boundary;
Boundaries(1).Y_boundary=Y_boundary;

% Load narrow boundaries:
load([inputfolder,inputfiles{2}],'X_boundary','Y_boundary')
Boundaries(2).X_boundary=X_boundary;
Boundaries(2).Y_boundary=Y_boundary;

% Load input grids:
load(gridfile)

% Make inset map:
figure(1)
contour(X/1000,Y/1000,Mask,.5*[1,1],'k')
hold on
contour(X/1000,Y/1000,isnan(Mask),.5*[1,1],'k')
plot([xlims,fliplr(xlims),xlims(1)],[ylims(1),ylims(1),ylims(2),ylims(2),ylims(1)],'k')
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'PlotBoxAspectRatio',[(X(2)-X(1))/(Y(2)-Y(1)),1,1])
set(gcf,'PaperSize',insetpagesize)
set(gcf,'PaperPosition',[0,0,insetpagesize])
print('-dpng',insetfigname,['-r',num2str(resolution)])

% Find truncating indices:
ind1=find(X>=xlims(1)*1000,1,'first'); % assuming that xlims are in km and X in m
ind2=find(X<=xlims(2)*1000,1,'last');
ind3=find(Y>=ylims(1)*1000,1,'first');
ind4=find(Y<=ylims(2)*1000,1,'last');

% Truncate input grids to region of interest:
X=X(ind1:ind2);
Y=Y(ind3:ind4);
BedElev=BedElev(ind3:ind4,ind1:ind2);
SurfElev=SurfElev(ind3:ind4,ind1:ind2);
IceThick=IceThick(ind3:ind4,ind1:ind2);
U=U(ind3:ind4,ind1:ind2);
V=V(ind3:ind4,ind1:ind2);
VelError=VelError(ind3:ind4,ind1:ind2);
SurfTemp=SurfTemp(ind3:ind4,ind1:ind2);
Accum=Accum(ind3:ind4,ind1:ind2);
Gflux=Gflux(ind3:ind4,ind1:ind2);
Mask=Mask(ind3:ind4,ind1:ind2);
if exist('U_bal','var') && exist('V_bal','var')
    U_bal=U_bal(ind3:ind4,ind1:ind2);
    V_bal=V_bal(ind3:ind4,ind1:ind2);
end

% Define grid size:
dx=X(2)-X(1); % assume dx=dy
[ysize,xsize]=size(SurfElev);

% Create gridded x/y vectors:
X_grid=repmat(X,[ysize,1]);
Y_grid=repmat(Y,[1,xsize]);

% Replace NaN's in surface and ice thickness data:
SurfElev(isnan(SurfElev))=0;
IceThick(isnan(IceThick))=0;

% Produce smooth surface data:
SurfElev_smooth=intuitive_lowpass2(SurfElev,surfwavelength/dx,surfwavelength/dx);

% Define downhill surface gradient:
[SurfGradX,SurfGradY]=gradient(-SurfElev_smooth,dx,dx);
SurfGradMag=sqrt(SurfGradX.^2+SurfGradY.^2);

% Compute surface gradient direction: (real=x, imag=y)
SurfGradDir=SurfGradX./SurfGradMag+1i*SurfGradY./SurfGradMag;

% Compute velocity direction:
VelMag=sqrt(U.^2+V.^2);
VelDir=U./VelMag+1i*V./VelMag;

% Adjust velocity error for baseline errors:
VelError=sqrt(VelError.^2+minvelerror^2);

%% Compute balance velocity:

% Check if input file contained balance velocities:
if exist('U_bal','var')==0 || exist('V_bal','var')==0
    
    % Communicate:
    disp('Computing balance velocities...')
    
    % Pre-allocate velocity on grid edges:
    U_bal_ew=zeros(ysize,xsize+1);
    V_bal_ns=zeros(ysize+1,xsize);
    
    % Compute coordinates on grid edges:
    X_ew=[X(1)-.5*dx,.5*(X(1:end-1)+X(2:end)),X(end)+.5*dx];
    Y_ns=[Y(1)-.5*dx;.5*(Y(1:end-1)+Y(2:end));Y(end)+.5*dx];
    
    % Interpolate ice thickness to grid edges:
    IceThick_ew=interp2(X,Y,IceThick,repmat(X_ew,[ysize,1]),repmat(Y,[1,xsize+1]));
    IceThick_ns=interp2(X,Y,IceThick,repmat(X,[ysize+1,1]),repmat(Y_ns,[1,xsize]));
    
    % Sort by surface elevation:
    [~,SortedInds]=sort(SurfElev_smooth(:),'descend');
    
    % Loop through grid cells:
    for ii=1:ysize*xsize
        
        % Get subscript indices:
        [d1,d2]=ind2sub([ysize,xsize],SortedInds(ii));
        
        % Check whether there is ice here and that we are not on the edge:
        if isnan(Mask(d1,d2)) || d1==1 || d1==ysize || d2==1 || d2==xsize
            continue
        end
        
        % Compute influx:
        influx=U_bal_ew(d1,d2)*IceThick_ew(d1,d2)*dx-U_bal_ew(d1,d2+1)*IceThick_ew(d1,d2+1)*dx...
            +V_bal_ns(d1,d2)*IceThick_ns(d1,d2)*dx-V_bal_ns(d1+1,d2)*IceThick_ns(d1+1,d2)*dx;
        
        % Compute outflux:
        outflux=influx++Accum(d1,d2)*(dx^2);
        
        % Compute total down-slope:
        totaldownslope=max(0,(IceThick_ns(d1,d2)~=0)*(SurfElev_smooth(d1,d2)-SurfElev_smooth(d1-1,d2)))...
            +max(0,(IceThick_ns(d1+1,d2)~=0)*(SurfElev_smooth(d1,d2)-SurfElev_smooth(d1+1,d2)))...
            +max(0,(IceThick_ew(d1,d2)~=0)*(SurfElev_smooth(d1,d2)-SurfElev_smooth(d1,d2-1)))...
            +max(0,(IceThick_ew(d1,d2+1)~=0)*(SurfElev_smooth(d1,d2)-SurfElev_smooth(d1,d2+1)));
        
        % Partition outflux:
        if totaldownslope>0
            if SurfElev_smooth(d1-1,d2)<SurfElev_smooth(d1,d2) && IceThick_ns(d1,d2)~=0
                V_bal_ns(d1,d2)=-((SurfElev_smooth(d1,d2)-SurfElev_smooth(d1-1,d2))/totaldownslope)*outflux/(IceThick_ns(d1,d2)*dx);
            end
            if SurfElev_smooth(d1+1,d2)<SurfElev_smooth(d1,d2) && IceThick_ns(d1+1,d2)~=0
                V_bal_ns(d1+1,d2)=((SurfElev_smooth(d1,d2)-SurfElev_smooth(d1+1,d2))/totaldownslope)*outflux/(IceThick_ns(d1+1,d2)*dx);
            end
            if SurfElev_smooth(d1,d2-1)<SurfElev_smooth(d1,d2) && IceThick_ew(d1,d2)~=0
                U_bal_ew(d1,d2)=-((SurfElev_smooth(d1,d2)-SurfElev_smooth(d1,d2-1))/totaldownslope)*outflux/(IceThick_ew(d1,d2)*dx);
            end
            if SurfElev_smooth(d1,d2+1)<SurfElev_smooth(d1,d2) && IceThick_ew(d1,d2+1)~=0
                U_bal_ew(d1,d2+1)=((SurfElev_smooth(d1,d2)-SurfElev_smooth(d1,d2+1))/totaldownslope)*outflux/(IceThick_ew(d1,d2+1)*dx);
            end
        end
        
    end
    
    % Interpolate balance velocities onto grid centers and compute magnitude:
    U_bal=.5*(U_bal_ew(:,1:end-1)+U_bal_ew(:,2:end));
    V_bal=.5*(V_bal_ns(1:end-1,:)+V_bal_ns(2:end,:));
    VelMag_bal=sqrt(U_bal.^2+V_bal.^2);
    
    % Save balance velocity to the input file:
    %save(inputfile,'U_bal','V_bal','surfwavelength','-append')
    
else
    % Compute balance velocity magnitude:
    VelMag_bal=sqrt(U_bal.^2+V_bal.^2);
end

%% Merge Velocity and Direction Fields:

% Communicate:
disp('Merging velocity data with balance velocity...')

% Compute weighting between velocity observations and balance velocity:
VelWeight=1-exp(-VelMag./(velerrorfactor*VelError));
VelWeight(isnan(VelWeight))=0;

% Assign direction weighting:
DirWeight=VelWeight;

% Adjust velocity weighting for areas where the observations indicate much
% slower flow than balance model:
VelWeight(VelMag<.1*VelMag_bal&VelMag_bal>velerrorfactor*VelError&VelMag~=0)=1;

% Compute combined velocity grid:
VelMagMerged=VelWeight.*VelMag+(1-VelWeight).*VelMag_bal;
VelMagMerged(isnan(VelMag))=VelMag_bal(isnan(VelMag));

% Compute combined direction field:
DirMerged=DirWeight.*VelDir+(1-DirWeight).*SurfGradDir;
DirMerged=DirMerged./abs(DirMerged);
DirMerged(isnan(DirMerged))=SurfGradDir(isnan(DirMerged));

%% Compute Along-Flow Distance:

% This only computes distance for the wide flowband boundaries.

% Check if we're going to be plotting distance contours:
if distcontour~=0
    
    % Communicate:
    disp('Computing along-flow distance...')
    
    % Load all boundary data for the wide flowband:
    load([inputfolder,inputfiles{1}],'*boundary*','*_gate','pathlength')
    
    % Define combined flowband boundary: (including upstream BC zone)
    X_bigboundary=[flipud(X_boundary4);X_boundary1;X_boundary3;flipud(X_boundary2);X_boundary5];
    Y_bigboundary=[flipud(Y_boundary4);Y_boundary1;Y_boundary3;flipud(Y_boundary2);Y_boundary5];
    
    % Identify grid cells within the flowband:
    InFlowBand=inpolygon(X_grid,Y_grid,X_bigboundary,Y_bigboundary);
    numtotalgridpoints=sum(sum(InFlowBand));
    
    % Identify grid cells within the extension:
    InExtension=inpolygon(X_grid,Y_grid,[X_boundary4;flipud(X_boundary5)],[Y_boundary4;flipud(Y_boundary5)]);
    numextensiongridpoints=sum(sum(InExtension));
    numnonextensiongridpoints=numtotalgridpoints-numextensiongridpoints;
    
    % Compute direction along the extension boundaries: (points downstream)
    % Compute gradient of boundaries:
    GradX4=[X_boundary4(1)-X_boundary1(1);X_boundary4(2:end)-X_boundary4(1:end-1)];
    GradY4=[Y_boundary4(1)-Y_boundary1(1);Y_boundary4(2:end)-Y_boundary4(1:end-1)];
    GradX5=[X_boundary5(1)-X_boundary2(1);X_boundary5(2:end)-X_boundary5(1:end-1)];
    GradY5=[Y_boundary5(1)-Y_boundary2(1);Y_boundary5(2:end)-Y_boundary5(1:end-1)];
    % Compute gradient magnitude of boundaries:
    GradMag4=sqrt(GradX4.^2+GradY4.^2);
    GradMag5=sqrt(GradX5.^2+GradY5.^2);
    % Compute direction of boundaries:
    Dir4=GradX4./GradMag4+1i*GradY4./GradMag4;
    Dir5=GradX5./GradMag5+1i*GradY5./GradMag5;
    
    % Interpolate direction from the boundaries into the extension:
    Interpolant=scatteredInterpolant([X_boundary4;flipud(X_boundary5)],...
        [Y_boundary4;flipud(Y_boundary5)],...
        [Dir4;flipud(Dir5)]);
    DirMerged(InExtension)=Interpolant(X_grid(InExtension),Y_grid(InExtension));
    
    % Create interpolated points along the flux gate:
    gatelength=sqrt((x_gate(2)-x_gate(1)).^2+(y_gate(2)-y_gate(1)).^2);
    numinterppoints=ceil(3*gatelength/dx);
    X_gateinterp=linspace(x_gate(1),x_gate(2),numinterppoints)';
    Y_gateinterp=linspace(y_gate(1),y_gate(2),numinterppoints)';
    
    % Compute straight-line distance from the middle of the flux gate:
    StraightLineDistance=sqrt((X_grid-mean(x_gate)).^2+(Y_grid-mean(y_gate)).^2);
    
    % Identify cells close to the flux gate:
    CloseCells=InFlowBand&StraightLineDistance<.75*gatelength;
    CloseInds=find(CloseCells);
    
    % Pre-allocate:
    Distance=zeros(ysize,xsize);
    ZeroInput=false(ysize,xsize);
    DoubleZeroInput=false(ysize,xsize);
    
    % Assign distance for cells on the flux gate:
    % Loop through close cells:
    for ii=1:length(CloseInds)
        % Get subscript indices:
        [d1,d2]=ind2sub([ysize,xsize],CloseInds(ii));
        % Compute distance from the flux gate:
        thisdist=min(sqrt((X(d2)-X_gateinterp).^2+(Y(d1)-Y_gateinterp).^2));
        % Check whether we are within one grid cell of the flux gate:
        if thisdist<dx
            % Assign distance:
            if InExtension(d1,d2)==0
                Distance(d1,d2)=thisdist;
            else
                Distance(d1,d2)=-thisdist;
            end
        end
    end
    
    % Iterate the distance integration:
    done=0;
    iteration=1;
    ActiveCells_last=false(ysize,xsize);
    while done==0
        
        % Identify cells we are working on in this iteration:
        ActiveCells=[false(ysize,1),[false(1,xsize-2);...
            Distance(2:end-1,2:end-1)==0&InFlowBand(2:end-1,2:end-1)&...
            (Distance(1:end-2,2:end-1)~=0|Distance(2:end-1,1:end-2)~=0|Distance(2:end-1,3:end)~=0|Distance(3:end,2:end-1)~=0);...
            false(1,xsize-2)],false(ysize,1)];
        ActiveCells=ActiveCells|ZeroInput;
        numactivecells=sum(sum(ActiveCells));
        ActiveInds=find(ActiveCells);
        
        % Loop through active cells:
        for ii=1:numactivecells
            
            % Get subscript indices:
            [d1,d2]=ind2sub([ysize,xsize],ActiveInds(ii));
            
            % Check whether we are in the extension:
            if InExtension(d1,d2)==0
                % Integrate upstream:
                % Check flow direction to assign downstream distance values:
                if real(DirMerged(d1,d2))>0
                    downstreamdistx=Distance(d1,d2+1);
                else
                    downstreamdistx=Distance(d1,d2-1);
                end
                if imag(DirMerged(d1,d2))>0
                    downstreamdisty=Distance(d1+1,d2);
                else
                    downstreamdisty=Distance(d1-1,d2);
                end
                % Check for zero inputs:
                ZeroInput(d1,d2)=0;
                if downstreamdistx==0 || downstreamdisty==0
                    ZeroInput(d1,d2)=1;
                end
                DoubleZeroInput(d1,d2)=0;
                if downstreamdistx==0 && downstreamdisty==0
                    DoubleZeroInput(d1,d2)=1;
                end
                % Integrate distance:
                if downstreamdistx==0 && downstreamdisty~=0
                    Distance(d1,d2)=downstreamdisty+dx;
                elseif downstreamdistx~=0 && downstreamdisty==0
                    Distance(d1,d2)=downstreamdistx+dx;
                elseif downstreamdistx~=0 && downstreamdisty~=0
                    Distance(d1,d2)=(dx+abs(real(DirMerged(d1,d2)))*downstreamdistx+abs(imag(DirMerged(d1,d2)))*downstreamdisty)/(abs(real(DirMerged(d1,d2)))+abs(imag(DirMerged(d1,d2))));
                end
            else
                % Integrate downstream:
                % Check flow direction to assign upstream distance values:
                if real(DirMerged(d1,d2))>0
                    upstreamdistx=Distance(d1,d2-1);
                else
                    upstreamdistx=Distance(d1,d2+1);
                end
                if imag(DirMerged(d1,d2))>0
                    upstreamdisty=Distance(d1-1,d2);
                else
                    upstreamdisty=Distance(d1+1,d2);
                end
                % Check for zero inputs:
                ZeroInput(d1,d2)=0;
                if upstreamdistx==0 || upstreamdisty==0
                    ZeroInput(d1,d2)=1;
                end
                DoubleZeroInput(d1,d2)=0;
                if upstreamdistx==0 && upstreamdisty==0
                    DoubleZeroInput(d1,d2)=1;
                end
                % Integrate distance:
                if upstreamdistx==0 && upstreamdisty~=0
                    Distance(d1,d2)=upstreamdisty-dx;
                elseif upstreamdistx~=0 && upstreamdisty==0
                    Distance(d1,d2)=upstreamdistx-dx;
                elseif upstreamdistx~=0 && upstreamdisty~=0
                    Distance(d1,d2)=(-dx+abs(real(DirMerged(d1,d2)))*upstreamdistx+abs(imag(DirMerged(d1,d2)))*upstreamdisty)/(abs(real(DirMerged(d1,d2)))+abs(imag(DirMerged(d1,d2))));
                end
            end
            
        end
        
        % Check whether to end the integration:
        if isequal(ActiveCells,ActiveCells_last)
            done=1;
        else
            ActiveCells_last=ActiveCells;
            iteration=iteration+1;
        end
        
    end
    
    % Iteratively fill in remaining zero distance values:
    % Identify zeros inside the flowband:
    ZeroCells=InFlowBand&Distance==0;
    numzerocells=sum(sum(ZeroCells));
    ZeroCellInds=find(ZeroCells);
    while numzerocells~=0
        % Loop through zero cells:
        for ii=1:numzerocells
            % Get subscript indices:
            [d1,d2]=ind2sub([ysize,xsize],ZeroCellInds(ii));
            % Isolate local 3x3 patch:
            TheseDistance=Distance(d1-1:d1+1,d2-1:d2+1);
            TheseUsefulCells=InFlowBand(d1-1:d1+1,d2-1:d2+1)&ZeroCells(d1-1:d1+1,d2-1:d2+1)==0;
            % Average non-zero neighbors that are inside the flowband:
            if sum(sum(TheseUsefulCells))~=0
                Distance(d1,d2)=mean(TheseDistance(TheseUsefulCells));
            end
        end
        % Identify zeros inside the flowband:
        ZeroCells=InFlowBand&Distance==0;
        numzerocells=sum(sum(ZeroCells));
        ZeroCellInds=find(ZeroCells);
    end
    
    % Identify cells just outside the border of the flowband:
    JustOutsideFlowBand=[false(ysize,1),[false(1,xsize-2);...
        InFlowBand(2:end-1,2:end-1)==0&(InFlowBand(1:end-2,1:end-2)|InFlowBand(1:end-2,2:end-1)|InFlowBand(1:end-2,3:end)|...
        InFlowBand(2:end-1,1:end-2)|InFlowBand(2:end-1,3:end)|...
        InFlowBand(3:end,1:end-2)|InFlowBand(3:end,2:end-1)|InFlowBand(3:end,3:end));...
        false(1,xsize-2)],false(ysize,1)];
    numjustoutsideflowband=sum(sum(JustOutsideFlowBand));
    JustOutsideInds=find(JustOutsideFlowBand);
    
    % Smooth gridded distance:
    if numsmoothingiterations>0
        % Iteratively apply smoothing:
        for iteration=1:numsmoothingiterations
            % Replace zeros around the border of the flowband:
            for ii=1:numjustoutsideflowband
                % Get subscript indices:
                [d1,d2]=ind2sub([ysize,xsize],JustOutsideInds(ii));
                % Isolate local 3x3 patch:
                TheseDistance=Distance(d1-1:d1+1,d2-1:d2+1);
                TheseUsefulCells=InFlowBand(d1-1:d1+1,d2-1:d2+1);
                % Average neighbors that are inside the flowband:
                Distance(d1,d2)=mean(TheseDistance(TheseUsefulCells));
            end
            % Smooth distance with a 3x3 filter:
            Distance=[zeros(ysize,1),[zeros(1,xsize-2);...
                (1/9)*(Distance(1:end-2,1:end-2)+Distance(1:end-2,2:end-1)+Distance(1:end-2,3:end)+...
                Distance(2:end-1,1:end-2)+Distance(2:end-1,2:end-1)+Distance(2:end-1,3:end)+...
                Distance(3:end,1:end-2)+Distance(3:end,2:end-1)+Distance(3:end,3:end));...
                zeros(1,xsize-2)],zeros(ysize,1)];
        end
    end
    
    % Replace smoothed values around the border of the flowband:
    for ii=1:numjustoutsideflowband
        % Get subscript indices:
        [d1,d2]=ind2sub([ysize,xsize],JustOutsideInds(ii));
        % Isolate local 3x3 patch:
        TheseDistance=Distance(d1-1:d1+1,d2-1:d2+1);
        TheseUsefulCells=InFlowBand(d1-1:d1+1,d2-1:d2+1);
        % Average neighbors that are inside the flowband:
        Distance(d1,d2)=mean(TheseDistance(TheseUsefulCells));
    end
    
    % Replace smoothed values one step further out:
    % Identify cells just outside the just outside:
    JustOutsideJustOutside=[false(ysize,1),[false(1,xsize-2);...
        InFlowBand(2:end-1,2:end-1)==0&JustOutsideFlowBand(2:end-1,2:end-1)==0&(JustOutsideFlowBand(1:end-2,1:end-2)|JustOutsideFlowBand(1:end-2,2:end-1)|JustOutsideFlowBand(1:end-2,3:end)|...
        JustOutsideFlowBand(2:end-1,1:end-2)|JustOutsideFlowBand(2:end-1,3:end)|...
        JustOutsideFlowBand(3:end,1:end-2)|JustOutsideFlowBand(3:end,2:end-1)|JustOutsideFlowBand(3:end,3:end));...
        false(1,xsize-2)],false(ysize,1)];
    numjustoutsidejustoutside=sum(sum(JustOutsideJustOutside));
    JustOutsideInds=find(JustOutsideJustOutside);
    % Loop through just outside just outside:
    for ii=1:numjustoutsidejustoutside
        % Get subscript indices:
        [d1,d2]=ind2sub([ysize,xsize],JustOutsideInds(ii));
        % Isolate local 3x3 patch:
        TheseDistance=Distance(d1-1:d1+1,d2-1:d2+1);
        TheseUsefulCells=JustOutsideFlowBand(d1-1:d1+1,d2-1:d2+1)|InFlowBand(d1-1:d1+1,d2-1:d2+1);
        % Average neighbors that are inside the flowband:
        Distance(d1,d2)=mean(TheseDistance(TheseUsefulCells));
    end
    
    % Ensure no balance velocity is used downstream:
    % (this is a fix for the problem I had in PIG where the velocity data and
    % the ice geometry data had a different front position, resulting in the
    % balance velocity reappearing once the real velocity data ran out)
    VelWeight(Distance<0)=1;
    VelMagMerged=VelWeight.*VelMag+(1-VelWeight).*VelMag_bal;
    VelMagMerged(isnan(VelMag))=VelMag_bal(isnan(VelMag));
    % Fix the mask (for the same reason):
    Mask(Distance<0&VelMag==0)=NaN;
    
end

%% Create Figure:

% Create a figure:
figure(2)

% Get colormaps:
Vcmap=colormap(vcmap);
if strcmp(bedcmap,'file')
    load(colormapfile,'thiscolormap')
    Bedcmap=thiscolormap;
else
    Bedcmap=colormap(bedcmap);
end

% Make first subplot:
subplot('Position',Boxes{1})

% Interpolate RGB values:
RGB=zeros(ysize,xsize,3);
for ii=1:3
    RGB(:,:,ii)=interp1(linspace(0,1,size(Vcmap,1))',Vcmap(:,ii),max(0,min(1,log(VelMagMerged/vlims(1))/log(vlims(2)/vlims(1)))),'linear');
end
RGB(repmat(isnan(Mask),[1,1,3]))=1;

% Plot velocity:
image(X/1000,Y/1000,RGB)
set(gca,'ydir','normal')
hold on
% Plot grounding line:
contour(X/1000,Y/1000,Mask,.5*[1,1],'LineWidth',1.5*linewidth,'Color','k')
% Plot ice margin:
contour(X/1000,Y/1000,isnan(Mask),.5*[1,1],'LineWidth',1.5*linewidth,'Color','k')
% Plot contours of distance:
if distcontour~=0
    Distance(InFlowBand==0)=NaN;
    contour(X/1000,Y/1000,Distance,[0:distcontour*1000:pathlength],'k','LineWidth',.75*linewidth)
    contour(X/1000,Y/1000,Distance,[0:-distcontour*1000:min(Distance(:))],'k','LineWidth',.75*linewidth)
end
% Plot wide boundaries:
plot([Boundaries(1).X_boundary;Boundaries(1).X_boundary(1)]/1000,...
    [Boundaries(1).Y_boundary;Boundaries(1).Y_boundary(1)]/1000,'k','LineWidth',1.5*linewidth)
%plot([Boundaries(1).X_boundary;Boundaries(1).X_boundary(1)]/1000,...
%    [Boundaries(1).Y_boundary;Boundaries(1).Y_boundary(1)]/1000,'k','LineWidth',.75*linewidth)
% Plot narrow boundaries:
% plot([Boundaries(2).X_boundary;Boundaries(2).X_boundary(1)]/1000,...
%     [Boundaries(2).Y_boundary;Boundaries(2).Y_boundary(1)]/1000,'w','LineWidth',1.5*linewidth)
% plot([Boundaries(2).X_boundary;Boundaries(2).X_boundary(1)]/1000,...
%     [Boundaries(2).Y_boundary;Boundaries(2).Y_boundary(1)]/1000,'k','LineWidth',.75*linewidth)
% Set title:
title('a) Surface Velocity','FontSize',fontsize)
% Fix up axes and aesthetic settings:
xlim(xlims)
ylim(ylims)
set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
xlabel('Relative Easting (km)','FontSize',fontsize)
ylabel('Relative Northing (km)','FontSize',fontsize)
set(gca,'XTick',xlims(1):xtick:xlims(2))
set(gca,'YTick',ylims(1):ytick:ylims(2))
set(gca,'XTickLabel',0:xtick:xlims(2)-xlims(1))
set(gca,'YTickLabel',0:ytick:ylims(2)-ylims(1))
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out')
set(gca,'TickLength',[ticklength,ticklength])
% Draw figure:
drawnow


% Make second subplot:
subplot('Position',Boxes{2})

% Interpolate RGB values:
RGB=zeros(ysize,xsize,3);
for ii=1:3
    RGB(:,:,ii)=interp1(linspace(0,1,size(Bedcmap,1))',Bedcmap(:,ii),max(0,min(1,(BedElev-bedlims(1))/(bedlims(2)-bedlims(1)))),'linear');
end

% Plot bed topography as a 3D surface:
surf(X/1000,Y/1000,BedElev,RGB,'FaceColor','interp','EdgeColor','none',...
    'AmbientStrength',totallight*(1-diffusefraction),'DiffuseStrength',totallight*diffusefraction,'SpecularStrength',0)
view([0,90])
hold on
% Make lighting objects:
light('Position',[0,1,sind(lightangle)])
light('Position',[1,0,sind(lightangle)])
% Plot grounding line:
contour(X/1000,Y/1000,Mask,.5*[1,1],'LineWidth',1.5*linewidth,'Color','k')
% Plot ice margin:
contour(X/1000,Y/1000,isnan(Mask),.5*[1,1],'LineWidth',1.5*linewidth,'Color','k')
% Plot contours of distance:
if distcontour~=0
    Distance(InFlowBand==0)=NaN;
    contour(X/1000,Y/1000,Distance,[0:distcontour*1000:pathlength],'k','LineWidth',.75*linewidth)
    contour(X/1000,Y/1000,Distance,[0:-distcontour*1000:min(Distance(:))],'k','LineWidth',.75*linewidth)
end
% Plot wide boundaries:
plot([Boundaries(1).X_boundary;Boundaries(1).X_boundary(1)]/1000,...
    [Boundaries(1).Y_boundary;Boundaries(1).Y_boundary(1)]/1000,'k','LineWidth',1.5*linewidth)
%plot([Boundaries(1).X_boundary;Boundaries(1).X_boundary(1)]/1000,...
%    [Boundaries(1).Y_boundary;Boundaries(1).Y_boundary(1)]/1000,'k','LineWidth',.75*linewidth)
% Plot narrow boundaries:
%plot([Boundaries(2).X_boundary;Boundaries(2).X_boundary(1)]/1000,...
%    [Boundaries(2).Y_boundary;Boundaries(2).Y_boundary(1)]/1000,'w','LineWidth',1.5*linewidth)
% plot([Boundaries(2).X_boundary;Boundaries(2).X_boundary(1)]/1000,...
%     [Boundaries(2).Y_boundary;Boundaries(2).Y_boundary(1)]/1000,'k','LineWidth',.75*linewidth)
% Set title:
title('b) Bed Topography','FontSize',fontsize)
% Fix up axes and aesthetic settings:
xlim(xlims)
ylim(ylims)
set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
xlabel('Relative Easting (km)','FontSize',fontsize)
set(gca,'XTick',xlims(1):xtick:xlims(2))
set(gca,'YTick',ylims(1):ytick:ylims(2))
set(gca,'XTickLabel',0:xtick:xlims(2)-xlims(1))
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out')
set(gca,'TickLength',[ticklength,ticklength])
% Draw figure:
drawnow


% Make third subplot:
subplot('Position',Boxes{3})

% Pre-allocate:
bedhandles=zeros(3,1);

% Loop through flowbands:
for flowband=1:3
    
    % Load this input file:
    load([inputfolder,inputfiles{flowband}],'X_input','BedElev_input','Icethick_input')
    
    % Determine this line color:
    if flowband==1
        thiscolor=widelincolor;
    elseif flowband==2
        thiscolor=narrowcolor;
    else
        thiscolor=widefluxcolor;
    end
    
    % Compute surface geometry:
    SurfElev_input=max(BedElev_input+Icethick_input,(1-rho_i/rho_sw)*Icethick_input);
    
    % Locate grounding line and calving front:
    lastgroundedind=find(SurfElev_input==BedElev_input+Icethick_input,1,'last');
    x_gl=X_input(lastgroundedind);
    lasticeind=find(isnan(Icethick_input)==0,1,'last');
    x_c=X_input(lasticeind);
    
    % Determine xlims:
    if flowband==1
        thesexlims=[min(-(X_input-x_gl)),max(-(X_input-x_gl))]/1000;
    else
        thesexlims=[min(thesexlims(1),min(-(X_input-x_gl))),max(thesexlims(2),max(-(X_input-x_gl)))]/1000;
    end
    
    % Compute elevations at grounding line and calving front:
    bedelev_gl=interp1(X_input,BedElev_input,x_gl,'linear',BedElev_input(end));
    surfelev_gl=interp1(X_input(1:lasticeind),SurfElev_input(1:lasticeind),x_gl,'linear',SurfElev_input(lasticeind));
    icebottom_c=interp1(X_input(1:lasticeind),SurfElev_input(1:lasticeind)-Icethick_input(1:lasticeind),x_c,'linear',SurfElev_input(lasticeind)-Icethick_input(lasticeind));
    surfelev_c=interp1(X_input(1:lasticeind),SurfElev_input(1:lasticeind),x_c,'linear',SurfElev_input(lasticeind));
    
    % Plot geometry:
    bedhandles(flowband)=plot(-(X_input-x_gl)/1000,BedElev_input,'Color',thiscolor,'LineWidth',linewidth);
    hold on
    plot(-(X_input-x_gl)/1000,SurfElev_input,'Color',thiscolor,'LineWidth',linewidth)
    plot(-(X_input-x_gl)/1000,SurfElev_input-Icethick_input,'Color',thiscolor,'LineWidth',linewidth)
    plot([0,0]/1000,[bedelev_gl,surfelev_gl],'Color',thiscolor,'LineWidth',linewidth)
    plot(-(x_c-x_gl)*[1,1]/1000,[icebottom_c,surfelev_c],'Color',thiscolor,'LineWidth',linewidth)
    
end

% Plot sea level:
plot(thesexlims,[0,0],'--k','LineWidth',linewidth)

% Make a legend:
legend(bedhandles,'A','B','C','location','NorthEast')

% Set up axes, title, etc:
xlim(thesexlims)
set(gca,'XTick',xtick*[ceil(thesexlims(1)/xtick):1:floor(thesexlims(2)/xtick)])
set(gca,'XDir','reverse')
ylim(zlims)
set(gca,'YTick',ztick*[ceil(zlims(1)/ztick):1:floor(zlims(2)/ztick)])
if zkm
    set(gca,'YTickLabel',ztick*[ceil(zlims(1)/ztick):1:floor(zlims(2)/ztick)]/1000)
    ylabel('Elevation (km)','FontSize',fontsize)
else
    ylabel('Elevation (m)','FontSize',fontsize)
end
xlabel({'Distance from';'Grounding Line (km)'},'FontSize',fontsize)
set(gca,'PlotBoxAspectRatio',[(1000/exag)*(thesexlims(2)-thesexlims(1))/(zlims(2)-zlims(1)),1,1])
set(gca,'FontSize',fontsize)
set(gca,'YAxisLocation','right')
set(gca,'TickDir','out')
set(gca,'TickLength',[ticklength,ticklength])
title('c) Width-Averaged Topography','FontSize',fontsize)
    

% Save figure:
set(gcf,'PaperSize',pagesize)
set(gcf,'PaperPosition',[0,0,pagesize])
print('-dpng',figname,['-r',num2str(resolution)])

% Final display:
disp('Done!')
toc