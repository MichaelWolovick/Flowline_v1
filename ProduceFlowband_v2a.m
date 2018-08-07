% ProduceFlowband_v2

% Mike Wolovick, 3/3/2016

% This script helps me produce the width-averaged profiles for a flowband
% model based on input datasets.  

% A merged velocity product with balance velocity and observations is used.

% The process starts by manually picking a flux gate near the calving
% front, and inland boundaries of the domain.  After that, the
% script automatically computes along-flow distance upstream from that flux
% gate using the merged velocity field (this is equivalent to solving the
% age equation when the velocity field has been normalized to unit
% magnitude).  I also pick downstream extensions to continue the model 
% domain into the fjord.

% The along-flow distance is used to compute width-averaged properties
% as a function of distance within the flowband.  Flowband width is
% computed in a way that conserves the area of the model domain.

% The model domain does not necessarily extend as far upstream as the
% upstream margin that you pick.  The parameter "pathlength" determines how
% far upstream it goes.

% The interpolation from the grid to the flowband is done in segments to
% save system memory (I crashed my computer the first time I tried to do
% Thwaites).  

% Picking lateral boundaries is done in two installments.  I found it was
% necessary to pick the boundaries close to the calving front in a close
% zoom, then zoom out to pick the inland boundaries.

% NOTE:  this script involves manual picking, but it does not have a
% polished UI.  You may have to quit and fiddle with things manually during
% the picking process.

% The script doesn't know how to handle the situation where the manually
% chosen upstream boundary intersects the middle of the upstream distance
% contour.  The workaround is to choose a smaller pathlength.

% The script has the option to do flux-weighted ice surface and bed
% topography in addition to simple width-averaging.  Flux-weighted
% averaging biases the profile towards faster-flowing areas with deep
% troughs.

clear all
close all
tic

%% Parameters:

% File names:
inputfile='/net/mjw/CombinedGrids/AntarcticaCombinedGrids_v1.mat';
pickfile='/home/mjw/Documents/FjordSillPaper/Flowbands/PineIslandFlowband_v8.mat'; % ignored if loadpicks=0
outputfile='/home/mjw/Documents/FjordSillPaper/Flowbands/PineIslandFlowband_v9.mat';
%inputfile='/home/mike/Documents/Research/CombinedGrids/GreenlandCombinedGrids_v1.mat';
%outputfile='/home/mike/Documents/Research/FjordSillPaper/KangerFlowband_v1.mat';

% Flowband parameters:
xwavelength=5e3;                % m
steplength=250;                 % m
pathlength=4.5e5;               % m (can be commented if you are loading picks)

% Truncate region: (comment to disable truncation)
% xlims=[-1900,-1000]*1e3;        % m (Thwaites)
% ylims=[-800,0]*1e3;             % m (Thwaites)
xlims=[-1900,-900]*1e3;         % m (PIG)
ylims=[-500,150]*1e3;           % m (PIG)

% How to weight ice properties: (simple width-averaging or flux-weighted
% width-averaging)
iceweight='flux';                % 'avg' or 'flux'

% Surface smoothing wavelength:  (for balance velocity and direction)
surfwavelength=2e4;             % m

% Minimum velocity error: (reported errors do not include baseline errors)
minvelerror=3;                  % m/yr

% Velocity error factor: (for weighting function)
velerrorfactor=5;               % unitless (>1 downweights obs)

% Number of distance smoothing iterations:
numsmoothingiterations=5;       % integer (0 deactivates)

% Interpolation segmentation:
segmentlength=1e5;              % m
overlapfactor=3;                % unitless (multiples of xwavelength)

% Ice fraction cutoffs:
icefraccutoff=.95;              % unitless [0,1]
groundedfraccutoff=.5;          % unitless [0,1]

% Work on existing picks?
loadpicks=1;                    % logical

% Produce flowlines to help guide boundary picking?
doflowlines=1;                  % logical
numflowlines=10;                % integer

%% Preparation:

% Communicate:
disp('Loading data, doing initial processing')

% Load input data:
load(inputfile)

% Load picks from existing flowband file:
if loadpicks
    load(pickfile,'*boundary*','*_gate','pathlength')
    %     if exist('pathlength','var')==0
    %         load(pickfile,'pathlength')
    %     end
end

% Find truncating indices:
if exist('xlims','var') && exist('ylims','var')
    ind1=find(X>=xlims(1),1,'first');
    ind2=find(X<=xlims(2),1,'last');
    ind3=find(Y>=ylims(1),1,'first');
    ind4=find(Y<=ylims(2),1,'last');
else
    ind1=1;
    ind2=length(X);
    ind3=1;
    ind4=length(Y);
end

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

% Mask kludge:
if strcmp(inputfile,'/net/mjw/CombinedGrids/AntarcticaCombinedGrids_v1.mat')
    % Make a new mask in the same format as Greenland:
    NewMask=zeros(size(Mask));
    NewMask(Mask==0)=2;
    NewMask(Mask==1)=3;
    Mask=NewMask;
    clear NewMask
end

% Define grid size:
dx=X(2)-X(1); % assume dx=dy
[ysize,xsize]=size(SurfElev);

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
    disp('Computing balance velocities')
    
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
        if (Mask(d1,d2)~=2 && Mask(d1,d2)~=3) || d1==1 || d1==ysize || d2==1 || d2==xsize
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
disp('Merging velocity data with balance velocity')

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

%% Pick Flowband Boundaries:

% Check whether we're picking at all:
if loadpicks==0
    
    % Communicate:
    disp('Commencing picking')
    
    % Create an initial display:
    figure(1)
    hold off
    %     imagesc(X/1000,Y/1000,IceThick)
    %     set(gca,'ydir','normal')
    %     caxis([0,500])
    %     xlims=get(gca,'Xlim');
    %     ylims=get(gca,'Ylim');
    %     set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    %     hold on
    %     contour(X/1000,Y/1000,Mask==2|Mask==3,[.5,.5],'w')
    imagesc(X/1000,Y/1000,log10(VelMagMerged))
    set(gca,'ydir','normal')
    colormap(jet)
    caxis([-1,4])
    xlims=get(gca,'Xlim');
    ylims=get(gca,'Ylim');
    set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    hold on
    contour(X/1000,Y/1000,Mask==2|Mask==3,[.5,.5],'w')
    %contour(X/1000,Y/1000,isnan(Mask),[.5,.5],'w')
    
    % Stop to fiddle with zoom:
    title('Fiddle with the zoom, please.  You are picking the flux gate next.')
    pause
    
    % Fix the aspect ratio:
    xlims=get(gca,'Xlim');
    ylims=get(gca,'Ylim');
    set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    
    % Pick a flux gate near the ice front:
    title('Pick the flux gate')
    [x_gate,y_gate]=ginput(2);
    x_gate=1000*x_gate;
    y_gate=1000*y_gate;
    hold on
    plot(x_gate/1000,y_gate/1000,'w','LineWidth',2)
    
    % Produce flowlines to help guide lateral boundary interpretation:
    if doflowlines
        % determine number of steps: (upstream buffer to deal with edge effects)
        numsteps=ceil((pathlength+overlapfactor*xwavelength)/steplength);
        % pre-allocate memory:
        X_flowlines=zeros(numsteps,numflowlines);
        Y_flowlines=zeros(numsteps,numflowlines);
        % assign starting locations:
        X_flowlines(1,:)=linspace(x_gate(1),x_gate(2),numflowlines);
        Y_flowlines(1,:)=linspace(y_gate(1),y_gate(2),numflowlines);
        % loop through steps:
        for ii=1:numsteps-1
            % Interpolate direction:
            TheseDir=interp2(X,Y,DirMerged,X_flowlines(ii,:),Y_flowlines(ii,:));
            % Normalize direction:
            TheseDir=TheseDir./abs(TheseDir);
            % Step upstream:
            X_flowlines(ii+1,:)=X_flowlines(ii,:)-steplength*real(TheseDir);
            Y_flowlines(ii+1,:)=Y_flowlines(ii,:)-steplength*imag(TheseDir);
        end
        % Plot flowlines:
        plot(X_flowlines/1000,Y_flowlines/1000,'w','LineWidth',1)
    else
        % Quiver plot of direction to help interpretation:
        %quiver(X(1:5:end)/1000,Y(1:5:end)/1000,real(DirMerged(1:5:end,1:5:end)),imag(DirMerged(1:5:end,1:5:end)),.8,'ShowArrowHead','off','Color','w')
        % Contour plot of bed elevation:
        %contour(X/1000,Y/1000,BedElev,[-500:100:500],'w')
    end
    
    % Stop to fiddle with zoom:
    title('Fiddle with the zoom, please.  You are picking the downstream lateral margins next.')
    pause
    
    % Fix the aspect ratio:
    xlims=get(gca,'Xlim');
    ylims=get(gca,'Ylim');
    set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    
    % Pick first lateral boundary:
    % mark which boundary is first:
    h1=plot(x_gate(1)/1000,y_gate(1)/1000,'ok');
    % Go pick:
    title('Pick the first lateral boundary')
    [X_boundary1_downstream,Y_boundary1_downstream]=ginput;
    X_boundary1_downstream=X_boundary1_downstream*1000;
    Y_boundary1_downstream=Y_boundary1_downstream*1000;
    plot(X_boundary1_downstream/1000,Y_boundary1_downstream/1000,'k','LineWidth',2)
    
    % Pick second lateral boundary:
    % mark which boundary is second:
    set(h1,'XData',x_gate(2)/1000)
    set(h1,'YData',y_gate(2)/1000)
    % Go pick:
    title('Pick the second lateral boundary')
    [X_boundary2_downstream,Y_boundary2_downstream]=ginput;
    X_boundary2_downstream=X_boundary2_downstream*1000;
    Y_boundary2_downstream=Y_boundary2_downstream*1000;
    plot(X_boundary2_downstream/1000,Y_boundary2_downstream/1000,'k','LineWidth',2)
    
    %     % Change plot to velocity:
    %     hold off
    %     imagesc(X/1000,Y/1000,log10(VelMagMerged))
    %     set(gca,'ydir','normal')
    %     caxis([-1,4])
    %     xlim(xlims)
    %     ylim(ylims)
    %     set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    %     hold on
    %     contour(X/1000,Y/1000,Mask==2|Mask==3,[.5,.5],'w')
    %     plot(x_gate/1000,y_gate/1000,'w','LineWidth',2)
    %     plot(X_boundary1_downstream/1000,Y_boundary1_downstream/1000,'k','LineWidth',2)
    %     plot(X_boundary2_downstream/1000,Y_boundary2_downstream/1000,'k','LineWidth',2)
    %     h1=plot(x_gate(2)/1000,y_gate(2)/1000,'ok');
    
    % Stop to fiddle with the zoom:
    title('Fiddle with the zoom, please.  You are picking the inland lateral margins next.')
    pause
    
    % Fix the aspect ratio:
    xlims=get(gca,'Xlim');
    ylims=get(gca,'Ylim');
    set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    
    % Mark which boundary is first:
    set(h1,'XData',X_boundary1_downstream(end)/1000)
    set(h1,'YData',Y_boundary1_downstream(end)/1000)
    % Go pick:
    title('Pick the first lateral boundary')
    [X_boundary1_upstream,Y_boundary1_upstream]=ginput;
    X_boundary1_upstream=X_boundary1_upstream*1000;
    Y_boundary1_upstream=Y_boundary1_upstream*1000;
    % Merge downstream and upstream picks:
    X_boundary1=[X_boundary1_downstream;X_boundary1_upstream];
    Y_boundary1=[Y_boundary1_downstream;Y_boundary1_upstream];
    % Plot combined boundary:
    plot(X_boundary1/1000,Y_boundary1/1000,'k','LineWidth',2)
    
    % Mark which boundary is second:
    set(h1,'XData',X_boundary2_downstream(end)/1000)
    set(h1,'YData',Y_boundary2_downstream(end)/1000)
    % Go pick:
    title('Pick the second lateral boundary')
    [X_boundary2_upstream,Y_boundary2_upstream]=ginput;
    X_boundary2_upstream=X_boundary2_upstream*1000;
    Y_boundary2_upstream=Y_boundary2_upstream*1000;
    % Merge downstream and upstream picks:
    X_boundary2=[X_boundary2_downstream;X_boundary2_upstream];
    Y_boundary2=[Y_boundary2_downstream;Y_boundary2_upstream];
    % Plot combined boundary:
    plot(X_boundary2/1000,Y_boundary2/1000,'k','LineWidth',2)
    
    % Clear extraneous variables:
    clear *_upstream *_downstream
    
    % Pick upstream boundary:
    % Mark where the pick starts:
    set(h1,'XData',X_boundary1(end)/1000)
    set(h1,'YData',Y_boundary1(end)/1000)
    % Go pick:
    title('Pick the upstream boundary')
    [X_boundary3,Y_boundary3]=ginput;
    X_boundary3=X_boundary3*1000;
    Y_boundary3=Y_boundary3*1000;
    plot(X_boundary3/1000,Y_boundary3/1000,'w','LineWidth',2)
    
    % Change display to surface:
    hold off
    imagesc(X/1000,Y/1000,SurfElev)
    set(gca,'ydir','normal')
    xlim(xlims)
    ylim(ylims)
    set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    caxis([0,500])
    hold on
    h1=plot(x_gate(1)/1000,y_gate(1)/1000,'ok');
    plot(x_gate/1000,y_gate/1000,'w','LineWidth',2)
    plot(X_boundary1/1000,Y_boundary1/1000,'w','LineWidth',2)
    plot(X_boundary2/1000,Y_boundary2/1000,'w','LineWidth',2)
    plot(X_boundary3/1000,Y_boundary3/1000,'w','LineWidth',2)
    
    % Stop to fiddle with zoom:
    title('Fiddle with the zoom, please.  You are picking the fjord/offshore extension next.')
    pause
    
    % Fix the aspect ratio:
    xlims=get(gca,'Xlim');
    ylims=get(gca,'Ylim');
    set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    
    % Mark which boundary is first:
    set(h1,'XData',x_gate(1)/1000)
    set(h1,'YData',y_gate(1)/1000)
    % Go pick:
    title('Pick the first extension boundary')
    [X_boundary4,Y_boundary4]=ginput;
    X_boundary4=[x_gate(1);X_boundary4*1000];
    Y_boundary4=[y_gate(1);Y_boundary4*1000];
    plot(X_boundary4/1000,Y_boundary4/1000,'w','LineWidth',2)
    
    % Mark which boundary is second:
    set(h1,'XData',x_gate(2)/1000)
    set(h1,'YData',y_gate(2)/1000)
    % Go pick:
    title('Pick the second extension boundary')
    [X_boundary5,Y_boundary5]=ginput;
    X_boundary5=[x_gate(2);X_boundary5*1000];
    Y_boundary5=[y_gate(2);Y_boundary5*1000];
    plot(X_boundary5/1000,Y_boundary5/1000,'w','LineWidth',2)
    delete(h1)
    
    % Draw now:
    drawnow
    
end


%% Create the Flowband and Compute Along-Flow Distance:

% Communicate:
disp('Computing along-flow distance')

% Define combined flowband boundary: (including upstream BC zone)
X_bigboundary=[flipud(X_boundary4);X_boundary1;X_boundary3;flipud(X_boundary2);X_boundary5];
Y_bigboundary=[flipud(Y_boundary4);Y_boundary1;Y_boundary3;flipud(Y_boundary2);Y_boundary5];

% Create gridded x/y vectors:
X_grid=repmat(X,[ysize,1]);
Y_grid=repmat(Y,[1,xsize]);

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
    
    % Display:
%     if rem(iteration,5)==0
%         figure(1)
%         hold off
%         imagesc(X/1000,Y/1000,Distance/1000,[-100,1.5*pathlength/1000])
%         set(gca,'ydir','normal')
%         xlims=get(gca,'Xlim');
%         ylims=get(gca,'Ylim');
%         set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
%         hold on
%         plot(x_gate/1000,y_gate/1000,'w')
%         title(['Iteration=',num2str(iteration)])
%         drawnow
%         pause(.1)
%     end
    
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

% Compute distance gradient direction vector: (points downstream)
[DistGradX,DistGradY]=gradient(-Distance,dx,dx);
DistGradMag=sqrt(DistGradX.^2+DistGradY.^2);
DistDir=DistGradX./DistGradMag+1i*DistGradY./DistGradMag;

% Compute dot product of distance direction with velocity direction:
DotProduct=real(DistDir).*real(DirMerged)+imag(DistDir).*imag(DirMerged);

% Ensure no balance velocity is used downstream:
% (this is a fix for the problem I had in PIG where the velocity data and
% the ice geometry data had a different front position, resulting in the
% balance velocity reappearing once the real velocity data ran out)
VelWeight(Distance<0)=1;
VelMagMerged=VelWeight.*VelMag+(1-VelWeight).*VelMag_bal;
VelMagMerged(isnan(VelMag))=VelMag_bal(isnan(VelMag));
% Fix the mask (for the same reason):
Mask(Distance<0&VelMag==0)=0;

%% Define Final Flowband Boundaries:

% Communicate:
disp('Defining flowband boundaries')

% Interpolate distance onto the boundaries:
Dist_boundary1=interp2(X,Y,Distance,X_boundary1,Y_boundary1);
Dist_boundary2=interp2(X,Y,Distance,X_boundary2,Y_boundary2);
Dist_boundary3=interp2(X,Y,Distance,X_boundary3,Y_boundary3);
Dist_boundary4=interp2(X,Y,Distance,X_boundary4,Y_boundary4);
Dist_boundary5=interp2(X,Y,Distance,X_boundary5,Y_boundary5);

% Compute most negative distance to include in the extension:
mindist=max([Dist_boundary4(end),Dist_boundary5(end)])+overlapfactor*xwavelength;

% Start building the final flowband boundaries with extension boundary:
X_boundary=[interp1(Dist_boundary4,X_boundary4,mindist);flipud(X_boundary4(Dist_boundary4>mindist))];
Y_boundary=[interp1(Dist_boundary4,Y_boundary4,mindist);flipud(Y_boundary4(Dist_boundary4>mindist))];

% Find where distance is less than the pathlength along lateral margins:
ind1=find(Dist_boundary1<=pathlength,1,'last');
ind2=find(Dist_boundary2<=pathlength,1,'last');

% Concatenate first lateral margin onto the final flowband boundary:
X_boundary=[X_boundary;X_boundary1(1:ind1)];
Y_boundary=[Y_boundary;Y_boundary1(1:ind1)];

% Set distance to NaN outside of the flowband:
Distance(InFlowBand==0)=NaN;

% Find distance contour corresponding to pathlength:
ThisContour=contourc(X,Y,Distance,pathlength*[1,1]);

% Identify longest continuous contour segment:
startind=1;
done=0;
while done==0
    % Find stopind:
    stopind=startind+ThisContour(2,startind);
    % Check that this contour segment is a majority of the total:
    if stopind-startind-1>size(ThisContour,2)/2
        % We're good:
        done=1;
        break
    else
        % Check whether to keep trying or quit:
        if stopind>=size(ThisContour,2)
            error('Unable to find a continuous distance contour.')
        else
            startind=stopind+1;
        end
    end
end

% Check whether to add some of the upstream boundary:
if Dist_boundary3(1)<pathlength
    ind3=find(Dist_boundary3>pathlength,1,'first')-1;
    % Concatenate some of the upstream boundary onto the final boundary:
    X_boundary=[X_boundary;X_boundary3(1:ind3)];
    Y_boundary=[Y_boundary;Y_boundary3(1:ind3)];
end

% Check orientation of distance contour:
if sqrt((ThisContour(1,startind+1)-X_boundary(end))^2+(ThisContour(2,startind+1)-Y_boundary(end))^2)<sqrt((ThisContour(1,stopind)-X_boundary(end))^2+(ThisContour(2,stopind)-Y_boundary(end))^2)
    % Concatenate distance contour onto the final flowband boundary:
    X_boundary=[X_boundary;ThisContour(1,startind+1:stopind)'];
    Y_boundary=[Y_boundary;ThisContour(2,startind+1:stopind)'];
else
    % Concatenate distance contour onto the final flowband boundary:
    X_boundary=[X_boundary;flipud(ThisContour(1,startind+1:stopind)')];
    Y_boundary=[Y_boundary;flipud(ThisContour(2,startind+1:stopind)')];
end

% Check whether to add some more of the upstream boundary:
if Dist_boundary3(end)<pathlength
    ind3=find(Dist_boundary3>pathlength,1,'last')+1;
    % Concatenate some of the upstream boundary onto the final boundary:
    X_boundary=[X_boundary;X_boundary3(ind3:end)];
    Y_boundary=[Y_boundary;Y_boundary3(ind3:end)];
end

% Concatenate second lateral margin onto the final flowband boundary:
X_boundary=[X_boundary;flipud(X_boundary2(1:ind2))];
Y_boundary=[Y_boundary;flipud(Y_boundary2(1:ind2))];

% Concatenate second extension onto the final flowband boundary:
X_boundary=[X_boundary;X_boundary5(Dist_boundary5>mindist);interp1(Dist_boundary5,X_boundary5,mindist)];
Y_boundary=[Y_boundary;Y_boundary5(Dist_boundary5>mindist);interp1(Dist_boundary5,Y_boundary5,mindist)];

%% Interpolate Properties from Grid to Flowband:

% Communicate:
disp('Interpolating from grid to flowband')

% Create the (long) flowband distance vector:
mindist=max([Dist_boundary4(end),Dist_boundary5(end)]);
flowbandlength=pathlength+overlapfactor*xwavelength-mindist;
numsteps_flowband=ceil(flowbandlength/steplength);
Distance_flowband=steplength*linspace(ceil(mindist/steplength),floor((pathlength+overlapfactor*xwavelength)/steplength),numsteps_flowband)';

% Pre-allocate memory:
BedElev_flowband=zeros(numsteps_flowband,1);
if strcmp(iceweight,'flux')
    FluxWeightedBedElev_flowband=zeros(numsteps_flowband,1);
end
SurfElev_flowband=zeros(numsteps_flowband,1);
IceThick_flowband=zeros(numsteps_flowband,1);
VelMag_flowband=zeros(numsteps_flowband,1);
VelWeight_flowband=zeros(numsteps_flowband,1);
SurfTemp_flowband=zeros(numsteps_flowband,1);
Accum_flowband=zeros(numsteps_flowband,1);
Gflux_flowband=zeros(numsteps_flowband,1);
AreaConservingWidth_flowband=zeros(numsteps_flowband,1);
FluxConservingWidth_flowband=zeros(numsteps_flowband,1);
IceFrac_flowband=zeros(numsteps_flowband,1);
GroundedFrac_flowband=zeros(numsteps_flowband,1);

% Determine number of interpolation segments:
numsegments=ceil(flowbandlength/segmentlength);

% Compute segment boundaries in along-flow distance:
bigsegmentboundaries=[mindist+segmentlength*[0:1:numsegments-1]'-overlapfactor*xwavelength,mindist+segmentlength*[1:numsegments]'+overlapfactor*xwavelength]; % includes overlap
segmentboundaries=[mindist+segmentlength*[0:1:numsegments-1]',mindist+segmentlength*[1:numsegments]']; % does not include overlap

% Loop through interpolation segments:
for ii=1:numsegments
    
    % Identify points in this (overlapping) segment:
    InSegment=InFlowBand&Distance>=bigsegmentboundaries(ii,1)&Distance<=bigsegmentboundaries(ii,2);
    InSegment_flowband=Distance_flowband>=bigsegmentboundaries(ii,1)&Distance_flowband<=bigsegmentboundaries(ii,2);
    ThisDistance=Distance_flowband(InSegment_flowband);
    numinsegment=sum(InSegment_flowband);
    
    % Create transfer matrices from grid to flowband:
    % Compute number of useful grid points:
    numgridpoints=sum(sum(InSegment));
    numicegridpoints=sum(sum(InSegment&(Mask==2|Mask==3)));
    % Compute Gaussian distance weighting:
    Weighting=exp(-.5*((repmat(Distance(InSegment)',[numinsegment,1])-repmat(ThisDistance,[1,numgridpoints]))/(.5*xwavelength)).^2);
    % Normalized transfer matrix:
    TransferMatrix=Weighting./repmat(sum(Weighting,2),[1,numgridpoints]);
    % Area-conserving tranfer matrix:
    AreaConservingTransferMatrix=Weighting./repmat(sum(Weighting,1),[numinsegment,1]); % prone to edge effects near upstream/downstream boundaries (hence the overlapping buffer zones)
    % Check if we can compute anything requiring ice cells only:
    if numicegridpoints>0
        % Compute Gaussian distance weighting:
        Weighting=exp(-.5*((repmat(Distance(InSegment&(Mask==2|Mask==3))',[numinsegment,1])-repmat(ThisDistance,[1,numicegridpoints]))/(.5*xwavelength)).^2);
        % Normalized transfer matrix:
        IceTransferMatrix=Weighting./repmat(sum(Weighting,2),[1,numicegridpoints]);
        % Check if we're doing flux-weighted ice geometry:
        if strcmp(iceweight,'flux')
            % Scale weighting by ice flux:
            Weighting=Weighting.*repmat(VelMagMerged(InSegment&(Mask==2|Mask==3))'.*IceThick(InSegment&(Mask==2|Mask==3))',[numinsegment,1]);
            % Normalized transfer matrix:
            FluxWeightedTransferMatrix=Weighting./repmat(sum(Weighting,2),[1,numicegridpoints]);
        elseif strcmp(iceweight,'avg')==0
            error('Parameter "iceweight" must be a string equal to "flux" or "avg".')
        end
    end
    
    % Compute flowband profiles: (simple averaging, all cells)
    ThisBedElev=TransferMatrix*BedElev(InSegment);
    ThisSurfTemp=TransferMatrix*SurfTemp(InSegment);
    ThisAccum=TransferMatrix*Accum(InSegment);
    ThisGflux=TransferMatrix*Gflux(InSegment);
    ThisIceFrac=TransferMatrix*(Mask(InSegment)==2|Mask(InSegment)==3);
    ThisGroundedFrac=TransferMatrix*(Mask(InSegment)==2);
    
    % Compute area-conserving flowband width:
    ThisWidth1=(AreaConservingTransferMatrix*(ones(numgridpoints,1)*(dx^2)))/steplength;
    
    % Compute flowband profiles: (ice cells only)
    if numicegridpoints>0
        % Check if we're doing flux-weighted ice properties:
        if strcmp(iceweight,'flux')
            % Flux-weighted, width-averaged geometry:
            ThisSurfElev=FluxWeightedTransferMatrix*SurfElev(InSegment&(Mask==2|Mask==3));
            ThisIceThick=FluxWeightedTransferMatrix*IceThick(InSegment&(Mask==2|Mask==3));
            ThisFluxWeightedBedElev=FluxWeightedTransferMatrix*BedElev(InSegment&(Mask==2|Mask==3));
            % Flux-weighted, width-averaged velocity:
            ThisVelMag=FluxWeightedTransferMatrix*(VelMagMerged(InSegment&(Mask==2|Mask==3)).*DotProduct(InSegment&(Mask==2|Mask==3)));
            ThisVelWeight=FluxWeightedTransferMatrix*VelWeight(InSegment&(Mask==2|Mask==3));
        else
            % Width-averaged geometry:
            ThisSurfElev=IceTransferMatrix*SurfElev(InSegment&(Mask==2|Mask==3));
            ThisIceThick=IceTransferMatrix*IceThick(InSegment&(Mask==2|Mask==3));
            % Velocity from width-averaged flux: 
            ThisVelMag=(IceTransferMatrix*(VelMagMerged(InSegment&(Mask==2|Mask==3)).*DotProduct(InSegment&(Mask==2|Mask==3)).*IceThick(InSegment&(Mask==2|Mask==3))))./ThisIceThick; % interpolates flux, then divides by icethick
            ThisVelWeight=(IceTransferMatrix*(VelWeight(InSegment&(Mask==2|Mask==3)).*IceThick(InSegment&(Mask==2|Mask==3))))./ThisIceThick; % weighting has same interpolation scheme as velocity
        end
        % Compute width-averaged specific flux:
        ThisSpecificFlux=IceTransferMatrix*(VelMagMerged(InSegment&(Mask==2|Mask==3)).*DotProduct(InSegment&(Mask==2|Mask==3)).*IceThick(InSegment&(Mask==2|Mask==3)));
        % Compute flux-conserving flowband width:
        ThisWidth2=ThisWidth1.*ThisSpecificFlux./(ThisVelMag.*ThisIceThick);
        % note that Width2=Width1 if we are doing simple width-averaged
        % properties, but not if we're doing flux-weighted properties.
    else
        % Assign flux-conserving width the same as area-conserving width:
        ThisWidth2=ThisWidth1;
    end
    
    % Determine how much of the big segment vector is inside the small
    % (non-overlapping) segment:
    ThisInSegment=ThisDistance>=segmentboundaries(ii,1)&ThisDistance<=segmentboundaries(ii,2);
    
    % Identify points in the small (non-overlapping) segment:
    InSegment_flowband=Distance_flowband>=segmentboundaries(ii,1)&Distance_flowband<=segmentboundaries(ii,2);
    
    % Assign the non-overlapping parts of the segment vectors:
    BedElev_flowband(InSegment_flowband)=ThisBedElev(ThisInSegment);
    SurfTemp_flowband(InSegment_flowband)=ThisSurfTemp(ThisInSegment);
    Accum_flowband(InSegment_flowband)=ThisAccum(ThisInSegment);
    Gflux_flowband(InSegment_flowband)=ThisGflux(ThisInSegment);
    AreaConservingWidth_flowband(InSegment_flowband)=ThisWidth1(ThisInSegment);
    FluxConservingWidth_flowband(InSegment_flowband)=ThisWidth2(ThisInSegment);
    IceFrac_flowband(InSegment_flowband)=ThisIceFrac(ThisInSegment);
    GroundedFrac_flowband(InSegment_flowband)=ThisGroundedFrac(ThisInSegment);
    if numicegridpoints>0
        SurfElev_flowband(InSegment_flowband)=ThisSurfElev(ThisInSegment);
        IceThick_flowband(InSegment_flowband)=ThisIceThick(ThisInSegment);
        if strcmp(iceweight,'flux')
            FluxWeightedBedElev_flowband(InSegment_flowband)=ThisFluxWeightedBedElev(ThisInSegment);
        end
        VelMag_flowband(InSegment_flowband)=ThisVelMag(ThisInSegment);
        VelWeight_flowband(InSegment_flowband)=ThisVelWeight(ThisInSegment);
    end
end

% Truncate ice sheet properties to where the ice sheet exists:
SurfElev_flowband(IceFrac_flowband<=icefraccutoff)=NaN;
IceThick_flowband(IceFrac_flowband<=icefraccutoff)=NaN;
VelMag_flowband(IceFrac_flowband<=icefraccutoff)=NaN;
VelWeight_flowband(IceFrac_flowband<=icefraccutoff)=NaN;
if strcmp(iceweight,'flux')
    FluxWeightedBedElev_flowband(IceFrac_flowband<=icefraccutoff)=NaN;
end

% Remove upstream and downstream edge effects from flowband profiles:
% recompute flowband length:
mindist=max([Dist_boundary4(end),Dist_boundary5(end)])+overlapfactor*xwavelength;
flowbandlength=pathlength-mindist;
numsteps_flowband=ceil(flowbandlength/steplength);
% remove edge effects:
BedElev_flowband=BedElev_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
if strcmp(iceweight,'flux')
    FluxWeightedBedElev_flowband=FluxWeightedBedElev_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
end
SurfElev_flowband=SurfElev_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
IceThick_flowband=IceThick_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
VelMag_flowband=VelMag_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
VelWeight_flowband=VelWeight_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
SurfTemp_flowband=SurfTemp_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
Accum_flowband=Accum_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
Gflux_flowband=Gflux_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
AreaConservingWidth_flowband=AreaConservingWidth_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
FluxConservingWidth_flowband=FluxConservingWidth_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
IceFrac_flowband=IceFrac_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
GroundedFrac_flowband=GroundedFrac_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
Distance_flowband=Distance_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);

% Clear transfer matrices:
clear *TransferMatrix* *Weighting InSegment*

% Fix flux-conserving width near the calving front: (starting two
% wavelengths in)
icefrontind=find(IceFrac_flowband>icefraccutoff,1,'first');
icefrontdist=Distance_flowband(icefrontind);
ind1=find(Distance_flowband>=icefrontdist+2*xwavelength,1,'first');
FluxConservingWidth_flowband(1:ind1)=FluxConservingWidth_flowband(ind1)+(1-exp(-.5*((Distance_flowband(1:ind1)-Distance_flowband(ind1))/(.5*xwavelength)).^2)).*(AreaConservingWidth_flowband(1:ind1)-FluxConservingWidth_flowband(ind1));

% % Kludge fix for the Jakobshavn flowband:
% ind1=find(Distance_flowband>=0,1,'first');
% ind2=find(Distance_flowband>=4*xwavelength,1,'first');
% FluxConservingWidth_flowband(1:ind1)=AreaConservingWidth_flowband(1:ind1);
% FluxConservingWidth_flowband(ind1:ind2)=linspace(AreaConservingWidth_flowband(ind1),FluxConservingWidth_flowband(ind2),ind2-ind1+1)';


%% Ensure Consistency Between Bed, Surface and Ice Thickness:

% There are multiple constraints on bed elevation, surface elevation, and
% ice thickness that must be respected in the flowband variables.  The
% definition of ice thickness (H=S-B) must hold where the ice is grounded.
% The flotation condition must be respected for floating ice, and grounded 
% ice must be thicker than the flotation thickness.  All three variables
% should, ideally, be close to the values derived from the flowband
% averaging process above.  And finally, we would like continuity of the
% final products across the grounding line.  Unfortunately, it is
% impossible to meet all of these constraints simultaneously.

% For grounded ice, the following equation must hold exactly:

% BedElev = SurfElev - IceThick

% In general, this equation does not hold after interpolation from the grid
% to the flowband.  The misfit is caused by two things: 1) bed
% interpolation includes grid cells without ice, and 2) ice thickness and
% surface elevation may be interpolated in a flux-weighted manner, while
% bed elevation is always interpolated as a simple across-flow average.
% (The bed cannot be interpolated in a flux-weighted manner because we need
% the bed profile to extend beyond the present-day ice front.)

% For grounded ice, it is simple enough to replace the interpolated bed
% elevation with one derived from surface and ice thickness.  In the case
% of flux-weighted interpolation, this has the effect of deepening the bed
% to account for fast-flowing troughs.  However, applying this correction
% for grounded ice has the effect of producing an unphysical cliff at the
% grounding line.  We would like a merged bed elevation product without a
% discontinuity in value or slope.


% Old Solution:

% We use a weighting scheme offshore of the grounding line to ensure a
% smooth bed product.  The weighting scheme gradually transitions between a
% linear extrapolation of the grounded bed and the regular width-averaged
% bed computed earlier.  This ensures both that bed, surface, and ice 
% thickness are exactly consistent for grounded ice, and also that the bed
% elevation is continuous in both value and slope across the grounding
% line.  The weighting scheme uses the same "xwavelength" as the
% interpolation from the grid to the flowband.
% This method assumes that there is only one grounding line in the
% flowband.

% Problem With Old Solution:
% The old scheme produced an unphysical sill just offshore of the grounding
% line when the bed just upstream is steeply overdeepened.  When I ran the 
% model forward in time, the grounding line advanced onto this sharp 
% unphysical sill and was artificially stabilized there.  


% New Solution:
% The flux-weighted bed is continued beyond the grounding line.  If a
% floating shelf is present, this will create a continuous bed product
% automatically.  If a floating shelf is not present then it doesn't
% matter, since I don't trust the fjord bathymetry grids in Greenland
% anyway and those are modified in a later script.


% Communicate:
disp('Ensuring consistency between bed, surface, and ice thickness.')

% Merge the weighted and the unweighted bed:
if strcmp(iceweight,'flux')
    MergedBed=IceFrac_flowband.*FluxWeightedBedElev_flowband+(1-IceFrac_flowband).*BedElev_flowband;
    MergedBed(isnan(FluxWeightedBedElev_flowband))=BedElev_flowband(isnan(FluxWeightedBedElev_flowband));
    BedElev_flowband=MergedBed;
    clear MergedBed
end

% Identify grounding line and ice front:
%glind=find(IceFrac_flowband>icefraccutoff&SurfElev_flowband-IceThick_flowband<=BedElev_flowband,1,'first');
frontind=find(IceFrac_flowband>=icefraccutoff,1,'first');
glind=max(find(GroundedFrac_flowband>=groundedfraccutoff|SurfElev_flowband-IceThick_flowband<=BedElev_flowband,1,'first'),frontind);

% Ensure exact consistency in grounded ice:
BedElev_flowband(glind:end)=SurfElev_flowband(glind:end)-IceThick_flowband(glind:end);

% Smooth over the merge point between weighted and unweighted beds:
if strcmp(iceweight,'flux')
    % Create smooth bed:
    SmoothBed=intuitive_lowpass(BedElev_flowband,xwavelength/steplength);
    % Create weighting:
    SmoothWeight=exp(-.5*((Distance_flowband-Distance_flowband(frontind))/(.5*xwavelength)).^2);
    % Create Merged Product:
    BedElev_flowband=SmoothWeight.*SmoothBed+(1-SmoothWeight).*BedElev_flowband;
end
error('end here')
%% Compute Flux Across the Lateral Boundaries of the Flowband:

% Communicate:
disp('Computing flux across the lateral boundaries')

% Identify cells just outside the border of the flowband:
JustOutsideFlowBand=[false(ysize,1),[false(1,xsize-2);...
    InFlowBand(2:end-1,2:end-1)==0&(InFlowBand(1:end-2,1:end-2)|InFlowBand(1:end-2,2:end-1)|InFlowBand(1:end-2,3:end)|...
    InFlowBand(2:end-1,1:end-2)|InFlowBand(2:end-1,3:end)|...
    InFlowBand(3:end,1:end-2)|InFlowBand(3:end,2:end-1)|InFlowBand(3:end,3:end));...
    false(1,xsize-2)],false(ysize,1)];
numjustoutsideflowband=sum(sum(JustOutsideFlowBand));
JustOutsideInds=find(JustOutsideFlowBand);

% Replace NaNs in distance around the border of the flowband:
for ii=1:numjustoutsideflowband
    % Get subscript indices:
    [d1,d2]=ind2sub([ysize,xsize],JustOutsideInds(ii));
    % Isolate local 3x3 patch:
    TheseDistance=Distance(d1-1:d1+1,d2-1:d2+1);
    TheseUsefulCells=InFlowBand(d1-1:d1+1,d2-1:d2+1);
    % Average neighbors that are inside the flowband:
    Distance(d1,d2)=mean(TheseDistance(TheseUsefulCells));
end

% Replace NaNs in distance one step further out:
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

% Interpolate distance onto the boundary:
Dist_boundary=interp2(X,Y,Distance,X_boundary,Y_boundary);

% Interpolate x/y coordinates of the lateral boundaries:
% (evenly spaced in distance)
% Find end of first lateral segment:
ind2=find(Dist_boundary>=pathlength-dx,1,'first');
% Interpolate x/y coordinates:
X1=interp1(Dist_boundary(1:ind2),X_boundary(1:ind2),Distance_flowband);
Y1=interp1(Dist_boundary(1:ind2),Y_boundary(1:ind2),Distance_flowband);
% Find start of second lateral segment:
ind1=ind2-1+find(Dist_boundary(ind2:end)>=pathlength-dx,1,'last');
% Interpolate x/y coordinates:
X2=interp1(Dist_boundary(ind1:end),X_boundary(ind1:end),Distance_flowband);
Y2=interp1(Dist_boundary(ind1:end),Y_boundary(ind1:end),Distance_flowband);

% Compute boundary normal vectors:
% Compute gradient of boundary coordinates:
GradX1=gradient(X1,steplength);
GradY1=gradient(Y1,steplength);
GradX2=-gradient(X2,steplength); % gradient w.r.t. increasing distance
GradY2=-gradient(Y2,steplength); % gradient w.r.t. increasing distance
% Compute boundary coordinate gradient magnitude:
GradMag1=sqrt(GradX1.^2+GradY1.^2); % normalized by nominal steplength
GradMag2=sqrt(GradX2.^2+GradY2.^2); % normalized by nominal steplength
% Compute along-boundary direction: (real=x, imag=y)
Dir1=GradX1./GradMag1+1i*GradY1./GradMag1;
Dir2=GradX2./GradMag2+1i*GradY2./GradMag2;
% Rotate 90 degrees to produce boundary normal vector:
BoundaryNormal1=imag(Dir1)-1i*real(Dir1);
BoundaryNormal2=-imag(Dir2)+1i*real(Dir2); % normal points inwards

% Compute direction of flux gate: (points from 1 to 2)
gatemag=sqrt((x_gate(2)-x_gate(1))^2+(y_gate(2)-y_gate(1))^2);
gatedir=(x_gate(2)-x_gate(1))/gatemag+1i*(y_gate(2)-y_gate(1))/gatemag;

% Find gate index:
ind1=find(Distance_flowband>=0,1,'first');

% Check boundary normals at the flux gate to ensure they point inward:
if real(BoundaryNormal1(ind1))*real(gatedir)+imag(BoundaryNormal1(ind1))*imag(gatedir)<0
    % Flip direction:
    BoundaryNormal1=-BoundaryNormal1;
end
if real(BoundaryNormal2(ind1))*real(gatedir)+imag(BoundaryNormal2(ind1))*imag(gatedir)>0
    % Flip direction:
    BoundaryNormal2=-BoundaryNormal2;
end

% Check that boundary normals point in opposite direction:
if real(BoundaryNormal1(ind1))*real(BoundaryNormal2(ind1))+imag(BoundaryNormal1(ind1))*imag(BoundaryNormal2(ind1))>0
    error('Lateral boundary normal vectors are facing in the same direction.')
end

% Interpolate flux magnitude and direction onto the boundary:
FluxMag1=interp2(X,Y,VelMagMerged.*IceThick,X1,Y1);
FluxMag2=interp2(X,Y,VelMagMerged.*IceThick,X2,Y2);
VelDir1=interp2(X,Y,DirMerged,X1,Y1);
VelDir2=interp2(X,Y,DirMerged,X2,Y2);

% Ensure interpolated velocity direction has unit amplitude:
VelDir1=VelDir1./abs(VelDir1);
VelDir2=VelDir2./abs(VelDir2);

% Compute flux across boundaries:
Influx_flowband=(real(VelDir1).*real(BoundaryNormal1)+imag(VelDir1).*imag(BoundaryNormal1)).*FluxMag1.*GradMag1...
    +(real(VelDir2).*real(BoundaryNormal2)+imag(VelDir2).*imag(BoundaryNormal2)).*FluxMag2.*GradMag2;

% Low-pass filter influx:
Influx_flowband(isnan(Influx_flowband))=0;
Influx_flowband=intuitive_lowpass(Influx_flowband,xwavelength/steplength);

% % Display normals to check them:
% figure(1)
% hold off
% plot(X1/1000,Y1/1000,'k')
% hold on
% plot(X2/1000,Y2/1000,'r')
% plot(x_gate(1)/1000,y_gate(1)/1000,'ok')
% plot(x_gate/1000,y_gate/1000,'k','LineWidth',2)
% xlims=get(gca,'Xlim');
% ylims=get(gca,'Ylim');
% set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
% figure(2)
% hold off
% plot(real(BoundaryNormal1),'k')
% hold on
% plot(real(BoundaryNormal2),'r')
% error('Check those mothafuckin normal vectors, homie!')

%% Display Results and Save Output:

% Communicate:
disp('Creating final figures')

% Create a final map figure:
% Basic figure:
figure(1)
hold off
imagesc(X/1000,Y/1000,log10(VelMagMerged))
hold on
set(gca,'ydir','normal')
caxis([-1,4])
xlims=get(gca,'Xlim');
ylims=get(gca,'Ylim');
set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
% Plot contours of distance:
Distance(InFlowBand==0)=NaN;
contour(X/1000,Y/1000,Distance,[0:2e4:pathlength],'k')
contour(X/1000,Y/1000,Distance,[0:-2e4:mindist],'w')
% Plot boundaries as closed shapes (white on black):
plot([X_boundary;X_boundary(1)]/1000,[Y_boundary;Y_boundary(1)]/1000,'k','LineWidth',2)
plot([X_boundary;X_boundary(1)]/1000,[Y_boundary;Y_boundary(1)]/1000,'w','LineWidth',1)
% Plot flux gate:
plot(x_gate/1000,y_gate/1000,'k')
% Pause to zoom properly:
title('Fiddle with the zoom, please (final map display)')
figure(1)
drawnow
pause
% Set up axes and title:
xlims=get(gca,'Xlim');
ylims=get(gca,'Ylim');
set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
xlabel('Easting (km)')
ylabel('Northing (km)')
title('Flowband Model Domain')
% Set up colorbar:
hc=colorbar;
set(get(hc,'YLabel'),'String','Velocity Magnitude, log_{10}[m/yr]')
% Save figure:
set(gcf,'PaperSize',[10,8])
set(gcf,'PaperPosition',[0,0,10,8])
figname=[outputfile(1:end-4),'_map.png'];
print('-dpng',figname,'-r300')

% figure(1)
% imagesc(X/1000,Y/1000,log10(VelMagMerged),[0,4])
% set(gca,'ydir','normal')
% hold on
% plot([X_boundary;X_boundary(1)]/1000,[Y_boundary;Y_boundary(1)]/1000,'w','LineWidth',1)
% contour(X/1000,Y/1000,Distance,Distance_flowband(frontind)*[1,1],'w')
% contour(X/1000,Y/1000,Distance,Distance_flowband(glind)*[1,1],'w')
% contour(X/1000,Y/1000,Mask,2.5*[1,1],'w')
% error('end here')

%% Create a profile figure:
figure(2)
hold off
% Plot bed, surface, and ice bottom:
subplot(2,4,1)
plot(Distance_flowband/1000,BedElev_flowband,'k')
hold on
plot(Distance_flowband/1000,SurfElev_flowband,'k')
plot(Distance_flowband/1000,SurfElev_flowband-IceThick_flowband,'k')
plot(Distance_flowband(frontind)*[1,1]/1000,[SurfElev_flowband(frontind),SurfElev_flowband(frontind)-IceThick_flowband(frontind)],'k')
plot([mindist,pathlength]/1000,[0,0],'--k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'TickDir','out')
set(gca,'XTickLabel',[])
ylabel('m')
title('Surface and Bed Profiles')
% Plot width:
subplot(2,4,2)
plot(Distance_flowband/1000,AreaConservingWidth_flowband/1000,'k')
if strcmp(iceweight,'flux')
    hold on
    plot(Distance_flowband/1000,FluxConservingWidth_flowband/1000,'--k')
    legend('Area-Conserving','Flux-Conserving','Location','NorthEast')
end
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
ylims=get(gca,'YLim');
ylim([0,ylims(2)])
set(gca,'TickDir','out')
set(gca,'XTickLabel',[])
ylabel('km')
title('Flowband Width')
% Plot velocity:
subplot(2,4,3)
semilogy(Distance_flowband/1000,VelMag_flowband,'k')
set(gca,'XDir','reverse')
set(gca,'YTick',10.^[-1:1:5])
xlim([mindist,pathlength]/1000)
set(gca,'TickDir','out')
set(gca,'XTickLabel',[])
ylabel('m/yr')
title('Ice Velocity')
% Plot velocity source:
subplot(2,4,4)
plot(Distance_flowband/1000,100*VelWeight_flowband,'k')
set(gca,'XDir','reverse')
ylim([0,100])
xlim([mindist,pathlength]/1000)
set(gca,'TickDir','out')
set(gca,'XTickLabel',[])
ylabel('Percent Data')
title('Velocity Source')
% Plot later influx:
subplot(2,4,5)
plot(Distance_flowband/1000,Influx_flowband./AreaConservingWidth_flowband,'k')
hold on
plot([mindist,pathlength]/1000,[0,0],'--k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'TickDir','out')
xlabel('Distance Upstream from Flux Gate (km)')
ylabel('m/yr')
title('Lateral Influx / Width (equivalent accum)')
% Plot surface temperature:
subplot(2,4,6)
plot(Distance_flowband/1000,SurfTemp_flowband,'k')
hold on
plot([mindist,pathlength]/1000,[0,0],'--k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'TickDir','out')
xlabel('Distance Upstream from Flux Gate (km)')
ylabel('\circ C')
title('Surface Temperature')
% Plot accumulation rate:
subplot(2,4,7)
plot(Distance_flowband/1000,Accum_flowband,'k')
hold on
plot([mindist,pathlength]/1000,[0,0],'--k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'TickDir','out')
xlabel('Distance Upstream from Flux Gate (km)')
ylabel('m/yr')
title('Surface Accumulation Rate')
% Plot geothermal flux:
subplot(2,4,8)
plot(Distance_flowband/1000,Gflux_flowband*1000,'k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'TickDir','out')
xlabel('Distance Upstream from Flux Gate (km)')
ylabel('mW/m^2')
title('Geothermal Flux')
% Save figure:
set(gcf,'PaperSize',[18,6])
set(gcf,'PaperPosition',[0,0,18,6])
figname=[outputfile(1:end-4),'_profiles.png'];
print('-dpng',figname,'-r300')

%% Save output:
save(outputfile,'*_flowband','*boundary*','num*','*_gate','*length','mindist','minvelerror','iceweight','*fraccutoff','velerrorfactor','glind','frontind','stereoprojection')

% Final display:
disp('Done!')
toc