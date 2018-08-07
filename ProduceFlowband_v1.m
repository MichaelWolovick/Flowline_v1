% ProduceFlowband_v1

% Mike Wolovick, 3/3/2016

% This script helps me produce the width-averaged profiles for a flowband
% model based on input datasets.  

% A merged velocity product with balance velocity and observations is used.

% The process starts by manually picking a flux gate near the calving
% front.  After that, the script propagates flowlines upstream from there.
% I choose a certain number of flowlines to discard from the edges,
% resulting in the final lateral boundaries of the flowband.  I also pick
% downstream extensions to continue the model domain into the fjord.

% Along-flow distance is interpolated from the flowlines onto the grid, and
% then smoothed.  The distance is used to compute width-averaged properties
% as a function of distance within the flowband.  Flowband width is
% computed in a way that conserves the area of the model domain.  

% NOTE:  this script does not have a polished UI.  Part of producing a
% flowband involves manual picking, so there are some places in the script
% where you have to stop and change parameters or run ginput.

clear all
%close all
tic

%% Parameters:

% File names:
inputfile='/home/mjw/Documents/CombinedGrids/GreenlandCombinedGrids_v1.mat';
outputfile='/home/mjw/Documents/FjordSillPaper/KangerFlowband_v1.mat';
%inputfile='/home/mike/Documents/Research/CombinedGrids/GreenlandCombinedGrids_v1.mat';
%outputfile='/home/mike/Documents/Research/FjordSillPaper/HelheimFlowband_v1.mat';

% Flowband lengthscales:
xwavelength=5e3;                % m
dx_flowband=1500;               % m

% Flowline tracking parameters: 
numflowlines=100;               % integer
steplength=250;                 % m
pathlength=1.5e5;               % m

% Surface smoothing wavelength:  (for balance velocity)
surfwavelength=2e4;             % m

% Minimum velocity error: (reported errors do not include baseline errors)
minvelerror=3;                  % m/yr

% Velocity error factor: (for weighting function)
velerrorfactor=1;               % unitless (>1 downweights obs)

% Number of distance smoothing iterations:
numsmoothingiterations=5;       % integer

% Work on existing picks?
loadpicks=0;                    % logical

%% Preparation:

% Load input data:
load(inputfile)

% Load picks from existing flowband file:
if loadpicks
    load(outputfile,'*_gate','numskipped*','numused*','*extension*')
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
    [~,SortedInd]=sort(SurfElev_smooth(:),'descend');
    
    % Loop through grid cells:
    for ii=1:ysize*xsize
        
        % Get subscript indices:
        [d1,d2]=ind2sub([ysize,xsize],SortedInd(ii));
        
        % Check whether there is ice here:
        if isnan(Mask(d1,d2))
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
    save(inputfile,'U_bal','V_bal','surfwavelength','-append')
    
else
    % Compute balance velocity magnitude:
    VelMag_bal=sqrt(U_bal.^2+V_bal.^2);
end

%% Merge Velocity and Direction Fields:

% Compute weighting between velocity observations and balance velocity:
VelWeight=1-exp(-VelMag./(velerrorfactor*VelError));
VelWeight(isnan(VelWeight))=0;

% Compute combined velocity grid:
VelMagMerged=VelWeight.*VelMag+(1-VelWeight).*VelMag_bal;
VelMagMerged(isnan(VelMag))=VelMag_bal(isnan(VelMag));

% Compute combined direction field:
DirMerged=VelWeight.*VelDir+(1-VelWeight).*SurfGradDir;
DirMerged(isnan(VelMag))=SurfGradDir(isnan(VelMag));
DirMerged=DirMerged./abs(DirMerged);

%% Pick Flux Gate and Compute Flowlines:

% Check whether we're picking a flux gate:
if loadpicks==0
    
    % Create an initial display:
    figure(1)
    hold off
    imagesc(X/1000,Y/1000,log10(VelMagMerged))
    set(gca,'ydir','normal')
    caxis([-1,4])
    xlims=get(gca,'Xlim');
    ylims=get(gca,'Ylim');
    set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    
    % Stop to fiddle with zoom:
    disp('Fiddle with the zoom, please')
    return
    
    % Create a flux gate near the ice front:
    [x_gate,y_gate]=ginput(2);
    x_gate=1000*x_gate;
    y_gate=1000*y_gate;
    hold on
    plot(x_gate/1000,y_gate/1000,'w')
    
end
    
%% Track flowlines upstream from the flux gate:
% determine number of steps: (upstream buffer to deal with edge effects)
numsteps=ceil((pathlength+3*xwavelength)/steplength);
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

% Check whether we're picking flowlines to skip at the edges:
if loadpicks==0
    
    % Plot flowlines:
    plot(X_flowlines/1000,Y_flowlines/1000,'w')
    
    % Plot lines of constant distance:
    plot(X_flowlines(1:100:end,:)'/1000,Y_flowlines(1:100:end,:)'/1000,'k')
    plot(X_flowlines(end,:)'/1000,Y_flowlines(end,:)'/1000,'k')
    
    % Plot first flowline:
    plot(X_flowlines(:,1)/1000,Y_flowlines(:,1)/1000,'k')
    
    % Pause to manually count flowlines in from the edge:
    disp('Tell me how many flowlines to kill from the edges')
    return
    
    % Change these values to indicate how far in to truncate the collection of
    % flowlines on each lateral edge.
    numskippedflowlines1=7;
    numskippedflowlines2=5;
    numusedflowlines=numflowlines-numskippedflowlines1-numskippedflowlines2;
    
    % Plot over the first flowline:
    plot(X_flowlines(:,1)/1000,Y_flowlines(:,1)/1000,'w')
    
    % Plot new lateral boundaries:
    plot(X_flowlines(:,numskippedflowlines1+1)/1000,Y_flowlines(:,numskippedflowlines1+1)/1000,'k')
    plot(X_flowlines(:,end-numskippedflowlines2)/1000,Y_flowlines(:,end-numskippedflowlines2)/1000,'k')
    
end

%% Pick lateral boundaries beyone the flux gate:

% Check if we're using existing picks:
if loadpicks==0
    
    % Display bed elevation:
    figure(1)
    hold off
    imagesc(X/1000,Y/1000,BedElev)
    hold on
    set(gca,'ydir','normal')
    caxis([-1000,2000])
    xlims=get(gca,'Xlim');
    ylims=get(gca,'Ylim');
    set(gca,'PlotBoxAspectRatio',[(xlims(2)-xlims(1))/(ylims(2)-ylims(1)),1,1])
    
    % Plot flowlines:
    plot(X_flowlines(:,numskippedflowlines1+1:end-numskippedflowlines2)/1000,Y_flowlines(:,numskippedflowlines1+1:end-numskippedflowlines2)/1000,'w')
    
    % Plot flux gate:
    plot(x_gate/1000,y_gate/1000,'k')
    
    % Stop to fiddle with zoom:
    disp('Fiddle with the zoom, please')
    return
    
    % Pick one boundary extension:
    [TheseX,TheseY]=ginput;
    X_extension1=[X_flowlines(1,numskippedflowlines1+1);TheseX*1000];
    Y_extension1=[Y_flowlines(1,numskippedflowlines1+1);TheseY*1000];
    plot(X_extension1/1000,Y_extension1/1000,'w')
    
    % Pick the other boundary extension:
    [TheseX,TheseY]=ginput;
    X_extension2=[X_flowlines(1,end-numskippedflowlines2);TheseX*1000];
    Y_extension2=[Y_flowlines(1,end-numskippedflowlines2);TheseY*1000];
    plot(X_extension2/1000,Y_extension2/1000,'w')
    
end

% Compute along-extension distance: (negative)
Dist_extension1=-cumsum([0;sqrt(diff(X_extension1).^2+diff(Y_extension1).^2)]);
Dist_extension2=-cumsum([0;sqrt(diff(X_extension2).^2+diff(Y_extension2).^2)]);

% Define flowband boundaries:
lastind=ceil(pathlength/steplength);
X_boundary=[flipud(X_extension1);X_flowlines(1:lastind,numskippedflowlines1+1);X_flowlines(lastind,numskippedflowlines1+1:end-numskippedflowlines2)';flipud(X_flowlines(1:lastind,end-numskippedflowlines2));X_extension2];
Y_boundary=[flipud(Y_extension1);Y_flowlines(1:lastind,numskippedflowlines1+1);Y_flowlines(lastind,numskippedflowlines1+1:end-numskippedflowlines2)';flipud(Y_flowlines(1:lastind,end-numskippedflowlines2));Y_extension2];

% Define expanded flowband boundaries: (including upstream BC zone)
X_bigboundary=[flipud(X_extension1);X_flowlines(:,numskippedflowlines1+1);X_flowlines(end,numskippedflowlines1+1:end-numskippedflowlines2)';flipud(X_flowlines(:,end-numskippedflowlines2));X_extension2];
Y_bigboundary=[flipud(Y_extension1);Y_flowlines(:,numskippedflowlines1+1);Y_flowlines(end,numskippedflowlines1+1:end-numskippedflowlines2)';flipud(Y_flowlines(:,end-numskippedflowlines2));Y_extension2];


%% Create the Flowband:

% Identify grid cells within the flowband:
InFlowBand=inpolygon(repmat(X,[ysize,1]),repmat(Y,[1,xsize]),X_bigboundary,Y_bigboundary);
numtotalgridpoints=sum(sum(InFlowBand));

% Identify grid cells within the extension:
InExtension=inpolygon(repmat(X,[ysize,1]),repmat(Y,[1,xsize]),[X_extension1;flipud(X_extension2)],[Y_extension1;flipud(Y_extension2)]);
numextensiongridpoints=sum(sum(InExtension));
numnonextensiongridpoints=numtotalgridpoints-numextensiongridpoints;

% Interpolate along-flow distance from the flowlines onto the grid:
% (along-flow distance measured from the flux gate)
% Create interpolant:
Interpolant1=scatteredInterpolant(reshape(X_flowlines(:,numskippedflowlines1+1:end-numskippedflowlines2),[],1),...
    reshape(Y_flowlines(:,numskippedflowlines1+1:end-numskippedflowlines2),[],1),...
    reshape(repmat(steplength*linspace(0,numsteps-1,numsteps)',[1,numusedflowlines]),[],1));
% Pre-allocate:
Distance=NaN*zeros(ysize,xsize);
X_grid=repmat(X,[ysize,1]);
Y_grid=repmat(Y,[1,xsize]);
% Interpolate:
Distance(InFlowBand)=Interpolant1(X_grid(InFlowBand),Y_grid(InFlowBand));

% Interpolate distance within the extension:
Interpolant2=scatteredInterpolant([X_extension1;flipud(X_extension2)],...
    [Y_extension1;flipud(Y_extension2)],...
    [Dist_extension1;flipud(Dist_extension2)]);
Distance(InExtension)=Interpolant2(X_grid(InExtension),Y_grid(InExtension));

% Smooth gridded distance:
if numsmoothingiterations>0
    % Identify cells just outside the border of the flowband:
    JustOutsideFlowBand=[false(ysize,1),[false(1,xsize-2);...
        InFlowBand(2:end-1,2:end-1)==0&(InFlowBand(1:end-2,1:end-2)|InFlowBand(1:end-2,2:end-1)|InFlowBand(1:end-2,3:end)|...
        InFlowBand(2:end-1,1:end-2)|InFlowBand(2:end-1,3:end)|...
        InFlowBand(3:end,1:end-2)|InFlowBand(3:end,2:end-1)|InFlowBand(3:end,3:end));...
        false(1,xsize-2)],false(ysize,1)];
    numjustoutsideflowband=sum(sum(JustOutsideFlowBand));
    JustOutsideInds=find(JustOutsideFlowBand);
    % Iteratively apply smoothing:
    for iteration=1:numsmoothingiterations
        % Replace NaN's around the border of the flowband:
        for ii=1:numjustoutsideflowband
            % Get subscript indices:
            [d1,d2]=ind2sub([ysize,xsize],JustOutsideInds(ii));
            % Isolate local 3x3 patch:
            TheseDistance=Distance(d1-1:d1+1,d2-1:d2+1);
            TheseInFlowBand=InFlowBand(d1-1:d1+1,d2-1:d2+1);
            % Average neighbors that are inside the flowband:
            Distance(d1,d2)=mean(TheseDistance(TheseInFlowBand));
        end
        % Smooth distance with a 3x3 filter:
        Distance=[NaN*zeros(ysize,1),[NaN*zeros(1,xsize-2);...
            (1/9)*(Distance(1:end-2,1:end-2)+Distance(1:end-2,2:end-1)+Distance(1:end-2,3:end)+...
            Distance(2:end-1,1:end-2)+Distance(2:end-1,2:end-1)+Distance(2:end-1,3:end)+...
            Distance(3:end,1:end-2)+Distance(3:end,2:end-1)+Distance(3:end,3:end));...
            NaN*zeros(1,xsize-2)],NaN*zeros(ysize,1)];
    end
end

% Create the (long) flowband distance vector:
mindist=max([Dist_extension1(end),Dist_extension2(end)]);
flowbandlength=pathlength+3*xwavelength-mindist;
numsteps_flowband=ceil(flowbandlength/steplength);
Distance_flowband=steplength*linspace(ceil(mindist/steplength),floor((pathlength+3*xwavelength)/steplength),numsteps_flowband)';

% Create transfer matrices from grid to flowband:
% Compute Gaussian weighting:
Weighting=exp(-.5*((repmat(Distance(InFlowBand)',[numsteps_flowband,1])-repmat(Distance_flowband,[1,numtotalgridpoints]))/(.5*xwavelength)).^2);
% Normalized transfer matrix:
NormalizedTransferMatrix=Weighting./repmat(sum(Weighting,2),[1,numtotalgridpoints]);
% Area-conserving tranfer matrix:
AreaConservingTransferMatrix=Weighting./repmat(sum(Weighting,1),[numsteps_flowband,1]); % prone to edge effects near upstream/downstream boundaries
% Compute Gaussian weighting, without extension:
Weighting=exp(-.5*((repmat(Distance(InFlowBand&InExtension==0)',[numsteps_flowband,1])-repmat(Distance_flowband,[1,numnonextensiongridpoints]))/(.5*xwavelength)).^2);
% Normalized transfer matrix:
NormalizedTransferMatrix_nonextension=Weighting./repmat(sum(Weighting,2),[1,numnonextensiongridpoints]);

% Compute flowband profiles: (simple averaging)
BedElev_flowband=NormalizedTransferMatrix*BedElev(InFlowBand);
SurfElev_flowband=NormalizedTransferMatrix*SurfElev(InFlowBand);
IceThick_flowband=SurfElev_flowband-BedElev_flowband;
SurfTemp_flowband=NormalizedTransferMatrix*SurfTemp(InFlowBand);
Accum_flowband=NormalizedTransferMatrix*Accum(InFlowBand);
Gflux_flowband=NormalizedTransferMatrix*Gflux(InFlowBand);

% Compute flowband velocity profile: (excludes extension cells)
VelMag_flowband=(NormalizedTransferMatrix_nonextension*(VelMagMerged(InFlowBand&InExtension==0).*IceThick(InFlowBand&InExtension==0)))./IceThick_flowband; % interpolates flux, then divides by icethick
VelWeight_flowband=(NormalizedTransferMatrix_nonextension*(VelWeight(InFlowBand&InExtension==0).*IceThick(InFlowBand&InExtension==0)))./IceThick_flowband; % weighting has same interpolation scheme as velocity

% Compute area-conserving flowband width:
Width_flowband=(AreaConservingTransferMatrix*(ones(numtotalgridpoints,1)*(dx^2)))/steplength;

% Truncate ice sheet properties to where the ice sheet exists:
SurfElev_flowband(Distance_flowband<0)=NaN;
IceThick_flowband(Distance_flowband<0)=NaN;
VelMag_flowband(Distance_flowband<0)=NaN;
VelWeight_flowband(Distance_flowband<0)=NaN;

% Remove upstream and downstream edge effects from flowband profiles:
% recompute flowband length:
mindist=max([Dist_extension1(end),Dist_extension2(end)])+3*xwavelength;
flowbandlength=pathlength-mindist;
numsteps_flowband=ceil(flowbandlength/steplength);
% remove edge effects:
BedElev_flowband=BedElev_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
SurfElev_flowband=SurfElev_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
IceThick_flowband=IceThick_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
VelMag_flowband=VelMag_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
VelWeight_flowband=VelWeight_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
SurfTemp_flowband=SurfTemp_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
Accum_flowband=Accum_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
Gflux_flowband=Gflux_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
Width_flowband=Width_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);
Distance_flowband=Distance_flowband(Distance_flowband>=mindist&Distance_flowband<=pathlength);


%% Display Results and Save Output:

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
% Plot flowlines:
plot(X_flowlines(1:lastind,numskippedflowlines1+1:end-numskippedflowlines2)/1000,Y_flowlines(1:lastind,numskippedflowlines1+1:end-numskippedflowlines2)/1000,'w')
% Plot boundaries:
plot(X_boundary/1000,Y_boundary/1000,'k')
% Plot lines of constant distance:
contour(X/1000,Y/1000,Distance,[0:1e4:pathlength],'k')
contour(X/1000,Y/1000,Distance,[0:-1e4:mindist],'w')
% Plot flux gate:
plot(x_gate/1000,y_gate/1000,'k')
% Plot extension:
plot(X_extension1/1000,Y_extension1/1000,'w')
plot(X_extension2/1000,Y_extension2/1000,'w')
% Pause to zoom properly:
return
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

% Create a profile figure:
figure(2)
hold off
% Plot bed and surface:
subplot(6,1,1)
plot(Distance_flowband/1000,BedElev_flowband,'k')
hold on
plot(Distance_flowband/1000,SurfElev_flowband,'k')
frontind=find(Distance_flowband>=0,1,'first');
plot([0,0],[BedElev_flowband(frontind),SurfElev_flowband(frontind)],'k')
plot([mindist,pathlength]/1000,[0,0],'--k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'XTickLabel',[])
ylabel('m')
title('Surface and Bed Profiles')
% Plot width:
subplot(6,1,2)
plot(Distance_flowband/1000,Width_flowband/1000,'k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'XTickLabel',[])
ylabel('km')
title('Flowband Width')
% Plot velocity:
subplot(6,1,3)
semilogy(Distance_flowband/1000,VelMag_flowband,'k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'XTickLabel',[])
ylabel('m/yr')
title('Ice Velocity')
% Plot surface temperature:
subplot(6,1,4)
plot(Distance_flowband/1000,SurfTemp_flowband,'k')
hold on
plot([mindist,pathlength]/1000,[0,0],'--k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'XTickLabel',[])
ylabel('\circ C')
title('Surface Temperature')
% Plot accumulation rate:
subplot(6,1,5)
plot(Distance_flowband/1000,Accum_flowband,'k')
hold on
plot([mindist,pathlength]/1000,[0,0],'--k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
set(gca,'XTickLabel',[])
ylabel('m/yr')
title('Surface Accumulation Rate')
% Plot geothermal flux:
subplot(6,1,6)
plot(Distance_flowband/1000,Gflux_flowband*1000,'k')
set(gca,'XDir','reverse')
xlim([mindist,pathlength]/1000)
xlabel('Distance Upstream from Calving Front (km)')
ylabel('mW/m^2')
title('Geothermal Flux')
% Save figure:
set(gcf,'PaperSize',[8,14])
set(gcf,'PaperPosition',[0,0,8,14])
figname=[outputfile(1:end-4),'_profiles.png'];
print('-dpng',figname,'-r300')

% Save output:
save(outputfile,'*_flowlines','*_flowband','num*','*_gate','*length','mindist','lastind','minvelerror','velerrorfactor','stereoprojection')

% Final display:
disp('Done!')
toc