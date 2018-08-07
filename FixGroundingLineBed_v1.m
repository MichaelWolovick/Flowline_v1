% FixGroundingLineBed_v1

% Mike Wolovick, 5/10/2016

% This script fixes the bed topography data for the Helheim, Kanger, and
% Jakobshavn flowbands.  In these flowbands, ice thickness and bed 
% elevation near the calving front are unreliable.  The width-averaged 
% input data imply a small floating tongue near the front.  If this tongue
% is allowed to float freely, it pops up into a position inconsistent with 
% the surface observations.  

% In addition, bed topography data from the input grids is unreliable in
% the fjords (contains obvious gridding artifacts).  Fjord profiles are 
% sometimes inconsistent with oceanographic results presented in Sutherland
% et al (2014) and several Straneo papers.  

% The ungrounded bed is therefore replaced by a flat bed at the fjord mean
% depth and a sill near the fjord mouth.  The fjord depth and sill depth
% are taken from Table 1 in Sutherland et al (2014).  The ice front is 
% moved back to the first grounded point.  The under-ice bed is smoothly
% joined with the fjord bed.

% The outer sill is gaussian shaped with a center 3 standard deviations in
% from the edge of the model domain.


% The erroneous floating tongues were caused by mixing the Morlighem bed
% with the Bamber surface.  When I use the Morlighem DEM for both bed and
% surface, they went away.  However, the fjord bathymetry still required
% correcting.

clear all
tic

%% Parameters:

% File names and paths:
inputfolder='/home/mjw/Documents/FjordSillPaper/Flowbands/';
inputfiles={'HelheimFlowband_v5.mat';...
    'HelheimFlowband_v6.mat';...
    'HelheimFlowband_v7.mat';...
    'KangerFlowband_v4.mat';...
    'KangerFlowband_v5.mat';...
    'KangerFlowband_v6.mat';...
    'JakobshavnFlowband_v8.mat';...
    'JakobshavnFlowband_v9.mat';...
    'JakobshavnFlowband_v10.mat'};

% Fjord and sill geometry:
fjorddepths=[780,780,780,840,840,840,770,770,770];  % [1xn] m >0
silldepths=[550,550,550,450,450,450,265,265,265];   % [1xn] m >0
sillwidth=1e4;                                      % m (2*stdev)

% Physical parameters:
rho_i=917;                 % kg/m^3
rho_sw=1028;               % kg/m^3
g=9.8;                     % m/s^2


%% Work:

% Loop through input files:
for thisfile=1:length(inputfiles)
    
    % Load input file:
    load([inputfolder,inputfiles{thisfile}],'Distance_flowband','BedElev_flowband','IceThick_flowband','SurfElev_flowband','steplength','xwavelength')
    
    % Record original variables before the change:
    OriginalBedElev_flowband=BedElev_flowband;
    OriginalIceThick_flowband=IceThick_flowband;
    OriginalSurfElev_flowband=SurfElev_flowband;
    
    % Compute hydraulic head:
    HydroHead_flowband=BedElev_flowband+(rho_i/rho_sw)*IceThick_flowband;
    
    % Locate first grounded cell:
    firsticeind=find(HydroHead_flowband>0,1,'first'); 
    
    % Locate old ice front:
    originalfirsticeind=find(isnan(IceThick_flowband)==0,1,'first');
    
    % Set erroneously floating cells to NaN:
    IceThick_flowband(1:firsticeind-1)=NaN;
    SurfElev_flowband(1:firsticeind-1)=NaN;
    
    % Compute bed value and slope at the grounding line:
    glbedvalue=BedElev_flowband(firsticeind);
    glbedslope=(BedElev_flowband(firsticeind)-BedElev_flowband(firsticeind+1))/steplength;
    
    % Create linear extrapolation of grounded bed:
    BedExtrap=glbedvalue+glbedslope*(Distance_flowband(firsticeind)-Distance_flowband(1:firsticeind-1));
    
    % Compute weighting function:
    BedExtrapWeight=exp(-.5*((Distance_flowband(1:firsticeind-1)-Distance_flowband(firsticeind))/(.5*xwavelength)).^2);
    
    % Create merged bed:
    BedElev_flowband(1:firsticeind-1)=-fjorddepths(thisfile)+BedExtrapWeight.*(BedExtrap+fjorddepths(thisfile));
    
    % Add sill:
    sillcentx=Distance_flowband(1)+1.5*sillwidth;
    BedElev_flowband(1:firsticeind-1)=BedElev_flowband(1:firsticeind-1)+(fjorddepths(thisfile)-silldepths(thisfile))*...
        exp(-.5*((Distance_flowband(1:firsticeind-1)-sillcentx)/(.5*sillwidth)).^2);
    
    % Make figure:
    figure(1)
    
    % Plot original geometry and geometry changes:
    hold off
    h1=plot(Distance_flowband/1000,OriginalBedElev_flowband,'k');
    hold on
    plot(Distance_flowband/1000,OriginalSurfElev_flowband,'k')
    plot(Distance_flowband(originalfirsticeind)*[1,1]/1000,[OriginalBedElev_flowband(originalfirsticeind),OriginalSurfElev_flowband(originalfirsticeind)],'k')
    if firsticeind~=originalfirsticeind
        plot(Distance_flowband(firsticeind:originalfirsticeind)/1000,(1-rho_i/rho_sw)*OriginalIcethick_flowband(firsticeind:originalfirsticeind),'r')
        plot(Distance_flowband(firsticeind:originalfirsticeind)/1000,(-rho_i/rho_sw)*OriginalIcethick_flowband(firsticeind:originalfirsticeind),'r')
        plot(Distance_flowband(originalfirsticeind)*[1,1]/1000,[(-rho_i/rho_sw)*OriginalIcethick_flowband(originalfirsticeind),(1-rho_i/rho_sw)*OriginalIcethick_flowband(originalfirsticeind)],'r')
        plot(Distance_flowband(firsticeind)*[1,1]/1000,[OriginalBedElev_flowband(firsticeind),OriginalSurfElev_flowband(firsticeind)],'r')
    end
    h2=plot(Distance_flowband(1:firsticeind-1)/1000,BedElev_flowband(1:firsticeind-1),'r');
    plot([Distance_flowband(1),Distance_flowband(end)]/1000,[0,0],'--k')
    set(gca,'XDir','reverse')
    xlim([Distance_flowband(1),Distance_flowband(end)]/1000)
    legend([h1;h2],'Original','Changes','location','NorthEast')
    xlabel('Distance (km)')
    ylabel('Elevation (m)')
    title(['Ice Geometry Changes in ',inputfiles{thisfile}(1:end-4)],'interpreter','none')
    
    % Pause to finish drawing:
    drawnow
    
    % Save figure:
    set(gcf,'PaperSize',[11,6])
    set(gcf,'PaperPosition',[0,0,11,6])
    figname=[inputfolder,inputfiles{thisfile}(1:end-4),'_GeometryModifications.png'];
    print('-dpng',figname,'-r300')
    
    % Save results:
    thissilldepth=silldepths(thisfile);
    thisfjorddepth=fjorddepths(thisfile);
    save([inputfolder,inputfiles{thisfile}],'*_flowband','Original*','sillcentx','sillwidth','thissilldepth','thisfjorddepth','-append')
    
end


% Final display:
disp('Done!')
toc