% FlowbandAblationCurves_v1

% Mike Wolovick, 4/29.2016

% This script converts the width- and time-averaged accumulation data from
% the flowbands into winter accumulation and summer melt as a function of
% elevation.  

% The script takes a linear fit to accumulation as a function of elevation
% at high elevation, then uses that fit to extrapolate accum to lower
% elevation.  The different between the actual accum and the extrapolated
% accum represents annually averaged melt.  

% The slope of the fit is used to define an exponential falloff for
% accumulation at elevations higher than the present-day maximum.

% Note: interpolation from the course climate model grid onto the fine ice
% grid (before flowband contruction) resulted in an unphysical dropoff of
% melt near the terminus at low elevation.  This is replaced by a linear
% extrapolation in the output melt profile.

% This method is only applicable to the Greenland flowbands.

clear all
tic

%% Parameters:

% File names and paths:
ablationinputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/ComplexInversions/';
ablationinputfiles={'HelheimFlowband_v5_complexinversion.mat';...
    'HelheimFlowband_v6_complexinversion.mat';...
    'HelheimFlowband_v7_complexinversion.mat';...
    'KangerFlowband_v4_complexinversion.mat';...
    'KangerFlowband_v5_complexinversion.mat';...
    'KangerFlowband_v6_complexinversion.mat';...
    'JakobshavnFlowband_v8_complexinversion.mat';...
    'JakobshavnFlowband_v9_complexinversion.mat';...
    'JakobshavnFlowband_v10_complexinversion.mat';...
    'PetermannFlowband_v8_complexinversion.mat';...
    'PetermannFlowband_v9_complexinversion.mat';...
    'PetermannFlowband_v10_complexinversion.mat'};

% Elevation scale for extrapolation slope:
meltextrapelevrange=500;          % m

% Maximum elevation of inputs:
zmax=5e3;                         % m

% Grid size of elevation-dependent inputs:
zsize_input=100;                  % integer

%% Work:

% Loop through input files:
for thisfile=1:length(ablationinputfiles)
    
    % Load input file:
    clear SurfElev*
    load([ablationinputfolder,ablationinputfiles{thisfile}],'X_input','Accum_yr','SurfElev*','SurfTemp','steplength','lasticeind','xwavelength','tmelt')
    
    % Standardize naming of surface elevation:
    if exist('SurfElev_input','var')==0
        SurfElev_input=SurfElev;
    end
    
    % Find zero-crossing of accumulation:
    ind1=find(Accum_yr<0,1,'first');
    
    % Take derivative of accumulation:
    AccumDiff=[0,(Accum_yr(2:end)-Accum_yr(1:end-1))/steplength,0];
    
    % Find all local maxima of accumulation:
    LocalAccumMaxima=AccumDiff(1:end-1)>0&AccumDiff(2:end)<0;
    
    % Find last local maximum of accum before the zero-crossing:
    ind2=find(LocalAccumMaxima&X_input<X_input(ind1),1,'last');
    
    % Compute linear best-fit to high-elevation accum:
    A=[ones(ind2,1),SurfElev_input(1:ind2)'];
    accumextrapequation=(A'*A)\A'*Accum_yr(1:ind2)';
    
    % Compute accumulation extrapolation:
    AccumExtrap_yr=accumextrapequation(1)+accumextrapequation(2)*SurfElev_input;
    
    % Compute annually averaged ablation:
    MeltMean_yr=AccumExtrap_yr-Accum_yr;
    
    % Set high-elevation ablation to zero:
    MeltMean_yr(1:find(MeltMean_yr<0,1,'last'))=0;
    
    % Find minimum of accumulation:
    ind4=find(Accum_yr==min(Accum_yr(isnan(SurfElev_input)==0)),1,'first');
    
    % Find area to be used for the ablation extrapolation:
    ind3=find(SurfElev_input<SurfElev_input(ind4)+meltextrapelevrange,1,'first');
    
    % Compute linear bestfit to ablation against elevation:
    A=[ones(ind4-ind3+1,1),SurfElev_input(ind3:ind4)'];
    meltextrapequation=(A'*A)\A'*MeltMean_yr(ind3:ind4)';
    
    % Compute ablation extrapolation:
    MeltMean_yr(ind4:end)=meltextrapequation(1)+meltextrapequation(2)*SurfElev_input(ind4:end);
    
    % Low-pass filter ablation rate:
    MeltMean_yr(1:lasticeind+1)=intuitive_lowpass(MeltMean_yr(1:lasticeind+1),xwavelength/steplength);
    
    % Ensure melt is always a positive number:
    MeltMean_yr=max(MeltMean_yr,0);
    
    % Create an elevation grid:
    Z_input=linspace(0,zmax,zsize_input)';
    
    % Create accumulation input:
    Accum_yr_input=accumextrapequation(1)+accumextrapequation(2)*Z_input;
    
    % Use an exponential falloff above the maximum observed elevation:
    ind5=find(Z_input<=max(SurfElev_input),1,'last');
    efoldingelev=-Accum_yr_input(ind5)/accumextrapequation(2); % positive e-folding elevation
    Accum_yr_input(ind5+1:end)=Accum_yr_input(ind5)*exp(-(Z_input(ind5+1:end)-Z_input(ind5))/efoldingelev); % minus sign for decay upwards
    
    % Create ablation input:
    AnnualMelt_yr_input=interp1(SurfElev_input(1:lasticeind+1),MeltMean_yr(1:lasticeind+1),Z_input,'linear');
    AnnualMelt_yr_input(Z_input>max(SurfElev_input))=0;
    
    % Extrapolate melt to zero elevation:
    AnnualMelt_yr_input(Z_input<min(SurfElev_input(1:lasticeind+1)))=MeltMean_yr(lasticeind+1)...
        +meltextrapequation(2)*(Z_input(Z_input<min(SurfElev_input(1:lasticeind+1)))-SurfElev_input(lasticeind+1));
    
    % Compute best-fit to surface temperature as a function of elevation:
    A=[ones(lasticeind+1,1),SurfElev_input(1:lasticeind+1)'];
    surftempequation=(A'*A)\A'*SurfTemp(1:lasticeind+1)';
    
    % Compute surface temperature input:
    SurfTemp_input=surftempequation(1)+surftempequation(2)*Z_input;
    
    % Make figure:
    figure(1)
    
    % Plot surface mass balance:
    subplot(1,2,1)
    hold off
    plot(Accum_yr_input,Z_input,'b')
    hold on
    plot(-AnnualMelt_yr_input,Z_input,'r')
    plot(Accum_yr_input-AnnualMelt_yr_input,Z_input,'g')
    plot(Accum_yr,SurfElev_input,'k')
    plot([0,0],[0,zmax],'--k')
    legend('Accumulation','Ablation','Net Annual SMB','Width-Averaged RACMO Data','location','NorthWest')
    xlabel('Surface Mass Balance (m/yr)')
    ylabel('Elevation (m)')
    title(['Surface Mass Balance Parameterization in ',ablationinputfiles{thisfile}(1:end-21)],'interpreter','none')
    
    % Plot surface temperature:
    subplot(1,2,2)
    hold off
    plot(SurfTemp_input,Z_input,'r')
    hold on
    plot(SurfTemp,SurfElev_input,'k')
    plot(tmelt*[1,1],[0,zmax],'--k')
    legend('Linear Fit','Width-Averaged RACMO Data','Melting Point','location','NorthEast')
    xlabel('Surface Temperature (K)')
    ylabel('Elevation (m)')
    set(gca,'YAxisLocation','right')
    title(['Surface Temperature Parameterization in ',ablationinputfiles{thisfile}(1:end-21)],'interpreter','none')
    
    % Pause to finish drawing:
    drawnow
    
    % Save figure:
    set(gcf,'PaperSize',[11,6])
    set(gcf,'PaperPosition',[0,0,11,6])
    figname=[ablationinputfolder,ablationinputfiles{thisfile}(1:end-21),'_climateparameterization.png'];
    print('-dpng',figname,'-r300')
    
    % Save results:
    save([ablationinputfolder,ablationinputfiles{thisfile}],'*_input','-append')
    
end


% Final display:
disp('Done!')
toc