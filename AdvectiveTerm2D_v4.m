% AdvectiveTerm2D_v1

% Mike Wolovick, 10/22/2011

% This script calculates the advective temperature term for FTM_2D.  It
% does so by integrating the snesible heat flux through each grid edge and
% dividing by the grid cell volume.  Temperature of the inflowing material
% on each edge is assigned in an upstream sense.  

% This script also calculates the advection of accreted ice, which can be
% later compared with the propagation of lagrangian tracers.

% v2:  the script is modified to have variable vertical grid spacing

% v3:  the advection of accreted ice cannot be disabled.  

% v4:  Flowline_v1.  accretion advection removed.  Variable width added.


%% Temperature Advection:

% Assign Temperature to Grid Edges:
% left/right edges:
Temp_lr_upwind=zeros(zsize,xsize+1);
CandidateTemp=[Temp_l,Temp_c];
Temp_lr_upwind(U_lr-Ugrid_lr>=0)=CandidateTemp(U_lr-Ugrid_lr>=0);
CandidateTemp=[Temp_c,Temp_r];
Temp_lr_upwind(U_lr-Ugrid_lr<0)=CandidateTemp(U_lr-Ugrid_lr<0);
% up/down edges:
Temp_ud_upwind=zeros(zsize+1,xsize);
CandidateTemp=[Temp_d;Temp_c];
Temp_ud_upwind(W_ud-Wgrid_ud>=0)=CandidateTemp(W_ud-Wgrid_ud>=0);
CandidateTemp=[Temp_c;SurfTemp_u];
Temp_ud_upwind(W_ud-Wgrid_ud<0)=CandidateTemp(W_ud-Wgrid_ud<0);

% Find Sensible Fluxes in Through Each Face:
SensibleFluxIn_l=rho_i*specheat_i*(U_lr(:,1:end-1)-Ugrid_lr(:,1:end-1)).*Temp_lr_upwind(:,1:end-1);
SensibleFluxIn_r=-rho_i*specheat_i*(U_lr(:,2:end)-Ugrid_lr(:,2:end)).*Temp_lr_upwind(:,2:end);
SensibleFluxIn_d=rho_i*specheat_i*(W_ud(1:end-1,:)-Wgrid_ud(1:end-1,:)).*Temp_ud_upwind(1:end-1,:);
SensibleFluxIn_u=-rho_i*specheat_i*(W_ud(2:end,:)-Wgrid_ud(2:end,:)).*Temp_ud_upwind(2:end,:);

% Integrate for Energy Into Each Cell:
EnergyIn_l=repmat(Width_lr(1:end-1),[zsize,1]).*DZ_lr_upwind(:,1:end-1).*SensibleFluxIn_l;
EnergyIn_r=repmat(Width_lr(2:end),[zsize,1]).*DZ_lr_upwind(:,2:end).*SensibleFluxIn_r;
EnergyIn_d=dx*repmat(Width_c,[zsize,1]).*SensibleFluxIn_d;
EnergyIn_u=dx*repmat(Width_c,[zsize,1]).*SensibleFluxIn_u;

% Compute correction for grid cell volume change:
NormalizedVolumeChangeRate_c=(Ugrid_lr(:,2:end).*DZ_lr_upwind(:,2:end).*repmat(Width_lr(2:end),[zsize,1])-Ugrid_lr(:,1:end-1).*DZ_lr_upwind(:,1:end-1).*repmat(Width_lr(1:end-1),[zsize,1]))./(DZ_c.*repmat(Width_c,[zsize,1])*dx)-...
    (Wgrid_ud(2:end,:)-Wgrid_ud(1:end-1,:))./DZ_c;
Correction_c=rho_i*specheat_i*Temp_c.*NormalizedVolumeChangeRate_c;

% Define Advective Term:
AdvectiveTerm_c=(EnergyIn_l+EnergyIn_r+EnergyIn_u+EnergyIn_d)./(dx*repmat(Width_c,[zsize,1]).*DZ_c)+Correction_c;


%% Clear Extraneous Variables:

clear SensibleFluxIn* EnergyIn* Correction_c 