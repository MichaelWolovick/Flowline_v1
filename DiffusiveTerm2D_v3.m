% DiffusiveTerm2D_v1

% Mike Wolovick, 10/22/2011

% This script calculates the diffusive temperature term for FTM_2D.  It
% does this by integrating the conductive heat flux through all grid edges
% and dividing by the grid cell volume.

% v2:  modified to have variable vertical grid spacing

% v3: Modified for Flowline_v1.  No more CosTheta!  Also can do variable
% width, and fixed distances at the edges (had been doing full dx instead 
% of half).


% Find Conductive Fluxes In through Each Face:
CondFluxIn_l=-cond_i*[2*(Temp_c(:,1)-Temp_l),Temp_c(:,2:end)-Temp_c(:,1:end-1)]/dx;
CondFluxIn_r=cond_i*[Temp_c(:,2:end)-Temp_c(:,1:end-1),2*(Temp_r-Temp_c(:,end))]/dx;
CondFluxIn_u=cond_i*[Temp_c(2:end,:)-Temp_c(1:end-1,:);SurfTemp_u-Temp_c(end,:)]./DZ_ud(2:end,:);
CondFluxIn_d=-cond_i*[Temp_c(1,:)-Temp_d;Temp_c(2:end,:)-Temp_c(1:end-1,:)]./DZ_ud(1:end-1,:);

% Integrate for Energy Into Each Cell:
EnergyIn_l=DZ_lr_center(:,1:end-1).*repmat(Width_lr(1:end-1),[zsize,1]).*CondFluxIn_l;
EnergyIn_r=DZ_lr_center(:,2:end).*repmat(Width_lr(2:end),[zsize,1]).*CondFluxIn_r;
EnergyIn_u=dx*repmat(Width_c,[zsize,1]).*CondFluxIn_u;
EnergyIn_d=dx*repmat(Width_c,[zsize,1]).*CondFluxIn_d;

% Define Diffusive Term:
DiffusiveTerm_c=(EnergyIn_l+EnergyIn_r+EnergyIn_u+EnergyIn_d)./(dx*repmat(Width_c,[zsize,1]).*DZ_c);

% Clear Extraneous Variables:
clear CondFluxIn* EnergyIn* 
