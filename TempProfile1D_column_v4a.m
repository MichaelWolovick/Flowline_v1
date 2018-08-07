function [Temp_c,basaltemp,meltrate]=TempProfile1D_column_v4a(T_l,SourceTerm_c,U_lr,W_ud,istied,meltrate,DZ_c,DZ_lr,dx,gflux,surftemp,meltpoint,MaterialParameters)

% TempProfile1D_column

% Mike Wolovick, 11/23/2011

% This function solves for a steady state temperature profile in a single
% vertical column, including vertical advection/diffusion and horizontal
% advection.

% This function does not parse its inputs.  It assumes you have called
% everything correctly.

% Changes on 12/22/2011:  the function is now compatable with variable
% vertical grid spacing.

% Changes for v2:  the function is now compatable with semi-lagrangian
% horizontal temperature advection.  Inputs must include a jump-off point
% temperature and a horizontal flux divergence.  The function also assumes
% dx=dy.

% Changes for v3:  the function is now compatable with variable material
% properties in FTM_3D_ss_v5.

% Changes for v4:  the function now expects vertical velocities to be
% defined on grid centers, and it does not require a flux divergence input.

% Changes for v4a:  the function now works with Flowline_v1 instead of
% FTM_2D.  The ability to do variable material parameters has been removed.
%  Velocity is now defined at grid edges, consistent with the rest of the
%  model.


%% Correspondance of inputs with variables in Flowline_v1:

% Assuming the script is run in the d2 column, flow is left->right:

% T_l=Temp_c_mid(:,d2-1);
% SourceTerm_c=StrainHeat_c(:,d2);
% U_lr=U_lr(:,d2:d2+1);
% W_ud=W_ud(:,d2); 
% DZ_c=DZ_c(:,d2);
% DZ_lr=DZ_lr_upwind(:,d2:d2+1);
% istied=IsTied(d2);
% meltrate=MeltRate_c(d2);
% dx=dx;
% gflux=gflux;
% surftemp=SurfTemp_c(d2);
% meltpoint=MeltPoint_d(d2);
% MaterialParameters=ThermalParameters;

%% Indexing key for EnergyBalanceMatrix:

% ind=(contributor-1)*n+equation
% where:
% ind=linear index of term in EnergyBalanceMatrix
% contributor= index (in sln vector) of contributor
% equation= index (in sln vector) of equation
% n= number of elements (=zsize+2)

% Equation/contributors:
% 1       = basal plane
% 2:n-1   = grid centers
% n       = surface plane

%% Preparation:

% Unpack MaterialParameters:
rho_i=MaterialParameters.rho_i;
specheat_i=MaterialParameters.specheat_i;
cond_i=MaterialParameters.cond_i;
latentheat=MaterialParameters.latentheat;

% Interpolate dz to top/bottom faces:
DZ_ud=[DZ_c(1);.5*(DZ_c(1:end-1)+DZ_c(2:end));DZ_c(end)];

% Pre-allocate memory:
zsize=length(T_l);
n=zsize+2;
EnergyBalanceMatrix=sparse([],[],[],n,n,3*n-2);
EnergyConstraint=zeros(n,1);

% Determine whether to do horizontal advection:
if isequal(isnan(T_l),zeros(size(T_l)))==1
    dohorzadv=1;
else
    dohorzadv=0;
end

% Evaluate net volume flux into grid cell:
if dohorzadv
    NormalizedVolumeIn=-(U_lr(:,2).*DZ_lr(:,2)-U_lr(:,1).*DZ_lr(:,1))./(DZ_c*dx)-(W_ud(2:end)-W_ud(1:end-1))./DZ_c;
else
    NormalizedVolumeIn=-(W_ud(2:end)-W_ud(1:end-1))./DZ_c;
end


%% Build LHS:

% Interior conduction:
TheseInds=([3:n-2]'-1)*n+[3:n-2]';
EnergyBalanceMatrix(TheseInds)=-(cond_i./(DZ_c(2:end-1).*DZ_ud(2:end-2))+cond_i./(DZ_c(2:end-1).*DZ_ud(3:end-1))); % contribution from T(z,x)
TheseInds=([4:n-1]'-1)*n+[3:n-2]';
EnergyBalanceMatrix(TheseInds)=cond_i./(DZ_c(2:end-1).*DZ_ud(3:end-1)); % contribution from T(z+1,x)
TheseInds=([2:n-3]'-1)*n+[3:n-2]';
EnergyBalanceMatrix(TheseInds)=cond_i./(DZ_c(2:end-1).*DZ_ud(2:end-2)); % contribution from T(z-1,x)

% Top full cell conduction:
thisind=(n-1-1)*n+n-1;
EnergyBalanceMatrix(thisind)=-(cond_i/(DZ_c(end)*DZ_ud(end-1))+2*cond_i/(DZ_c(end)*DZ_ud(end))); % contribution from T(z,x)
thisind=(n-1)*n+n-1;
EnergyBalanceMatrix(thisind)=2*cond_i/(DZ_c(end)*DZ_ud(end)); % contribution from T(z+1,x)
thisind=((n-2)-1)*n+n-1;
EnergyBalanceMatrix(thisind)=cond_i/(DZ_c(end)*DZ_ud(end-1)); % contribution from T(z-1,x)

% Bottom full cell conduction:
thisind=n+2;
EnergyBalanceMatrix(thisind)=-(2*cond_i/(DZ_c(1)*DZ_ud(1))+cond_i/(DZ_c(1)*DZ_ud(2))); % contribution from T(z,x)
thisind=2*n+2;
EnergyBalanceMatrix(thisind)=cond_i/(DZ_c(1)*DZ_ud(2)); % contribution from T(z+1,x)
if istied==0
    thisind=2;
    EnergyBalanceMatrix(thisind)=2*cond_i/(DZ_c(1)*DZ_ud(1)); % contribution from T(z-1,x)
end

% Bottom edge conduction:
if istied==0
    thisind=1;
    EnergyBalanceMatrix(thisind)=-2*cond_i/DZ_ud(1); % contribution from T(z,x)
else
    thisind=1;
    EnergyBalanceMatrix(thisind)=-rho_i*latentheat; % contribution from m(y,x)
end
thisind=n+1;
EnergyBalanceMatrix(thisind)=2*cond_i/DZ_ud(1); % contribution from T(z+1,x)

% Interior vertical advection:
TheseInds=([2:n-1]'-1)*n+[2:n-1]';
TheseW=W_ud(2:end);
EnergyBalanceMatrix(TheseInds(TheseW>0))=EnergyBalanceMatrix(TheseInds(TheseW>0))-(rho_i*specheat_i./DZ_c(TheseW>0)).*TheseW(TheseW>0);  % contribution from T(z,x), up-moving, top edge
TheseInds=([2:n-2]'-1)*n+[3:n-1]';
TheseW=W_ud(2:end-1);
TheseDZ=DZ_c(2:end);
EnergyBalanceMatrix(TheseInds(TheseW>0))=EnergyBalanceMatrix(TheseInds(TheseW>0))+(rho_i*specheat_i./TheseDZ(TheseW>0)).*TheseW(TheseW>0);  % contribution from T(z-1,x), up-moving, bottom edge
TheseInds=([2:n-1]'-1)*n+[2:n-1]';
TheseW=W_ud(1:end-1);
EnergyBalanceMatrix(TheseInds(TheseW<0))=EnergyBalanceMatrix(TheseInds(TheseW<0))+(rho_i*specheat_i./DZ_c(TheseW<0)).*TheseW(TheseW<0);  % contribution from T(z,x), down-moving, bottom edge
TheseInds=([3:n-1]'-1)*n+[2:n-2]';
TheseW=W_ud(2:end-1);
TheseDZ=DZ_c(1:end-1);
EnergyBalanceMatrix(TheseInds(TheseW<0))=EnergyBalanceMatrix(TheseInds(TheseW<0))-(rho_i*specheat_i./TheseDZ(TheseW<0)).*TheseW(TheseW<0);  % contribution from T(z+1,x), down-moving, top edge

% Top vertical advection:  (contribution from Tsurf, down-moving, top edge
% of top cell)
if W_ud(end)<0
    thisind=(n-1)*n+n-1;
    EnergyBalanceMatrix(thisind)=EnergyBalanceMatrix(thisind)-(rho_i*specheat_i/DZ_c(end))*W_ud(end); % contribution from T(z+1,x), down-moving
end

% Bottom vertical advection: (contribution from Tbed, up-moving, bottom
% edge of bottom cell)
if W_ud(1)>0 && istied==0
    thisind=2;
    EnergyBalanceMatrix(thisind)=EnergyBalanceMatrix(thisind)+(rho_i*specheat_i/DZ_c(1))*W_ud(1); % contribution from T(z-1,x), up-moving
end

% Interior horizontal advection:
if dohorzadv
    TheseInds=([2:n-1]'-1)*n+[2:n-1]';
    EnergyBalanceMatrix(TheseInds)=EnergyBalanceMatrix(TheseInds)-rho_i*specheat_i.*U_lr(:,2).*DZ_lr(:,2)./(DZ_c*dx);   % contribution from T(z,x)
end

% Correction term for volume changes:
TheseInds=([2:n-1]'-1)*n+[2:n-1]';
EnergyBalanceMatrix(TheseInds)=EnergyBalanceMatrix(TheseInds)-rho_i*specheat_i*NormalizedVolumeIn;

% Surface Dirichlet "equation":
thisind=n^2;
EnergyBalanceMatrix(thisind)=1;

%% Build RHS:

% Top Dirichlet constraint:
thisind=n;
EnergyConstraint(thisind)=surftemp; % hard constraint

% Bottom full cell conduction:
if istied==1
    thisind=2;
    EnergyConstraint(thisind)=-2*meltpoint*cond_i/(DZ_c(1)*DZ_ud(1)); % contribution from T(pmp) on RHS
end

% Bottom edge conduction:
if istied==1
    thisind=1;
    EnergyConstraint(thisind)=2*meltpoint*cond_i/DZ_ud(1); % contribution from T(pmp) on RHS
else
    thisind=1;
    EnergyConstraint(thisind)=rho_i*latentheat*meltrate; % contribution from m(y,x) on RHS
end
EnergyConstraint(thisind)=EnergyConstraint(thisind)-gflux; % contribution from G(y,x) on RHS

% Bottom vertical advection: (upward-moving, contribution from Tbed to the
% bottom edge of the bottom cell)
if istied==1 && W_ud(1)>0
    thisind=2;
    EnergyConstraint(thisind)=EnergyConstraint(thisind)-(rho_i*specheat_i/DZ_c(1))*meltpoint*W_ud(1); % contribution from T(pmp) on RHS
end

% Interior horizontal advection:
if dohorzadv
    TheseInds=[2:n-1]';
    EnergyConstraint(TheseInds)=EnergyConstraint(TheseInds)-rho_i*specheat_i*U_lr(:,1).*T_l.*DZ_lr(:,1)./(DZ_c*dx);  % contribution from T(z,x-1) on RHS
end

% Strain Heating:
TheseInds=[2:n-1]';
EnergyConstraint(TheseInds)=EnergyConstraint(TheseInds)-SourceTerm_c;

%% Solve Equations and Parse Solution Vector:

% Solve Ax=b:
Temp=EnergyBalanceMatrix\EnergyConstraint;

% Separate the melt rate and basal temperature:
Temp_c=Temp(2:end-1);
if istied==1
    meltrate=Temp(1);
    basaltemp=meltpoint;
else
    basaltemp=Temp(1);
end


