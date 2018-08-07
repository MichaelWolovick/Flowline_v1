function [Depth,Temp,meltrate,gflux]=TempProfile1D_simple_v2(BoundaryConditions,MaterialParameters,mode,gridsize)

% TempProfile1D_simple_v2

% Mike Wolovick, 2/4/2016

% Based on TempProfile1D_column_v4a

% This function solves for a steady state temperature profile in a single
% vertical column, including vertical advection/diffusion and latent heat
% of melting/freezing at the basal plane.

% The vertical velocity profile assumes a constant strain rate between the
% melting rate at the bed and the surface accumulation rate.

% The function has 3 basic modes.  These modes describe how the goethermal
% flux and melt rate inputs/outputs behave.  These modes also determine
% whether basal temperature is constrained to equal the melting point.

% 1.  Untied: basal temperature free to be anything, above or below the
%     melting point.  Input geothermal flux is used and input melt rate is
%     ignored.  Output melt rate is zero and output geothermal flux is
%     equal to the input.

% 2.  Tied, melt rate unknown:  basal temperature fixed to the melting
%     point.  Input geothermal flux is used and input melt rate is ignored.
%     Output melt rate is computed based on basal heat flow and output
%     geothermal flux is equal to the input.

% 3.  Tied, geothermal flux unknown:  basal temperature fixed to the
%     melting point.  Input melt rate is used and input geothermal flux is
%     ignored.  Output geothermal flux is computed based on basal heat flow
%     and output melt rate is equal to the input.

% Note that in modes 1 and 2, vertical velocity is set to zero at the bed.


% Inputs:

% 1.  BoundaryConditions.  A structure with the following fields:
%        surftemp     surface temperature                        deg C
%        accum        surface accumulation rate                  m/yr
%        icethick     ice thickness                              m
%        gflux        geothermal flux                            W/m^2
%        meltrate     meltrate                                   m/yr
%        g            gravitational acceleration                 m/s^2

% 2.  MaterialParameters.  A structure with the following fields:
%        cond_i       thermal conductivity of ice                W/(m*K)
%        rho_i        density of ice                             kg/m^3
%        specheat_i   specific heat of ice                       J/(kg*K)
%        latentheat   latent heat of fusion                      J/kg
%        beta         pressure coefficient of melting point (>0) K/Pa

% 3.  mode.  A scalar equal to 1, 2, or 3.  Meaning described above.

% 4.  gridsize.  A scalar indicating the number of points in the output.


% Outputs:
% 1.  Depth          Depth, in meters  (size=[gridsize,1])
% 2.  Temp           Temperature, in degrees C  (size=[gridsize,1])
% 3.  meltrate       Melt rate, in m/yr  (size=[1,1])
% 4.  gflux          geothermal flux, in W/m^2  (size=[1,1])


% Units note:  inputs and outputs use units of years (ie, m/yr for
% velocities).  Internal operations within this function are in units of
% seconds.  The function automatically converts between them.

% Gridding note:  The output grid is interpolated onto grid top/bottom
% edges relative to the internal grid.  This is so that the output can be
% on an evenly spaced depth grid including the bed and the surface, while
% the internal operations can conserve energy within whole grid cells.  The
% internal operations follow the grid center/edges conventions of my other 
% modeling scripts.  Suffixes: "_c" indicates grid cell centers, "_ud" 
% indicates grid cell up/down edges.  In addition, the output is flipped
% with respect to the internal operations (internal operations index up
% from the bed, output depth and temperature vectors index down from the
% surface).


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

%% Check Inputs:

% Check BoundaryConditions:
if isfield(BoundaryConditions,'surftemp')==0
    error('Input structure "BoundaryConditions" must contain a parameter called "surftemp".')
end
if isfield(BoundaryConditions,'accum')==0
    error('Input structure "BoundaryConditions" must contain a parameter called "accum".')
end
if isfield(BoundaryConditions,'icethick')==0
    error('Input structure "BoundaryConditions" must contain a parameter called "icethick".')
end
if isfield(BoundaryConditions,'gflux')==0
    error('Input structure "BoundaryConditions" must contain a parameter called "gflux".')
end
if isfield(BoundaryConditions,'meltrate')==0
    error('Input structure "BoundaryConditions" must contain a parameter called "meltrate".')
end
if isfield(BoundaryConditions,'g')==0
    error('Input structure "BoundaryConditions" must contain a parameter called "g".')
end

% Check MaterialParameters:
if isfield(MaterialParameters,'rho_i')==0
    error('Input structure "MaterialParameters" must contain a parameter called "rho_i".')
end
if isfield(MaterialParameters,'specheat_i')==0
    error('Input structure "MaterialParameters" must contain a parameter called "specheat_i".')
end
if isfield(MaterialParameters,'cond_i')==0
    error('Input structure "MaterialParameters" must contain a parameter called "cond_i".')
end
if isfield(MaterialParameters,'latentheat')==0
    error('Input structure "MaterialParameters" must contain a parameter called "latentheat".')
end
if isfield(MaterialParameters,'beta')==0
    error('Input structure "MaterialParameters" must contain a parameter called "beta".')
end

% Check mode:
if mode~=1 && mode~=2 && mode~=3
    error('Input "mode" must be a scalar equal to 1, 2, or 3.')
end

% Check grid size:
if gridsize<4
    error('Input "gridsize" must be a scalar greater than or equal to 4.')
end

%% Preparation:

% Compute the number of seconds in a year:
secondsperyear=60*60*24*365.25;

% Unpack MaterialParameters:
rho_i=MaterialParameters.rho_i;
specheat_i=MaterialParameters.specheat_i;
cond_i=MaterialParameters.cond_i;
latentheat=MaterialParameters.latentheat;
beta=MaterialParameters.beta;

% Unpack BoundaryConditions:
surftemp=BoundaryConditions.surftemp;
accum=BoundaryConditions.accum/secondsperyear;
icethick=BoundaryConditions.icethick;
gflux=BoundaryConditions.gflux;
meltrate=BoundaryConditions.meltrate/secondsperyear;
g=BoundaryConditions.g;

% Compute meltpoint at the bed:
meltpoint=-rho_i*g*icethick*beta;

% Compute number of grid centers:
zsize=gridsize-1;

% Compute dz:
dz=icethick/zsize;

% Compute vertical velocity on grid edges:
if mode==3
    W_ud=-(meltrate+(accum-meltrate)*linspace(0,1,zsize+1)');
else
    W_ud=-accum*linspace(0,1,zsize+1)';
end

% Compute depth vector:
Depth=linspace(0,icethick,gridsize);

% Evaluate net volume flux into grid cell:
NormalizedVolumeIn=-(W_ud(2:end)-W_ud(1:end-1))/dz;

% Pre-allocate memory:
n=zsize+2;
EnergyBalanceMatrix=sparse([],[],[],n,n,3*n-2);
EnergyConstraint=zeros(n,1);


%% Build LHS:

% Interior conduction:
TheseInds=([3:n-2]'-1)*n+[3:n-2]';
EnergyBalanceMatrix(TheseInds)=-2*cond_i/(dz^2); % contribution from T(z,x)
TheseInds=([4:n-1]'-1)*n+[3:n-2]';
EnergyBalanceMatrix(TheseInds)=cond_i/(dz^2); % contribution from T(z+1,x)
TheseInds=([2:n-3]'-1)*n+[3:n-2]';
EnergyBalanceMatrix(TheseInds)=cond_i/(dz^2); % contribution from T(z-1,x)

% Top full cell conduction:
thisind=(n-1-1)*n+n-1;
EnergyBalanceMatrix(thisind)=-3*cond_i/(dz^2); % contribution from T(z,x)
thisind=(n-1)*n+n-1;
EnergyBalanceMatrix(thisind)=2*cond_i/(dz^2); % contribution from T(z+1,x)
thisind=((n-2)-1)*n+n-1;
EnergyBalanceMatrix(thisind)=cond_i/(dz^2); % contribution from T(z-1,x)

% Bottom full cell conduction:
thisind=n+2;
EnergyBalanceMatrix(thisind)=-3*cond_i/(dz^2); % contribution from T(z,x)
thisind=2*n+2;
EnergyBalanceMatrix(thisind)=cond_i/(dz^2); % contribution from T(z+1,x)
if mode==1
    thisind=2;
    EnergyBalanceMatrix(thisind)=2*cond_i/(dz^2); % contribution from T(z-1,x)
end

% Bottom edge conduction:
if mode==1
    thisind=1;
    EnergyBalanceMatrix(thisind)=-2*cond_i/dz; % contribution from basal temp
elseif mode==2
    thisind=1;
    EnergyBalanceMatrix(thisind)=-rho_i*latentheat; % contribution from melt rate
else
    thisind=1;
    EnergyBalanceMatrix(thisind)=1; % contribution from geothermal flux
end
thisind=n+1;
EnergyBalanceMatrix(thisind)=2*cond_i/dz; % contribution from T(z+1,x)

% Interior vertical advection:
TheseInds=([2:n-1]'-1)*n+[2:n-1]';
TheseW=W_ud(2:end);
EnergyBalanceMatrix(TheseInds(TheseW>0))=EnergyBalanceMatrix(TheseInds(TheseW>0))-(rho_i*specheat_i/dz)*TheseW(TheseW>0);  % contribution from T(z,x), up-moving, top edge
TheseInds=([2:n-2]'-1)*n+[3:n-1]';
TheseW=W_ud(2:end-1);
EnergyBalanceMatrix(TheseInds(TheseW>0))=EnergyBalanceMatrix(TheseInds(TheseW>0))+(rho_i*specheat_i/dz)*TheseW(TheseW>0);  % contribution from T(z-1,x), up-moving, bottom edge
TheseInds=([2:n-1]'-1)*n+[2:n-1]';
TheseW=W_ud(1:end-1);
EnergyBalanceMatrix(TheseInds(TheseW<0))=EnergyBalanceMatrix(TheseInds(TheseW<0))+(rho_i*specheat_i/dz)*TheseW(TheseW<0);  % contribution from T(z,x), down-moving, bottom edge
TheseInds=([3:n-1]'-1)*n+[2:n-2]';
TheseW=W_ud(2:end-1);
EnergyBalanceMatrix(TheseInds(TheseW<0))=EnergyBalanceMatrix(TheseInds(TheseW<0))-(rho_i*specheat_i/dz)*TheseW(TheseW<0);  % contribution from T(z+1,x), down-moving, top edge

% Top vertical advection:  (contribution from Tsurf, down-moving, top edge
% of top cell)
if W_ud(end)<0
    thisind=(n-1)*n+n-1;
    EnergyBalanceMatrix(thisind)=EnergyBalanceMatrix(thisind)-(rho_i*specheat_i/dz)*W_ud(end); % contribution from T(z+1,x), down-moving
end

% Bottom vertical advection: (contribution from Tbed, up-moving, bottom
% edge of bottom cell)
if W_ud(1)>0 && mode==1
    thisind=2;
    EnergyBalanceMatrix(thisind)=EnergyBalanceMatrix(thisind)+(rho_i*specheat_i/dz)*W_ud(1); % contribution from T(z-1,x), up-moving
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
if mode~=1
    thisind=2;
    EnergyConstraint(thisind)=-2*meltpoint*cond_i/(dz^2); % contribution from T(pmp) on RHS
end

% Bottom edge conduction:
thisind=1;
if mode==1
    EnergyConstraint(thisind)=-gflux; % contribution from gflux on RHS
elseif mode==2
    EnergyConstraint(thisind)=2*meltpoint*cond_i/dz-gflux; % contribution from T(pmp) and gflux on RHS
else
    EnergyConstraint(thisind)=2*meltpoint*cond_i/dz+rho_i*latentheat*meltrate; % contribution from T(pmp) and melt rate on RHS
end

% Bottom vertical advection: (upward-moving, contribution from Tbed to the
% bottom edge of the bottom cell)
if mode~=1 && W_ud(1)>0
    thisind=2;
    EnergyConstraint(thisind)=EnergyConstraint(thisind)-(rho_i*specheat_i/dz)*meltpoint*W_ud(1); % contribution from T(pmp) on RHS
end

%% Solve Equations and Parse Solution Vector:

% Solve Ax=b:
SolutionVector=EnergyBalanceMatrix\EnergyConstraint;

% Separate the outputs:
Temp_c=SolutionVector(2:end-1);
if mode==1
    meltrate=0;
    basaltemp=SolutionVector(1);
elseif mode==2
    meltrate=SolutionVector(1)*secondsperyear;
    basaltemp=meltpoint; 
else
    meltrate=meltrate*secondsperyear;
    gflux=SolutionVector(1);
    basaltemp=meltpoint; 
end

% Interpolate temperature to grid edges and flip vertically:
Temp=flipud([basaltemp;.5*(Temp_c(1:end-1)+Temp_c(2:end));surftemp]);


