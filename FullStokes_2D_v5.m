% FullStokes_2D_v4

% Mike Wolovick, 6/2/2012

% This script solves the Stokes equations on a skewed grid in 2D.

% The solution method is to linearize the problem for a given viscosity
% field and iterate for the nonlinear viscosity field.  The linear problem
% is solved with either Matlab's direct sparse matrix solver "\" or one of 
% the built-in iterative solvers, and the nonlinear iteration method is 
% damped fixed-point iteration.  The nonlinear iteration also solves for a 
% basal drag coefficient and a side drag coefficient.  The damping is
% performed on a logarithmic scale.

% v4:  runs inside Flowline_v1 instead of FTM_2D.  Basal drag coefficient
% is multiplied by grounded fraction.

% v5: includes side drag!  The side drag uses a linearized drag coefficient
% just like basal drag and effective viscosity.  A variable width has also
% been added to the conservation equations.  For now, I assume that width
% is constant vertically.  Therefore, width cancels for the terms related
% to stresses and fluxes on the top/bottom edges of grid cells.  In
% addition, I do the proper full stokes treatment of the ice shelf bottom
% (ie, I dont assume hydrostatic geometry).

% The majority of the script is spent setting up the terms in the linear
% system of equations.

% Boundary Conditions:

% Top BC:     stress_xz=0
%             stress_zz-P=0

% Bottom BC:  power-law sliding 
%             w=-meltrate(x)  (grounded)
%             normal stress = water pressure (ungrounded)
%   Note that the sliding law does not see the contribution to shear from 
%   dw/dx (angular momentum is not conserved at the bed)
%   The drag coefficient in the sliding law is set to zero in ungrounded
%   portions.
%   The ice bottom evolves dynamically in floating regions.  ie, we don't
%   make the assumption of hydrostatic equilibrium in full stokes.

% Side BC, 3 Options:
%   1.  flux condition:  u=u(z)
%   2.  fixed icethick:  w=w(z)
%   3.  moving front:    stress_xx-P=stress_xx(z) (only on right side)
%   stress_xz=0 for all options

% Viscosity BC: no gradient  (required for interpolation to grid cell corners)

% Dynamic pressure in the model is an anomaly relative to hydrostatic.  The
% body force is the driving stress gradient (ie, the static pressure 
% gradient plus the component of gravity acting along the skewed grid).

% Vertical velocity due to englacial melting is computed as a correction
% after the linear solve is performed.

% This script cannot function outside of the context of the variables
% created by FTM_2D_v2

% Description of the linear problem: (changes depending on side BC)
% Number of unknowns:            3*zsize*xsize+[-1,0]*zsize
% Number of equations:           3*zsize*xsize+[-1,0]*zsize
% Max number of terms possible:  26*zsize*xsize-5*zsize-8*xsize

% INDEXING KEYS:

% The vectors of equations and unknowns change their indexing depending on
% the side boundary conditions.  As good practice, I have the script 
% generate indexing maps at the begining of the model run.  The indexing
% maps are the same size as the grids they represent (ex: U_lr).  Elements
% that are not in play for a given set of side BC are set to NaN.  The 
% indexing maps are used in all subsequent commands to build the 
% components of the matrix that represents the linearized system of
% equations.  

% StokesMatrix indexing for a single term in the linear problem:

% MatrixInd=(VariableInd-1)*numvars+EquationInd

% where "VariableInd" and "EquationInd" are the indexing maps for the
% variable and equation under consideration

% example:

% StokesMatrix((UInds_lr(d1,d2)-1)*numvars+HorzForceInds_lr(d1,d2))= <insert linear coefficient here> ;

%% Produce Indexing Keys:

% MODIFY THIS SECTION AT YOUR OWN RISK!

% Only do this once:
if exist('Uinds','var')==0
    
    % Pre-allocate indexing maps:
    UInds_lr=NaN*zeros(zsize,xsize+1);
    WInds_ud=NaN*zeros(zsize+1,xsize);
    PInds_c=NaN*zeros(zsize,xsize);
    HorzForceInds_lr=NaN*zeros(zsize,xsize+1);
    VertForceInds_ud=NaN*zeros(zsize+1,xsize);
    MassBalInds_c=NaN*zeros(zsize,xsize);
    FreeSurfInds_u=NaN*zeros(1,xsize);
    if strcmp(rightbctype,'front') 
        SideStressInds_r=NaN*zeros(zsize,1);
    end
    
    % Loop through grid cells:
    lasteqnind=0;
    lastvarind=0;
    for ii=1:zsize*xsize
        
        % Convert to subscript index:
        [d1,d2]=ind2sub([zsize,xsize],ii);
        
        % Assign horizontal force indexing map (interior cells):
        if d2~=xsize
            HorzForceInds_lr(d1,d2+1)=lasteqnind+1;
            lasteqnind=lasteqnind+1;
        end
        
        % Assign vertical force indexing map (interior cells):
        if d1~=zsize
            VertForceInds_ud(d1+1,d2)=lasteqnind+1;
            lasteqnind=lasteqnind+1;
        end
        
        % Assign mass balance indexing map (all cells):
        MassBalInds_c(d1,d2)=lasteqnind+1;
        lasteqnind=lasteqnind+1;
        
        % Assign top free surface indexing map:
        if d1==zsize
            FreeSurfInds_u(d2)=lasteqnind+1;
            lasteqnind=lasteqnind+1;
        end
        
        % Assign right side stress indexing map:
        if d2==xsize && strcmp(rightbctype,'front') 
            SideStressInds_r(d1)=lasteqnind+1;
            lasteqnind=lasteqnind+1;
        end
        
        % Assign horizontal velocity indexing map, internal and left edge:
        if d2~=1 || strcmp(leftbctype,'flux')==0
            UInds_lr(d1,d2)=lastvarind+1;
            lastvarind=lastvarind+1;
        end
        
        % Assign horizontal velocity indexing map, right edge:
        if d2==xsize && strcmp(rightbctype,'flux')==0
            UInds_lr(d1,d2+1)=lastvarind+1;
            lastvarind=lastvarind+1;
        end
        
        % Assign vertical velocity indexing map:
        if (d2~=1 || strcmp(leftbctype,'fixedicethick')==0) && (d2~=xsize || strcmp(rightbctype,'fixedicethick')==0)
            WInds_ud(d1+1,d2)=lastvarind+1;
            lastvarind=lastvarind+1;
        end
        
        % Assign pressure indexing map:
        PInds_c(d1,d2)=lastvarind+1;
        lastvarind=lastvarind+1;
        
    end
    
    % Check indexing maps:
    if lasteqnind~=lastvarind
        error('Number of variables did not match number of equations.')
    end
    
    % Compute number of variables:
    numvars=lastvarind;
    
end

%% Solve Stokes Equation:

% Create solution vector for the first time:
if timestep==1 && iteration_major==1 && strcmp(solvertype,'direct')==0
    SolutionVector=zeros(numvars,1);
end

% Evaulate rheological constant on grid centers:
A_c(Temp_c_mid>=t0)=a0*exp(-(q_big/r)*((1./Temp_c_mid(Temp_c_mid>=t0))-(1/t0)));
A_c(Temp_c_mid<t0)=a0*exp(-(q_small/r)*((1./Temp_c_mid(Temp_c_mid<t0))-(1/t0)));

% Interpolate rheological constant to grid edges: (no gradient BC)
A_lr=[A_c(:,1),.5*(A_c(:,1:end-1)+A_c(:,2:end)),A_c(:,end)];

% Nonlinear iteration:
converged_stokes=0;
iteration_stokes=1;
forcerepeat=0;
triedreset=0;
Viscosity_c_last=Viscosity_c;
DragCoefficient_lrd_last=DragCoefficient_lrd;
if dosidedrag
    DragCoefficient_lr_last=DragCoefficient_lr;
end
FlowDir_lr_last=FlowDir_lr;
while converged_stokes==0
    
    % Reset coefficients, indices, and constraint:
    StokesIndices=zeros(26*xsize*zsize-8*xsize-5*zsize,1);
    StokesCoefficients=zeros(26*xsize*zsize-8*xsize-5*zsize,1);
    StokesConstraint=zeros(numvars,1);
    lastterm=0;
    
    % Linearly interpolate viscosity to grid corners:  (no gradient BC)
    Viscosity_lrud=[Viscosity_c(1,1),.5*(Viscosity_c(1,1:end-1)+Viscosity_c(1,2:end)),Viscosity_c(1,end);... % bot row all corners
        Viscosity_c(1:end-1,1).*(1-.5*DZhat_c(1:end-1)./DZhat_ud(2:end-1))+Viscosity_c(2:end,1)*.5.*DZhat_c(1:end-1)./DZhat_ud(2:end-1),...  % left edge interior corners
        .5*(Viscosity_c(1:end-1,1:end-1).*repmat(1-.5*DZhat_c(1:end-1)./DZhat_ud(2:end-1),[1,xsize-1])+Viscosity_c(2:end,1:end-1).*repmat(.5*DZhat_c(1:end-1)./DZhat_ud(2:end-1),[1,xsize-1])... % internal linear interpolation, part 1
        +Viscosity_c(1:end-1,2:end).*repmat(1-.5*DZhat_c(1:end-1)./DZhat_ud(2:end-1),[1,xsize-1])+Viscosity_c(2:end,2:end).*repmat(.5*DZhat_c(1:end-1)./DZhat_ud(2:end-1),[1,xsize-1])),...  % internal linear interpolation, part 2
        Viscosity_c(1:end-1,end).*(1-.5*DZhat_c(1:end-1)./DZhat_ud(2:end-1))+Viscosity_c(2:end,end)*.5.*DZhat_c(1:end-1)./DZhat_ud(2:end-1);... % right edge interior corners
        Viscosity_c(end,1),.5*(Viscosity_c(end,1:end-1)+Viscosity_c(end,2:end)),Viscosity_c(end,end)]; % top row all corners
        
    % Horizontal Force Balance, internal cells:  (z,x+.5)
    % Contribution from u(z,x+.5):  (applies to left/right columns as well)
    StokesIndices(lastterm+1:lastterm+(zsize-2)*(xsize-1))=reshape((UInds_lr(2:end-1,2:end-1)-1)*numvars+HorzForceInds_lr(2:end-1,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize-2)*(xsize-1))=reshape(-.5*Viscosity_lrud(3:end-1,2:end-1)./(DZ_lr_center(2:end-1,2:end-1).*DZ_lrud(3:end-1,2:end-1))-.5*Viscosity_lrud(2:end-2,2:end-1)./(DZ_lr_center(2:end-1,2:end-1).*DZ_lrud(2:end-2,2:end-1))...% shear stress on top/bot
        -(Viscosity_c(2:end-1,2:end).*DZ_c(2:end-1,2:end).*Width_c(2:end-1,2:end)+Viscosity_c(2:end-1,1:end-1).*DZ_c(2:end-1,1:end-1).*Width_c(2:end-1,1:end-1))./((dx^2)*DZ_lr_center(2:end-1,2:end-1).*Width_lr(2:end-1,2:end-1))...  % long stress on left/right
        -2*DragCoefficient_lr(2:end-1,2:end-1)./Width_lr(2:end-1,2:end-1),[],1); % side drag
    lastterm=lastterm+(zsize-2)*(xsize-1);
    % Contribution from u(z+1,x+.5):  (applies to left/right columns and bottom row as well)
    StokesIndices(lastterm+1:lastterm+(zsize-1)*(xsize-1))=reshape((UInds_lr(2:end,2:end-1)-1)*numvars+HorzForceInds_lr(1:end-1,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize-1)*(xsize-1))=reshape(.5*Viscosity_lrud(2:end-1,2:end-1)./(DZ_lr_center(1:end-1,2:end-1).*DZ_lrud(2:end-1,2:end-1)),[],1); % shear stress on top
    lastterm=lastterm+(zsize-1)*(xsize-1);
    % Contribution from u(z-1,x+.5):  (applies to left/right columns and top row as well)
    StokesIndices(lastterm+1:lastterm+(zsize-1)*(xsize-1))=reshape((UInds_lr(1:end-1,2:end-1)-1)*numvars+HorzForceInds_lr(2:end,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize-1)*(xsize-1))=reshape(.5*Viscosity_lrud(2:end-1,2:end-1)./(DZ_lr_center(2:end,2:end-1).*DZ_lrud(2:end-1,2:end-1)),[],1); % shear stress on bot
    lastterm=lastterm+(zsize-1)*(xsize-1);
    % Contribution from u(z,x+1.5): (applies to top/bot rows and left column as well)
    StokesIndices(lastterm+1:lastterm+(zsize)*(xsize-2))=reshape((UInds_lr(:,3:end-1)-1)*numvars+HorzForceInds_lr(:,2:end-2),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize)*(xsize-2))=reshape(Viscosity_c(:,2:end-1).*DZ_c(:,2:end-1).*Width_c(:,2:end-1)./((dx^2)*DZ_lr_center(:,2:end-2).*Width_lr(:,2:end-2)),[],1); % long stress on right
    lastterm=lastterm+(zsize)*(xsize-2);
    % Contribution from u(z,x-.5): (applies to top/bot rows and right column as well)
    StokesIndices(lastterm+1:lastterm+(zsize)*(xsize-2))=reshape((UInds_lr(:,2:end-2)-1)*numvars+HorzForceInds_lr(:,3:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize)*(xsize-2))=reshape(Viscosity_c(:,2:end-1).*DZ_c(:,2:end-1).*Width_c(:,2:end-1)./((dx^2)*DZ_lr_center(:,3:end-1).*Width_lr(:,3:end-1)),[],1); % long stress on left
    lastterm=lastterm+(zsize)*(xsize-2);
    % Contribution from w(z+.5,x+1):  (applies to bot row and left column as well)
    StokesIndices(lastterm+1:lastterm+(zsize-1)*(xsize-2))=reshape((WInds_ud(2:end-1,2:end-1)-1)*numvars+HorzForceInds_lr(1:end-1,2:end-2),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize-1)*(xsize-2))=reshape(.5*Viscosity_lrud(2:end-1,2:end-2)./(dx*DZ_lr_center(1:end-1,2:end-2)),[],1); % shear stress on top
    lastterm=lastterm+(zsize-1)*(xsize-2);
    % Contribution from w(z+.5,x):  (applies to bot row and right column as well)
    StokesIndices(lastterm+1:lastterm+(zsize-1)*(xsize-2))=reshape((WInds_ud(2:end-1,2:end-1)-1)*numvars+HorzForceInds_lr(1:end-1,3:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize-1)*(xsize-2))=reshape(-.5*Viscosity_lrud(2:end-1,3:end-1)./(dx*DZ_lr_center(1:end-1,3:end-1)),[],1); % shear stress on top
    lastterm=lastterm+(zsize-1)*(xsize-2);
    % Contribution from w(z-.5,x+1):  (applies to top row and left column as well)
    StokesIndices(lastterm+1:lastterm+(zsize-1)*(xsize-2))=reshape((WInds_ud(2:end-1,2:end-1)-1)*numvars+HorzForceInds_lr(2:end,2:end-2),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize-1)*(xsize-2))=reshape(-.5*Viscosity_lrud(2:end-1,2:end-2)./(dx*DZ_lr_center(2:end,2:end-2)),[],1); % shear stress on bot
    lastterm=lastterm+(zsize-1)*(xsize-2);
    % Contribution from w(z-.5,x):  (applies to top row and right column as well)
    StokesIndices(lastterm+1:lastterm+(zsize-1)*(xsize-2))=reshape((WInds_ud(2:end-1,2:end-1)-1)*numvars+HorzForceInds_lr(2:end,3:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize-1)*(xsize-2))=reshape(.5*Viscosity_lrud(2:end-1,3:end-1)./(dx*DZ_lr_center(2:end,3:end-1)),[],1); % shear stress on bot
    lastterm=lastterm+(zsize-1)*(xsize-2);
    % Contribution from P(z,x+1):  (applies everywhere)
    StokesIndices(lastterm+1:lastterm+(zsize)*(xsize-1))=reshape((PInds_c(:,2:end)-1)*numvars+HorzForceInds_lr(:,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize)*(xsize-1))=reshape(-DZ_c(:,2:end).*Width_c(:,2:end)./(dx*DZ_lr_center(:,2:end-1).*Width_lr(:,2:end-1)),[],1); % pressure on right
    lastterm=lastterm+(zsize)*(xsize-1);
    % Contribution from P(z,x):  (applies everywhere)
    StokesIndices(lastterm+1:lastterm+(zsize)*(xsize-1))=reshape((PInds_c(:,1:end-1)-1)*numvars+HorzForceInds_lr(:,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(zsize)*(xsize-1))=reshape(DZ_c(:,1:end-1).*Width_c(:,1:end-1)./(dx*DZ_lr_center(:,2:end-1).*Width_lr(:,2:end-1)),[],1); % pressure on left
    lastterm=lastterm+(zsize)*(xsize-1);
    
    % Horizontal Force Balance, top/bot boundaries:  (z,x+.5)
    % Contribution from u(z,x+.5):  (bottom row, excludes left/right corners)
    StokesIndices(lastterm+1:lastterm+(xsize-1))=reshape((UInds_lr(1,2:end-1)-1)*numvars+HorzForceInds_lr(1,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-1))=-.5*Viscosity_lrud(2,2:end-1)./(DZ_lr_center(1,2:end-1).*DZ_lrud(2,2:end-1))... % shear stress on top
        -DragCoefficient_lrd(2:end-1)./((2*DragCoefficient_lrd(2:end-1).*DZ_lrud(1,2:end-1)./Viscosity_lrud(1,2:end-1)+1).*DZ_lr_center(1,2:end-1))... % shear stress on bot
        -(Viscosity_c(1,2:end).*DZ_c(1,2:end).*Width_c(1,2:end)+Viscosity_c(1,1:end-1).*DZ_c(1,1:end-1).*Width_c(1,1:end-1))./((dx^2)*DZ_lr_center(1,2:end-1).*Width_lr(1,2:end-1))... % long stress on left/right
        -2*DragCoefficient_lr(1,2:end-1)./Width_lr(1,2:end-1); % side drag
    lastterm=lastterm+(xsize-1);
    % Contribution from u(z,x+.5):  (top row, excludes left/right corners)
    StokesIndices(lastterm+1:lastterm+(xsize-1))=reshape((UInds_lr(end,2:end-1)-1)*numvars+HorzForceInds_lr(end,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-1))=-.5*Viscosity_lrud(end-1,2:end-1)./(DZ_lr_center(end,2:end-1).*DZ_lrud(end-1,2:end-1))... % shear stress on bot
        -(Viscosity_c(end,2:end).*DZ_c(end,2:end).*Width_c(end,2:end)+Viscosity_c(end,1:end-1).*DZ_c(end,1:end-1).*Width_c(end,1:end-1))./((dx^2)*DZ_lr_center(end,2:end-1).*Width_lr(end,2:end-1))... % long stress on left/right
        -2*DragCoefficient_lr(end,2:end-1)./Width_lr(end,2:end-1); % side drag
    lastterm=lastterm+(xsize-1);
    
    % Horizontal Force Balance, left/right boundaries:  (z,x+.5)
    % (fills in terms neglected in the last internal cell for fixed u case)
    % Left Column:
    if strcmp(leftbctype,'flux')==1
        % Contribution from w(z+.5,x):  (applies to bot corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,1)-1)*numvars+HorzForceInds_lr(1:end-1,2);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=-.5*Viscosity_lrud(2:end-1,2)./(dx*DZ_lr_center(1:end-1,2)); % shear stress on top
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z-.5,x):  (applies to top corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,1)-1)*numvars+HorzForceInds_lr(2:end,2);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=.5*Viscosity_lrud(2:end-1,2)./(dx*DZ_lr_center(2:end,2)); % shear stress on bot
        lastterm=lastterm+(zsize-1);
    elseif strcmp(leftbctype,'fixedicethick')==1
        % Contribution from u(z,x-.5): (left column, applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+zsize)=(UInds_lr(:,1)-1)*numvars+HorzForceInds_lr(:,2);
        StokesCoefficients(lastterm+1:lastterm+zsize)=Viscosity_c(:,1).*DZ_c(:,1).*Width_c(:,1)./((dx^2)*DZ_lr_center(:,2).*Width_lr(:,2)); % long stress on left edge
        lastterm=lastterm+zsize;
    end
    % Right Column:
    if strcmp(rightbctype,'flux')==1
        % Contribution from w(z+.5,x+1):  (right column, applies to bot corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,end)-1)*numvars+HorzForceInds_lr(1:end-1,end-1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=.5*Viscosity_lrud(2:end-1,end-1)./(dx*DZ_lr_center(1:end-1,end-1)); % shear stress on top
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z-.5,x+1):  (right column, applies to top corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,end)-1)*numvars+HorzForceInds_lr(2:end,end-1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=-.5*Viscosity_lrud(2:end-1,end-1)./(dx*DZ_lr_center(2:end,end-1)); % shear stress on bot
        lastterm=lastterm+(zsize-1);
    elseif strcmp(rightbctype,'fixedicethick')==1
        % Contribution from u(z,x+1.5): (right column, applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+zsize)=(UInds_lr(:,end)-1)*numvars+HorzForceInds_lr(:,end-1);
        StokesCoefficients(lastterm+1:lastterm+zsize)=Viscosity_c(:,end).*DZ_c(:,end).*Width_c(:,end)./((dx^2)*DZ_lr_center(:,end-1).*Width_lr(:,end-1)); % long stress on right edge
        lastterm=lastterm+zsize;
    else % front
        % Contribution from u(z,x+1.5): (right column, applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+zsize)=(UInds_lr(:,end)-1)*numvars+HorzForceInds_lr(:,end-1);
        StokesCoefficients(lastterm+1:lastterm+zsize)=Viscosity_c(:,end).*DZ_c(:,end).*Width_c(:,end)./((dx^2)*DZ_lr_center(:,end-1).*Width_lr(:,end-1)); % long stress on right edge
        lastterm=lastterm+zsize;
        % Contribution from w(z+.5,x+1):  (right column, applies to bot corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,end)-1)*numvars+HorzForceInds_lr(1:end-1,end-1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=.5*Viscosity_lrud(2:end-1,end-1)./(dx*DZ_lr_center(1:end-1,end-1)); % shear stress on top
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z-.5,x+1):  (right column, applies to top corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,end)-1)*numvars+HorzForceInds_lr(2:end,end-1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=-.5*Viscosity_lrud(2:end-1,end-1)./(dx*DZ_lr_center(2:end,end-1)); % shear stress on bot
        lastterm=lastterm+(zsize-1);
    end
    
    % Horizontal Force Balance Constraints:
    % All cells:
    StokesConstraint(reshape(HorzForceInds_lr(:,2:end-1),[],1))...
        =reshape(-DrivingStressGrad_lr(:,2:end-1),[],1);
    % Left column:
    if strcmp(leftbctype,'flux')==1
        StokesConstraint(reshape(HorzForceInds_lr(:,2),[],1))=StokesConstraint(reshape(HorzForceInds_lr(:,2),[],1))...
            -Viscosity_c(:,1).*DZ_c(:,1).*Width_c(:,1).*U_l./((dx^2)*DZ_lr_center(:,2).*Width_lr(:,2)); % long stress on left edge
    elseif strcmp(leftbctype,'fixedicethick')==1
        % Left column, internal only:
        StokesConstraint(reshape(HorzForceInds_lr(2:end-1,2),[],1))=StokesConstraint(reshape(HorzForceInds_lr(2:end-1,2),[],1))...
            +.5*(Viscosity_lrud(3:end-1,2).*W_ud(3:end-1,1)-Viscosity_lrud(2:end-2,2).*W_ud(2:end-2,1))./(dx*DZ_lr_center(2:end-1,2)); % shear stress on top/bot edges
        % Top left corner:
        StokesConstraint(reshape(HorzForceInds_lr(end,2),[],1))=StokesConstraint(reshape(HorzForceInds_lr(end,2),[],1))...
            -.5*Viscosity_lrud(end-1,2).*W_ud(end-1,1)./(dx*DZ_lr_center(end,2));  % shear stress on bot edge of grid cell
        % Bottom left corner:
        StokesConstraint(reshape(HorzForceInds_lr(1,2),[],1))=StokesConstraint(reshape(HorzForceInds_lr(1,2),[],1))...
            +.5*Viscosity_lrud(2,2).*W_ud(2,1)./(dx*DZ_lr_center(1,2));  % shear stress on top edge of grid cell
    end
    % Right column:
    if strcmp(rightbctype,'flux')==1
        StokesConstraint(reshape(HorzForceInds_lr(:,end-1),[],1))=StokesConstraint(reshape(HorzForceInds_lr(:,end-1),[],1))...
            -Viscosity_c(:,end).*DZ_c(:,end).*Width_c(:,end).*U_r./((dx^2)*DZ_lr_center(:,end-1).*Width_lr(:,end-1)); % long stress on right edge
    elseif strcmp(rightbctype,'fixedicethick')==1
        % Right column, internal only:
        StokesConstraint(reshape(HorzForceInds_lr(2:end-1,end-1),[],1))=StokesConstraint(reshape(HorzForceInds_lr(2:end-1,end-1),[],1))...
            -.5*(Viscosity_lrud(3:end-1,end-1).*W_ud(3:end-1,end)-Viscosity_lrud(2:end-2,end-1).*W_ud(2:end-2,end))./(dx*DZ_lr_center(2:end-1,end-1));  % shear stress on top/bot edges
        % Top right corner:
        StokesConstraint(reshape(HorzForceInds_lr(end,end-1),[],1))=StokesConstraint(reshape(HorzForceInds_lr(end,end-1),[],1))...
            +.5*Viscosity_lrud(end-1,end-1).*W_ud(end-1,end)./(dx*DZ_lr_center(end,end-1)); % shear stress on bot edge of grid cell
        % Bottom right corner:
        StokesConstraint(reshape(HorzForceInds_lr(1,end-1),[],1))=StokesConstraint(reshape(HorzForceInds_lr(1,end-1),[],1))...
            -.5*Viscosity_lrud(2,end-1).*W_ud(2,end)./(dx*DZ_lr_center(1,end-1)); % shear stress on top edge of grid cell
    end
    
    % Vertical Force Balance, internal cells:  (z+.5,x)
    % Contribution from w(z+.5,x):  (applies to top/bot rows as well)
    StokesIndices(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape((WInds_ud(2:end-1,2:end-1)-1)*numvars+VertForceInds_ud(2:end-1,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape(-.5*(Viscosity_lrud(2:end-1,3:end-1).*DZ_lrud(2:end-1,3:end-1).*Width_lrud(2:end-1,3:end-1)+Viscosity_lrud(2:end-1,2:end-2).*DZ_lrud(2:end-1,2:end-2).*Width_lrud(2:end-1,2:end-2))./((dx^2)*DZ_ud(2:end-1,2:end-1).*Width_ud(2:end-1,2:end-1))... % shear stress on left/right edges
        -(Viscosity_c(2:end,2:end-1)./(DZ_c(2:end,2:end-1).*DZ_ud(2:end-1,2:end-1))+Viscosity_c(1:end-1,2:end-1)./(DZ_c(1:end-1,2:end-1).*DZ_ud(2:end-1,2:end-1))),[],1); % long stress on top/bot edges
    lastterm=lastterm+(xsize-2)*(zsize-1);
    % Contribution from w(z+.5,x+1):  (applies to top/bot rows and left column as well, does not apply to right two columns)
    StokesIndices(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape((WInds_ud(2:end-1,2:end-1)-1)*numvars+VertForceInds_ud(2:end-1,1:end-2),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape(.5*Viscosity_lrud(2:end-1,2:end-2).*DZ_lrud(2:end-1,2:end-2).*Width_lrud(2:end-1,2:end-2)./((dx^2)*DZ_ud(2:end-1,1:end-2).*Width_ud(2:end-1,1:end-2)),[],1); % shear stress on right edge
    lastterm=lastterm+(xsize-2)*(zsize-1);
    % Contribution from w(z+.5,x-1):  (applies to top/bot rows and right column as well, does not apply to left two columns)
    StokesIndices(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape((WInds_ud(2:end-1,2:end-1)-1)*numvars+VertForceInds_ud(2:end-1,3:end),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape(.5*Viscosity_lrud(2:end-1,3:end-1).*DZ_lrud(2:end-1,3:end-1).*Width_lrud(2:end-1,3:end-1)./((dx^2)*DZ_ud(2:end-1,3:end).*Width_ud(2:end-1,3:end)),[],1); % shear stress on left edge
    lastterm=lastterm+(xsize-2)*(zsize-1);
    % Contribution from w(z+1.5,x):  (applies to top/bot rows as well)
    StokesIndices(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape((WInds_ud(3:end,2:end-1)-1)*numvars+VertForceInds_ud(2:end-1,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape(Viscosity_c(2:end,2:end-1)./(DZ_c(2:end,2:end-1).*DZ_ud(2:end-1,2:end-1)),[],1); % long stress on top edge
    lastterm=lastterm+(xsize-2)*(zsize-1);
    % Contribution from w(z-.5,x):  (applies to top row as well)
    StokesIndices(lastterm+1:lastterm+(xsize-2)*(zsize-2))=reshape((WInds_ud(2:end-2,2:end-1)-1)*numvars+VertForceInds_ud(3:end-1,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2)*(zsize-2))=reshape(Viscosity_c(2:end-1,2:end-1)./(DZ_c(2:end-1,2:end-1).*DZ_ud(3:end-1,2:end-1)),[],1); % long stress on bot edge
    lastterm=lastterm+(xsize-2)*(zsize-2);
    % Contribution from u(z+1,x+.5):  (applies to top/bot rows and left column as well)
    StokesIndices(lastterm+1:lastterm+(xsize-1)*(zsize-1))=reshape((UInds_lr(2:end,2:end-1)-1)*numvars+VertForceInds_ud(2:end-1,1:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-1)*(zsize-1))=reshape(.5*Viscosity_lrud(2:end-1,2:end-1).*Width_lrud(2:end-1,2:end-1)./(dx*DZ_ud(2:end-1,1:end-1).*Width_ud(2:end-1,1:end-1)),[],1); % shear stress on right edge
    lastterm=lastterm+(xsize-1)*(zsize-1);
    % Contribution from u(z,x+.5):  (applies to top/bot rows and left column as well)
    StokesIndices(lastterm+1:lastterm+(xsize-1)*(zsize-1))=reshape((UInds_lr(1:end-1,2:end-1)-1)*numvars+VertForceInds_ud(2:end-1,1:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-1)*(zsize-1))=reshape(-.5*Viscosity_lrud(2:end-1,2:end-1).*Width_lrud(2:end-1,2:end-1)./(dx*DZ_ud(2:end-1,1:end-1).*Width_ud(2:end-1,1:end-1)),[],1); % shear stress on right edge
    lastterm=lastterm+(xsize-1)*(zsize-1);
    % Contribution from u(z+1,x-.5):  (applies to top/bot rows and right column as well)
    StokesIndices(lastterm+1:lastterm+(xsize-1)*(zsize-1))=reshape((UInds_lr(2:end,2:end-1)-1)*numvars+VertForceInds_ud(2:end-1,2:end),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-1)*(zsize-1))=reshape(-.5*Viscosity_lrud(2:end-1,2:end-1).*Width_lrud(2:end-1,2:end-1)./(dx*DZ_ud(2:end-1,2:end).*Width_ud(2:end-1,2:end)),[],1); % shear stress on left edge
    lastterm=lastterm+(xsize-1)*(zsize-1);
    % Contribution from u(z,x-.5):  (applies to top/bot rows and right column as well)
    StokesIndices(lastterm+1:lastterm+(xsize-1)*(zsize-1))=reshape((UInds_lr(1:end-1,2:end-1)-1)*numvars+VertForceInds_ud(2:end-1,2:end),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-1)*(zsize-1))=reshape(.5*Viscosity_lrud(2:end-1,2:end-1).*Width_lrud(2:end-1,2:end-1)./(dx*DZ_ud(2:end-1,2:end).*Width_ud(2:end-1,2:end)),[],1);  % shear stress on left edge
    lastterm=lastterm+(xsize-1)*(zsize-1);
    % Contribution from P(z+1,x):  (applies everywhere)
    StokesIndices(lastterm+1:lastterm+(xsize)*(zsize-1))=reshape((PInds_c(2:end,:)-1)*numvars+VertForceInds_ud(2:end-1,:),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize)*(zsize-1))=reshape(-1./DZ_ud(2:end-1,:),[],1); % pressure on top edge
    lastterm=lastterm+(xsize)*(zsize-1);
    % Contribution from P(z,x):  (applies everywhere)
    StokesIndices(lastterm+1:lastterm+(xsize)*(zsize-1))=reshape((PInds_c(1:end-1,:)-1)*numvars+VertForceInds_ud(2:end-1,:),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize)*(zsize-1))=reshape(1./DZ_ud(2:end-1,:),[],1); % pressure on bottom edge
    lastterm=lastterm+(xsize)*(zsize-1);
    
    % Vertical Force Balance, boundaries:
    % Left Column:
    if strcmp(leftbctype,'flux')==1
        % Contribution from w(z+.5,x):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,1)-1)*numvars+VertForceInds_ud(2:end-1,1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=-.5*Viscosity_lrud(2:end-1,2).*DZ_lrud(2:end-1,2).*Width_lrud(2:end-1,2)./((dx^2)*DZ_ud(2:end-1,1).*Width_ud(2:end-1,1))... % shear stress on right (inner) edge
            -(Viscosity_c(2:end,1)./(DZ_c(2:end,1).*DZ_ud(2:end-1,1))+Viscosity_c(1:end-1,1)./(DZ_c(1:end-1,1).*DZ_ud(2:end-1,1))); % long stress on top/bot edges
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z+.5,x-1):  (one column in from left)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,1)-1)*numvars+VertForceInds_ud(2:end-1,2);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=.5*Viscosity_lrud(2:end-1,2).*DZ_lrud(2:end-1,2).*Width_lrud(2:end-1,2)./((dx^2)*DZ_ud(2:end-1,2).*Width_ud(2:end-1,2)); % shear stress on left edge
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z+1.5,x):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(3:end,1)-1)*numvars+VertForceInds_ud(2:end-1,1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=Viscosity_c(2:end,1)./(DZ_c(2:end,1).*DZ_ud(2:end-1,1)); % long stress on top edge
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z-.5,x):  (applies to top corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-2))=(WInds_ud(2:end-2,1)-1)*numvars+VertForceInds_ud(3:end-1,1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-2))=Viscosity_c(2:end-1,1)./(DZ_c(2:end-1,1).*DZ_ud(3:end-1,1)); % long stress on bot edge
        lastterm=lastterm+(zsize-2);
    end
    % Right Column:
    if strcmp(rightbctype,'flux') || strcmp(rightbctype,'front') 
        % Contribution from w(z+.5,x):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,end)-1)*numvars+VertForceInds_ud(2:end-1,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=-.5*Viscosity_lrud(2:end-1,end-1).*DZ_lrud(2:end-1,end-1).*Width_lrud(2:end-1,end-1)./((dx^2)*DZ_ud(2:end-1,end).*Width_ud(2:end-1,end))... % shear stress on left (inner) edge
            -(Viscosity_c(2:end,end)./(DZ_c(2:end,end).*DZ_ud(2:end-1,end))+Viscosity_c(1:end-1,end)./(DZ_c(1:end-1,end).*DZ_ud(2:end-1,end))); % long stress on top/bot edges
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z+.5,x+1):  (one column in from right)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,end)-1)*numvars+VertForceInds_ud(2:end-1,end-1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=.5*Viscosity_lrud(2:end-1,end-1).*DZ_lrud(2:end-1,end-1).*Width_lrud(2:end-1,end-1)./((dx^2)*DZ_ud(2:end-1,end-1).*Width_ud(2:end-1,end-1)); % shear stress on right edge
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z+1.5,x):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(3:end,end)-1)*numvars+VertForceInds_ud(2:end-1,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=Viscosity_c(2:end,end)./(DZ_c(2:end,end).*DZ_ud(2:end-1,end)); % long stress on top
        lastterm=lastterm+(zsize-1);
        % Contribution from w(z-.5,x):  (applies to top corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-2))=(WInds_ud(2:end-2,end)-1)*numvars+VertForceInds_ud(3:end-1,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize-2))=Viscosity_c(2:end-1,end)./(DZ_c(2:end-1,end).*DZ_ud(3:end-1,end)); % long stress on bot
        lastterm=lastterm+(zsize-2);
    end
    
    % Vertical Force Balance Constraints:
    % Bottom row:
    StokesConstraint(reshape(VertForceInds_ud(2,:),[],1))...
        =reshape(Viscosity_c(1,:).*MeltRate_d./(DZ_ud(2,:).*DZ_c(1,:)),[],1);
    % Left Column:
    if strcmp(leftbctype,'fixedicethick')==1
        % Actual left edge:
        StokesConstraint(reshape(VertForceInds_ud(2:end-1,1),[],1))=...
            reshape(.5*Viscosity_lrud(2:end-1,2).*DZ_lrud(2:end-1,2).*Width_lrud(2:end-1,2).*W_ud(2:end-1,1)./((dx^2)*DZ_ud(2:end-1,1).*Width_ud(2:end-1,1))... % shear stress on right (inner) edge
            -Viscosity_c(2:end,1).*(W_ud(3:end,1)-W_ud(2:end-1,1))./(DZ_c(2:end,1).*DZ_ud(2:end-1,1))... % long stress on top
            +Viscosity_c(1:end-1,1).*(W_ud(2:end-1,1)-W_ud(1:end-2,1))./(DZ_c(1:end-1,1).*DZ_ud(2:end-1,1)),[],1); % long stress on bot
        % One column in from left edge:
        StokesConstraint(reshape(VertForceInds_ud(2:end-1,2),[],1))=StokesConstraint(reshape(VertForceInds_ud(2:end-1,2),[],1))...
            +reshape(-.5*Viscosity_lrud(2:end-1,2).*DZ_lrud(2:end-1,2).*Width_lrud(2:end-1,2).*W_ud(2:end-1,1)./((dx^2)*DZ_ud(2:end-1,2).*Width_ud(2:end-1,2)),[],1); % shear stress on left edge
    end
    % Right Column:
    if strcmp(rightbctype,'fixedicethick')==1
        % Actual right edge:
        StokesConstraint(reshape(VertForceInds_ud(2:end-1,end),[],1))=...
            reshape(.5*Viscosity_lrud(2:end-1,end-1).*DZ_lrud(2:end-1,end-1).*Width_lrud(2:end-1,end-1).*W_ud(2:end-1,end)./((dx^2)*DZ_ud(2:end-1,end).*Width_ud(2:end-1,end))... % shear stress on left (inner) edge
            -Viscosity_c(2:end,end).*(W_ud(3:end,end)-W_ud(2:end-1,end))./(DZ_c(2:end,end).*DZ_ud(2:end-1,end))... % long stress on top
            +Viscosity_c(1:end-1,end).*(W_ud(2:end-1,end)-W_ud(1:end-2,end))./(DZ_c(1:end-1,end).*DZ_ud(2:end-1,end)),[],1); % long stress on bot
        % One column in from right edge:
        StokesConstraint(reshape(VertForceInds_ud(2:end-1,end-1),[],1))=StokesConstraint(reshape(VertForceInds_ud(2:end-1,end-1),[],1))...
            +reshape(-.5*Viscosity_lrud(2:end-1,end-1).*DZ_lrud(2:end-1,end-1).*W_ud(2:end-1,end)./((dx^2)*DZ_ud(2:end-1,end-1)),[],1);
    end
    
    % Mass Balance, internal cells:  (z,x)
    % Contribution from u(z,x+.5):  (applies to top/bot rows and left column as well)
    StokesIndices(lastterm+1:lastterm+(xsize-1)*(zsize))=reshape((UInds_lr(:,2:end-1)-1)*numvars+MassBalInds_c(:,1:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-1)*(zsize))=reshape(-DZ_lr_upwind(:,2:end-1).*Width_lr(:,2:end-1)./(dx*DZ_c(:,1:end-1).*Width_c(:,1:end-1)),[],1); % flux through right edge
    lastterm=lastterm+(xsize-1)*(zsize);
    % Contribution from u(z,x-.5):  (applies to top/bot rows and right column as well)
    StokesIndices(lastterm+1:lastterm+(xsize-1)*(zsize))=reshape((UInds_lr(:,2:end-1)-1)*numvars+MassBalInds_c(:,2:end),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-1)*(zsize))=reshape(DZ_lr_upwind(:,2:end-1).*Width_lr(:,2:end-1)./(dx*DZ_c(:,2:end).*Width_c(:,2:end)),[],1); % flux through left edge
    lastterm=lastterm+(xsize-1)*(zsize);
    % Contribution from w(z+.5,x):  (applies to top/bot rows as well)
    StokesIndices(lastterm+1:lastterm+(xsize-2)*(zsize))=reshape((WInds_ud(2:end,2:end-1)-1)*numvars+MassBalInds_c(:,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2)*(zsize))=reshape(-1./DZ_c(:,2:end-1),[],1); % flux through top
    lastterm=lastterm+(xsize-2)*(zsize);
    % Contribution from w(z-.5,x):  (applies to top row as well)
    StokesIndices(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape((WInds_ud(2:end-1,2:end-1)-1)*numvars+MassBalInds_c(2:end,2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2)*(zsize-1))=reshape(1./DZ_c(2:end,2:end-1),[],1); % flux through bottom
    lastterm=lastterm+(xsize-2)*(zsize-1);
    
    % Mass Balance, boundaries:  (z,x)
    % Left Column:
    if strcmp(leftbctype,'flux')==1
        % Contribution from w(z+.5,x):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize))=(WInds_ud(2:end,1)-1)*numvars+MassBalInds_c(:,1);
        StokesCoefficients(lastterm+1:lastterm+(zsize))=-1./DZ_c(:,1); % flux through top
        lastterm=lastterm+(zsize);
        % Contribution from w(z-.5,x):  (applies to top corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,1)-1)*numvars+MassBalInds_c(2:end,1);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=1./DZ_c(2:end,1); % flux through bot
        lastterm=lastterm+(zsize-1);
    else % fixed icethick
        % Contribution from u(z,x-.5):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize))=(UInds_lr(:,1)-1)*numvars+MassBalInds_c(:,1);
        StokesCoefficients(lastterm+1:lastterm+(zsize))=DZ_lr_upwind(:,1).*Width_lr(:,1)./(dx*DZ_c(:,1).*Width_c(:,1)); % flux through left edge
        lastterm=lastterm+(zsize);
    end
    % Right Column:
    if strcmp(rightbctype,'flux')==1
        % Contribution from w(z+.5,x):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize))=(WInds_ud(2:end,end)-1)*numvars+MassBalInds_c(:,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize))=-1./DZ_c(:,end); % flux through top
        lastterm=lastterm+(zsize);
        % Contribution from w(z-.5,x):  (applies to top corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,end)-1)*numvars+MassBalInds_c(2:end,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=1./DZ_c(2:end,end); % flux through bot
        lastterm=lastterm+(zsize-1);
    elseif strcmp(rightbctype,'fixedicethick')==1
        % Contribution from u(z,x+.5):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize))=(UInds_lr(:,end)-1)*numvars+MassBalInds_c(:,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize))=-DZ_lr_upwind(:,end).*Width_lr(:,end)./(dx*DZ_c(:,end).*Width_c(:,end)); % flux through right edge
        lastterm=lastterm+(zsize);
    else % front
        % Contribution from u(z,x+.5):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize))=(UInds_lr(:,end)-1)*numvars+MassBalInds_c(:,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize))=-DZ_lr_upwind(:,end).*Width_lr(:,end)./(dx*DZ_c(:,end).*Width_c(:,end)); % flux through right edge
        lastterm=lastterm+(zsize);
        % Contribution from w(z+.5,x):  (applies to top/bot corners as well)
        StokesIndices(lastterm+1:lastterm+(zsize))=(WInds_ud(2:end,end)-1)*numvars+MassBalInds_c(:,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize))=-1./DZ_c(:,end); % flux through top
        lastterm=lastterm+(zsize);
        % Contribution from w(z-.5,x):  (applies to top corner as well)
        StokesIndices(lastterm+1:lastterm+(zsize-1))=(WInds_ud(2:end-1,end)-1)*numvars+MassBalInds_c(2:end,end);
        StokesCoefficients(lastterm+1:lastterm+(zsize-1))=1./DZ_c(2:end,end); % flux through bot
        lastterm=lastterm+(zsize-1);
    end
    
    % Mass Balance Constraints:
    % Bottom row:
    StokesConstraint(reshape(MassBalInds_c(1,:),[],1))...
        =MeltRate_d(1,:)./DZ_c(1,:); % flux through bottom edge
    % Left column:
    if strcmp(leftbctype,'flux')==1
        StokesConstraint(reshape(MassBalInds_c(:,1),[],1))=StokesConstraint(reshape(MassBalInds_c(:,1),[],1))...
            -DZ_lr_upwind(:,1).*Width_lr(:,1).*U_l./(dx*DZ_c(:,1).*Width_c(:,1)); % flux through left edge
    elseif strcmp(leftbctype,'fixedicethick')==1
        StokesConstraint(reshape(MassBalInds_c(:,1),[],1))=StokesConstraint(reshape(MassBalInds_c(:,1),[],1))...
            +(W_ud(2:end,1)-W_ud(1:end-1,1))./DZ_c(:,1); % flux through top/bot edges (prescribed vertical strain rate)
    end
    % Right column:
    if strcmp(rightbctype,'flux')==1
        StokesConstraint(reshape(MassBalInds_c(:,end),[],1))=StokesConstraint(reshape(MassBalInds_c(:,end),[],1))...
            +DZ_lr_upwind(:,end).*Width_lr(:,end).*U_r./(dx*DZ_c(:,end).*Width_c(:,end)); % flux through right edge
    elseif strcmp(rightbctype,'fixedicethick')==1
        StokesConstraint(reshape(MassBalInds_c(:,end),[],1))=StokesConstraint(reshape(MassBalInds_c(:,end),[],1))...
            +(W_ud(2:end,end)-W_ud(1:end-1,end))./DZ_c(:,end); % flux through top/bot edges (prescribed vertical strain rate)
    end
    
    % Top Free Surface Equation: (z,x)
    % Contribution from w(z+.5,x): (internal only)
    StokesIndices(lastterm+1:lastterm+(xsize-2))=reshape((WInds_ud(end,2:end-1)-1)*numvars+FreeSurfInds_u(2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2))=reshape(-Viscosity_c(end,2:end-1)./DZ_c(end,2:end-1),[],1); % long stress on surface
    lastterm=lastterm+(xsize-2);
    % Contribution from w(z-.5,x): (internal only)
    StokesIndices(lastterm+1:lastterm+(xsize-2))=reshape((WInds_ud(end-1,2:end-1)-1)*numvars+FreeSurfInds_u(2:end-1),[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize-2))=reshape(Viscosity_c(end,2:end-1)./DZ_c(end,2:end-1),[],1); % long stress on surface
    lastterm=lastterm+(xsize-2);
    % Top left corner:
    if strcmp(leftbctype,'fixedicethick')==0
        % Contribution from w(z+.5,x): (top left corner)
        StokesIndices(lastterm+1)=(WInds_ud(end,1)-1)*numvars+FreeSurfInds_u(1);
        StokesCoefficients(lastterm+1)=-Viscosity_c(end,1)./DZ_c(end,1); % long stress on surface
        lastterm=lastterm+1;
        % Contribution from w(z-.5,x): (top left corner)
        StokesIndices(lastterm+1)=(WInds_ud(end-1,1)-1)*numvars+FreeSurfInds_u(1);
        StokesCoefficients(lastterm+1)=Viscosity_c(end,1)./DZ_c(end,1); % long stress on surface
        lastterm=lastterm+1;
    end
    % Top right corner:
    if strcmp(rightbctype,'fixedicethick')==0
        % Contribution from w(z+.5,x): (top right corner)
        StokesIndices(lastterm+1)=(WInds_ud(end,end)-1)*numvars+FreeSurfInds_u(end);
        StokesCoefficients(lastterm+1)=-Viscosity_c(end,end)./DZ_c(end,end); % long stress on surface
        lastterm=lastterm+1;
        % Contribution from w(z-.5,x): (top right corner)
        StokesIndices(lastterm+1)=(WInds_ud(end-1,end)-1)*numvars+FreeSurfInds_u(end);
        StokesCoefficients(lastterm+1)=Viscosity_c(end,end)./DZ_c(end,end); % long stress on surface
        lastterm=lastterm+1;
    end
    % Contribution from P(z,x):
    StokesIndices(lastterm+1:lastterm+(xsize))=reshape((PInds_c(end,:)-1)*numvars+FreeSurfInds_u,[],1);
    StokesCoefficients(lastterm+1:lastterm+(xsize))=1; % pressure on surface
    lastterm=lastterm+(xsize);
    
    % Top Free Surface Constraint:
    if strcmp(leftbctype,'fixedicethick')==1
        StokesConstraint(FreeSurfInds_u(1))=...
            Viscosity_c(end,1)*(W_ud(end,1)-W_ud(end-1,1))/DZ_c(end,1); % long stress on surface
    end
    if strcmp(rightbctype,'fixedicethick')==1
        StokesConstraint(FreeSurfInds_u(end))=...
            Viscosity_c(end,end)*(W_ud(end,end)-W_ud(end-1,end))/DZ_c(end,end); % long stress on surface
    end
    
    % Right Side Stress Equation: (z,x)
    if strcmp(rightbctype,'front')==1 
        % Contribution from u(z,x+.5):
        StokesIndices(lastterm+1:lastterm+(zsize))=(UInds_lr(:,end)-1)*numvars+SideStressInds_r;
        StokesCoefficients(lastterm+1:lastterm+(zsize))=-Viscosity_c(:,end)/dx;
        lastterm=lastterm+(zsize);
        % Contribution from u(z,x-.5):
        StokesIndices(lastterm+1:lastterm+(zsize))=(UInds_lr(:,end-1)-1)*numvars+SideStressInds_r;
        StokesCoefficients(lastterm+1:lastterm+(zsize))=Viscosity_c(:,end)/dx;
        lastterm=lastterm+(zsize);
        % Contribution from P(z,x):
        StokesIndices(lastterm+1:lastterm+(zsize))=(PInds_c(:,end)-1)*numvars+SideStressInds_r;
        StokesCoefficients(lastterm+1:lastterm+(zsize))=1;
        lastterm=lastterm+(zsize);
    end
    
    % Right side stress constraint:
    if strcmp(rightbctype,'front')==1 
        StokesConstraint(reshape(SideStressInds_r,[],1))=NormalStress_r;
    end
    
    % Create the sparse matrix:
    [StokesIndices1,StokesIndices2]=ind2sub([numvars,numvars],StokesIndices(1:lastterm));
    StokesMatrix=sparse(StokesIndices1,StokesIndices2,StokesCoefficients(1:lastterm),numvars,numvars,lastterm);
    
    % Solve The Linear Problem:
    if strcmp(solvertype,'pcg')==1
        [SolutionVector,flag,linerror,numiterations]=pcg(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'bicg')==1
        [SolutionVector,flag,linerror,numiterations]=bicg(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'bicgstab')==1
        [SolutionVector,flag,linerror,numiterations]=bicgstab(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'bicgstabl')==1
        [SolutionVector,flag,linerror,numiterations]=bicgstabl(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'cgs')==1
        [SolutionVector,flag,linerror,numiterations]=cgs(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'gmres')==1
        [SolutionVector,flag,linerror,numiterations]=gmres(StokesMatrix,StokesConstraint,[],tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'lsqr')==1
        [SolutionVector,flag,linerror,numiterations]=lsqr(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'minres')==1
        [SolutionVector,flag,linerror,numiterations]=minres(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'qmr')==1
        [SolutionVector,flag,linerror,numiterations]=qmr(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'symmlq')==1
        [SolutionVector,flag,linerror,numiterations]=symmlq(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    elseif strcmp(solvertype,'tfqmr')==1
        [SolutionVector,flag,linerror,numiterations]=tfqmr(StokesMatrix,StokesConstraint,tolerance,maxiterations,[],[],SolutionVector);
    else % direct solve
        SolutionVector=full(StokesMatrix\StokesConstraint);
    end
    
    % Check that the iterative linear solve was successful:
    if exist('flag','var')
        if flag~=0
            disp(['Solver Failure Flag= ',num2str(flag)])
            error('Linear Velocity Solver Failed.')
        end
    end
    
    % Assign horizontal velocity from SolutionVector to grid:
    U_lr(:,2:end-1)=SolutionVector(UInds_lr(:,2:end-1));
    if strcmp(leftbctype,'flux')==1
        U_lr(:,1)=U_l;
    else
        U_lr(:,1)=SolutionVector(UInds_lr(:,1));
    end
    if strcmp(rightbctype,'flux')==1
        U_lr(:,end)=U_r;
    else
        U_lr(:,end)=SolutionVector(UInds_lr(:,end));
    end
    
    % Assign vertical velocity from SolutionVector to grid:
    if strcmp(leftbctype,'fixedicethick')==1 && strcmp(rightbctype,'fixedicethick')==1
        W_ud(2:end,2:end-1)=SolutionVector(WInds_ud(2:end,2:end-1));
    elseif strcmp(leftbctype,'fixedicethick')==1
        W_ud(2:end,2:end)=SolutionVector(WInds_ud(2:end,2:end));
    elseif strcmp(rightbctype,'fixedicethick')==1
        W_ud(2:end,1:end-1)=SolutionVector(WInds_ud(2:end,1:end-1));
    else
        W_ud(2:end,:)=SolutionVector(WInds_ud(2:end,:));
    end
    W_ud(1,:)=-MeltRate_d;
    
    % Assign dynamic pressure from SolutionVector to grid:
    PressureDynamic_c=SolutionVector(PInds_c);
    
    % Deal with horizontal regriding in moving front case:
    if strcmp(rightbctype,'front') && strcmp(calvingtype,'fixed')==0
        % Compute front velocity:
        u_r=sum(U_lr(:,end).*DZhat_c)-calvingrate_r-meltrate_r;  
        % Calculate grid velocity:
        Ugrid_lr=repmat(u_r*X_lr/domainwidth_mid,[zsize,1]);
    else
        % Grid is not moving horizontally:
        u_r=0;
        Ugrid_lr=zeros(zsize,xsize+1);
    end
    
    % Set lower corners horizontal velocity:
    u_ld=leftbcslidefrac*sum(U_lr(:,1).*DZhat_c);
    u_rd=rightbcslidefrac*sum(U_lr(:,end).*DZhat_c);
    
    % Compute vertically averaged flow direction: (relative to grid)
    FlowDir_lr=sign(sum((U_lr-Ugrid_lr).*repmat(DZhat_c,[1,xsize+1]),1));
    
    % Check if flow direction has changed:
    if isequal(FlowDir_lr,FlowDir_lr_last)==0
        
        % Reassign upwind icethick:
        Icethick_lr_upwind=zeros(1,xsize+1);
        CandidateIcethick=[icethick_l,Icethick_c_mid];
        Icethick_lr_upwind(FlowDir_lr>=0)=CandidateIcethick(FlowDir_lr>=0);
        CandidateIcethick=[Icethick_c_mid,icethick_r];
        Icethick_lr_upwind(FlowDir_lr<0)=CandidateIcethick(FlowDir_lr<0);
        clear CandidateIcethick
        
        % Recompute upwind DZ:
        DZ_lr_upwind=repmat(DZhat_c,[1,xsize+1]).*repmat(Icethick_lr_upwind,[zsize,1]);
        
        % Force a repeat iteration:
        forcerepeat=1;
        
    end
    
    % Compute vertical velocity from englacial melting:
    W_ud=W_ud+cumsum([zeros(1,xsize);-MeltRate_c.*DZ_c],1);
    
    % Compute bottom edge horizontal velocities:
    U_lrd(IsGrounded_lr~=0&IsInternal_lrd)=U_lr(1,IsGrounded_lr~=0&IsInternal_lrd)./(2*DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd).*DZ_lrud(1,IsGrounded_lr~=0&IsInternal_lrd)./Viscosity_lrud(1,IsGrounded_lr~=0&IsInternal_lrd)+1);
    U_lrd(IsGrounded_lr==0&IsInternal_lrd)=U_lr(1,IsGrounded_lr==0&IsInternal_lrd);
    U_lrd(1)=u_ld;
    U_lrd(end)=u_rd;
    
    % Compute basal drag:
    Drag_lrd(2:end-1)=DragCoefficient_lrd(2:end-1).*GroundedFraction_lr(2:end-1).*U_lrd(2:end-1);
    
    % Recompute basal drag coefficient:
    DragCoefficient_lrd(2:end-1)=GroundedFraction_lr(2:end-1).*(SlidingConstant_lrd(2:end-1).^(-1/m)).*(U_lrd(2:end-1).^((1-m)/m));
    
    % Check drag coefficient:
    if sum(isnan(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd))|isinf(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)))~=0
        error('Drag Coefficient has impossible values.')
    end
    
    % Recompute side drag coefficient:
    if dosidedrag
        DragCoefficient_lr=((A_lr.*Width_lr).^(-1/n)).*(U_lr.^((1-n)/n));
    end
    
    % Check side drag coefficient:
    if dosidedrag
        if sum(isnan(DragCoefficient_lr)|isinf(DragCoefficient_lr))~=0
            error('Side drag coefficient has impossible values.')
        end
    end
    
    % Compute shear strain rates:
    StrainRate_xz_lrud=zeros(zsize+1,xsize+1);
    StrainRate_xz_lrud(1:end-1,2:end-1)=.5*([U_lr(1,2:end-1)-U_lrd(2:end-1);(U_lr(2:end,2:end-1)-U_lr(1:end-1,2:end-1))]./DZ_lrud(1:end-1,2:end-1)...     % du/dz
        +(W_ud(1:end-1,2:end)-W_ud(1:end-1,1:end-1))/dx);  % dw/dx
    % Interpolate shear strain rate to grid centers:
    StrainRate_xz_c=.25*(StrainRate_xz_lrud(1:end-1,1:end-1)+StrainRate_xz_lrud(2:end,1:end-1)+StrainRate_xz_lrud(1:end-1,2:end)+StrainRate_xz_lrud(2:end,2:end)); 
    
    % Compute longitudinal strain rates:
    StrainRate_xx_c=(U_lr(:,2:end)-U_lr(:,1:end-1))/dx;      % du/dx
    StrainRate_zz_c=(W_ud(2:end,:)-W_ud(1:end-1,:))./DZ_c;   % dw/dz
    
    % Interpolate horizontal velocity to grid centers:
    U_c=.5*(U_lr(:,1:end-1)+U_lr(:,2:end));
    
    % Compute cross-flow strain rate:
    StrainRate_yy_c=-U_c.*(Width_lr(:,2:end)-Width_lr(:,1:end-1))./(dx*Width_c);
    
    % Correct vertical strain rate for grid geometry:
    StrainRate_zz_c=StrainRate_zz_c+U_c.*repmat((Icethick_lr_upwind(2:end)-Icethick_lr_upwind(1:end-1))./(Icethick_c_mid*dx),[zsize,1]);
    
    % Note: the grid geometry correction has an error related to velocity
    % curvature.  In a test I ran while coding this, the error produced 
    % up to a 10% mismatch between the xx and zz strain rates where 
    % velocity had strong curvature (~0.3% elsewhere).  However, this error
    % only effects the strain rate used for the effective viscosity 
    % computation.  Mass conservation is unaffected.
    
    % Recompute effective viscosity:
    Viscosity_c=(A_c.^(-1/n)).*((2*StrainRate_xz_c.^2+StrainRate_xx_c.^2+StrainRate_zz_c.^2+strainratestabilizer^2).^((1-n)/(2*n)));
    
    % Check viscosity:
    if sum(sum(isnan(Viscosity_c)|isinf(Viscosity_c)|Viscosity_c==0))~=0
        error('Viscosity has impossible values.')
    end
    
    % Compute misfit:
    if dosidedrag
        misfit_stokes=max([sqrt(mean(mean(log(Viscosity_c./Viscosity_c_last).^2)))/viscositylogscale,... % viscosity misfit
            sqrt(mean(log(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)./DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd)).^2))/dragcoefflogscale,... % basal drag coefficient misfit (grounded only)
            sqrt(mean(mean(log(DragCoefficient_lr./DragCoefficient_lr_last).^2)))/sidedragcoefflogscale]); % side drag coefficient misfit
    else
        misfit_stokes=max([sqrt(mean(mean(log(Viscosity_c./Viscosity_c_last).^2)))/viscositylogscale,... % viscosity misfit
            sqrt(mean(log(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)./DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd)).^2))/dragcoefflogscale]); % basal drag coefficient misfit (grounded only)
    end
    
    % Apply logarithmic damping_visc to nonlinear iteration:
    if damping_visc~=0
        % Viscosity:
        Viscosity_c=Viscosity_c_last.*((Viscosity_c./Viscosity_c_last).^(1-damping_visc));
        % Side drag:
        if dosidedrag
            DragCoefficient_lr=DragCoefficient_lr_last.*((DragCoefficient_lr./DragCoefficient_lr_last).^(1-damping_visc));
        end
        % Basal drag:
        if iteration_stokes==1 
            % Omit damping for newly grounded cells:
            DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd&NewlyGrounded_lr==0)=DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd&NewlyGrounded_lr==0).*((DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd&NewlyGrounded_lr==0)./DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd&NewlyGrounded_lr==0)).^(1-damping_drag));
        else
            % Damp all grounded cells:
            DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)=DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd).*((DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)./DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd)).^(1-damping_drag));
        end
    end
    
    % Break from loop:
    if ((misfit_stokes<tolerance_visc && iteration_stokes>=miniterations_visc) || (n==1 && m==1)) && forcerepeat==0
        converged_stokes=1;
    elseif iteration_stokes>=maxiterations_visc
        % Try resetting the iteration:
        if triedreset==0
            % Reset everything to initial guess:
            Viscosity_c=viscosityguess*ones(zsize,xsize);
            Viscosity_c_last=Viscosity_c;
            DragCoefficient_lrd=dragcoefficientguess*ones(1,xsize+1);
            DragCoefficient_lrd_last=DragCoefficient_lrd;
            if dosidedrag
                DragCoefficient_lr=sidedragcoefficientguess*ones(zsize,xsize+1);
                DragCoefficient_lr_last=DragCoefficient_lr;
            end
            FlowDir_lr_last=FlowDir_lr;
            % Start iterating from the top:
            iteration_stokes=1;
            triedreset=1;
        else
            % Throw the error message:
            error('Unable to converge on a viscosity solution.')
        end
    else
        % Record last guess:
        Viscosity_c_last=Viscosity_c;
        DragCoefficient_lrd_last=DragCoefficient_lrd;
        if dosidedrag
            DragCoefficient_lr_last=DragCoefficient_lr;
        end
        FlowDir_lr_last=FlowDir_lr;
        % Count iterations:
        iteration_stokes=iteration_stokes+1;
        % Reset forcerepeat variable:
        forcerepeat=0;
    end
    
end
