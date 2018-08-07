% ShallowShelfApproximation_2D_v1

% Mike Wolovick, 2/27/2016

% This script solves the shallow shelf approximation to the Stokes 
% equations on a skewed grid in 2D.  The shallow shelf approximation
% only considers longitudinal stresses (sigma_xx) and gets vertical 
% velocity from mass balance.  

% The solution method is to linearize the problem for a given viscosity
% field and iterate for the nonlinear viscosity field.  The linear problem
% is solved with either Matlab's direct sparse matrix solver "\" or one of 
% the built-in iterative solvers, and the nonlinear iteration method is 
% damped fixed-point iteration.  The nonlinear iteration also solves for a 
% basal drag coefficient and a side drag coefficient.  

% The shallow shelf approximation computes vertically averaged velocity
% based on vertically averaged force balance

% When computing vertically averaged viscosity, the vertical average of 
% A^(-1/n) is used for the rate factor

% Boundary Conditions:

% Top BC:     stress_xz=0

% Bottom BC:  power-law sliding
%             w=-meltrate(x)

% Side BC, 3 Options:
%   1.  flux condition:  u=u(z)
%   2.  fixed icethick:  w=w(z)
%   3.  moving front:    stress_xx=stress_xx(z)

% Viscosity BC: no gradient  (required for interpolation to all grid cell corners)

% This script cannot function outside of the context of the variables
% created by Flowline_v1

% Description of the linear problem:  
% Number of unknowns:            xsize+1
% Number of equations:           xsize+1
% Max number of terms possible:  3*xsize+1

% This script does not use indexing keys because the SSA is
% one-dimensional, and indexing in the Stokes Matrix is therefore
% straightforward.

% StokesMatrix indexing for a single term in the linear problem:

% MatrixInd=(VariableInd-1)*numvars+EquationInd

% where "VariableInd" and "EquationInd" are the x-indices for the variable
% and equation under consideration.

%% Solve Stokes Equation:

% Create solution vector for the first time:
numvars=xsize+1;
if timestep==1 && strcmp(solvertype,'direct')==0
    SolutionVector=zeros(numvars,1);
end

% Evaulate rheological constant on grid centers:
if useinputa
    A_c=a_input*ones(zsize,xsize);
else
    A_c(Temp_c>=t0)=a0*exp(-(q_big/r)*((1./Temp_c(Temp_c>=t0))-(1/t0)));
    A_c(Temp_c<t0)=a0*exp(-(q_small/r)*((1./Temp_c(Temp_c<t0))-(1/t0)));
end

% Take vertically average of A^(-1/n) to compute viscosity factor:
ViscosityFactor_c=sum((A_c.^(-1/n)).*repmat(DZhat_c,[1,xsize]),1);

% Compute side drag factor: (interpolate viscosity factor and incorporate
% width)
if useinputsidea
    SideDragFactor_lr=(a_side_input*Width_lr).^(-1/n);
else
    SideDragFactor_lr=[ViscosityFactor_c(1),.5*(ViscosityFactor_c(1:end-1)+ViscosityFactor_c(2:end)),ViscosityFactor_c(end)].*(Width_lr.^(-1/n));
end

% Nonlinear iteration:
converged_stokes=0;
iteration_stokes=1;
triedreset=0;
Viscosity_c_last=Viscosity_c;
DragCoefficient_lrd_last=DragCoefficient_lrd;
if dosidedrag
    DragCoefficient_lr_last=DragCoefficient_lr;
end
FlowDir_lr_last=FlowDir_lr;
while converged_stokes==0
    
    %% Reset stokes matrix:
    StokesMatrix=sparse([],[],[],xsize+1,xsize+1,3*xsize+1);
    StokesConstraint=zeros(xsize+1,1);
    
    % Horizontal Force Balance, internal cells:  (x+.5)
    % Contribution from u(x+.5): (central velocity) 
    StokesMatrix((xsize+1)+2:xsize+2:((xsize)-1)*(xsize+1)+xsize)=-(Width_c(1:end-1).*Icethick_c(1:end-1).*Viscosity_c(1:end-1)+Width_c(2:end).*Icethick_c(2:end).*Viscosity_c(2:end))/(dx^2)...
        -(DragCoefficient_lrd(2:end-1).*Width_lr(2:end-1)+2*DragCoefficient_lr(2:end-1).*Icethick_lr_center(2:end-1));
    % Contribution from u(x-.5): (upstream velocity)
    StokesMatrix(2:xsize+2:((xsize-1)-1)*(xsize+1)+xsize)=Width_c(1:end-1).*Icethick_c(1:end-1).*Viscosity_c(1:end-1)/(dx^2);
    % Contribution from u(x+1.5): (downstream velocity)
    StokesMatrix(2*(xsize+1)+2:xsize+2:((xsize+1)-1)*(xsize+1)+xsize)=Width_c(2:end).*Icethick_c(2:end).*Viscosity_c(2:end)/(dx^2);
    
    % Add the correction for implicit ice surface:
    if doimplicitsurface
        % Contribution from u(x+.5): (central velocity)
        StokesMatrix((xsize+1)+2:xsize+2:((xsize)-1)*(xsize+1)+xsize)=StokesMatrix((xsize+1)+2:xsize+2:((xsize)-1)*(xsize+1)+xsize)-(rho_i*g*dt*Icethick_lr_center(2:end-1).*Icethick_lr_upwind(2:end-1).*(Width_lr(2:end-1).^2)/(dx^2))...
            .*((1-(rho_i/rho_sw)*(1-IsGrounded_c(2:end)))./Width_c(2:end)+(1-(rho_i/rho_sw)*(1-IsGrounded_c(1:end-1)))./Width_c(1:end-1));
        % Contribution from u(x-.5): (upstream velocity)
        StokesMatrix(2:xsize+2:((xsize-1)-1)*(xsize+1)+xsize)=StokesMatrix(2:xsize+2:((xsize-1)-1)*(xsize+1)+xsize)+rho_i*g*dt*Icethick_lr_center(2:end-1).*Icethick_lr_upwind(1:end-2).*Width_lr(2:end-1).*Width_lr(1:end-2).*(1-(rho_i/rho_sw)*(1-IsGrounded_c(1:end-1)))./((dx^2)*Width_c(1:end-1));
        % Contribution from u(x+1.5): (downstream velocity)
        StokesMatrix(2*(xsize+1)+2:xsize+2:((xsize+1)-1)*(xsize+1)+xsize)=StokesMatrix(2*(xsize+1)+2:xsize+2:((xsize+1)-1)*(xsize+1)+xsize)+rho_i*g*dt*Icethick_lr_center(2:end-1).*Icethick_lr_upwind(3:end).*Width_lr(2:end-1).*Width_lr(3:end).*(1-(rho_i/rho_sw)*(1-IsGrounded_c(2:end)))./((dx^2)*Width_c(2:end));
    end
    
    % Left boundary:
    if strcmp(leftbctype,'flux')
        % Contribution from u(boundary):
        StokesMatrix(1)=Viscosity_c(1)*Icethick_c(1)*Width_c(1)/(dx^2); % Dirichlet "equation", scaled up
    else % fixedicethick
        % Contribution from u(x-.5):
        StokesMatrix(1)=Icethick_lr_upwind(1)*Width_lr(1)./(dx*Icethick_c(1).*Width_c(1));
        % Contribution from u(x+.5):
        StokesMatrix(xsize+1+1)=-Icethick_lr_upwind(2)*Width_lr(2)./(dx*Icethick_c(1).*Width_c(1));
    end
    
    % Right boundary:
    if strcmp(rightbctype,'flux')
        % Contribution from u(boundary):
        StokesMatrix((xsize)*(xsize+1)+xsize+1)=Viscosity_c(end)*Icethick_c(end)*Width_c(end)/(dx^2); % Dirichlet "equation", scaled up
    elseif strcmp(rightbctype,'fixedicethick')
        % Contribution from u(x-.5):
        StokesMatrix((xsize-1)*(xsize+1)+xsize+1)=Icethick_lr_upwind(end)*Width_lr(end)./(dx*Icethick_c(end).*Width_c(end));
        % Contribution from u(x+.5):
        StokesMatrix((xsize)*(xsize+1)+xsize+1)=-Icethick_lr_upwind(end)*Width_lr(end)./(dx*Icethick_c(end).*Width_c(end));
    else % front
        % Contribution from u(x-.5):
        StokesMatrix((xsize-1)*(xsize+1)+xsize+1)=Viscosity_c(end)/dx;
        % Contribution from u(x+.5):
        StokesMatrix((xsize)*(xsize+1)+xsize+1)=-Viscosity_c(end)/dx;
    end
    
    % Constraints:
    % Basic driving stress:
    StokesConstraint(2:end-1)=sum(-DrivingStressGrad_lr(:,2:end-1).*DZ_lr_center(:,2:end-1),1).*Width_lr(2:end-1);
    % Implicit surface correction term:
    if doimplicitsurface
        StokesConstraint(2:end-1)=StokesConstraint(2:end-1)+(rho_i*g*dt*Icethick_lr_center(2:end-1).*...
            ((1-(rho_i/rho_sw)*(1-IsGrounded_c(2:end))).*(Accum_c(2:end)-MeltRate_d(2:end)-sum(MeltRate_c(:,2:end).*DZ_c(:,2:end),1))...
            -(1-(rho_i/rho_sw)*(1-IsGrounded_c(1:end-1))).*(Accum_c(1:end-1)-MeltRate_d(1:end-1)-sum(MeltRate_c(:,1:end-1).*DZ_c(:,1:end-1),1)))/dx)';
    end
    % Left boundary:
    if strcmp(leftbctype,'flux')
        % Set velocity:
        StokesConstraint(1)=sum(U_l.*DZhat_c)*Viscosity_c(1)*Icethick_c(1)*Width_c(1)/(dx^2); 
    else % fixedicethick
        % Mass removed:
        StokesConstraint(1)=-Width_c(1)*dx*(Accum_c(1)-MeltRate_d(1)-sum(MeltRate_c(:,1).*DZ_c(:,1))-icethickchangerate_l);
    end
    % Right boundary:
    if strcmp(rightbctype,'flux')
        % Set velocity:
        StokesConstraint(end)=sum(U_r.*DZhat_c)*Viscosity_c(end)*Icethick_c(end)*Width_c(end)/(dx^2); 
    elseif strcmp(rightbctype,'fixedicethick')
        % Mass removed:
        StokesConstraint(end)=-Width_c(end)*dx*(Accum_c(end)-MeltRate_d(end)-sum(MeltRate_c(:,end).*DZ_c(:,end))-icethickchangerate_r);
    else % front
        % Unbalanced pressure:
        StokesConstraint(end)=sum(NormalStress_r.*DZhat_c);
    end
    
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
        
    %% Assign horizontal velocity from SolutionVector to grid:
    U_lr(:,2:end-1)=repmat(SolutionVector(2:end-1)',[zsize,1]);
    if strcmp(leftbctype,'flux')==1
        U_lr(:,1)=U_l;
    else
        U_lr(:,1)=SolutionVector(1);
    end
    if strcmp(rightbctype,'flux')==1
        U_lr(:,end)=U_r;
    else
        U_lr(:,end)=SolutionVector(end);
    end
    
    % Set lower corners horizontal velocity:
    u_ld=U_lr(1,1);
    u_rd=U_lr(1,end);
    
    % Compute vertically averaged flow direction: 
    FlowDir_lr=sign(U_lr(1,:));
    
    % Check if flow direction has changed:
    if isequal(FlowDir_lr,FlowDir_lr_last)==0
        
        % Reassign upwind icethick:
        Icethick_lr_upwind=zeros(1,xsize+1);
        CandidateIcethick=[icethick_l,Icethick_c];
        Icethick_lr_upwind(FlowDir_lr>=0)=CandidateIcethick(FlowDir_lr>=0);
        CandidateIcethick=[Icethick_c,icethick_r];
        Icethick_lr_upwind(FlowDir_lr<0)=CandidateIcethick(FlowDir_lr<0);
        clear CandidateIcethick
        
        % Recompute upwind DZ:
        DZ_lr_upwind=repmat(DZhat_c,[1,xsize+1]).*repmat(Icethick_lr_upwind,[zsize,1]);
        
    end
    
    % Compute vertical velocity from mass balance:
    W_ud=cumsum([-MeltRate_d;-(U_lr(:,2:end).*DZ_lr_upwind(:,2:end).*repmat(Width_lr(:,2:end),[zsize,1])-U_lr(:,1:end-1).*DZ_lr_upwind(:,1:end-1).*repmat(Width_lr(:,1:end-1),[zsize,1]))./repmat(dx*Width_c,[zsize,1])],1);
    
    % Assign bottom edge horizontal velocities:
    U_lrd=U_lr(1,:);
    
    % Compute basal drag:
    Drag_lrd(2:end-1)=DragCoefficient_lrd(2:end-1).*GroundedFraction_lr(2:end-1).*U_lrd(2:end-1);
    
    % Recompute basal drag coefficient:
    DragCoefficient_lrd(2:end-1)=GroundedFraction_lr(2:end-1).*(SlidingConstant_lrd(2:end-1).^(-1./M_lrd(2:end-1))).*((abs(U_lrd(2:end-1))+slidingstabilizer).^((1-M_lrd(2:end-1))./M_lrd(2:end-1)));
    
    % Check drag coefficient:
    if sum(isnan(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd))|isinf(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd))|DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)<0|imag(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd))~=0)~=0
        error('Basal drag coefficient has impossible values.')
    end
    
    % Recompute side drag coefficient:
    if dosidedrag
        DragCoefficient_lr=SideDragFactor_lr.*((abs(U_lr(1,:))+slidingstabilizer).^((1-n)/n));
    end
    
    % Check side drag coefficient:
    if dosidedrag
        if sum(isnan(DragCoefficient_lr(2:end-1))|isinf(DragCoefficient_lr(2:end-1))|DragCoefficient_lr(2:end-1)<0|imag(DragCoefficient_lr(2:end-1))~=0)~=0
            error('Side drag coefficient has impossible values.')
        end
    end
    
    % Compute longitudinal strain rates:
    StrainRate_xx_c=(U_lr(1,2:end)-U_lr(1,1:end-1))/dx;      % du/dx
    StrainRate_zz_c=(W_ud(2,:)-W_ud(1,:))./DZ_c(1,:);   % dw/dz
    
    % Interpolate horizontal velocity to grid centers:
    U_c=.5*(U_lr(1,1:end-1)+U_lr(1,2:end));
    
    % Compute cross-flow strain rate:
    StrainRate_yy_c=-U_c.*(Width_lr(2:end)-Width_lr(1:end-1))./(dx*Width_c);
    
    % Compute transverse strain rate:
    if dosidedrag
        StrainRate_xy_c=U_c./Width_c; % .5 cancels in numerator and denominator
    else
        StrainRate_xy_c=zeros(1,xsize);
    end
    
    % Correct vertical strain rate for grid geometry:
    StrainRate_zz_c=StrainRate_zz_c+U_c.*(Icethick_lr_upwind(2:end)-Icethick_lr_upwind(1:end-1))./(Icethick_c*dx);
    
    % Note: the grid geometry correction has an error related to velocity
    % curvature.  In a test I ran while coding this, the error produced 
    % up to a 10% mismatch between the xx and zz strain rates where 
    % velocity had strong curvature (~0.3% elsewhere).  However, this error
    % only effects the strain rate used for the effective viscosity 
    % computation.  Mass conservation is unaffected.
    
    % Recompute effective viscosity:
    %Viscosity_c=ViscosityFactor_c.*((StrainRate_xx_c.^2+StrainRate_yy_c.^2+StrainRate_zz_c.^2+2*StrainRate_xy_c.^2+strainratestabilizer^2).^((1-n)/(2*n)));
    Viscosity_c=2*ViscosityFactor_c.*((StrainRate_xx_c.^2+StrainRate_xy_c.^2+strainratestabilizer^2).^((1-n)/(2*n)));
    
    
    % Check viscosity:
    if sum(sum(isnan(Viscosity_c)|isinf(Viscosity_c)|Viscosity_c<=0|imag(Viscosity_c)~=0))~=0
        error('Viscosity has impossible values.')
    end
    
    % Compute misfit:
    if dosidedrag
        misfit_stokes=max([sqrt(mean(mean(log(Viscosity_c./Viscosity_c_last).^2)))/viscositylogscale,... % viscosity misfit
            sqrt(mean(log(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)./DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd)).^2))/dragcoefflogscale,... % basal drag coefficient misfit (grounded only)
            sqrt(mean(mean(log(DragCoefficient_lr(2:end-1)./DragCoefficient_lr_last(2:end-1)).^2)))/sidedragcoefflogscale]); % side drag coefficient misfit
    else
        misfit_stokes=max([sqrt(mean(mean(log(Viscosity_c./Viscosity_c_last).^2)))/viscositylogscale,... % viscosity misfit
            sqrt(mean(log(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)./DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd)).^2))/dragcoefflogscale]); % basal drag coefficient misfit (grounded only)
    end
    
    % Apply damping to nonlinear iteration:
    if damping_visc~=0
        % Viscosity:
        Viscosity_c=Viscosity_c_last+(Viscosity_c-Viscosity_c_last)*(1-damping_visc);
        % Side drag:
        if dosidedrag
            DragCoefficient_lr=DragCoefficient_lr_last+(DragCoefficient_lr-DragCoefficient_lr_last)*(1-damping_visc);
        end
        % Basal drag:
        if iteration_stokes==1 
            % Omit damping for newly grounded cells:
            DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd&NewlyGrounded_lr==0)=DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd&NewlyGrounded_lr==0)+(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd&NewlyGrounded_lr==0)-DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd&NewlyGrounded_lr==0))*(1-damping_drag);
        else
            % Damp all grounded cells:
            DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)=DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd)+(DragCoefficient_lrd(IsGrounded_lr~=0&IsInternal_lrd)-DragCoefficient_lrd_last(IsGrounded_lr~=0&IsInternal_lrd))*(1-damping_drag);
        end
    end
    
    % Break from loop:
    if (misfit_stokes<tolerance && iteration_stokes>=miniterations) || (n==1 && m==1)  
        converged_stokes=1;
        %disp(num2str(iteration_stokes))
    elseif iteration_stokes>=maxiterations
        % Try resetting the iteration:
        if triedreset==0
            % Reset everything to initial guess:
            Viscosity_c=viscosityguess*ones(1,xsize);
            Viscosity_c_last=Viscosity_c;
            DragCoefficient_lrd=dragcoefficientguess*IsGrounded_lr;
            DragCoefficient_lrd_last=DragCoefficient_lrd;
            if dosidedrag
                DragCoefficient_lr=sidedragcoefficientguess*ones(1,xsize+1);
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
