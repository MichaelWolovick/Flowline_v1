% SteadyStateIC_2D

% Mike Wolovick, 10/21/2011

% This script makes the temperature and water thickness initial conditions
% for FTM_2D.  The temperature IC is a 1D steady-state solution assuming
% constant vertical strain rate and neglecting horizontal processes 
% (Horizontal advection added in v3).  The water thickness IC is a balance 
% flux solution in equilibrium with the temperature IC. 

% This algorithm fills closed basins that have meltwater input, and it puts
% zero water thickness at peaks by default.

% This script will not produce the true steady state of the coupled model,
% but it is designed to give a decent approximation to it (at least on the 
% thermal side).

% Changes for v2:  meltpoint is defined on the bottom surface of the grid,
% instead of at grid centers, and the script has been modified for inflow
% water BC on the left.

% Changes for v3:  the script is now modified for variable vertical grid
% spacing and includes horizontal temperature advection in the IC.
% Horizontal velocity is assumed to be constant with depth, consistent with
% a Nye vertical velocity profile.  In addition, flow of ice is assumed to
% always be from left to right.  The same assumption is made for water if
% quasistatic hydrology is active.

% Changes for v4:  the script now runs with Flowline_v1 instead of FTM_2D.
% Dynamic hydrology has been removed.  Ungrounded cells are automatically
% tied to the melting point.  Ungrounded cells still receive geothermal
% flux for simplicity.  Velocity is now defined on grid edges, consistent
% with the rest of the model.


% Define horizontal velocity in first grid edge:
if strcmp(leftbctype,'flux')
    U_lr(:,1)=influx_l/icethick_l;
else
    U_lr(:,1)=0;
end

% Pre-allocate water flux and basal logical state:
WaterFlux_lrd=zeros(1,xsize+1);
WaterFlux_lrd(1)=waterinflux_l;
IsTied_d=true(1,xsize);

% Loop through columns:
for d2=1:xsize
    
    % Firct guess velocity (no melting):
    W_ud(:,d2)=-Zhat_ud*(Accum_u(d2)-AnnualMelt_u(d2));
    U_lr(:,d2+1)=U_lr(:,d2)+(dx/Icethick_lr_upwind(d2+1))*(Accum_u(d2)-AnnualMelt_u(d2));
    
    % Compute first guess of basal shear heating:
    if dostrainheat==1
        SlideHeat_d(d2)=abs(.5*(DrivingStressGrad_lr(1,d2)*Icethick_lr_center(d2)*U_lr(1,d2)+DrivingStressGrad_lr(1,d2+1)*Icethick_lr_center(d2+1)*U_lr(1,d2+1)));
    end
    
    % Compute hydraulically limited freeze rate:
    meltrate_hydro=-(rho_fw/rho_i)*WaterFlux_lrd(d2)/dx;  % m/s of ice
    
    % First guess temperature:
    if d2>1 && U_lr(1,d2)>0 && U_lr(1,d2+1)>0 && doichorzadv==1
        [Temp_c(:,d2),Temp_d(d2),MeltRate_d(d2)]=TempProfile1D_column_v4a(Temp_c(:,d2-1),zeros(zsize,1),U_lr(:,d2:d2+1),W_ud(:,d2),IsTied_d(d2),MeltRate_d(d2),DZ_c(:,d2),DZ_lr_upwind(:,d2:d2+1),dx,Gflux_c(d2)+SlideHeat_d(d2),SurfTemp_u(d2),MeltPoint_d(d2),ThermalParameters);
    else
        [Temp_c(:,d2),Temp_d(d2),MeltRate_d(d2)]=TempProfile1D_column_v4a(NaN*zeros(zsize,1),zeros(zsize,1),U_lr(:,d2:d2+1),W_ud(:,d2),IsTied_d(d2),MeltRate_d(d2),DZ_c(:,d2),DZ_lr_upwind(:,d2:d2+1),dx,Gflux_c(d2)+SlideHeat_d(d2),SurfTemp_u(d2),MeltPoint_d(d2),ThermalParameters);
    end
    
    % Iterate for basal thermal state (logical state and melt rate):
    meltratelast=0;
    istiedlast=1;
    converged=0;
    iteration=1;
    while converged==0
        
        % Solve for steady-state temperature profile:
        if d2>1 && U_lr(1,d2)>0 && U_lr(1,d2+1)>0 && doichorzadv==1
            [Temp_c(:,d2),Temp_d(d2),meltrate_thermal]=TempProfile1D_column_v4a(Temp_c(:,d2-1),zeros(zsize,1),U_lr(:,d2:d2+1),W_ud(:,d2),IsTied_d(d2),MeltRate_d(d2),DZ_c(:,d2),DZ_lr_upwind(:,d2:d2+1),dx,Gflux_c(d2)+SlideHeat_d(d2),SurfTemp_u(d2),MeltPoint_d(d2),ThermalParameters);
        else
            [Temp_c(:,d2),Temp_d(d2),meltrate_thermal]=TempProfile1D_column_v4a(NaN*zeros(zsize,1),zeros(zsize,1),U_lr(:,d2:d2+1),W_ud(:,d2),IsTied_d(d2),MeltRate_d(d2),DZ_c(:,d2),DZ_lr_upwind(:,d2:d2+1),dx,Gflux_c(d2)+SlideHeat_d(d2),SurfTemp_u(d2),MeltPoint_d(d2),ThermalParameters);
        end
        
        % Change basal logical state:
        if (IsTied_d(d2)==0 && Temp_d(d2)>=MeltPoint_d(d2)) || (IsTied_d(d2)==1 && meltrate_thermal>0) || (IsTied_d(d2)==1 && meltrate_hydro<meltrate_thermal) || IsGrounded_c(d2)==0
            IsTied_d(d2)=1;
        elseif (IsTied_d(d2)==1 && meltrate_thermal<0 && meltrate_hydro>meltrate_thermal) || (IsTied_d(d2)==0 && Temp_d(d2)<MeltPoint_d(d2))
            IsTied_d(d2)=0;
        end
        
        % Check melt rate against water supply limitation:
        if IsGrounded_c(d2)
            MeltRate_d(d2)=max([meltrate_hydro,meltrate_thermal]);
        else
            MeltRate_d(d2)=meltrate_thermal;
        end
        
        % Adjust velocity for melt/freeze:
        W_ud(:,d2)=-(MeltRate_d(d2)+Zhat_ud*(Accum_u(d2)-AnnualMelt_u(d2)-MeltRate_d(d2)));
        U_lr(:,d2+1)=U_lr(:,d2)+(dx/Icethick_lr_upwind(d2+1))*(Accum_u(d2)-AnnualMelt_u(d2)-MeltRate_d(d2));
        
        % Compute basal shear heating:
        if dostrainheat==1
            SlideHeat_d(d2)=abs(.5*(DrivingStressGrad_lr(1,d2)*Icethick_lr_center(d2)*U_lr(1,d2)+DrivingStressGrad_lr(1,d2+1)*Icethick_lr_center(d2+1)*U_lr(1,d2+1)));
        end
        
        % Check for melt rate convergance:
        denom=max([abs(meltratelast),abs(MeltRate_d(d2))]);
        misfit=abs(meltratelast-MeltRate_d(d2))/denom;
        
        % Break from loop:
        if (misfit<tolerance || denom==0) && iteration>=miniterations && istiedlast==IsTied_d(d2)
            converged=1;
        elseif iteration>=maxiterations
            error('Unable to converge on a 1D temperature distribution.')
        else
            iteration=iteration+1;
            meltratelast=MeltRate_d(d2);
            istiedlast=IsTied_d(d2);
        end
        
    end
    
    % Compute water outflux (m^2/s of water):
    WaterFlux_lrd(d2+1)=WaterFlux_lrd(d2)+(rho_i/rho_fw)*MeltRate_d(d2)*dx;
    
end


