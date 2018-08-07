% BasalHydrology_2D_quasistatic_v1

% Mike Wolovick, 10/23/2011

% This script does the basal hydrology for FTM_2D.  The "2D" in the name
% refers to the parent model, but the basal hydrology is a one-dimensional
% problem.  

% This script is completely different from BasalHydrology_2D_v1 and
% BasalHydrology_2D_v2

% Water flow is now quasi-static.  Water flow is balance flux in
% equilibrium with the ice model at any given timestep.  

% NOTE:  water is assumed to flow continuously from left to right.  Closed
% basins in hydraulic potential are ignored.  

% v2:  Now runs with Flowline_v1 instead of FTM_2D.  Water flux in 
% ungrounded grid cells is set to zero.  Ungrounded cells are tied to the
% melting point, and use a thermally limited melt rate without geothermal
% flux, hydro heating, or slide heating.  This melt rate can later be
% overruled by a plume model.

% The land-hydrology setting extends into the first ungrounded grid cell,
% in case the sub-grid grounding line position overlaps the grid edge.
% This is important when interpolating plume model melt rates onto the ice
% grid.

% Note: "WaterFlux_lrd" has units of m^2/s (ie, flux/width)


% Assign left edge influx:
WaterFlux_lrd(1)=waterinflux_l; % m^2/s of water

% Define hydraulic potential and its gradient:
HydroPot_c=rho_i*g*Icethick_c+rho_fw*g*IceBottom_c;
hydropot_l=rho_i*g*icethick_l+rho_fw*g*icebottom_l;
hydropot_r=rho_i*g*icethick_r+rho_fw*g*icebottom_r;
HydroPotGrad_lr=-[2*(HydroPot_c(1)-hydropot_l),(HydroPot_c(2:end)-HydroPot_c(1:end-1)),2*(hydropot_r-HydroPot_c(end))]/dx;  % this gradient points downhill

% Check for hydraulic basins:
if isequal(sign(HydroPotGrad_lr),ones(1,xsize+1))==0
    hasbasins=1;
end

% Define pressure gradient:
PressureGrad_lr=[2*(PressureStatic_d(1)-rho_i*g*icethick_l),PressureStatic_d(2:end)-PressureStatic_d(1:end-1),2*(rho_i*g*icethick_r-PressureStatic_d(end))]/dx;  % i need to be specific about side BC for water model

% Compute source term from englacial melting and surface ablation:
SourceTerm_d=(rho_i/rho_fw)*(sum(MeltRate_c.*DZ_c,1)+MeltRate_u); % m/s of water

% Integrate downstream:
for d2=1:xsize
    
    % Check if grid cell is grounded:
    if IsGrounded_c(d2) || d2==lastgroundedind+1
        % Iterate for hydraulic heating:
        converged_minor=0;
        iteration_minor=1;
        hydroheat_c_last=HydroHeat_d(d2);
        while converged_minor==0
            % Compute hydraulic heating:
            HydroHeat_lrd(d2:d2+1)=WaterFlux_lrd(d2:d2+1).*(HydroPotGrad_lr(d2:d2+1)+meltingpointslope*rho_fw*specheat_w*PressureGrad_lr(d2:d2+1)); % W/m^2
            HydroHeat_d(d2)=.5*(HydroHeat_lrd(d2)+HydroHeat_lrd(d2+1));
            % Calculate basal temperature independent of phase changes:
            Temp_d(d2)=(DZ_ud(1,d2)/cond_i)*(Gflux_c(d2)+SlideHeat_d(d2)+HydroHeat_d(d2)+cond_i*Temp_c(1,d2)/DZ_ud(1,d2));
            % Compute inflow:
            influx=WaterFlux_lrd(d2)*Width_lr(d2)+dx*Width_c(d2)*SourceTerm_d(d2); % m^3/s of water
            % Compute thermally limited melt rate:
            meltratethermal=(Gflux_c(d2)+SlideHeat_d(d2)+HydroHeat_d(d2)+cond_i*(Temp_c(1,d2)-MeltPoint_d(d2))/DZ_ud(1,d2))/(rho_i*latentheat); % m/s of ice
            % Compute hydraulically limited melt rate:
            meltratehydro=-(rho_fw/rho_i)*influx/(Width_c(d2)*dx); % m/s of ice
            % Compute actual melt rate and basal temperature:
            if Temp_d(d2)>=MeltPoint_d(d2) || (meltratehydro<=meltratethermal && meltratehydro<0 && meltratethermal<0) % thermally limited melt/freeze
                MeltRate_d(d2)=meltratethermal;
                Temp_d(d2)=MeltPoint_d(d2);
                IsTied(d2)=1;
            elseif meltratehydro>meltratethermal && meltratehydro<0 && meltratethermal<0  % hydraulically limited freezing
                MeltRate_d(d2)=meltratehydro;
                Temp_d(d2)=(DZ_ud(1,d2)/cond_i)*(Gflux_c(d2)+SlideHeat_d(d2)+HydroHeat_d(d2)-rho_i*latentheat*MeltRate_d(d2)+cond_i*Temp_c(1,d2)/DZ_ud(1,d2));
                IsTied(d2)=0;
            else
                MeltRate_d(d2)=0;
                IsTied(d2)=0;
            end
            % Compute outflow:
            WaterFlux_lrd(d2+1)=(influx+(rho_i/rho_fw)*MeltRate_d(d2)*dx*Width_c(d2))/Width_lr(d2+1); % m^2/s of water
            % Compute misfit:
            misfit_minor=abs(hydroheat_c_last-HydroHeat_d(d2))/Gflux_c(d2);
            % Break from loop:
            if misfit_minor<tolerance && iteration_minor>=miniterations
                converged_minor=1;
            elseif iteration_minor>=maxiterations
                error('MATLAB:FTM_2D:hydroheatnonconvergent','Unable to converge on a hydraulic heating solution.')
            else
                hydroheat_c_last=HydroHeat_d(d2);
                iteration_minor=iteration_minor+1;
            end
        end
        
    else % If ungrounded:
        % Tied to melting point:
        IsTied(d2)=1;
        Temp_d(d2)=MeltPoint_d(d2);
        % Thermally limited melt rate (conduction only):
        MeltRate_d(d2)=(cond_i*(Temp_c(1,d2)-MeltPoint_d(d2))/DZ_ud(1,d2))/(rho_i*latentheat); % m/s of ice
        % Outflow is zero:
        WaterFlux_lrd(d2+1)=0;
    end
    
end

