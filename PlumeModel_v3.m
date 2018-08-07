% PlumeModel_v1

% Mike Wolovick, 1/27/2016

% This script solves a 1D buoyant plume model at the ice-ocean interface.
% It runs inside of Flowline_v1

% The model is based on Jenkins, [1991], and Jenkins, [2011].  The buoyant
% plume is described by 4 coupled ODEs.  

% Once the plume detaches from the ice, the remaining oceanic melt rates
% are calculated based on the difference between the ambient temperature
% and the freezing rate.  An assumed ocean velocity is used to get the
% constant of proportionality.

% If the ice shelf is less than 2 grid cells of the big model, the plume
% model acts as if the shelf doesn't exist and applies all of the grounding
% line flux to the vertical face.  This is because an ephemeral shelf in
% the last grid cell sometimes pops into and out of existence, and the
% ephemeral one-grid-cell-shelf has a poorly defined geometry.  

% The integration procedure can be either iterated midpoint, iterated
% endpoint, or startpoint (Euler).  However, iterated endpoint is the most
% stable and I usually have it set to that.

% This script requires ambient water temperature/salinity, ice geometry,
% conductive heat flux into the ice, and grounding line water flux.  If
% conductive heat flux is not available (if the ice sheet thermal model is
% deactivated) the plume model uses a different formulation for the melt
% rate.

% This script converts grounding line water flux into velocity with an
% assumed thickness. 

% This script operates on its own grid, and then interpolates melt rates
% back onto the ice model grid.  One grid node is always located at the
% corner between the underside of the ice shelf and the vertical ice front.
% When the plume grid size is less than the ice grid size, interpolated 
% melt rates represent the mean of all plume grid cells within the ice grid
% cell, with allowance made for fractional inclusion on the boundaries.
% When the plume grid size is larger than the ice grid size, a simple
% linear interpolation is used.

% Note that the "lr" suffix in the plume grid describes left/right edges in
% the semi-horizontal sub-shelf portion and also up/down edges in the 
% vertical cliff face portion.

% v2: The plume grid can be bunched near the grounding line.

% v3:  The model uses matlab's built-in ODE solvers instead of my own
% solver.  I've left flowband width out of the calculations.

%% Preparation:
tic
% Determine if a floating shelf is present:
if domainwidth-x_gl>2*dx
    hasshelf=1;
else
    hasshelf=0;
    % Determine if any ocean front is present at all:
    if icebottom_r>=sealevel
        oceanmeltrate_r=0;
        return
    end
end

% Create plume model horizontal grid:
X_lr_plume(1:xsize_plume+1)=linspace(x_gl,domainwidth,xsize_plume+1);
X_lr_plume(xsize_plume+2:xsize_plume+zsize_plume+1)=domainwidth;
X_c_plume=.5*(X_lr_plume(1:end-1)+X_lr_plume(2:end));
dx_plume=(domainwidth-x_gl)/xsize_plume;

% Compute plume model elevation:
if hasshelf
    % Interpolate elevation under the floating shelf:
    Z_lr_plume(1:xsize_plume+1)=interp1([x_gl,X_c(IsGrounded_c==0),domainwidth],...
        [bedelev_gl,MonotonicIceBottom_c(IsGrounded_c==0),icebottom_r],...
        X_lr_plume(1:xsize_plume+1),'linear');
else
    Z_lr_plume(1:xsize_plume+1)=icebottom_r;
end
% Compute elevation along the vertical face:
Z_lr_plume(xsize_plume+1:xsize_plume+zsize_plume+1)=linspace(icebottom_r,sealevel,zsize_plume+1);
% Interpolate to grid centers:
Z_c_plume=.5*(Z_lr_plume(1:end-1)+Z_lr_plume(2:end));

% Compute along-track distance (S) of plume model:
DS_c_plume=sqrt((X_lr_plume(2:end)-X_lr_plume(1:end-1)).^2+(Z_lr_plume(2:end)-Z_lr_plume(1:end-1)).^2);
S_lr_plume=[0,cumsum(DS_c_plume)];

% Compute pressure and pressure interpolant:
Pressure_lr_plume=rho_sw*g*(sealevel-Z_lr_plume);
PressureInterpolant=griddedInterpolant(S_lr_plume,Pressure_lr_plume,'linear');

% Compute orientation of plume model:
% Compute bottom gradient on ice model grid:
if hasshelf
    % Bottom gradient on internal cells:
    if lastgroundedind<=xsize-2
        BottomGradient_lr(lastgroundedind+2:end-1)=(MonotonicIceBottom_c(lastgroundedind+2:end)-MonotonicIceBottom_c(lastgroundedind+1:end-1))/dx;
    end
    % Special expression for bottom gradient at the grounding line:
    BottomGradient_lr(lastgroundedind+1)=max(0,(MonotonicIceBottom_c(lastgroundedind+1)-bedelev_gl)/((1-GroundedFraction_lr(lastgroundedind+1))*dx));
    % Bottom gradient at the last grid cell: (no curvature BC)
    BottomGradient_lr(end)=BottomGradient_lr(end-1);
    % Interpolate to plume model grid, convert to sintheta:
    SinTheta_lr_plume(1:xsize_plume+1)=interp1([x_gl,X_lr(lastgroundedind+2:end)],sin(atan(BottomGradient_lr(lastgroundedind+1:end))),X_lr_plume(1:xsize_plume+1));
end
% Vertical orientation:
SinTheta_lr_plume(xsize_plume+2:end)=1;
% Create sintheta interpolant:
SinThetaInterpolant=griddedInterpolant(S_lr_plume,SinTheta_lr_plume,'linear');

% Compute conductive flux into ice sheet:
CondFlux_d=-2*cond_i*(Temp_c(1,:)-Temp_d)./DZ_c(1,:);
CondFlux_r=-2*cond_i*(Temp_c(:,end)-Temp_r)/dx;

% Interpolate conductive flux onto plume model:
if dotemp
    CondFlux_lr_plume(1:xsize_plume+1)=interp1(X_c,CondFlux_d,X_lr_plume(1:xsize_plume+1),'linear',CondFlux_d(end));
    if zsize>1
        CondFlux_lr_plume(xsize_plume+2:end)=interp1(GridElev_lr(:,end),CondFlux_r,Z_lr_plume(xsize_plume+2:end),'linear',CondFlux_r(1));
    else
        CondFlux_lr_plume(xsize_plume+2:end)=CondFlux_r;
    end
    % Create interpolant for conductive flux:
    CondFluxInterpolant=griddedInterpolant(S_lr_plume,CondFlux_lr_plume,'linear');
end

% Find sill depth:
silldepth=sealevel-max(BedElev_input(X_input>x_gl)+SillThick_input(X_input>x_gl));

% Compute sill-blocked ambient water properties:
% Check if any inputs are time-dependent:
if exist('Time_input','var') && oceantimedependence
    % Check if ocean inputs are time-dependent:
    if size(Temperature_input,2)==length(Time_input)
        % Interpolate time-dependent ocean properties:
        ThisTemperature_input=interp1(Time_input,Temperature_input',time,'linear');
        ThisSalinity_input=interp1(Time_input,Salinity_input',time,'linear');
        % Extrapolate last profile forward in time:
        if time>Time_input(end)
            ThisTemperature_input=Temperature_input(:,end);
            ThisSalinity_input=Salinity_input(:,end);
        end
    else
        % Ocean properties are time-independent:
        ThisTemperature_input=Temperature_input;
        ThisSalinity_input=Salinity_input;
    end
else
    % Ocean properties are time-independent:
    ThisTemperature_input=Temperature_input(:,1);
    ThisSalinity_input=Salinity_input(:,1);
end
% Compute sill blocking:
abovesillind=find(Depth_input<=silldepth,1,'last');
ThisTemperature_input(Depth_input>silldepth)=ThisTemperature_input(abovesillind);
ThisSalinity_input(Depth_input>silldepth)=ThisSalinity_input(abovesillind);

% Interpolate ambient T/S:
AmbientTemp_lr_plume=interp1(Depth_input,ThisTemperature_input,sealevel-Z_lr_plume,'linear',ThisTemperature_input(end));
AmbientSal_lr_plume=interp1(Depth_input,ThisSalinity_input,sealevel-Z_lr_plume,'linear',ThisSalinity_input(end));

% Construct an interpolant for ambient T/S:
TempInterpolant=griddedInterpolant(S_lr_plume,AmbientTemp_lr_plume,'linear');
SalInterpolant=griddedInterpolant(S_lr_plume,AmbientSal_lr_plume,'linear');

%% Integrate ODEs:

% Assemble OtherInputs structure:
OtherInputs=struct('PlumeParameters',PlumeParameters,'ThermalParameters',ThermalParameters,'dotemp',dotemp,'g',g,'sealevel',sealevel,...
    'PressureInterpolant',PressureInterpolant,'SinThetaInterpolant',SinThetaInterpolant,'TempInterpolant',TempInterpolant,'SalInterpolant',SalInterpolant);
if dotemp
    OtherInputs.CondFluxInterpolant=CondFluxInterpolant;
end

% Initial Conditions:
% Compute subglacial inputs:
if hasshelf
    % Compute water flux at grounding line:
    waterflux_gl=WaterFlux_lrd(lastgroundedind)+(rho_i/rho_fw)*dx*GroundedFraction_lr(lastgroundedind+1)*...
        (MeltRate_d(lastgroundedind)+sum(MeltRate_c(:,lastgroundedind).*DZ_c(:,lastgroundedind))+MeltRate_u(lastgroundedind));
    % Compute melt point at grounding line:
    initialplumetemp=tmelt-meltingpointslope*Pressure_lr_plume(1);
    % Compute initial density contrast:
    initialdeltarho=-beta*AmbientSal_lr_plume(1)-alpha*(initialplumetemp-AmbientTemp_lr_plume(1));
    % Convert flux into velocity:
    initialplumeu=(-initialdeltarho*g*SinTheta_lr_plume(1)*waterflux_gl/plumedragcoeff)^(1/3);
else
    % Compute water flux at grounding line:
    waterflux_gl=WaterFlux_lrd(lastgroundedind);
    % Compute melt point at grounding line:
    initialplumetemp=tmelt-meltingpointslope*Pressure_lr_plume(xsize_plume+1);
    % Compute initial density contrast:
    initialdeltarho=-beta*AmbientSal_lr_plume(xsize_plume+1)-alpha*(initialplumetemp-AmbientTemp_lr_plume(xsize_plume+1));
    % Convert flux into velocity:
    initialplumeu=(-initialdeltarho*g*SinTheta_lr_plume(xsize_plume+1)*waterflux_gl/plumedragcoeff)^(1/3);
end
% Assign IC:
InitialValues=[waterflux_gl/initialplumeu;initialplumeu;initialplumetemp-tmelt;0];

% Create ODE options structure:
%solveroptions=odeset('Events',@(s,Values)TerminatePlumeModel_v1(s,Values,OtherInputs),'Mass',@EvaluatePlumeMatrix_v1);
solveroptions=odeset('Events',@(s,Values)TerminatePlumeModel_v1(s,Values,OtherInputs),'RelTol',1e-3,'InitialStep',100);

% Olga's settings:
% odeset('RelTol',1e-9,'AbsTol',1e-9,'InitialStep',1,'MaxStep',5e3);

% Solve ODE under floating shelf:
if hasshelf
    % Run solver:
    [~,ODE_output,plumeseparations,~,~]=ode15s(@(s,Values)EvaluatePlumeRHS_v1(s,Values,OtherInputs),S_lr_plume(1:xsize_plume+1),InitialValues,solveroptions);
    % Parse solver output:
    Thick_lr_plume(1:xsize_plume+1)=ODE_output(:,1)';
    U_lr_plume(1:xsize_plume+1)=ODE_output(:,2)';
    Temp_lr_plume(1:xsize_plume+1)=ODE_output(:,3)'+tmelt;
    Salinity_lr_plume(1:xsize_plume+1)=ODE_output(:,4)';
    % Assign new IC:
    InitialValues=[Thick_lr_plume(xsize_plume+1);U_lr_plume(xsize_plume+1);Temp_lr_plume(xsize_plume+1)-tmelt;Salinity_lr_plume(xsize_plume+1)];
end

% Solve ODE on the vertical calving front:
if hasshelf==0 || isempty(plumeseparations)
    % Run solver:
    [~,ODE_output,plumeseparations,~,~]=ode15s(@(s,Values)EvaluatePlumeRHS_v1(s,Values,OtherInputs),S_lr_plume(xsize_plume+1:end),InitialValues,solveroptions); 
    % Parse solver output:
    Thick_lr_plume(xsize_plume+1:end)=ODE_output(:,1)';
    U_lr_plume(xsize_plume+1:end)=ODE_output(:,2)';
    Temp_lr_plume(xsize_plume+1:end)=ODE_output(:,3)'+tmelt;
    Salinity_lr_plume(xsize_plume+1:end)=ODE_output(:,4)';
end

% Interpolate separation location:
if isempty(plumeseparations)==0
    plumeseparationind=find(S_lr_plume<=plumeseparations,1,'last'); % last index covered by model
    plumeseparationx=interp1(S_lr_plume,X_lr_plume,plumeseparations,'linear');
    plumeseparationz=interp1(S_lr_plume,Z_lr_plume,plumeseparations,'linear');
else
    plumeseparationind=xsize_plume+zsize_plume+1;
    plumeseparationx=domainwidth;
    plumeseparationz=sealevel;
end

% Fill in values beyond the separation point:
if isempty(plumeseparations)==0
    Thick_lr_plume(plumeseparationind+1:end)=0;
    U_lr_plume(plumeseparationind+1:end)=oceanu;
    Temp_lr_plume(plumeseparationind+1:end)=AmbientTemp_lr_plume(plumeseparationind+1:end);
    Salinity_lr_plume(plumeseparationind+1:end)=AmbientSal_lr_plume(plumeseparationind+1:end);
end

% Compute melt rates:
MeltPoint_lr_plume=tmelt-meltingpointslope*Pressure_lr_plume-salmeltcoeff*Salinity_lr_plume;
if dotemp
    MeltRate_lr_plume=(rho_sw*specheat_sw*stantonnumber.*U_lr_plume.*(Temp_lr_plume-MeltPoint_lr_plume)-CondFlux_lr_plume)/(rho_sw*latentheat);
else
    MeltRate_lr_plume=specheat_sw*stantonnumber.*U_lr_plume.*(Temp_lr_plume-MeltPoint_lr_plume)./(latentheat+specheat_i*(MeltPoint_lr_plume-consttemp));
end

% Interpolate melt rates to grid cell centers:
MeltRate_c_plume=.5*(MeltRate_lr_plume(1:end-1)+MeltRate_lr_plume(2:end));

%% Interpolate melt rate onto ice model grid:

% In order to conserve mass, the interpolation takes an average of all
% plume grid cells within the ice grid cell, including a fractional value
% for plume cells that straddle the border of the ice cell.

% The interpolation also converts between seawater-equivalent and
% ice-equivalent melt rates.

% If grounded fraction on the last grid edge is less than 0.5, then the
% last grounded grid center experiences ocean melting.  
% If grounded fraction on the last grid edge is greater than 0.5, then the
% first floating grid center experiences grounded melting.

% Interpolate floating cells:
if hasshelf
    % Check grid size:
    if dx<=dx_plume
        % Linear interpolation for fully floating cells:
        if lastgroundedind<=xsize-2
            MeltRate_d(lastgroundedind+2:end)=interp1([x_gl,X_c_plume(1:xsize_plume),domainwidth],...
                [MeltRate_c_plume(1),MeltRate_c_plume(1:xsize_plume),MeltRate_c_plume(xsize_plume)]*rho_sw/rho_i,...
                X_c(lastgroundedind+2:end),'linear');
        end
        % Weighting between grounded and ocean melt rates near the grounding line:
        if GroundedFraction_lr(lastgroundedind+1)<.5
            % Last grounded cell gets a bit of both melt rates:
            MeltRate_d(lastgroundedind)=(.5+GroundedFraction_lr(lastgroundedind+1))*MeltRate_d(lastgroundedind)+... % grounded melt rate
                (.5-GroundedFraction_lr(lastgroundedind+1))*MeltRate_c_plume(1)*rho_sw/rho_i; % floating melt rate
            % First floating cell is only the ocean melt rate:
            MeltRate_d(lastgroundedind+1)=interp1([x_gl,X_c_plume(1:xsize_plume),domainwidth],...
                [MeltRate_c_plume(1),MeltRate_c_plume(1:xsize_plume),MeltRate_c_plume(xsize_plume)]*rho_sw/rho_i,...
                X_c(lastgroundedind+1),'linear');
        else
            % First floating cell gets a bit of both melt rates:
            MeltRate_d(lastgroundedind+1)=(GroundedFraction_lr(lastgroundedind+1)-.5)*MeltRate_d(lastgroundedind+1)+... % grounded melt rate
                (1.5-GroundedFraction_lr(lastgroundedind+1))*MeltRate_c_plume(1)*rho_sw/rho_i; % floating melt rate
        end
    else
        % Down sample in a way that conserves mass:
        % Downweight grounded melt rate and figure out where to start in the ice grid:
        if GroundedFraction_lr(lastgroundedind+1)<.5
            % Downweight grounded melt rate:
            startind=lastgroundedind;
            MeltRate_d(lastgroundedind)=(.5+GroundedFraction_lr(lastgroundedind+1))*MeltRate_d(lastgroundedind);
            % Zero out ungrounded melt rates:
            MeltRate_d(lastgroundedind+1:end)=0;
        else
            % Downweight grounded melt rate:
            startind=lastgroundedind+1;
            MeltRate_d(lastgroundedind+1)=(GroundedFraction_lr(lastgroundedind+1)-.5)*MeltRate_d(lastgroundedind+1);
            % Zero out ungrounded melt rates:
            MeltRate_d(lastgroundedind+2:end)=0;
        end
        % Loop through ice grid cells:
        thisind=1;
        for d2=startind:xsize
            % Loop through plume grid cells that are within this ice cell:
            donedownsample=0;
            while donedownsample==0
                % Compute weighting for this plume cell:
                % Weighting is fraction of ice cell overlapping the plume
                % cell:
                thisweight=(min([X_lr_plume(thisind+1),X_lr(d2+1)])-max([X_lr_plume(thisind),X_lr(d2)]))/dx;
                % Add this plume cell's contribution:
                MeltRate_d(d2)=MeltRate_d(d2)+thisweight*MeltRate_c_plume(thisind)*rho_sw/rho_i;
                % Check if this ice cell is done:
                if X_lr_plume(thisind+1)==X_lr(d2+1) 
                    donedownsample=1;
                    thisind=thisind+1;
                elseif X_lr_plume(thisind+1)>X_lr(d2+1)
                    donedownsample=1;
                else
                    thisind=thisind+1;
                end
            end
        end
    end
end

% Average the cliff face cells and convert to units of ice volume:
oceanmeltrate_r=(1/(sealevel-icebottom_r))*sum(MeltRate_c_plume(xsize_plume+1:end).*DS_c_plume(xsize_plume+1:end))*rho_sw/rho_i;

toc
