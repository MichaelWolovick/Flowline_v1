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
% endpoint, or startpoint (Euler).

% This script requires ambient water temperature/salinity, ice geometry,
% conductive heat flux into the ice, and grounding line water flux.

% This script converts grounding line water flux into velocity with an
% assumed thickness. 

% This script operates on its own grid, and then interpolates melt rates
% back onto the ice model grid.  One grid node is always located at the
% corner between the underside of the ice shelf and the vertical ice front.
% Interpolated melt rates represent the mean of all plume grid cells within
% the ice grid cell, with allowance made for fractional inclusion on the
% boundaries.

% Note that the "lr" suffix in the plume grid describes left/right edges in
% the semi-horizontal sub-shelf portion and also up/down edges in the 
% vertical cliff face portion.

%% Preparation:

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

% Compute plume model grid size:
dx_plume=(domainwidth-x_gl)/xsize_plume;

% Create plume model horizontal grid:
X_lr_plume(1:xsize_plume+1)=linspace(x_gl,domainwidth,xsize_plume+1);
X_lr_plume(xsize_plume+2:xsize_plume+zsize_plume+1)=domainwidth;
X_c_plume=.5*(X_lr_plume(1:end-1)+X_lr_plume(2:end));

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

% Compute pressure:
Pressure_c_plume=rho_sw*g*(sealevel-Z_c_plume);

% Compute along-track distance of plume model:
DS_c_plume=sqrt((X_lr_plume(2:end)-X_lr_plume(1:end-1)).^2+(Z_lr_plume(2:end)-Z_lr_plume(1:end-1)).^2);

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
    SinTheta_c_plume(1:xsize_plume)=interp1([x_gl,X_lr(lastgroundedind+2:end)],sin(atan(BottomGradient_lr(lastgroundedind+1:end))),X_c_plume(1:xsize_plume));
end
% Vertical orientation:
SinTheta_c_plume(xsize_plume+1:end)=1;

% Compute conductive flux into ice sheet:
CondFlux_d=-2*cond_i*(Temp_c(1,:)-Temp_d)./DZ_c(1,:);
CondFlux_r=-2*cond_i*(Temp_c(:,end)-Temp_r)/dx;

% Interpolate conductive flux onto plume model:
CondFlux_c_plume(1:xsize_plume)=interp1(X_c,CondFlux_d,X_c_plume(1:xsize_plume),'linear',CondFlux_d(end));
if zsize>1
    CondFlux_c_plume(xsize_plume+1:end)=interp1(GridElev_lr(:,end),CondFlux_r,Z_c_plume(xsize_plume+1:end),'linear',CondFlux_r(1));
else
    CondFlux_c_plume(xsize_plume+1:end)=CondFlux_r;
end

% Interpolate width onto plume model:
Width_lr_plume=interp1(X_lr,Width_lr,X_lr_plume);
Width_c_plume=.5*(Width_lr_plume(1:end-1)+Width_lr_plume(2:end));

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
AmbientTemp_c_plume=interp1(Depth_input,ThisTemperature_input,sealevel-Z_c_plume,'linear',ThisTemperature_input(end));
AmbientSal_c_plume=interp1(Depth_input,ThisSalinity_input,sealevel-Z_c_plume,'linear',ThisSalinity_input(end));

% Set initial conditions:
if hasshelf
    % Compute water flux at grounding line:
    waterflux_gl=WaterFlux_lrd(lastgroundedind)+(rho_i/rho_fw)*dx*GroundedFraction_lr(lastgroundedind+1)*...
        (MeltRate_d(lastgroundedind)+sum(MeltRate_c(:,lastgroundedind).*DZ_c(:,lastgroundedind))+MeltRate_u(lastgroundedind));
    % Compute melt point at grounding line:
    initialplumetemp=tmelt-meltingpointslope*Pressure_c_plume(1);
    % Compute initial density contrast:
    initialdeltarho=-beta*AmbientSal_c_plume(1)-alpha*(initialplumetemp-AmbientTemp_c_plume(1));
    % Convert flux into velocity:
    initialplumeu=(-initialdeltarho*g*SinTheta_c_plume(1)*waterflux_gl/plumedragcoeff)^(1/3);
    % Assign IC:
    Thick_lr_plume(1)=waterflux_gl/initialplumeu;
    U_lr_plume(1)=initialplumeu;
    Temp_lr_plume(1)=tmelt-meltingpointslope*Pressure_c_plume(1);
    Salinity_lr_plume(1)=0;
    % Start model at beginning:
    startind=1;
else
    % Compute water flux at grounding line:
    waterflux_gl=WaterFlux_lrd(lastgroundedind);
    % Compute melt point at grounding line:
    initialplumetemp=tmelt-meltingpointslope*Pressure_c_plume(xsize_plume+1);
    % Compute initial density contrast:
    initialdeltarho=-beta*AmbientSal_c_plume(xsize_plume+1)-alpha*(initialplumetemp-AmbientTemp_c_plume(xsize_plume+1));
    % Convert flux into velocity:
    initialplumeu=(-initialdeltarho*g*SinTheta_c_plume(xsize_plume+1)*waterflux_gl/plumedragcoeff)^(1/3);
    % Assign IC:
    Thick_lr_plume(xsize_plume+1)=waterflux_gl/initialplumeu;
    U_lr_plume(xsize_plume+1)=initialplumeu;
    Temp_lr_plume(xsize_plume+1)=tmelt-meltingpointslope*Pressure_c_plume(xsize_plume+1);
    Salinity_lr_plume(xsize_plume+1)=0;
    % Start model at vertical face:
    startind=xsize_plume+1;
end

%% Integrate ODEs:

% Loop through plume model grid cells:
doneplume=0;
meanplumeiterations=0;
for thisind=startind:xsize_plume+zsize_plume
    
    % First guess of properties: (unchanging)
    Thick_lr_plume(thisind+1)=Thick_lr_plume(thisind);
    U_lr_plume(thisind+1)=U_lr_plume(thisind);
    Temp_lr_plume(thisind+1)=Temp_lr_plume(thisind);
    Salinity_lr_plume(thisind+1)=Salinity_lr_plume(thisind);
    
    % Record first guess:
    lastthick=Thick_lr_plume(thisind);
    lastu=U_lr_plume(thisind);
    lasttemp=Temp_lr_plume(thisind);
    lastsal=Salinity_lr_plume(thisind);
    
    %     % Pre-allocate iteration record:
    %     Thick_plume_iteration=zeros(maxiterations,1);
    %     U_plume_iteration=zeros(maxiterations,1);
    %     Temp_plume_iteration=zeros(maxiterations,1);
    %     Salinity_plume_iteration=zeros(maxiterations,1);
    %     Meltrate_plume_iteration=zeros(maxiterations,1);
    
    % Iterate for convergence on properties:
    converged_plume=0;
    iteration_plume=1;
    forcerepeat=0;
    while converged_plume==0
        
        % Assign grid-center properties:
        if strcmp(plumesolvertype,'endpoint')
            Thick_c_plume(thisind)=Thick_lr_plume(thisind+1);
            U_c_plume(thisind)=U_lr_plume(thisind+1);
            Temp_c_plume(thisind)=Temp_lr_plume(thisind+1);
            Salinity_c_plume(thisind)=Salinity_lr_plume(thisind+1);
        elseif strcmp(plumesolvertype,'midpoint')
            Thick_c_plume(thisind)=.5*(Thick_lr_plume(thisind)+Thick_lr_plume(thisind+1));
            U_c_plume(thisind)=.5*(U_lr_plume(thisind)+U_lr_plume(thisind+1));
            Temp_c_plume(thisind)=.5*(Temp_lr_plume(thisind)+Temp_lr_plume(thisind+1));
            Salinity_c_plume(thisind)=.5*(Salinity_lr_plume(thisind)+Salinity_lr_plume(thisind+1));
        else % 'startpoint'
            Thick_c_plume(thisind)=Thick_lr_plume(thisind);
            U_c_plume(thisind)=U_lr_plume(thisind);
            Temp_c_plume(thisind)=Temp_lr_plume(thisind);
            Salinity_c_plume(thisind)=Salinity_lr_plume(thisind);
        end
        
        % Compute freezing temperature of plume water:
        MeltPoint_c_plume(thisind)=tmelt-meltingpointslope*Pressure_c_plume(thisind)-salmeltcoeff*Salinity_c_plume(thisind);
        
        % Compute melt rate:
        if thisind==1 && iteration_plume==1
            % First guess melt rate from ambient temp:
            if dotemp
                MeltRate_c_plume(thisind)=(rho_sw*specheat_sw*U_c_plume(thisind)*stantonnumber*(AmbientTemp_c_plume(thisind)-MeltPoint_c_plume(thisind))-CondFlux_c_plume(thisind))/(rho_sw*latentheat);
            else
                MeltRate_c_plume(thisind)=specheat_sw*U_c_plume(thisind)*stantonnumber*(AmbientTemp_c_plume(thisind)-MeltPoint_c_plume(thisind))/(latentheat+specheat_i*(MeltPoint_c_plume(thisind)-consttemp));
            end
            %MeltRate_c_plume(thisind)=0;
        else
            % Normal melt rate calculation:
            if dotemp
                MeltRate_c_plume(thisind)=(rho_sw*specheat_sw*U_c_plume(thisind)*stantonnumber*(Temp_c_plume(thisind)-MeltPoint_c_plume(thisind))-CondFlux_c_plume(thisind))/(rho_sw*latentheat);
            else
                MeltRate_c_plume(thisind)=specheat_sw*U_c_plume(thisind)*stantonnumber*(Temp_c_plume(thisind)-MeltPoint_c_plume(thisind))/(latentheat+specheat_i*(MeltPoint_c_plume(thisind)-consttemp));
            end
        end
        
        % Compute entrainment rate:
        entrainment=e0*abs(U_c_plume(thisind)*SinTheta_c_plume(thisind));
        
        % Compute normalized density contrast: (negative for buoyant plume)
        deltarho=beta*(Salinity_c_plume(thisind)-AmbientSal_c_plume(thisind))-alpha*(Temp_c_plume(thisind)-AmbientTemp_c_plume(thisind));
        
        % Check that the plume is sufficiently buoyant:
        if deltarho>=maxdeltarho
            % Plume separation:
            plumeseparationx=X_c_plume(thisind);
            plumeseparationz=Z_c_plume(thisind);
            % Quit plume model integration:
            doneplume=1;
            break
        end
        
        % Compute in-mass and out-mass:
        inmass=Width_lr_plume(thisind)*U_lr_plume(thisind)*Thick_lr_plume(thisind);
        outmass=inmass+(entrainment+MeltRate_c_plume(thisind))*DS_c_plume(thisind)*Width_c_plume(thisind);
        
        % Compute in-momentum:
        inmomentum=(U_lr_plume(thisind)^2)*Thick_lr_plume(thisind)*Width_lr_plume(thisind);
        
        % Solve quadratic equation for out-momentum:
        % note: this is an endpoint-style solver regardless of what you
        % have chosen for the rest of the problem.
        % Assign coefficients and solve equation:
        thesecoeffs=[1,(outmass^2)/(DS_c_plume(thisind)*Width_c_plume(thisind)*plumedragcoeff),(outmass^2)*(DS_c_plume(thisind)*Thick_c_plume(thisind)*Width_c_plume(thisind)*deltarho*g*SinTheta_c_plume(thisind)-inmomentum)/(DS_c_plume(thisind)*Width_c_plume(thisind)*plumedragcoeff)];
        theseroots=roots(thesecoeffs);
        % Check for positive roots:
        if sum(theseroots>0)==1
            % Assign out-momentum:
            outmomentum=theseroots(theseroots>0);
        else
            % Absence of positive roots indicates plume separation:
            plumeseparationx=X_c_plume(thisind);
            plumeseparationz=Z_c_plume(thisind);
            % Quit plume model integration:
            doneplume=1;
            break
        end
        
        % Split up outmass into thickness and velocity:
        U_lr_plume(thisind+1)=outmomentum/outmass;
        Thick_lr_plume(thisind+1)=outmass/(Width_lr_plume(thisind+1)*U_lr_plume(thisind+1));
        
        % Compute in-heat and out-heat:
        inheat=Width_lr_plume(thisind)*U_lr_plume(thisind)*Thick_lr_plume(thisind)*Temp_lr_plume(thisind);
        outheat=inheat+DS_c_plume(thisind)*Width_c_plume(thisind)*(AmbientTemp_c_plume(thisind)*entrainment+MeltPoint_c_plume(thisind)*MeltRate_c_plume(thisind)+stantonnumber*U_c_plume(thisind)*(MeltPoint_c_plume(thisind)-Temp_c_plume(thisind)));
        
        %         % Compute out-heat: (always endpoint-style, solves for melt rate too)
        %         outheat=(ambienttemp*entrainment-((MeltPoint_c_plume(thisind)^2)*specheat_sw*U_c_plume(thisind)*stantonnumber)/(latentheat+specheat_i*(MeltPoint_c_plume(thisind)-consttemp))...
        %             +MeltPoint_c_plume(thisind)*stantonnumber*U_c_plume(thisind)+inheat/(DS_c_plume(thisind)*Width_c_plume(thisind)))/...
        %             (1/(DS_c_plume(thisind)*Width_c_plume(thisind))+stantonnumber*U_c_plume(thisind)/outmass-(MeltPoint_c_plume(thisind)*specheat_sw*stantonnumber*U_c_plume(thisind))/(outmass*(latentheat+specheat_i*(MeltPoint_c_plume(thisind)-consttemp))));
        
        % Compute temperature from out-heat:
        Temp_lr_plume(thisind+1)=outheat/outmass;
        
        % Truncate temperature at the melting point:
        if Temp_lr_plume(thisind+1)<MeltPoint_c_plume(thisind)
            Temp_lr_plume(thisind+1)=MeltPoint_c_plume(thisind);
            forcerepeat=1;
        end
        
        % Compute in-salt and out-salt:
        insalt=Width_lr_plume(thisind)*U_lr_plume(thisind)*Thick_lr_plume(thisind)*Salinity_lr_plume(thisind);
        outsalt=insalt+DS_c_plume(thisind)*Width_c_plume(thisind)*(AmbientSal_c_plume(thisind)*entrainment);
        
        % Compute salinity from out-salt:
        Salinity_lr_plume(thisind+1)=outsalt/outmass;
        
        %         % Record iteration values:
        %         Thick_plume_iteration(iteration_plume)=Thick_lr_plume(thisind+1);
        %         U_plume_iteration(iteration_plume)=U_lr_plume(thisind+1);
        %         Temp_plume_iteration(iteration_plume)=Temp_lr_plume(thisind+1);
        %         Salinity_plume_iteration(iteration_plume)=Salinity_lr_plume(thisind+1);
        %         Meltrate_plume_iteration(iteration_plume)=MeltRate_c_plume(thisind);
        
        % Compute misfit:
        misfit_plume=max([abs(Thick_lr_plume(thisind+1)-lastthick)/plumethickscale,abs(U_lr_plume(thisind+1)-lastu)/plumeuscale,...
            abs(Temp_lr_plume(thisind+1)-lasttemp)/plumetempscale,abs(Salinity_lr_plume(thisind+1)-lastsal)/plumesalscale]);
        
        % Break from loop:
        if (misfit_plume<tolerance && iteration_plume>=miniterations_plume && forcerepeat==0 && Thick_lr_plume(thisind+1)~=0 && U_lr_plume(thisind+1)~=0 && Temp_lr_plume(thisind+1)~=0 && Salinity_lr_plume(thisind+1)~=0) || strcmp(plumesolvertype,'startpoint')
            converged_plume=1;
            meanplumeiterations=meanplumeiterations+iteration_plume;
        elseif iteration_plume>=maxiterations_plume
            error(['Unable to converge on a plume model solution in plume grid cell number ',num2str(thisind),'.'])
        else
             % Compute this damping factor:
            if thisind<=xsize_plume 
                thisdamping=damping_plume1+(damping_plume0-damping_plume1)*exp(-(thisind-1)/efoldingiteration_plume);
            else
                thisdamping=damping_plume1+(damping_plume0-damping_plume1)*exp(-(thisind-xsize_plume-1)/efoldingiteration_plume);
            end
            % Apply iteration damping:
            Thick_lr_plume(thisind+1)=lastthick+(1-thisdamping)*(Thick_lr_plume(thisind+1)-lastthick);
            U_lr_plume(thisind+1)=lastu+(1-thisdamping)*(U_lr_plume(thisind+1)-lastu);
            Temp_lr_plume(thisind+1)=lasttemp+(1-thisdamping)*(Temp_lr_plume(thisind+1)-lasttemp);
            Salinity_lr_plume(thisind+1)=lastsal+(1-thisdamping)*(Salinity_lr_plume(thisind+1)-lastsal);
            % Count iterations:
            iteration_plume=iteration_plume+1;
            % Record last guess:
            lastthick=Thick_lr_plume(thisind+1);
            lastu=U_lr_plume(thisind+1);
            lasttemp=Temp_lr_plume(thisind+1);
            lastsal=Salinity_lr_plume(thisind+1);
            % Reset forcerepeat variable:
            forcerepeat=0;
        end
    end
    
    % Check whether to quit plume integration:
    if doneplume
        break
    end
    
end

% Finish computing mean number of plume iterations:
meanplumeiterations=meanplumeiterations/(thisind-startind);

% Assign separation location if it hasn't been assigned already:
if thisind==xsize_plume+zsize_plume
    plumeseparationx=domainwidth;
    plumeseparationz=sealevel;
end

% Do simple melt rates where the plume model wasn't used:
% Reset grid-center thickness to zero:
Thick_c_plume(thisind:end)=0;
% Set velocity to oceanic background value:
U_c_plume(thisind:end)=oceanu;
% Set temperature and salinity to ambient:
Temp_c_plume(thisind:end)=interp1(Depth_input,Temperature_input,sealevel-Z_c_plume(thisind:end),'linear',Temperature_input(end));
Salinity_c_plume(thisind:end)=interp1(Depth_input,Salinity_input,sealevel-Z_c_plume(thisind:end),'linear',Temperature_input(end));
% Compute melt point:
MeltPoint_c_plume(thisind:end)=tmelt-meltingpointslope*Pressure_c_plume(thisind:end)-salmeltcoeff*Salinity_c_plume(thisind:end);
% Compute melt rate:
MeltRate_c_plume(thisind:end)=(rho_sw*specheat_sw*stantonnumber.*U_c_plume(thisind:end).*abs(SinTheta_c_plume(thisind:end)).*(Temp_c_plume(thisind:end)-MeltPoint_c_plume(thisind:end))-CondFlux_c_plume(thisind:end))/(rho_sw*latentheat);

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
oceanmeltrate_r=mean(MeltRate_c_plume(xsize_plume+1:end))*rho_sw/rho_i;
