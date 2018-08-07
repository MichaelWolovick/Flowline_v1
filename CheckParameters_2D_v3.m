% CheckParameters_2D_v1

% Mike Wolovick, 8/8/2012

% This script checks the input parameters for FTM_2D_v2 for
% inconsistencies.  It also displays many of the dependencies on the
% command line.

% v2:  If climate is time-variable, the input file must have these two
% variables:  Time_yr_input, ElevAnomaly_input.  The elevation anomalies
% are conceptualized as a vertical shift of the climate.  ie, positive
% numbers cause warm temperatures to move upwards (as though the surface
% moved down).  There is also the option for time-variability in all the
% side bc parameters.

% v3:  runs with Flowline_v1 instead of FTM_2D.  Time-variable surface
% climate can also be constant shifts in addition to elevation-dependent
% changes.

% All string inputs are tested to make sure they have one of the allowed
% values.

% I make no guarantees that all inconsistencies will be found by this
% script.

% Start command-line display:
if verbose
    disp('...')
    disp('Checking Input Parameters')
end

%% Geometry:

% Display:
if verbose
    disp('...')
    disp('Checking domain geometry:')
end

% Check basic variables loaded from input file:
if exist('X_input','var')==0 || exist('BedElev_input','var')==0 || exist('Icethick_input','var')==0
    error('Input file must contain variables named "X_input", "BedElev_input", and "Icethick_input".')
else
    if verbose
        disp('...bed geometry and ice thickness located')
    end
end

% Check for width input:
if exist('Width_input','var')==0
    if dosidedrag
        error('If input parameter "dosidedrag" is set to 1, the input file must contain a variable named "Width_input".')
    else
        if verbose
            disp('...flowband width: constant')
        end
        variablewidth=0;
        Width_input=ones(1,length(X_input));
    end
else
    % Check the size of the width input:
    if isequal(size(Width_input),size(X_input))
        if verbose
            disp('...flowband width: variable')
        end
        variablewidth=1;
    else
        error('Input variable "Width_input" must be the same size as input variable "X_input".')
    end
end

% Check for side mass input:
if usesidemassinput==1
    % Check that the input vector exists:
    if exist('Influx_yr_input','var')
        % Check that the input vector is the right size:
        if isequal(size(Influx_yr_input),size(X_input))
            % Check that the input had a variable width:
            if variablewidth==1
                if verbose
                    disp('...side mass input: enabled')
                end
                SideInflux_input=Influx_yr_input/secondsperyear;
            else
                % Explanation for this error message:
                % Side input is a flux (units: m^2/yr), with typical values
                % on the order of 10^4 or 10^5.  The side input flux is
                % divided by flowband width in the mass conservation
                % equation, producing an "equivalent SMB" (units: m/yr).
                % Real widths are of the order 10^3-10^5 m, producing
                % reasonable "equivalent SMBs".  However, when width has
                % not been specified, the model assumes that width=1m
                % everywhere.  As a result "the equivalent SMB" values will
                % be outlandishly large.
                error('If input parameter "usesidemassinput" is set to 1, the input file must contain a variable "Width_input".')
            end
        else
            error('Input variable "Influx_yr_input" must be the same size as input variable "X_input".')
        end
    else
        error('If input parameter "usesidemassinput" is set to 1, the input file must contain a variable named "Influx_yr_input".')
    end
elseif usesidemassinput==0
    if verbose
        disp('...side mass input: disabled')
    end
else
    error('Input parameter "usesidemassinput" must be equal to 0 or 1.')
end

%% Flow Model:

% Display:
if verbose
    disp('...')
    disp('Checking ice flow model:')
end

% Check velocity equation type:
if strcmp(velocitytype,'SIA')==1
    if verbose
        disp('...velocity equations: Shallow Ice Approximation')
    end
elseif strcmp(velocitytype,'SSA')==1
    if verbose
        disp('...velocity equations: Shallow Shelf Approximation')
    end
    % Check boundary slide fractions:
    if leftbcslidefrac~=1
        leftbcslidefrac=1;
        if verbose
            disp('......reset "leftbcslidefrac" to 1 for consistency with SSA equations')
        end
    end
    if rightbcslidefrac~=1
        rightbcslidefrac=1;
        if verbose
            disp('......reset "rightbcslidefrac" to 1 for consistency with SSA equations')
        end
    end
elseif strcmp(velocitytype,'horzforce')==1
    if verbose
        disp('...velocity equations: Blatter-Pattyn')
    end
elseif strcmp(velocitytype,'fullstokes')==1
    if verbose
        disp('...velocity equations: Full Stokes')
    end
else
    error('Parameter "velocitytype" must be a string equal to "SIA", "SSA", "horzforce", or "fullstokes".')
end

% Check velocity solver type:
if strcmp(solvertype,'pcg')==1
    if verbose
        disp('...velocity solver: "pcg"')
    end
elseif strcmp(solvertype,'bicg')==1
    if verbose
        disp('...velocity solver: "bicg"')
    end
elseif strcmp(solvertype,'bicgstab')==1
    if verbose
        disp('...velocity solver: "bicgstab"')
    end
elseif strcmp(solvertype,'bicgstabl')==1
    if verbose
        disp('...velocity solver: "bicgstabl"')
    end
elseif strcmp(solvertype,'cgs')==1
    if verbose
        disp('...velocity solver: "cgs"')
    end
elseif strcmp(solvertype,'gmres')==1
    if verbose
        disp('...velocity solver: "gmres"')
    end
elseif strcmp(solvertype,'lsqr')==1
    if verbose
        disp('...velocity solver: "lsqr"')
    end
elseif strcmp(solvertype,'minres')==1
    if verbose
        disp('...velocity solver: "minres"')
    end
elseif strcmp(solvertype,'qmr')==1
    if verbose
        disp('...velocity solver: "qmr"')
    end
elseif strcmp(solvertype,'symmlq')==1
    if verbose
        disp('...velocity solver: "symmlq"')
    end
elseif strcmp(solvertype,'tfqmr')==1
    if verbose
        disp('...velocity solver: "tfqmr"')
    end
elseif strcmp(solvertype,'direct')==1
    if verbose
        disp('...velocity solver: "direct"')
    end
else
    error('Sorry, Parameter "solvertype" must be one of these very specific strings: "direct", "pcg", "bicg", "bicgstab", "bicgstabl", "cgs", "gmres", "lsqr", "minres", "qmr", "symmlq", or "tfqmr"')
end

% Check side drag setting:
if dosidedrag==1
    if verbose
        disp('...side-drag: enabled')
    end
elseif dosidedrag==0
    if verbose
        disp('...side-drag: disabled')
    end
else
    error('Input parameter "dosidedrag" must be equal to 0 or 1.')
end

% Check for input rate factor:
if useinputa==1
    if exist('a_input','var')
        if verbose
            disp('...rate factor: loaded from file')
        end
    else
        error('If parameter "useinputa" is set to 1, the input file must contain a variable named "a_input".')
    end
elseif useinputa==0
    if verbose
        disp('...rate factor: constant or modeled')
    end
else
    error('Input parameter "useinputa" must be equal to 0 or 1.')
end

% Check for separate rate factor for side drag:
if useinputsidea==1
    if exist('a_side_input','var')
        if verbose
            disp('...side-drag rate factor: separate factor loaded from file')
        end
    else
        error('If parameter "useinputsidea" is set to 1, the input file must contain a variable named "a_side_input".')
    end
elseif useinputsidea==0
    if verbose
        disp('...side-drag rate factor: same factor')
    end
else
    error('Input parameter "useinputsidea" must be equal to 0 or 1.')
end

% Check grounding line interpolation style:
if strcmp(glinterp,'linear')
    if verbose
        disp('...grounding line interpolation: linear')
    end
elseif strcmp(glinterp,'cubic')
    if verbose
        disp('...grounding line interpolation: cubic')
    end
else
    error('Input parameter "glinterp" must be a string equal to "linear" or "cubic".')
end

% Check implicit ice surface:
if doimplicitsurface==1
    % Check velocity type:
    if strcmp(velocitytype,'horzforce')
        error('Parameter "doimplicitsurface" cannot equal 1 when parameter "velocitytype" equals "horzforce".')
    end
    if verbose
        disp('...ice surface evolution: implicit')
    end
elseif doimplicitsurface==0
    if verbose
        disp('...ice surface evolution: explicit')
    end
else
    error('Parameter "doimplicitsurface" must be equal to 0 or 1.')
end

% Check dynamic ice bottom:
if dodynamicbottom==1
    % Check velocity type:
    if strcmp(velocitytype,'fullstokes')==0
        error('Parameter "dodynamicbottom" can only equal 1 when parameter "velocitytype" equals "fullstokes".')
    end
    if verbose
        disp('...ice bottom evolution: dynamic contact problem')
    end
elseif dodynamicbottom==0
    if verbose
        disp('...ice bottom evolution: hydrostatic approximation')
    end
else
    error('Parameter "dodynamicbottom" must be equal to 0 or 1.')
end

%% Thermal Model:

% Check thermal model setting:
if verbose
    disp('...')
    disp('Checking thermal model:')
end
if dotemp==1
    if verbose
        disp('...thermal model: enabled')
    end
elseif dotemp==0
    if verbose
        disp('...thermal model: disabled')
    end
else
    error('Input parameter "dotemp" must be equal to 0 or 1.')
end

%% Boundary Conditions:

% Display:
if verbose
    disp('...')
    disp('Checking boundary conditions:')
end

% Check input time vector and convert to seconds:
if exist('Time_yr_input','var')
    if size(Time_yr_input,2)==1
        Time_yr_input=Time_yr_input';
    end
    Time_input=Time_yr_input*secondsperyear;
    if verbose
        disp('...time-variable forcing: enabled')
    end
else
    if verbose
        disp('...time-variable forcing: disabled')
    end
end

% Check surface climate temporal dependence:
if strcmp(climatetimedependence,'elev')
    % Check that the necessary variables exist:
    if exist('ElevAnomaly_input','var')==0 || exist('Time_yr_input','var')==0 
        error('If parameter "climatetimedependence" is set to "elev", then the input file must contain variables named "Time_yr_input", and "ElevAnomaly_input".')
    end
    % Check that the variables are the right size:
    if isequal(size(Time_yr_input),size(ElevAnomaly_input))==0 
        error('The input file variable "ElevAnomaly_input" must be the same size as "Time_yr_input".')
    end
    if verbose
        disp('...surface climate time dependence: vertical shifts')
    end
elseif strcmp(climatetimedependence,'direct')
    % Check that the necessary variables exist:
    if exist('AccumAnomaly_yr_input','var')==0 || exist('AnnualMeltAnomaly_yr_input','var')==0 || exist('Time_yr_input','var')==0 || exist('SurfTempAnomaly_input','var')==0
        error('If parameter "climatetimedependence" is set to "direct", then the input file must contain variables named "Time_yr_input", "AccumAnomaly_yr_input", "AnnualMeltAnomaly_yr_input", and "SurfTempAnomaly_input".')
    end
    % Check that the variables are the right size:
    if isequal(size(Time_yr_input),size(AccumAnomaly_yr_input))==0 || isequal(size(Time_yr_input),size(SurfTempAnomaly_input))==0 
        error('The input file variables "AccumAnomaly_yr_input" and "SurfTempAnomaly_input" must be the same size as "Time_yr_input".')
    end
    % Convert years to seconds:
    AccumAnomaly_input=AccumAnomaly_yr_input/secondsperyear;
    AnnualMeltAnomaly_input=AnnualMeltAnomaly_yr_input/secondsperyear;
    % Communicate:
    if verbose
        disp('...surface climate time dependence: uniform shifts')
    end
elseif strcmp(climatetimedependence,'mix')
    % Check that the necessary variables exist:
    if exist('AccumMultiplier_input','var')==0 || exist('ElevAnomaly_input','var')==0 || exist('Time_yr_input','var')==0 
        error('If parameter "climatetimedependence" is set to "mix", then the input file must contain variables named "Time_yr_input", "AccumMultiplier_input", and "ElevAnomaly_input".')
    end
    % Check that the variables are the right size:
    if isequal(size(Time_yr_input),size(ElevAnomaly_input))==0 || isequal(size(Time_yr_input),size(AccumMultiplier_input))==0 
        error('The input file variables "ElevAnomaly_input" and "AccumMultiplier_input" must be the same size as "Time_yr_input".')
    end
    if verbose
        disp('...surface climate time dependence: custom mixture')
    end
elseif strcmp(climatetimedependence,'none')
    if verbose
        disp('...surface climate time dependence: none')
    end
else
    error('Input parameter "climatetimedependence" must be a string equal to "elev", "direct", or "none".')
end

% Check surface climate spatial dependence:
if strcmp(climatespatialdependence,'x')
    % Check that climate variables exist:
    if exist('Accum_yr_input','var')==0 || exist('AnnualMelt_yr_input','var')==0 || exist('SurfTemp_input','var')==0
        error('If parameter "climatespatialdependence" is set to "x", then the input file must contain variables named "Accum_yr_input", "AnnualMelt_yr_input", and "SurfTemp_input".')
    end
    % Check that climate variables are the right size:
    if isequal(size(X_input),size(Accum_yr_input))==0 || isequal(size(X_input),size(AnnualMelt_yr_input))==0 || isequal(size(X_input),size(SurfTemp_input))==0
        error('If parameter "climatespatialdependence" is set to "x", then the input file variables "Accum_yr_input", "AnnualMelt_yr_input", and "SurfTemp_input" all must be the same size as "X_input".')
    end
    % Convert years to seconds:
    Accum_input=Accum_yr_input/secondsperyear;
    AnnualMelt_input=AnnualMelt_yr_input/secondsperyear;
    % Compute summer peak melt:
    PeakMelt_input=AnnualMelt_input*pi/(2*meltseasonduration);
    % Communicate:
    if verbose
        disp('...surface climate spatial dependence: position-dependent.')
    end
elseif strcmp(climatespatialdependence,'z')
    % Check that climate variables exist:
    if exist('Z_input','var')==0 || exist('Accum_yr_input','var')==0 || exist('AnnualMelt_yr_input','var')==0 || exist('SurfTemp_input','var')==0
        error('If parameter "climatespatialdependence" is set to "z", then the input file must contain variables named "Z_input", "Accum_yr_input", "AnnualMelt_yr_input", and "SurfTemp_input".')
    end
    % Check that climate variables are the right size:
    if isequal(size(Z_input),size(Accum_yr_input))==0 || isequal(size(Z_input),size(AnnualMelt_yr_input))==0 || isequal(size(Z_input),size(SurfTemp_input))==0
        error('If parameter "climatespatialdependence" is set to "z", then the input file variables "Accum_yr_input", "AnnualMelt_yr_input", and "SurfTemp_input" all must be the same size as "Z_input".')
    end
    % Convert years to seconds:
    Accum_input=Accum_yr_input/secondsperyear;
    AnnualMelt_input=AnnualMelt_yr_input/secondsperyear;
    % Compute extrapolation slopes:
    accumextrapslope=(Accum_input(2)-Accum_input(1))/(Z_input(2)-Z_input(1));
    annualmeltextrapslope=(AnnualMelt_input(2)-AnnualMelt_input(1))/(Z_input(2)-Z_input(1));
    surftempextrapslope=(SurfTemp_input(2)-SurfTemp_input(1))/(Z_input(2)-Z_input(1));
    % Communicate:
    if verbose
        disp('...surface climate spatial dependence: elevation-dependent.')
    end
else
    error('Input parameter "climatespatialdependence" must be a string equal to "x" or "z".')
end

% Left side ice BC:
% Check type:
if strcmp(leftbctype,'flux')==1
    if verbose
        disp('...upstream boundary type: prescribed flux')
    end
    influx_l=leftbcparam/secondsperyear;
elseif strcmp(leftbctype,'fixedicethick')==1
    if verbose
        disp('...upstream boundary type: prescribed thickness')
    end
elseif strcmp(leftbctype,'front')==1
    error('Sorry, ice BC type "front" is not available on the left side.')
else
    error('Parameter "leftbctype" must be a string equal to "flux" or "fixedicethick".')
end
% Check if time-variable:
if strcmp(leftbcparam,'file')
    if exist('Time_yr_input','var')==0 || exist('LeftBCParam_input','var')==0
        error('If parameter "leftbcparam" is set to "file", the input file must contain variables named "Time_yr_input" and "LeftBCParam_input".')
    end
    if size(LeftBCParam_input,1)==1
        LeftBCParam_input=LeftBCParam_input';
    end
    if verbose
        disp('......time-variable parameters')
    end
else
    if verbose
        disp('......constant parameters')
    end
end

% Right side ice BC:
if strcmp(rightbctype,'flux')==1
    if verbose
        disp('...downstream boundary type: prescribed flux')
    end
    influx_r=rightbcparam/secondsperyear;
elseif strcmp(rightbctype,'fixedicethick')==1
    if verbose
        disp('...downstream boundary type: prescribed thickness')
    end
elseif strcmp(rightbctype,'front')==1
    if verbose
        disp('...downstream boundary type: calving front condition')
    end
    % Check calving type:
    if strcmp(calvingtype,'fixed')
        if verbose
            disp('......calving law: fixed position')
        end
    elseif strcmp(calvingtype,'thick')
        if verbose
            disp('......calving law: inverse ice thickness')
        end
        refcalvingicethick=calvingparam1;
        refcalvingrate=calvingparam2/secondsperyear;
    elseif strcmp(calvingtype,'uthick')
        if verbose
            disp('......calving law: velocity/thickness')
        end
        refcalvingicethick=calvingparam1;
    elseif strcmp(calvingtype,'meltthick')
        if verbose
            disp('......calving law: melt/thickness')
        end
        refcalvingicethick=calvingparam1;
        refcalvingrate=calvingparam2/secondsperyear;
        refcalvingmeltrate=calvingparam3/secondsperyear;
    elseif strcmp(calvingtype,'vonmises')
        if verbose
            disp('......calving law: Von Mises')
        end
        icetensilestrength=calvingparam1;
    elseif strcmp(calvingtype,'meltmultiplier')
        if verbose
            disp('......calving law: melt multiplier')
        end
        meltmultiplier=calvingparam1;
    elseif strcmp(calvingtype,'file') && exist('CalvingRate_yr_input','var')
        % Check length:
        if length(CalvingRate_yr_input)~=length(X_input)
            error('Input file variables "X_input" and "CalvingRate_yr_input" must be the same length.')
        end
        % Communicate:
        if verbose
            disp('......calving law: position-dependent')
        end
        % Convert years to seconds:
        CalvingRate_input=CalvingRate_yr_input/secondsperyear;
    elseif strcmp(calvingtype,'file') && exist('CalvingRate_yr_input','var')==0
        error('If parameter "calvingtype" is set to "file", the input file must contain a variable named "CalvingRate_yr_input".')
    else
        error('Parameter "calvingtype" must be a string equal to "fixed", "thick", "uthick", "meltthick", "vonmises", "meltmultiplier", or "file".')
    end
else
    error('Parameter "rightbctype" must be a string equal to "flux", "fixedicethick", or "front".')
end
% Check if time-variable:
if strcmp(rightbcparam,'file')
    if exist('Time_yr_input','var')==0 || exist('RightBCParam_input','var')==0
        error('If parameter "rightbcparam" is set to "file", the input file must contain variables named "Time_yr_input" and "RightBCParam_input".')
    end
    if size(RightBCParam_input,1)<=3
        RightBCParam_input=RightBCParam_input';
    end
    if verbose
        disp('......time-variable parameters')
    end
else
    if verbose
        disp('......constant parameters')
    end
end

% Check for time-variable left side water BC:
if strcmp(waterinflux_l_yr,'file')
    if exist('Time_yr_input','var')==0 || exist('WaterInflux_l_yr_input','var')==0
        error('If parameter "waterinflux_l_yr" is set to "file", the input file must contain variables named "Time_yr_input" and "WaterInflux_l_yr_input".')
    end
    if size(WaterInflux_l_yr_input,1)==1
        WaterInflux_l_yr_input=WaterInflux_l_yr_input';
    end
    if verbose
        disp('...upstream water influx: time-variable')
    end
else
    waterinflux_l=waterinflux_l_yr/secondsperyear;
    if verbose
        disp('...upstream water influx: constant')
    end
end

% Check sliding velocity scale:
if ischar(slidingvelocityscale_yr) && strcmp(slidingvelocityscale_yr,'file')
    if exist('SlidingVelocityScale_yr_input','var')==0
        error('If input parameter "slidingvelocityscale_yr" is set to "file", input file must contain a variable named "SlidingVelocityScale_yr_input".')
    end
    if verbose
        disp('...sliding velocity scale: spatially variable')
    end
elseif ischar(slidingvelocityscale_yr) && strcmp(slidingvelocityscale_yr,'file')==0
    error('If input parameter "slidingvelocityscale_yr" is type "char" it must be a string equal to "file".')
else
    if verbose
        disp('...sliding velocity scale: constant')
    end
end

% Check sliding stress scale:
if ischar(slidingstressscale) && strcmp(slidingstressscale,'file')
    if exist('SlidingStressScale_input','var')==0
        error('If input parameter "slidingstressscale" is set to "file", input file must contain a variable named "SlidingStressScale_input".')
    end
    if verbose
        disp('...sliding stress scale: spatially variable')
    end
elseif ischar(slidingstressscale) && strcmp(slidingstressscale,'file')==0
    error('If input parameter "slidingstressscale" is type "char" it must be a string equal to "file".')
else
    if verbose
        disp('...sliding stress scale: constant')
    end
end

% Check sliding exponent:
if ischar(m) && strcmp(m,'file')
    if exist('M_input','var')==0
        error('If input parameter "m" is set to "file", input file must contain a variable named "M_input".')
    end
    if verbose
        disp('...sliding exponent: spatially variable')
    end
elseif ischar(m) && strcmp(m,'file')==0
    error('If input parameter "m" is type "char" it must be a string equal to "file".')
else
    if verbose
        disp('...sliding exponent: constant')
    end
end

% Check geothermal flux:
if ischar(gflux) && strcmp(gflux,'file')
    if exist('Gflux_input','var')==0
        error('If input parameter "gflux" is set to "file", input file must contain a variable named "Gflux_input".')
    end
    if verbose
        disp('...geothermal flux: spatially variable')
    end
elseif ischar(gflux) && strcmp(gflux,'file')==0
    error('If input parameter "gflux" is type "char" it must be a string equal to "file".')
else
    if verbose
        disp('...geothermal flux: constant')
    end
end


%% Plume Model:

% Display:
if verbose
    disp('...')
    disp('Checking plume model:')
end

% Check plume model setting:
if doplume==1
    % Check that the temperature/salinity inputs are present:
    if exist('Depth_input','var')==0 || exist('Temperature_input','var')==0 || exist('Salinity_input','var')==0
        error('If input parameter "doplume" is set to 1, then the input file must contain variables named "Depth_input", "Temperature_input", and "Salinity_input".')
    end
    % Check that the temperature/salinity inputs are the right size:
    if length(Depth_input)~=size(Temperature_input,1) || length(Depth_input)~=size(Salinity_input,1)
        error('Input variables "Temperature_input", and "Salinity_input" must have the same number of rows as "Depth_input".')
    end
    % Communicate:
    if verbose
        disp('...plume melt model: enabled')
    end
    % Check units:
    if max(max(Temperature_input))<100
        % Convert from celsius to kelvin and communicate:
        Temperature_input=Temperature_input+tmelt;
        if verbose
            disp('......ocean temperature converted from Celsius to Kelvin')
        end
    end
    if max(max(Salinity_input))>1
        % Convert from psu to concentration and communicate:
        Salinity_input=Salinity_input/1000;
        if verbose
            disp('......ocean salinity converted from ppt/psu to concentration')
        end
    end
    % Check plume solver type:
    if strcmp(plumesolvertype,'custom')
        if verbose
            disp('...plume model ODE solver: custom-made')
        end
    elseif strcmp(plumesolvertype,'builtin')
        if verbose
            disp('...plume model ODE solver: Matlab built-in')
        end
    else
        error('Input parameter "plumesolvertype" must be a string equal to "custom" or "builtin".')
    end
elseif doplume==0
    if verbose
        disp('...plume melt model: disabled')
    end
else
    error('Input parameter "doplume" must be equal to 0 or 1.')
end

%% Particle Tracking:

% Display:
if verbose
    disp('...')
    disp('Checking particle tracking model:')
end

% Check whether to do tracers:
if dotracers==0
    if verbose
        disp('...particle tracking: deactivated')
    end
elseif dotracers==1
    if verbose
        disp('...particle tracking: activated')
    end
    % Check tracer side BC:
    if strcmp(tracersidebc,'throughput')
        if verbose
            disp('......tracer side boundary conditions:  "throughput"')
        end
    elseif strcmp(tracersidebc,'periodic')
        error('Sorry, Tracer Side Boundary Condition "periodic" Is Not Available Yet.')
    else
        error('Parameter "tracersidebc" must be a string equal to "throughput" or "periodic".')
    end
else
    error('Parameter "dotracers" must be equal to 1 or 0.')
end

%% Other Things:

% Control the ill-conditioned warnings:
if strcmp(nearsingularwarning,'on')==1 || strcmp(nearsingularwarning,'off')==1
    warning(nearsingularwarning,'MATLAB:nearlySingularMatrix')
else
    error('Parameter "nearsingularwarning" must be a string equal to "off" or "on".')
end
if strcmp(singularwarning,'on')==1 || strcmp(singularwarning,'off')==1
    warning(singularwarning,'MATLAB:SingularMatrix')
else
    error('Parameter "singularwarning" must be a string equal to "off" or "on".')
end

% Check that the outputfile does not already exist:
test1=dir(outputfile);
if isempty(test1)==0 && overwriteoutput==0
    error('"outputfile" already exists on disc and parameter "overwriteoutput" is disabled.')
else
    if verbose
        disp('...')
        disp(['Writing output to: ',outputfile])
    end
end
clear test1


