% Flowline_v1

% Mike Wolovick, 1/19/2016

% Based on FTM_2D_v2b

%% Functional Forms:

% There are various ways to use this model as a function.  Uncomment the
% appropriate expression.

% Basic function:
%function Flowline_v1(inputfile,outputsuffix,climatespatialdependence,climatetimedependence,oceantimedependence,leftbcparam,calvingtype,calvingparam1,calvingparam2,calvingparam3,stantonnumber,dosill,silltopz,usesidemassinput,m,sillstressscale,sillerosionparam,sillblockingfraction)

% Plume model calibration:
%function [initialslope,finalslope]=Flowline_v1(inputfile,climatespatialdependence,leftbcparam,stantonnumber,usesidemassinput,plumecalibrationtime_yr)

%% Basic Description:

% This script does thermomechanical ice flow on a 2D (x,z) flowband domain
% with prescribed cross-flow width.  The model can be coupled to a basal 
% hydrology model where the ice is grounded and to a buoyant plume model 
% where the ice is floating.

% The buoyant plume model solves the Jenkins [1991,2011] ODEs.

% Ice flow can be the shallow shelf approximation (SSA), the Blatter-Pattyn
% approximation, or full stokes.

% In full stokes mode, you have the option to either solve the full contact
% problem for grounding line motion, or to use the hydrostatic
% approximation.  The other two flow modes use the hydrostatic
% approximation.

% Grounding line motion includes an interpolated subgrid position
% for partially grounded cells.  The subgrid position is used to compute a
% grounding fraction for use with the basal drag coefficient.
% Interpolation can be linear or cubic.

% A moving calving front is accomodated by dynamic regridding, so that the
% ice front is always located at the edge of the model domain.

% Ice temperature is the standard advection-diffusion equation with upwind
% differencing for advection.  If temperature exceeds the melt point an
% enthalpy formulation is used (ie, englacial melt is produced).  All
% englacial melt immediately drains to the bed; the vertical subsidence due
% to englacial melting does not affect temperature advection.

% Basal hydrology is quasistatic.  Quasistatic neglects storage changes and
% assumes that water flux is in equilibrium with the thermal forcing
% provided by the ice sheet in any timestep.  The dynamic option that was
% present in FTM_2D has been removed and the WaterThick variable now
% represents the ocean cavity.

% The model can include a buoyant plume underneath the floating shelf or at
% the ice front.

% The model has the option for inflow from the sides (across the lateral
% boundaries of the flowband).  The side influx is prescribed as a function
% of X and adjusted in response to surface elevation changes.

% Time integration is explicit, with the option for implicit ice dynamics.
% If you choose implicit ice dynamics, then the linear velocity solve 
% includes a correction for the change in driving stress at the end of the 
% time step.  Note that this option is not available in the Blatter-Pattyn
% ("horzforce") approximation, because it destroys the computational
% benefit of Blatter-Pattyn.

% Spatial gridding within the ice is a skewed grid (along lines of constant
% dimensionless elevation) and horizontal velocity is assumed parallel to
% this grid.  In other words, the horizontal velocity is contructed to
% satisfy a no-penetration condition at the bed and the surface.  All
% violations of the no-penetration condition (surface uplift/subsidence,
% floating shelf uplift/subsidence, basal melt/freeze, and surface 
% accumulation/ablation) act through the vertical velocity only.

% Transport velocity is defined on grid edges.  This means that the 
% velocity solvers operate on a staggered grid.

% Stratigraphy is modeled with particle tracking.  The model keeps track
% of the connectivity of the tracers in coherent layers.

% Model output is saved at regular time intervals, even if the timestep can
% be variable within the model.  Output records represent the mean of a
% variable within that time period, except for stratigraphy, which is a 
% snapshot at the record middle.  Note that this makes it impossible to
% exactly restart an old run at an arbitrary point in time.  The model
% saves the exact final state, but restarting from any arbitrary model
% record means restarting from a temporal average.

% The above limitation notwithstanding, the model has the ability to
% restart from any arbitrary time in a previous model run.  It also has the
% ability to change any time-dependent inputs after restarting; this is
% useful for multiple model experiments with a common spinup.

% Many of the features described above can be enabled or disabled depending
% on the problem you want to solve.

%% Model Conservation Equations:

% In ice:
% mass
% energy 
% momentum

% In water: (basal hydrology)
% mass
% energy

% In water:  (buoyant plume)
% mass
% momentum
% heat
% salt

%% Model Conditions:

% Shape of model domain:
% bed elevation and initial ice thickness loaded from file
% surface (and possibly bottom) elevation evolves as a free surface

% Top Boundary Conditions:
% prescribed temperature 
% prescribed accumulation rate 
% no stress (shear or normal)

% Annual cycle:  the model starts in the middle of winter.  The peak of
% summer melt and the peak of surface temperature are assumed to coincide
% with one another, at the halfway point of each year.  Surface temperature
% varies on an annual sinusoid and ablation rate has a specified melt
% season duration.

% Bottom boundary conditions:
% geothermal flux 
% power-law sliding
% temperature<=melting point
% ocean pressure
% basal hydrology or ocean plume melt (1D models coupled to the main model)

% Side Boundary Conditions:
% Temperature: 
%      1D steady state solution when boundary is ice, melt point when
%      boundary is calving front
% Velocity and Ice Thickness: (multiple options)
%   1. prescribed flux (dirichlet for velocity, neumann for icethick)
%   2. prescribed icethick (neumann for velocity, dirichlet for icethick)
%   3. ice front (neumann for velocity, plus calving law)
% Basal Hydrology:
%      prescribed influx (left side only)

% Calving Law Options:
%   1. fixed front (calving rate matches ice velocity)
%   2. thickness-dependent (rate inversely proportional to ice thickness)
%   3. thickness and velocity dependent (rate scales with ice velocity,
%      inversely with thickness)
%   4. thickness and melt dependent (rate scales with frontal melt rate,
%      inversely with thickness)
%   5. von mises (rate scales with tensile stress and ice velocity)
%   6. melt multiplier (rate scales with frontal melt rate)
%   7. position-dependent (calving rate function loaded from input file) 

% Global inequality condition:
% temperature <= pressure-dependent melting point (PMP)
%   If temperature rises above the PMP the extra thermal energy is devoted
%   to volumetric melting instead (ie, we resort to an implicit enthalpy 
%   formulation).  All volumetric melting immediately drains to the bed.

% Initial conditions:
% Temperature from an approximate steady state solution
% Ice thickness from file
% Stratigraphy flat (parallel to dimensionless elevation)

%% Variable Nomenclature:

% Scalars are lowercase, while vectors, matrices, and structures are uppercase.

% Variables defined on the grid have the following suffix scheme:
% Name_[c/lr/ud/lrud]
% where  c    =grid centers 
%        lr   =grid left/right edges
%        ud   =grid up/down edges
%        lrud =grid corners 

% Variables that have a distinction between upwind and centered values
% have the additional suffix:
% _lr_[center/upwind]
% where center and upwind have the obvious meanings.  Note that this
% distinction only applies in assignments to grid cell edges.

% Variables that have upwind values:
% Icethick    for ice thickness advection
% DZ          for consistency with ice thickness
% Temp        for temperature advection

% Boundary conditions have the following suffix scheme:
% Name_[l/r/u/d]
% where   l =left boundary
%         r =right boundary
%         u =up boundary
%         d =down boundary
% note that "Name" will be lowercase if the boundary condition is
% scalar-valued.

% Tracers have the suffix:
% _t

% Variables that are loaded from the input file have the suffix:
% _input

% Variables that are saved from the last iteration have the suffix:
% _last


%% Input File Format:

% The input file must contain at least the following variables: 
% X_input               m
% BedElev_input         m
% Icethick_input        m
% Accum_yr_input        m/yr
% AnnualMelt_yr_input   m/yr
% SurfTemp_input        K

% These define the model geometry and surface climate.  The model 
% automatically figures out which grid cells are floating based on the 
% hydrostatic condition.

% If you want to have a buoyant plume, the input file must also contain these
% variables:
% Depth_input         m
% Temp_input          K
% Salinity_input      unitless

% These variables describe the ambient ocean conditions that act as a
% forcing on the buoyant plume model.  

% In addition, the input file may also contain:
% Time_yr_input                yr
% Z_input                      m
% LeftBCParam_input            variable units
% RightBCParam_input           variable units
% WaterInflux_l_yr_input       m^2/yr
% ElevAnomaly_input            m
% AccumAnomaly_yr_input        m/yr
% AnnualMeltAnomaly_yr_input   m/yr
% SurfTempAnomaly_input        K
% Gflux_input                  W/m^2

% These define time-variable or spatially-variable boundary conditions.

% Icethick_input can have NaN's.  In that case, the last non-NaN element of
% Icethick_input defines the starting size of the model domain.  This is
% useful in the moving front side BC, so that bed topography can be
% specified in places the ice has not yet reached.  

% The ice front cannot advance beyond the last element of X_input.  I
% assume that the input flowbands stop where the ice would spread out into
% the ocean and form an unconfined shelf.

% Inputs are not filtered before interpolation onto the model grid is
% performed.  If they contain a lot of high-frequency content, aliasing 
% will occur and the model may produce unexpected errors.

% Input variables that are interpolated in time use the last value to
% extrapolate to out-of-range times (or the first value, if time=0 is out
% of range).  

% If resuming an old run, "inputfile{1}" should refer to a previous model
% output.  "inputfile" is a cell array of strings, so you can add an
% additional file with time-variable inputs different from the old run you
% are resuming.  inputfile{2} is only used for time-variable forcings.

%% Start Model:

% Clear the workspace, start the timer:
%clearvars
%t1=tic;

% Loop through model runs:
%for runnum=2

% Start another timer:
t2=tic;
    
% Run inside of a try/catch pair:
%try

%% Parameters:

% Are you starting a new run or resuming an old run?
resumeoldrun=0;                   % logical (all parameters except "inputfile", "outputfile", "resumetime_yr", "runtime_yr", "shutoffhour", and "restarthour" are ignored if this equals 1)
resumetime_yr=[];                 % yr (if this is empty, the model resumes at the end of the old run)

% Input file:
%inputfile={['/work/Michael.Wolovick/FjordSillPaper/ModelOutput/SchoofTest_v3_run',num2str(runnum),'.mat']};  % cell array of strings
%inputfile={'/net/mjw/FjordSillPaper/ModelOutput/ThwaitesFlowband_v1_Warming_v2.mat'};
%inputfile={'/work/Michael.Wolovick/FjordSillPaper/ModelOutput/PetermannFlowband_v2_Sill_v3.mat'};
%inputfile={'/home/mjw/Documents/FjordSillPaper/ModelInput/ThwaitesA_v2.mat'};
%inputfile={'/home/mike/Documents/Research/FjordSillPaper/ModelInput/TestInput_v1.mat'};
inputfile={'/home/wolovick/Dropbox/FjordSillPaper/ModelInput/ThwaitesA_v2.mat'};

% Output file:
%outputsuffix='_Warming_v3';
%outputfile=inputfile{1};
%outputfile=['/work/Michael.Wolovick/FjordSillPaper/ModelOutput/',inputfile{1}(47:end-4),outputsuffix,'.mat'];
%outputfile=['/work/Michael.Wolovick/FjordSillPaper/ModelOutput/ClimateExperiment/',inputfile{1}(71:end-4),outputsuffix,'.mat'];
%outputfile=['/work/Michael.Wolovick/FjordSillPaper/ModelOutput/SchoofTest_v3_run',num2str(runnum),'.mat'];
%outputfile=['/home/mjw/Documents/FjordSillPaper/ModelOutput/SchoofTest_v2_run',num2str(runnum),'.mat'];
%outputfile='/net/mjw/FjordSillPaper/ModelOutput/Test1.mat';
%outputfile='/work/Michael.Wolovick/FjordSillPaper/ModelOutput/Test1.mat';
%outputfile='/home/mike/Documents/Research/FjordSillPaper/ModelOutput/Test1.mat';
outputfile='/home/wolovick/Dropbox/FjordSillPaper/ModelOutput/Test1.mat';

% Overall Model Controls:
dotemp=0;                         % 0 or 1
dosidedrag=1;                     % 0 or 1
useinputsidea=1;                  % 0 or 1
useinputa=1;                      % 0 or 1
usesidemassinput=0;               % 0 or 1
doplume=1;                        % 0 or 1
velocitytype='SSA';               % 'SSA', 'horzforce', or 'fullstokes' 
glinterp='cubic';                 % 'linear' or 'cubic'
doimplicitsurface=0;              % 0 or 1
dodynamicbottom=0;                % 0 or 1
dotracers=0;                      % 0 or 1
doenglacialmelting=1;             % 0 or 1
doichorzadv=1;                    % 0 or 1
dostrainheat=1;                   % 0 or 1
maxicethickcurvature=1e-3;        % 1/m
interpstyle='pchip';              % valid interp1 method
ModelParameters=struct('dotemp',dotemp,'dosidedrag',dosidedrag,'usesidea',useinputsidea,'useinputa',useinputa,'usesidemassinput',usesidemassinput,'doplume',doplume,'velocitytype',velocitytype,'glinterp',glinterp,'doimplicitsurface',doimplicitsurface,'dotracers',dotracers,'doenglacialmelting',doenglacialmelting,'doichorzadv',doichorzadv,'dostrainheat',dostrainheat,'maxicethickcurvature',maxicethickcurvature,'interpstyle',interpstyle);

% Sill building parameters:
dosill=0;                         % logical
sillblockingfraction=.5;          % unitless [0,1]
dosillerosion=0;                  % logical
dosilldrag=0;                     % logical
sillstressscale=1e4;              % Pa
sillerosionparam=1e-4;            % unitless
silldragthickscale=1;             % m
sillcentx=423500;                 % m or 'auto' 
sillwidth=7.5e3;                  % m
silltopz=-250;                    % m (<0)
sillstarttime_yr=100;             % yr
sillconstructiontime_yr=10;       % yr
silldistbuffer=1.5;               % unitless
SillParameters=struct('dosill',dosill,'sillblockingfraction',sillblockingfraction,'dosillerosion',dosillerosion,'dosilldrag',dosilldrag,'sillstressscale',sillstressscale,'sillerosionparam',sillerosionparam,'silldragthickscale',silldragthickscale,'sillcentx',sillcentx,'sillwidth',sillwidth,'silltopz',silltopz,'sillstarttime_yr',sillstarttime_yr,'sillconstructiontime_yr',sillconstructiontime_yr,'silldistbuffer',silldistbuffer,'sillamp',[]);

% Top Boundary Parameters:
climatespatialdependence='x';     % 'x' or 'z'
climatetimedependence='mix';      % 'elev', 'direct', 'mix', or 'none'
elevanomtemplapserate=7e-3;       % K/m >0 (used if climatetimedependence='elev' but climatespatialdependence='x')
elevanomaccumlapserate_yr=1e-4;   % (m/yr)/m >0 (used if climatetimedependence='elev' but climatespatialdependence='x')
elevanomablatelapserate_yr=1e-3;  % (m/yr)/m >0 (used if climatetimedependence='elev' but climatespatialdependence='x')
meltseasonduration=1/3;           % fraction of a year
seasonaltempcycle=15;             % K (peak-mean amplitude)
TopParameters=struct('climatespatialdependence',climatespatialdependence,'climatetimedependence',climatetimedependence,'elevanomtemplapserate',elevanomtemplapserate,'elevanomaccumlapserate_yr',elevanomaccumlapserate_yr,'elevanomablatelapserate_yr',elevanomablatelapserate_yr,'meltseasonduration',meltseasonduration,'seasonaltempcycle',seasonaltempcycle);

% Bottom Boundary Parameters:  
gflux='file';                     % W/m^2 or 'file'
m=1;                              % unitless>0 or 'file'
slidingstressscale='file';        % Pa, or 'file'
slidingvelocityscale_yr='file';   % m/yr, or 'file'
slidingtempscale=1;               % K (set to zero to deactivate subfreezing slip)
BottomParameters=struct('gflux',gflux,'m',m,'slidingstressscale',slidingstressscale,'slidingvelocityscale_yr',slidingvelocityscale_yr,'slidingtempscale',slidingtempscale);

% Side Boundary Parameters:
leftbctype='flux';                % 'flux', or 'fixedicethick'
leftbcparam=2.14e4;               % m^2/yr (pos for influx), m (icethick), or 'file' (same units in file)
leftbcslidefrac=1;                % unitless [0,1]
waterinflux_l_yr=0;               % m^2/yr, or 'file' (same units in file)
rightbctype='front';              % 'flux', 'fixedicethick', or 'front'
rightbcparam=0;                   % m^2/yr (pos for influx), m (icethick), m (sea level), or 'file' (same units in file)
rightbcslidefrac=1;               % unitless [0,1]
SideParameters=struct('leftbctype',leftbctype,'leftbcparam',leftbcparam,'leftbcslidefrac',leftbcslidefrac,'waterinflux_l_yr',waterinflux_l_yr,'icethickadjustment_l',[],'rightbctype',rightbctype,'rightbcparam',rightbcparam,'rightbcslidefrac',rightbcslidefrac,'icethickadjustment_r',[]);

% Calving Parameters:
calvingtype='fixed';              % 'fixed', 'thick', 'uthick', 'meltthick', 'vonmises', 'meltmultiplier', or 'file'
calvingparam1=300;                % refthick, m (thick, uthick, meltthick); yieldstress, Pa (vonmises); unitless (meltmultiplier)
calvingparam2=2082;               % refcalvingrate, m/yr (thick, meltthick)
calvingparam3=52;                 % refmeltrate, m/yr (meltthick)
CalvingParameters=struct('calvingtype',calvingtype,'calvingparam1',calvingparam1,'calvingparam2',calvingparam2,'calvingparam3',calvingparam3);

% Plume Model Parameters:
oceantimedependence=1;            % 0 or 1
rho_sw=1028;                      % kg/m^3
e0=.036;                          % unitless (Jenkins 2011 value: .036)
plumedragcoeff=2.5e-3;            % unitless
stantonnumber=3.25e-4;            % unitless (Jenkins 2011 value: 5.9e-4)
salmeltcoeff=57.3;                % K (>0)
specheat_sw=3.974e3;              % J/(kg*K)
beta=7.86e-1;                     % unitless
alpha=3.87e-5;                    % 1/K
u_tidal=.1;                       % m/s 
plumethickscale=10;               % m
plumeuscale=1;                    % m/s
plumetempscale=4;                 % K
plumesalscale=.03;                % unitless
plumesolvertype='custom';         % 'custom' or 'builtin'
maxdeltarho=1e-6;                 % unitless (pos indicates plume denser than ambient)
minsintheta=-1e-4;                % unitless
minbotgrad=2e-3;                  % unitless>0
PlumeParameters=struct('oceantimedependence',oceantimedependence,'rho_sw',rho_sw,'e0',e0,'plumedragcoeff',plumedragcoeff,'stantonnumber',stantonnumber,'salmeltcoeff',salmeltcoeff,'specheat_sw',specheat_sw,'beta',beta,'alpha',alpha,'u_tidal',u_tidal,'plumethickscale',plumethickscale,'plumeuscale',plumeuscale,'plumetempscale',plumetempscale,'plumesalscale',plumesalscale,'maxdeltarho',maxdeltarho,'minsintheta',minsintheta,'minbotgrad',minbotgrad);

% Thermal Parameters:
consttemp=268;                    % K (used if dotemp=0 and useinputa=0)
constbasalmelt_yr=1e-3;           % m/yr (used if dotemp=0 or if waterflux_gl=0 in plume model)
cond_i=2.4;                       % W/(K*m)
rho_i=917;                        % kg/m^3
rho_fw=1000;                      % kg/m^3
latentheat=3.335e5;               % J/kg
specheat_i=1.9e3;                 % J/(kg*K)
specheat_w=4.2e3;                 % J/(kg*K)
tmelt=273.15;                     % K
meltingpointslope=7.4e-8;         % K/Pa
ThermalParameters=struct('consttemp',consttemp,'cond_i',cond_i,'rho_i',rho_i,'rho_fw',rho_fw,'specheat_i',specheat_i,'specheat_w',specheat_w,'latentheat',latentheat','tmelt',tmelt,'meltingpointslope',meltingpointslope);

% Rheological Parameters:
n=3;                              % unitless
a0=4.9e-25;                       % Pa^-n*s^-1
t0=263;                           % K
q_big=1.39e5;                     % J/mol
q_small=6e4;                      % J/mol
RheologyParameters=struct('n',n,'a0',a0,'t0',t0,'q_big',q_big,'q_small',q_small);

% Grid Parameters:
zsize=1;                          % integer (number of grid centers)
targetdx=500;                     % m
densifypower=1;                   % unitless
maxdensify=1;                     % unitless (1,inf) (too much can produce instability in SteadyStateIC_2D)
minicethick=1;                    % m
targetdx_plume=500;               % m
targetdz_plume=3;                 % m
minxsize_plume=50;                % integer
minzsize_plume=10;                % integer
plumemaxdensify=1;                % unitless
plumedensifypower=1;              % unitless
GridParameters=struct('zsize',zsize,'xsize',[],'targetdx',targetdx,'densifypower',densifypower,'maxdensify',maxdensify,'minicethick',minicethick,'targetdx_plume',targetdx_plume,'targetdz_plume',targetdz_plume,'minxsize_plume',minxsize_plume,'minzsize_plume',minzsize_plume,'xsize_plume',[],'zsize_plume',[],'plumemaxdensify',plumemaxdensify,'plumedensifypower',plumedensifypower);

% Timing Parameters:
dt_yr=.02;                        % yr
runtime_yr=1e3;                   % yr
recordinterval_yr=1;              % yr
spinuptime_yr=0;                  % yr
rapidfactor=10;                   % unitless
TimingParameters=struct('dt_yr',dt_yr,'runtime_yr',runtime_yr,'recordinterval_yr',recordinterval_yr,'spinuptime_yr',spinuptime_yr,'rapidfactor',rapidfactor);

% Iteration Parameters:
tolerance=1e-3;                   % unitless
damping_visc=.75;                 % unitless [0,1)
damping_drag=.5;                  % unitless [0,1)
damping_plume0=.9;                % unitless [0,1) (first grid cell)
damping_plume1=.25;               % unitless [0,1) (late grid cells)
efoldingiteration_plume=7.5;      % # iterations (can be noninteger) (was 7.5)
miniterations=2;                  % integer
miniterations_plume=1;            % integer
maxiterations=100;                % integer
maxiterations_plume=300;          % integer
IterationParameters=struct('tolerance',tolerance,'damping_visc',damping_visc,'damping_drag',damping_drag,'damping_plume0',damping_plume0,'damping_plume1',damping_plume1,'efoldingiteration_plume',efoldingiteration_plume,'miniterations',miniterations,'miniterations_plume',miniterations_plume,'maxiterations',maxiterations,'maxiterations_plume',maxiterations_plume);

% Velocity Solver Parameters:
solvertype='direct';              % 'direct', 'pcg', 'bicg', 'bicgstab', 'bicgstabl', 'cgs', 'gmres', 'lsqr', 'minres', 'qmr', 'symmlq', 'tfqmr'  (names of matlab iterative solvers)
viscosityguess=1e15;              % Pa*s
dragcoefficientguess=1e10;        % Pa/(m/s)
sidedragcoefficientguess=1e10;    % Pa/(m/s)
strainratestabilizer_yr=1e-6;     % 1/yr
slidingstabilizer_yr=1e-3;        % m/yr
singularwarning='on';             % 'on' or 'off'
nearsingularwarning='off';        % 'on' or 'off'
VelocitySolverParameters=struct('solvertype',solvertype,'viscosityguess',viscosityguess,'dragcoefficientguess',dragcoefficientguess,'sidedragcoefficientguess',sidedragcoefficientguess,'strainratestabilizer_yr',strainratestabilizer_yr,'slidingstabilizer_yr',slidingstabilizer_yr,'singularwarning',singularwarning,'nearsingularwarning',nearsingularwarning);

% Tracer Parameters:
tracerzsize=5;                    % integer
tracerxsize=2000;                 % integer
spawninterval_u_yr=4000;          % yr
spawninterval_d_yr=2000;          % yr
spawninterval_l_yr=25;            % yr
tracersidebc='throughput';        % 'throughput' or 'periodic'
efzheight=.2;                     % unitless (0,1) (ignored in periodic case)
icacczhat=.1;                     % unitless (0,1)
interpstyle_t='linear';           % valid interp2 method
storagebuffer=10;                 % integer
TracerParameters=struct('tracerzsize',tracerzsize,'tracerxsize',tracerxsize,'spawninterval_u_yr',spawninterval_u_yr,'spawninterval_d_yr',spawninterval_d_yr,'spawninterval_l_yr',spawninterval_l_yr,'tracersizebc',tracersidebc,'efzheight',efzheight,'icacczhat',icacczhat,'interpstyle_t',interpstyle_t,'storagebuffer',storagebuffer);

% Characteristic Scales (computed automatically):
ScaleParameters=struct('tempscale',[],'icethickscale',[],'viscositylogscale',[],'dragcoefflogscale',[],'sidedragcoefflogscale',[]);

% Other Parameters:
g=9.8;                            % m/s^2     (gravitational acceleration)
r=8.314;                          % J/(mol*K) (ideal gas constant)
secondsperyear=60*60*24*365.25;   % s/yr
verysmallnumber=1e-12;            % unitless
displayinterval=5000;             % integer 



overwriteoutput=1;                % logical




shutoffhour=NaN;                  % time of day (24 hr clock, decimal hours) (NaN deactivates)
restarthour=17;                   % time of day (24 hr clock, decimal hours) (NaN deactivates)
verbose=1;                        % 0 or 1
OtherParameters=struct('g',g,'r',r,'secondsperyear',secondsperyear,'verysmallnumber',verysmallnumber,'displayinterval',displayinterval,'overwriteoutput',overwriteoutput,'shutoffhour',shutoffhour,'restarthour',restarthour,'verbose',verbose);

%% Preparation:

% Check if resuming an old run:
if resumeoldrun==0
    
    % Load input file:
    load(inputfile{1},'*_input*')
    
    % Set initial computer time:
    computertime=0;
    
    % Check input parameters:
    CheckParameters_2D_v3;
    
    % Disable 1D interpolation warning:
    warning('off','MATLAB:interp1:NaNinY')
    
    % Convert years to seconds:
    constbasalmelt=constbasalmelt_yr/secondsperyear;
    dt=dt_yr*secondsperyear;
    dt_basic=dt;
    recordinterval=recordinterval_yr*secondsperyear;
    spinuptime=spinuptime_yr*secondsperyear;
    runtime=runtime_yr*secondsperyear;
    strainratestabilizer=strainratestabilizer_yr/secondsperyear;
    slidingstabilizer=slidingstabilizer_yr/secondsperyear;
    spawninterval_u=spawninterval_u_yr*secondsperyear;
    spawninterval_d=spawninterval_d_yr*secondsperyear;
    spawninterval_l=spawninterval_l_yr*secondsperyear;
    elevanomaccumlapserate=elevanomaccumlapserate_yr/secondsperyear;
    elevanomablatelapserate=elevanomablatelapserate_yr/secondsperyear;
    if dosill
        sillstarttime=sillstarttime_yr*secondsperyear;
        sillconstructiontime=sillconstructiontime_yr*secondsperyear;
    end
    
    % Define dimensionless vertical grid:
    if maxdensify<1
        error('Parameter "maxdensify" must be greater than or equal to 1.')
    end
    Zhat_ud=linspace(0,1,zsize+1)'/maxdensify+(1-1/maxdensify)*linspace(0,1,zsize+1)'.^densifypower;
    Zhat_c=.5*(Zhat_ud(1:end-1)+Zhat_ud(2:end));
    DZhat_c=Zhat_ud(2:end)-Zhat_ud(1:end-1);
    DZhat_ud=[Zhat_c(1);Zhat_c(2:end)-Zhat_c(1:end-1);1-Zhat_c(end)];  % first and last are half-cells
    
    % Ensure that input position vector starts at zero:
    X_input=X_input-X_input(1);
    
    % Define domain width:
    if sum(isnan(Icethick_input))~=0
        lasticeind=find(isnan(Icethick_input)==0,1,'last');
        domainwidth=X_input(lasticeind)*(1-1e-9); % arbitrary small number
    else
        lasticeind=length(X_input);
        domainwidth=X_input(end);
    end
    
    % Compute horizontal grid size for ice:
    inputxsize=length(X_input);
    xsize=round(domainwidth/targetdx);
    GridParameters.xsize=xsize;
    
    % Compute grid size for the plume model:
    if doplume
        % Locate grounding line in input profile:
        lastgroundedind=find(BedElev_input+(rho_i/rho_sw)*Icethick_input<0|isnan(Icethick_input)==0,1,'last');
        % Compute plume model horizontal grid size:
        xsize_plume=max(minxsize_plume,round((X_input(lasticeind)-X_input(lastgroundedind))/targetdx_plume));
        % Compute plume model vertical grid size:
        zsize_plume=max(minzsize_plume,round(min((rho_i/rho_sw)*Icethick_input(lasticeind),-BedElev_input(lasticeind))/targetdz_plume));
        % Record plume grid sizes:
        GridParameters.xsize_plume=xsize_plume;
        GridParameters.zsize_plume=zsize_plume;
    end
    
    % Define maximum domain width:
    maxdomainwidth=X_input(end);
    
    % Define sill profile:
    if dosill
        % Compute central location:
        if strcmp(sillcentx,'auto')
            sillcentx=maxdomainwidth-silldistbuffer*sillwidth;
            SillParameters.sillcentx=sillcentx;
        end
        % Compute sill amplitude:
        sillamp=silltopz-interp1(X_input,BedElev_input,sillcentx,interpstyle,BedElev_input(end));
        SillParameters.sillamp=sillamp;
        % Compute sill thickness profile:
        if sillstarttime==0 && sillconstructiontime==0
            SillThick_input=sillamp*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2));
        else
            SillThick_input=zeros(size(X_input));
        end
    else
        SillThick_input=zeros(size(X_input));
    end
    
    % Define initial horizontal grid:
    dx=domainwidth/xsize;
    X_lr=linspace(0,domainwidth,xsize+1);
    X_c=linspace(.5*dx,domainwidth-.5*dx,xsize);
    
    % Interpolate bed:
    BedElev_c=interp1(X_input,BedElev_input+SillThick_input,X_c,interpstyle,BedElev_input(end)+SillThick_input(end));
    bedelev_l=interp1(X_input,BedElev_input+SillThick_input,0,interpstyle,BedElev_input(1)+SillThick_input(1));
    bedelev_r=interp1(X_input,BedElev_input+SillThick_input,domainwidth,interpstyle,BedElev_input(end)+SillThick_input(end));
    BedElev_lr=[bedelev_l,.5*(BedElev_c(1:end-1)+BedElev_c(2:end)),bedelev_r];
    
    % Interpolate ice thickness:
    Icethick_c=interp1(X_input(1:lasticeind),Icethick_input(1:lasticeind),X_c,interpstyle,Icethick_input(end));
    
    % Check area-conserving vs flux-conserving width:
    if exist('Width_ac_input','var')
        if isequal(Width_ac_input,Width_input)==0
            % Inform command line:
            if verbose
                disp('...')
                disp('Swapping area-conserving flowband width for flux-conserving flowband width')
            end
            % Rename input width variables:
            Width_fc_input=Width_input;
            Width_input=Width_ac_input;
        end
    end
    
    % Interpolate width: (assumes width is constant vertically)
    Width_lr=interp1(X_input,Width_input,X_lr,interpstyle,Width_input(end));
    Width_c=.5*(Width_lr(1:end-1)+Width_lr(2:end));
    
    % Adjust ice thickness to match side BC, left side:
    if strcmp(leftbctype,'fixedicethick') && isequal(Icethick_c(1),leftbcparam)==0
        if strcmp(leftbcparam,'file')==1
            icethickadjustment_l=interp1(Time_input,LeftBCParam_input,0,interpstyle,LeftBCParam_input(1))-Icethick_c(1);
        else
            icethickadjustment_l=leftbcparam-Icethick_c(1);
        end
        Icethick_c=Icethick_c+icethickadjustment_l*(1-(X_c-.5*dx)/(domainwidth-dx));
        if verbose
            disp('Modified ice thickness to match left side BC')
        end
    else
        icethickadjustment_l=0;
    end
    SideParameters.icethickadjustment_l=icethickadjustment_l;
    % Adjust ice thickness to match side BC, right side:
    if strcmp(rightbctype,'fixedicethick') && isequal(Icethick_c(end),rightbcparam)==0
        if strcmp(rightbcparam,'file')==1
            icethickadjustment_r=interp1(Time_input,RightBCParam_input,0,interpstyle,RightBCParam_input(1))-Icethick_c(end);
        else
            icethickadjustment_r=rightbcparam-Icethick_c(end);
        end
        Icethick_c=Icethick_c+icethickadjustment_r*(X_c-.5*dx)/(domainwidth-dx);
        if verbose
            disp('Modified ice surface to match right side BC')
        end
    else
        icethickadjustment_r=0;
    end
    SideParameters.icethickadjustment_r=icethickadjustment_r;
    
    % Regularize time intervals:
    spawninterval_l_yr=spawninterval_u_yr/ceil(spawninterval_u_yr/spawninterval_l_yr);
    
    % Record regularized time intervals:
    TracerParameters.spawninterval_l_yr=spawninterval_l_yr;
    
    % Define sliding parameters:
    if strcmp(slidingvelocityscale_yr,'file')
        SlidingVelocityScale_lrd=interp1(X_input,SlidingVelocityScale_yr_input,X_lr,interpstyle)/secondsperyear;
    else
        SlidingVelocityScale_lrd=(slidingvelocityscale_yr/secondsperyear)*ones(1,xsize+1);
    end
    if strcmp(slidingstressscale,'file')
        SlidingStressScale_lrd=interp1(X_input,SlidingStressScale_input,X_lr,interpstyle);
    else
        SlidingStressScale_lrd=slidingstressscale*ones(1,xsize+1);
    end
    if strcmp(m,'file')
        M_lrd=interp1(X_input,M_input,X_lr,interpstyle);
    else
        M_lrd=m*ones(1,xsize+1);
    end
    if dosill && dosilldrag
        SillThick_lrd=interp1(X_input,SillThick_input,X_lr,interpstyle,0);
        SlidingStressScale_lrd=(1-exp(-SillThick_lrd/silldragthickscale))*sillstressscale+...
            exp(-SillThick_lrd/silldragthickscale).*SlidingStressScale_lrd;
    end
    
    % Define geothermal flux:
    if strcmp(gflux,'file')
        Gflux_c=interp1(X_input,Gflux_input,X_c,interpstyle,Gflux_input(end));
    else
        Gflux_c=gflux*ones(xsize,1);
    end
    
    % Define side influx:
    if usesidemassinput==1
        SideInflux_c=interp1(X_input,SideInflux_input,X_c,interpstyle,SideInflux_input(end));
    else
        SideInflux_c=zeros(1,xsize);
    end
    
    % Shorten initial timesteps:
    if spinuptime_yr>0
        dt=dt/rapidfactor;
        dt_yr=dt/secondsperyear;
    end
    
    % Compute number of tracers:
    if dotracers
        numlivetracers=tracerxsize*(tracerzsize+1);
        numliveconnectors=(tracerxsize-1)*(tracerzsize+1);
        numtracers=storagebuffer*numlivetracers;
        numconnectors=storagebuffer*numliveconnectors;
    end
    
    % Label the internal edges:
    IsInternal_lrd=[false(1),true(1,xsize-1),false(1)];
    
    % Pre-allocate model variables:
    GridElev_lr=zeros(zsize,xsize+1);
    Temp_c=zeros(zsize,xsize);
    Temp_d=zeros(1,xsize);
    IceBottom_c=zeros(1,xsize);
    MeltRate_d=zeros(1,xsize);
    MeltRate_c=zeros(zsize,xsize);
    A_c=zeros(zsize,xsize);
    if strcmp(velocitytype,'SSA')
        Viscosity_c=viscosityguess*ones(1,xsize);
    else
        Viscosity_c=viscosityguess*ones(zsize,xsize);
    end
    U_lr=zeros(zsize,xsize+1);
    W_ud=zeros(zsize+1,xsize);
    U_lrd=zeros(1,xsize+1);
    Drag_lrd=zeros(1,xsize+1);
    SlideHeat_d=zeros(1,xsize);
    if strcmp(velocitytype,'fullstokes')
        PressureDynamic_c=zeros(zsize,xsize);
    end
    HydroHeat_lrd=zeros(1,xsize+1);
    HydroHeat_d=zeros(1,xsize);
    if dotracers
        U_t=zeros(numtracers,1);
        What_t=zeros(numtracers,1);
    end
    if doplume
        MonotonicIceBottom_c=zeros(1,xsize);
        BottomGradient_lr=zeros(1,xsize+1);
        X_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        Z_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        if dotemp
            if strcmp(plumesolvertype,'custom')
                CondFlux_c_plume=zeros(1,xsize_plume+zsize_plume);
            else
                CondFlux_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
            end
        end
        SinTheta_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        Thick_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        U_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        Temp_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        Salinity_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        MeltPoint_c_plume=zeros(1,xsize_plume+zsize_plume);
        MeltRate_c_plume=zeros(1,xsize_plume+zsize_plume);
    end
    
    % Interpolate time-variable side parameters to start time:
    % Left side ice BC:
    if strcmp(leftbcparam,'file') && strcmp(leftbctype,'flux')
        influx_l=interp1(Time_input,LeftBCParam_input,0,interpstyle,LeftBCParam_input(1))/secondsperyear;
    elseif strcmp(leftbcparam,'file')
        icethickchangerate_l=(interp1(Time_input,LeftBCParam_input,dt,interpstyle,LeftBCParam_input(1))-interp1(Time_input,LeftBCParam_input,0,interpstyle,LeftBCParam_input(1)))/dt;
    end
    % Right side ice BC:
    if strcmp(rightbcparam,'file') && strcmp(rightbctype,'flux')
        influx_r=interp1(Time_input,RightBCParam_input,0,interpstyle,RightBCParam_input(1))/secondsperyear;
    elseif strcmp(rightbcparam,'file') && strcmp(rightbctype,'fixedicethick')
        icethickchangerate_r=(interp1(Time_input,RightBCParam_input,dt,interpstyle,RightBCParam_input(1))-interp1(Time_input,RightBCParam_input,0,interpstyle,RightBCParam_input(1)))/dt;
    elseif strcmp(rightbcparam,'file') && strcmp(rightbctype,'front')
        sealevel=interp1(Time_input,RightBCParam_input,0,interpstyle,RightBCParam_input(1));
    end
    % Left side water BC:
    if strcmp(waterinflux_l_yr,'file') 
        waterinflux_l=interp1(Time_input,WaterInflux_l_yr_input,0,interpstyle,WaterInflux_l_yr_input(1))/secondsperyear;
    end
    
    % Set side ice thickness change rates to zero if they are not changing:
    if strcmp(leftbctype,'fixedicethick') && strcmp(leftbcparam,'file')==0
        icethickchangerate_l=0;
    end
    if strcmp(rightbctype,'fixedicethick') && strcmp(rightbcparam,'file')==0
        icethickchangerate_r=0;
    end
    
    % Assign sea level:
    if strcmp(rightbctype,'front') && strcmp(rightbcparam,'file')==0
        sealevel=rightbcparam;
    elseif strcmp(rightbctype,'front')
        sealevel=interp1(Time_input,RightBCParam_input,0,interpstyle,RightBCParam_input(1));
    else
        sealevel=min(BedElev_input)-1;
    end
    
    % Define input ice surface:
    if exist('SurfElev_input','var')==0
        SurfElev_input=max(BedElev_input+Icethick_input,sealevel+(1-rho_i/rho_sw)*Icethick_input);
        SurfElev_input(isnan(SurfElev_input))=sealevel;
    end
    
    % Create first grid:
    firstgrid=1;
    GridSetUp_2D_v3;
    firstgrid=0;
    
    % Compute time-variable climate anomalies:
    if strcmp(climatetimedependence,'elev')
        % Interpolate elevation anomaly:
        elevanomaly=interp1(Time_input,ElevAnomaly_input,0,interpstyle,ElevAnomaly_input(1));
        % Check if spatial dependence also goes by elevation:
        if strcmp(climatespatialdependence,'z')
            % Direct anomalies are zero:
            accumanomaly=0;
            annualmeltanomaly=0;
            surftempanomaly=0;
        else
            % Special  behavior to be compatible with x-dependence:
            % Constant accum and surftemp anomalies:
            accumanomaly=elevanomaly*elevanomaccumlapserate;
            surftempanomaly=elevanomaly*elevanomtemplapserate;
            % Annual melt anomaly is spatially variable:
            AnnualMeltAnomaly_u=zeros(1,xsize);
            % Linear profile below elevation anomaly:
            AnnualMeltAnomaly_u(SurfElev_c<=elevanomaly)=elevanomablatelapserate*(elevanomaly-SurfElev_c(SurfElev_c<=elevanomaly));
        end
        % Accum multiplier is one:
        accummultiplier=1;
    elseif strcmp(climatetimedependence,'direct')
        % Interpolate direct anomalies:
        surftempanomaly=interp1(Time_input,SurfTempAnomaly_input,0,interpstyle,SurfTempAnomaly_input(1));
        accumanomaly=interp1(Time_input,AccumAnomaly_input,0,interpstyle,AccumAnomaly_input(1));
        annualmeltanomaly=interp1(Time_input,AnnualMeltAnomaly_input,0,interpstyle,AnnualMeltAnomaly_input(1));
        % Elevation anomaly is zero:
        elevanomaly=0;
        % Accum multiplier is one:
        accummultiplier=1;
    elseif strcmp(climatetimedependence,'mix')
        % Interpolate elevation anomaly:
        elevanomaly=interp1(Time_input,ElevAnomaly_input,0,interpstyle,ElevAnomaly_input(1));
        % Interpolate accumulation multiplier:
        accummultiplier=interp1(Time_input,AccumMultiplier_input,0,interpstyle,AccumMultiplier_input(1));
        % Check if spatial dependence also goes by elevation:
        if strcmp(climatespatialdependence,'z')
            % Direct anomalies are zero:
            accumanomaly=0;
            annualmeltanomaly=0;
            surftempanomaly=0;
        else
            % Special  behavior to be compatible with x-dependence:
            % Constant surftemp anomaly:
            surftempanomaly=elevanomaly*elevanomtemplapserate;
            % Accum anomaly is zero:
            accumanomaly=0;
            % Annual melt anomaly is spatially variable:
            AnnualMeltAnomaly_u=zeros(1,xsize);
            % Linear profile below elevation anomaly:
            AnnualMeltAnomaly_u(SurfElev_c<=elevanomaly)=elevanomablatelapserate*(elevanomaly-SurfElev_c(SurfElev_c<=elevanomaly));
        end
    else
        % All anomalies are zero:
        accumanomaly=0;
        annualmeltanomaly=0;
        surftempanomaly=0;
        elevanomaly=0;
        % Accum multiplier is one:
        accummultiplier=1;
    end
    
    % Adjust surface temperature anomaly for seasonal cycle:
    surftempanomaly=surftempanomaly-(seasonaltempcycle/(2*pi*dt_yr))*sin(2*pi*dt_yr); % averages a sinusoid across the timestep (starting at time=0)
    
    % Compute climate variables:
    if strcmp(climatespatialdependence,'x')
        % Interpolate w.r.t. x:
        SurfTemp_u=min(interp1(X_input,SurfTemp_input,X_c,interpstyle,SurfTemp_input(end))+surftempanomaly,tmelt);
        Accum_u=accummultiplier*(interp1(X_input,Accum_input,X_c,interpstyle,Accum_input(end))+accumanomaly);
        if strcmp(climatetimedependence,'elev') || strcmp(climatetimedependence,'mix')
            AnnualMelt_u=interp1(X_input,AnnualMelt_input,X_c,interpstyle,AnnualMelt_input(end))+AnnualMeltAnomaly_u;
            annualmelt_r=interp1(X_input,AnnualMelt_input,domainwidth,interpstyle,AnnualMelt_input(end))+elevanomablatelapserate*max(0,elevanomaly-.5*surfelev_r);
        else
            AnnualMelt_u=interp1(X_input,AnnualMelt_input,X_c,interpstyle,AnnualMelt_input(end))+annualmeltanomaly;
            annualmelt_r=interp1(X_input,AnnualMelt_input,domainwidth,interpstyle,AnnualMelt_input(end))+annualmeltanomaly;
        end
    else
        % Interpolate w.r.t. z:
        SurfTemp_u=interp1(Z_input,SurfTemp_input,SurfElev_c-elevanomaly,interpstyle)+surftempanomaly; % positive elevanomaly means warmer conditions.
        if strcmp(climatetimedependence,'mix')
            Accum_u=accummultiplier*(interp1(Z_input,Accum_input,SurfElev_c,interpstyle)+accumanomaly); % no vertical shift for accum in mixed version
        else
            Accum_u=accummultiplier*(interp1(Z_input,Accum_input,SurfElev_c-elevanomaly,interpstyle)+accumanomaly);
        end
        AnnualMelt_u=interp1(Z_input,AnnualMelt_input,SurfElev_c-elevanomaly,interpstyle)+annualmeltanomaly;
        annualmelt_r=interp1(Z_input,AnnualMelt_input,.5*(surfelev_r+sealevel),interpstyle)+annualmeltanomaly; % neglects curvature in melt/elevation profile
        % Use lapse rate for out-of-range low values:
        SurfTemp_u(SurfElev_c-elevanomaly<Z_input(1))=SurfTemp_input(1)+surftempanomaly+surftempextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1));   % assumes Z_input is monotonically increasing
        if strcmp(climatetimedependence,'mix')
            Accum_u(SurfElev_c-elevanomaly<Z_input(1))=accummultiplier*(Accum_input(1)+accumanomaly+accumextrapslope*(SurfElev_c(SurfElev_c<Z_input(1))-Z_input(1)));   % no vertical shift for accum in mixed version
        else
            Accum_u(SurfElev_c-elevanomaly<Z_input(1))=accummultiplier*(Accum_input(1)+accumanomaly+accumextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1)));   % assumes Z_input is monotonically increasing
        end
        AnnualMelt_u(SurfElev_c-elevanomaly<Z_input(1))=AnnualMelt_input(1)+annualmeltanomaly+annualmeltextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1));   % assumes Z_input is monotonically increasing
        if .5*(surfelev_r+sealevel)<Z_input(1)
            annualmelt_r=AnnualMelt_input(1)+annualmeltanomaly+annualmeltextrapslope*(.5*(surfelev_r+sealevel)-elevanomaly-Z_input(1));   % neglects curvature in melt/elevation profile, assumes Z_input is monotonically increasing
        end
        % Truncate temperature at the melting point:
        SurfTemp_u=min(SurfTemp_u,tmelt);
    end
    
    % Zero surface ablation when the model starts in midwinter:
    MeltRate_u=zeros(1,xsize);
    
    % Make IC:
    if dotemp
        % Make the steady-state thermal IC:
        SteadyStateIC_2D_v4;
    else
        % Set everything to a constant temperature:
        Temp_c(:)=consttemp;
        % Balance velocity:
        U_lr=repmat(cumsum([influx_l,dx*Accum_u])./Icethick_lr_upwind,[zsize,1]); 
        W_ud=repmat(Zhat_ud,[1,xsize]).*repmat(-Accum_u,[zsize+1,1]);
        % Zero water flux:
        WaterFlux_lrd=zeros(1,xsize+1);
        % Zero melting:
        MeltRate_d(:)=0;
        % Every grid cell at the melting point:
        IsTied_d=true(1,xsize);
        % Basal shear heating:
        SlideHeat_d=abs(.5*(DrivingStressGrad_lr(1,1:end-1).*Icethick_lr_center(1:end-1).*U_lr(1,1:end-1)+DrivingStressGrad_lr(1,2:end).*Icethick_lr_center(2:end).*U_lr(1,2:end)));
        % Basal temperature at the melting point:
        Temp_d=MeltPoint_d;
    end
    % Truncate temperature at the melting point:
    Temp_c=min(Temp_c,MeltPoint_c);
    
    % Set side temperature BC:
    temp_ld=Temp_d(1);
    Temp_l=Temp_c(:,1);
    if strcmp(rightbctype,'front')
        if bedelev_r>sealevel
            temp_rd=min(refsurftemp+templapserate*(bedelev_r-refelev_temp),MeltPoint_d(end));
            Temp_r=min(refsurftemp+templapserate*(bedelev_r+Zhat_c*icethick_r-refelev_temp),MeltPoint_r);
        else
            temp_rd=MeltPoint_d(end);
            Temp_r=zeros(zsize,1);
            if sum(GridElev_lr(:,end)>sealevel)>0
                Temp_r(GridElev_lr(:,end)>sealevel)=min(refsurftemp+templapserate*(GridElev_lr(GridElev_lr(:,end)>sealevel,end)-refelev_temp),MeltPoint_r(GridElev_lr(:,end)>sealevel));
            end
            Temp_r(GridElev_lr(:,end)<=sealevel)=MeltPoint_r(GridElev_lr(:,end)<=sealevel);
        end
    else
        temp_rd=Temp_d(end);
        Temp_r=Temp_c(:,end);
    end
    
    % Compute left side BC shape functions:
    % Interpolate temperature to grid edge corners:
    Temp_lud=zeros(zsize+1,1);
    Temp_lud(end)=SurfTemp_u(1);
    Temp_lud(2:end-1)=.5*(Temp_l(1:end-1)+Temp_l(2:end));
    Temp_lud(1)=temp_ld;
    % Evaulate rheological constant on grid edge corners:
    A_lud=zeros(zsize+1,1);
    A_lud(Temp_lud>=t0)=a0*exp(-(q_big/r)*((1./Temp_lud(Temp_lud>=t0))-(1/t0)));
    A_lud(Temp_lud<t0)=a0*exp(-(q_small/r)*((1./Temp_lud(Temp_lud<t0))-(1/t0)));
    % Get normalized strain rate on grid corners:
    StrainRate_lud=A_lud.*((1-Zhat_ud).^n);
    % Integrate for normalized horizontal velocity:
    U_l=cumsum(DZhat_ud(1:end-1).*StrainRate_lud(1:end-1));
    % Define deformation shape function for horizontal velocity:
    ShapeFunctionU_l=U_l/sum(U_l.*DZhat_c);
    % Compute horizontal shape function:
    ShapeFunctionU_l=leftbcslidefrac*ones(zsize,1)+(1-leftbcslidefrac)*ShapeFunctionU_l;
    % Compute vertical shape function:
    ShapeFunctionW_lud=[0;cumsum(ShapeFunctionU_l.*DZhat_c)];
    
    % Compute right side BC shape functions:
    % Interpolate temperature to grid edge corners:
    Temp_rud=zeros(zsize+1,1);
    Temp_rud(end)=SurfTemp_u(end);
    Temp_rud(2:end-1)=.5*(Temp_r(1:end-1)+Temp_r(2:end));
    Temp_rud(1)=temp_rd;
    % Evaulate rheological constant on grid edge corners:
    A_rud=zeros(zsize+1,1);
    A_rud(Temp_rud>=t0)=a0*exp(-(q_big/r)*((1./Temp_rud(Temp_rud>=t0))-(1/t0)));
    A_rud(Temp_rud<t0)=a0*exp(-(q_small/r)*((1./Temp_rud(Temp_rud<t0))-(1/t0)));
    % Get normalized strain rate on grid corners:
    StrainRate_rud=A_rud.*((1-Zhat_ud).^n);
    % Integrate for normalized horizontal velocity:
    U_r=cumsum(DZhat_ud(1:end-1).*StrainRate_rud(1:end-1));
    % Define deformation shape function for horizontal velocity:
    ShapeFunctionU_r=U_r/sum(U_r.*DZhat_c);
    % Compute horizontal shape function:
    ShapeFunctionU_r=rightbcslidefrac*ones(zsize,1)+(1-rightbcslidefrac)*ShapeFunctionU_r;
    % Compute vertical shape function:
    ShapeFunctionW_rud=[0;cumsum(ShapeFunctionU_r.*DZhat_c)];
    
    % Create second grid:
    GridSetUp_2D_v3;
    
    % Compute sliding constant:
    % Pre-allocate memory:
    SlidingConstant_lrd=zeros(1,xsize+1);
    % Compute temperature effect:
    Temp_lrd=[Temp_d(1),.5*(Temp_d(1:end-1)+Temp_d(2:end)),Temp_d(end)];
    MeltPoint_lrd=[MeltPoint_d(1),.5*(MeltPoint_d(1:end-1)+MeltPoint_d(2:end)),MeltPoint_d(end)];
    TempEffect_lrd=zeros(1,xsize+1);
    if slidingtempscale==0
        TempEffect_lrd(2:end-1)=Temp_lrd(2:end-1)==MeltPoint_lrd(2:end-1);
    else
        TempEffect_lrd(2:end-1)=exp((Temp_lrd(2:end-1)-MeltPoint_lrd(2:end-1))/slidingtempscale);
    end
    % Finish sliding constant:
    SlidingConstant_lrd(2:end-1)=max(slidingstabilizer,SlidingVelocityScale_lrd(2:end-1).*TempEffect_lrd(2:end-1))./(SlidingStressScale_lrd(2:end-1).^M_lrd(2:end-1)); 
    
    % Compute first guess of basal drag coefficient:
    DragCoefficient_lrd=zeros(1,xsize+1); 
    DragCoefficient_lrd(2:end-1)=GroundedFraction_lr(2:end-1).*abs((DrivingStressGrad_lr(1,2:end-1).*Icethick_lr_center(2:end-1)).^(1-M_lrd(2:end-1)))./SlidingConstant_lrd(2:end-1);
    DragCoefficient_lrd(GroundedFraction_lr==0)=0;
    
    % Compute characteristic scales (temp, icethick):
    tempscale=mean(abs(MeltPoint_d-SurfTemp_u));
    icethickscale=mean(Icethick_c);
    
    % Compute characteristic scales (viscosity):  (NOTE STABILIZER)
    amin=a0*exp(-(q_small/r)*((1./min(SurfTemp_u))-(1/t0)));
    amax=a0*exp(-(q_big/r)*((1./tmelt)-(1/t0)));
    viscositylogscale=sqrt(log(((amin^(-1/n))*((mean(abs(Accum_u))/icethickscale)^((1-n)/n)))/((amax^-1)*(mean(abs(DrivingStressGrad_lr(1,:).*Icethick_lr_center))^(1-n))))^2+1);
    
    % Compute characteristic scales (basal drag coefficient):  (NOTE STABILIZER)
    dragcoefflogscale=sqrt(log(max(DragCoefficient_lrd(IsGrounded_lr==1&IsInternal_lrd))/min(DragCoefficient_lrd(IsGrounded_lr==1&IsInternal_lrd)))^2+1);
    
    % Compute first guess of side drag coefficient:
    if dosidedrag
        if strcmp(velocitytype,'SSA')
            DragCoefficient_lr=((exp(mean(log([amin,amax])))*Width_lr).^(-1/n)).*(U_lr(1,:).^((1-n)/n));
        else
            DragCoefficient_lr=((exp(mean(log([amin,amax])))*repmat(Width_lr,[zsize,1])).^(-1/n)).*(U_lr.^((1-n)/n));
        end
    else
        if strcmp(velocitytype,'SSA')
            DragCoefficient_lr=zeros(1,xsize+1);
        else
            DragCoefficient_lr=zeros(zsize,xsize+1);
        end
    end
    
    % Compute characteristic scales (side drag coefficient):
    sidedragcoefflogscale=sqrt(log(max(max(DragCoefficient_lr(:,2:end-1)))/min(min(DragCoefficient_lr(:,2:end-1))))^2+1);
    
    % Check characteristic scales:
    if isnan(tempscale) || isinf(tempscale) || tempscale<=0
        error('Calculation of characteristic temperature scale failed.')
    end
    if isnan(icethickscale) || isinf(icethickscale) || icethickscale<=0
        error('Calculation of characteristic ice thickness scale failed.')
    end
    if isnan(viscositylogscale) || isinf(viscositylogscale) || viscositylogscale<=0
        error('Calculation of characteristic viscosity scale failed.')
    end
    if isnan(dragcoefflogscale) || isinf(dragcoefflogscale) || dragcoefflogscale<=0
        error('Calculation of characteristic basal drag coefficient scale failed.')
    end
    if (isnan(sidedragcoefflogscale) || isinf(sidedragcoefflogscale) || sidedragcoefflogscale<=0) && dosidedrag
        error('Calculation of characteristic side drag coefficient scale failed.')
    end
    
    % Record characteristic scales:
    ScaleParameters.tempscale=tempscale;
    ScaleParameters.icethickscale=icethickscale;
    ScaleParameters.viscositylogscale=viscositylogscale;
    ScaleParameters.dragcoefflogscale=dragcoefflogscale;
    ScaleParameters.sidedragcoefflogscale=sidedragcoefflogscale;
    
    % Initial time variables:
    time=0;
    time_yr=0;
    if spinuptime_yr>0
        spinuptimedone=0;
    else
        spinuptimedone=1;
    end
    nextrecordtime=recordinterval;
    recordcounter=1;
    
    % Spawn initial tracer distribution:
    if dotracers
        SpawnTracers_v3;
    end
    
    % Pre-allocate model time series:
    totalnumrecords=ceil(runtime_yr/recordinterval_yr);
    ModelTimeSeries=struct('filename',inputfile,'Time',zeros(totalnumrecords,1),...
        'X_gl',zeros(totalnumrecords,1),'Domainwidth',zeros(totalnumrecords,1),...
        'Volume',zeros(totalnumrecords,1),'VAF',zeros(totalnumrecords,1),...
        'Accumulation',zeros(totalnumrecords,1),'Ablation',zeros(totalnumrecords,1),'OceanMelt',zeros(totalnumrecords,1),'Calving',zeros(totalnumrecords,1),'SideInflux',zeros(totalnumrecords,1));
    
    % Create initial conditions structure:
    % Basic initial conditions:
    InitialConditions=struct('Icethick_c',Icethick_c,'IceBottom_c',IceBottom_c,...
            'Accum_u',Accum_u,'AnnualMelt_u',AnnualMelt_u,'domainwidth',domainwidth,'x_gl',x_gl,'dx',dx);
    % Add thermal conditions:
    if dotemp
        InitialConditions.Temp_c=Temp_c;
        InitialConditions.Temp_l=Temp_l;
        InitialConditions.Temp_r=Temp_r;
        InitialConditions.Temp_d=Temp_d;
        InitialConditions.MeltRate_d=MeltRate_d;
        InitialConditions.SurfTemp_u=SurfTemp_u;
    end
    % Add tracer conditions:
    if dotracers
        InitialConditions.ID_t=ID_t;
        InitialConditions.X_t=X_t;
        InitialConditions.Zhat_t=Zhat_t;
        InitialConditions.Age_t_yr=Age_t_yr;
        InitialConditions.InDomain_t=InDomain_t;
        InitialConditions.IsOriginal_t=IsOriginal_t;
        InitialConditions.IsAcc_t=IsAcc_t;
        InitialConditions.IsLooseEnd_t=IsLooseEnd_t;
        InitialConditions.Connectivity_t=Connectivity_t;
        InitialConditions.DeadConnector_t=DeadConnector_t;
    end
    % Add sill thickness:
    if dosill && dosillerosion
        InitialConditions.SillThick_input=SillThick_input;
    end
    
    % Initialize record file:
    numrecords=1;
    numdigits=floor(log10(floor(runtime/recordinterval)+1))+1;
    save(outputfile,'*Parameters','InitialConditions','*_input','Zhat_*','DZhat_*','numdigits','inputfile','resumeoldrun','-v7.3')
    
    
    
    
    
    % Resume an old run:
elseif resumeoldrun==1
    
    % Communicate with human:
    if verbose
        disp(['Resuming old model run:',inputfile{1}])
    end
    
    % Check if input file and output file are identical:
    if strcmp(outputfile,inputfile{1})
        if verbose
            disp('Appending model run to old output file')
        end
        branchrun=0;
        if isempty(resumetime_yr)==0
            if verbose
                disp('WARNING: overwriting old model output file past resume point!')
            end
        end
    else
        if verbose
            disp('Branching model run to new output file')
        end
        branchrun=1;
    end
    
    % Load basic parameters and time series:
    load(inputfile{1},'ModelTimeSeries','*Parameters','*_input*','Zhat_*','DZhat_*','numdigits','numrecords','computertime')
    
    % Determine size of inputs:
    inputxsize=length(X_input);
    
    % Determine last ice index in inputs:
    if sum(isnan(Icethick_input))~=0
        lasticeind=find(isnan(Icethick_input)==0,1,'last');
    else
        lasticeind=length(X_input);
    end
    
    % Kludge for the old runs that didn't record computertime:
    if exist('computertime','var')==0
        computertime=0;
    end
    
    % Load initial conditions:
    if branchrun==0
        load(inputfile{1},'InitialConditions')
    end
    
    % Load new time-variable inputs: 
    if length(inputfile)>1
        % Rename old time-variable inputs:
        oldinputs=who('*_input*');
        for ii=1:length(oldinputs)
            eval([oldinputs{ii},'_old=',oldinputs{ii},';'])
        end
        % Load new inputs:
        load(inputfile{2},'*_input*')
    end
    
    % Remember shutoffhour and restarthour:
    shutoffhour1=shutoffhour;
    restarthour1=restarthour;
    
    % Unpack all the parameter structures except timing:
    unpack(ModelParameters)
    unpack(SillParameters)
    unpack(TopParameters)
    unpack(BottomParameters)
    unpack(SideParameters)
    unpack(CalvingParameters)
    unpack(PlumeParameters)
    unpack(ThermalParameters)
    unpack(RheologyParameters)
    unpack(GridParameters)
    unpack(IterationParameters)
    unpack(VelocitySolverParameters)
    unpack(TracerParameters)
    unpack(ScaleParameters)
    unpack(OtherParameters)
    
    % Reset shutoffhour and restarthour:
    shutoffhour=shutoffhour1;
    restarthour=restarthour1;
    
    % Unpack timing parameters and change runtime:
    newruntime_yr=runtime_yr;
    unpack(TimingParameters);
    runtime_yr=newruntime_yr;
    TimingParameters.runtime_yr=runtime_yr;
    
    % Check input parameters:
    CheckParameters_2D_v3;
    
    % Load model record for the resume point:
    if isempty(resumetime_yr)
        load(inputfile{1},'FinalConditions')
        ModelRecord=FinalConditions;
        clear FinalConditions
        resumetime_yr=ModelRecord.time_yr;
    else
        prefix='0'*ones(1,numdigits(end)-floor(log10(numrecords(end)))-1);
        idnumber=num2str(max(1,floor(resumetime_yr/recordinterval_yr+.5)));  % assumes that record time is the middle of the record
        load(inputfile{1},['ModelRecord_',prefix,idnumber])
        eval(['ModelRecord=ModelRecord_',prefix,idnumber,';'])
        clear(['ModelRecord_',prefix,idnumber])
    end
    
    % Assign initial conditions for a branch run:
    if branchrun==1
        InitialConditions=ModelRecord;
    end
    
    % Disable 1D interpolation warning:
    warning('off','MATLAB:interp1:NaNinY')
    
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
    
    % Restart time:
    time_yr=ModelRecord.time_yr+.5*recordinterval_yr;
    
    % Convert years to seconds:
    constbasalmelt=constbasalmelt_yr/secondsperyear;
    dt=dt_yr*secondsperyear;
    dt_basic=dt;
    time=time_yr*secondsperyear;
    recordinterval=recordinterval_yr*secondsperyear;
    spinuptime=spinuptime_yr*secondsperyear;
    runtime=runtime_yr*secondsperyear;
    strainratestabilizer=strainratestabilizer_yr/secondsperyear;
    slidingstabilizer=slidingstabilizer_yr/secondsperyear;
    spawninterval_u=spawninterval_u_yr*secondsperyear;
    spawninterval_d=spawninterval_d_yr*secondsperyear;
    spawninterval_l=spawninterval_l_yr*secondsperyear;
    elevanomaccumlapserate=elevanomaccumlapserate_yr/secondsperyear;
    elevanomablatelapserate=elevanomablatelapserate_yr/secondsperyear;
    if dosill
        sillstarttime=sillstarttime_yr*secondsperyear;
        sillconstructiontime=sillconstructiontime_yr*secondsperyear;
    end
    
    % Compute sill thickness profile:
    if dosill && dosillerosion
        SillThick_input=ModelRecord.SillThick_input;
    elseif dosill
        SillThick_input=sillamp*max(0,min(1,(time-sillstarttime)/(sillconstructiontime+1)))*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2)); % note stabilizer
    else
        SillThick_input=zeros(size(X_input));
    end
    
    % Define domain width:
    domainwidth=ModelRecord.domainwidth;
    
    % Define maximum domain width:
    maxdomainwidth=X_input(end);
    
    % Define horizontal grid:
    dx=domainwidth/xsize;
    X_lr=linspace(0,domainwidth,xsize+1);
    X_c=linspace(.5*dx,domainwidth-.5*dx,xsize);
    
    % Interpolate bed:
    BedElev_c=interp1(X_input,BedElev_input+SillThick_input,X_c,interpstyle,BedElev_input(end)+SillThick_input(end));
    bedelev_l=interp1(X_input,BedElev_input+SillThick_input,0,interpstyle,BedElev_input(1)+SillThick_input(1));
    bedelev_r=interp1(X_input,BedElev_input+SillThick_input,domainwidth,interpstyle,BedElev_input(end)+SillThick_input(end));
    BedElev_lr=[bedelev_l,.5*(BedElev_c(1:end-1)+BedElev_c(2:end)),bedelev_r];
    
    % Define ice thickness:
    Icethick_c=ModelRecord.Icethick_c;
    
    % Interpolate width: (assumes width is constant vertically)
    Width_lr=interp1(X_input,Width_input,X_lr,interpstyle,Width_input(end));
    Width_c=.5*(Width_lr(1:end-1)+Width_lr(2:end));
    
    % Define Temperature:
    if dotemp
        % Assign from model record:
        Temp_c=ModelRecord.Temp_c;
        Temp_l=ModelRecord.Temp_l;
        Temp_r=ModelRecord.Temp_r;
        Temp_d=ModelRecord.Temp_d;
    else
        % Assign from constant temperature:
        Temp_c=consttemp*ones(zsize,xsize);
        Temp_l=consttemp*ones(zsize,1);
        Temp_r=consttemp*ones(zsize,1);
        MeltPoint_d=tmelt-meltingpointslope*rho_i*g*Icethick_c;
        Temp_d=MeltPoint_d;
    end
    
    % Set up tracers:
    if dotracers
        % Compute number of tracers:
        numtracers=length(ModelRecord.X_t);
        numconnectors=size(ModelRecord.Connectivity_t,1);
        numlivetracers=sum(ModelRecord.InDomain_t);
        numliveconnectors=sum(ModelRecord.DeadConnector_t==0);
        % Define tracer variables:
        Inds_t=linspace(1,numtracers,numtracers)';
        X_t=ModelRecord.X_t;
        Zhat_t=ModelRecord.Zhat_t;
        ID_t=ModelRecord.ID_t;
        Age_t_yr=ModelRecord.Age_t_yr;
        IsOriginal_t=ModelRecord.IsOriginal_t;
        IsAcc_t=ModelRecord.IsAcc_t;
        IsLooseEnd_t=ModelRecord.IsLooseEnd_t;
        InDomain_t=ModelRecord.InDomain_t;
        Connectivity_t=ModelRecord.Connectivity_t;
        DeadConnector_t=ModelRecord.DeadConnector_t;
        ConnectorAge_t_yr=[Age_t_yr(Connectivity_t(1:numliveconnectors,1));NaN*ones(numconnectors-numliveconnectors,1)];
        nextid=max(ID_t(InDomain_t))+1;
        % Define left edge tracer spawn points:
        Zhat_t_l=sort(Zhat_t(IsLooseEnd_t&IsAcc_t==0),'descend');
        Age_t_l_yr=sort(Age_t_yr(IsLooseEnd_t&IsAcc_t==0),'ascend');
        dzhat_u_l=Zhat_t_l(1)-Zhat_t_l(2);
        % Compute next spawn times:
        nextspawntime_u=min(Age_t_yr(InDomain_t&IsAcc_t==0))*secondsperyear+spawninterval_u;
        nextspawntime_d=min(Age_t_yr(InDomain_t&IsAcc_t))*secondsperyear+spawninterval_d;
        nextspawntime_l=time+spawninterval_l;
    end
    
    % Unpack velocity-related records:
    U_lr=ModelRecord.U_lr;
    U_lrd=ModelRecord.U_lrd;
    u_r=ModelRecord.u_r;
    calvingrate_r=ModelRecord.calvingrate_r;
    W_ud=ModelRecord.W_ud;
    Drag_lrd=ModelRecord.Drag_lrd;
    Viscosity_c=ModelRecord.Viscosity_c;
    if strcmp(velocitytype,'fullstokes')
        PressureDynamic_c=ModelRecord.PressureDynamic_c;
    end
    
    % Make a guess of the drag coefficient:
    DragCoefficient_lrd=Drag_lrd./U_lrd;
    
    % Label the internal edges:
    IsInternal_lrd=[false(1),true(1,xsize-1),false(1)];
    
    % Define ice bottom:
    if dodynamicbottom
        IceBottom_c=ModelRecord.IceBottom_c;
    else
        IceBottom_c=zeros(1,xsize);
    end
    
    % Pre-allocate model variables:
    MeltRate_d=zeros(1,xsize);
    MeltRate_c=zeros(zsize,xsize);
    A_c=zeros(zsize,xsize);
    SlideHeat_d=zeros(1,xsize);
    HydroHeat_lrd=zeros(1,xsize+1);
    HydroHeat_d=zeros(1,xsize);
    if dotracers
        U_t=zeros(numtracers,1);
        What_t=zeros(numtracers,1);
    end
    if doplume
        MonotonicIceBottom_c=zeros(1,xsize);
        BottomGradient_lr=zeros(1,xsize+1);
        X_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        Z_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        if dotemp
            if strcmp(plumesolvertype,'custom')
                CondFlux_c_plume=zeros(1,xsize_plume+zsize_plume);
            else
                CondFlux_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
            end
        end
        SinTheta_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        Thick_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        U_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        Temp_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        Salinity_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
        MeltPoint_c_plume=zeros(1,xsize_plume+zsize_plume);
        MeltRate_c_plume=zeros(1,xsize_plume+zsize_plume);
    end
    
    % Define sliding parameters:
    if strcmp(slidingvelocityscale_yr,'file')
        SlidingVelocityScale_lrd=interp1(X_input,SlidingVelocityScale_yr_input,X_lr,interpstyle)/secondsperyear;
    else
        SlidingVelocityScale_lrd=(slidingvelocityscale_yr/secondsperyear)*ones(1,xsize+1);
    end
    if strcmp(slidingstressscale,'file')
        SlidingStressScale_lrd=interp1(X_input,SlidingStressScale_input,X_lr,interpstyle);
    else
        SlidingStressScale_lrd=slidingstressscale*ones(1,xsize+1);
    end
    if strcmp(m,'file')
        M_lrd=interp1(X_input,M_input,X_lr,interpstyle);
    else
        M_lrd=m*ones(1,xsize+1);
    end
    if dosill && dosilldrag
        SillThick_lrd=interp1(X_input,SillThick_input,X_lr,interpstyle,0);
        SlidingStressScale_lrd=(1-exp(-SillThick_lrd/silldragthickscale))*sillstressscale+...
            exp(-SillThick_lrd/silldragthickscale).*SlidingStressScale_lrd;
    end
    
    % Define geothermal flux:
    if strcmp(gflux,'file')
        Gflux_c=interp1(X_input,Gflux_input,X_c,interpstyle,Gflux_input(end));
    else
        Gflux_c=gflux*ones(xsize,1);
    end
    
    % Define side input:
    if usesidemassinput==1
        SideInflux_c=interp1(X_input,SideInflux_input,X_c,interpstyle,SideInflux_input(end));
    else
        SideInflux_c=zeros(1,xsize);
    end
    
    % Compute sliding constant:
    % Pre-allocate memory:
    SlidingConstant_lrd=zeros(1,xsize+1);
    % Compute temperature effect:
    Temp_lrd=[Temp_d(1),.5*(Temp_d(1:end-1)+Temp_d(2:end)),Temp_d(end)];
    MeltPoint_lrd=[MeltPoint_d(1),.5*(MeltPoint_d(1:end-1)+MeltPoint_d(2:end)),MeltPoint_d(end)];
    TempEffect_lrd=zeros(1,xsize+1);
    if slidingtempscale==0
        TempEffect_lrd(2:end-1)=Temp_lrd(2:end-1)==MeltPoint_lrd(2:end-1);
    else
        TempEffect_lrd(2:end-1)=exp((Temp_lrd(2:end-1)-MeltPoint_lrd(2:end-1))/slidingtempscale);
    end
    % Finish sliding constant:
    SlidingConstant_lrd(2:end-1)=max(slidingstabilizer,SlidingVelocityScale_lrd(2:end-1).*TempEffect_lrd(2:end-1))./(SlidingStressScale_lrd(2:end-1).^M_lrd(2:end-1)); 
    
    % Sort out side BC:
    % Left side ice BC:
    if strcmp(leftbcparam,'file') && strcmp(leftbctype,'flux')
        influx_l=interp1(Time_input,LeftBCParam_input,time,interpstyle,LeftBCParam_input(1))/secondsperyear;
    elseif strcmp(leftbcparam,'file')
        icethickchangerate_l=(interp1(Time_input,LeftBCParam_input,time+dt,interpstyle,LeftBCParam_input(1))-interp1(Time_input,LeftBCParam_input,time,interpstyle,LeftBCParam_input(1)))/dt;
    elseif strcmp(leftbctype,'flux')
        influx_l=leftbcparam/secondsperyear;
    end
    % Right side ice BC:
    if strcmp(rightbcparam,'file') && strcmp(rightbctype,'flux')
        influx_r=interp1(Time_input,RightBCParam_input,time,interpstyle,RightBCParam_input(1))/secondsperyear;
    elseif strcmp(rightbcparam,'file') && strcmp(rightbctype,'fixedicethick')
        icethickchangerate_r=(interp1(Time_input,RightBCParam_input,time+dt,interpstyle,RightBCParam_input(1))-interp1(Time_input,RightBCParam_input,time,interpstyle,RightBCParam_input(1)))/dt;
    elseif strcmp(rightbcparam,'file') && strcmp(rightbctype,'front')
        sealevel=interp1(Time_input,RightBCParam_input,time,interpstyle,RightBCParam_input(1));
    elseif strcmp(rightbctype,'flux')
        influx_r=rightbcparam/secondsperyear;
    end
    % Left side water BC:
    if strcmp(waterinflux_l_yr,'file') 
        waterinflux_l=interp1(Time_input,WaterInflux_l_yr_input,time,interpstyle,WaterInflux_l_yr_input(1))/secondsperyear;
    else
        waterinflux_l=waterinflux_l_yr/secondsperyear;
    end
    
    % Set side ice thickness change rates to zero if they are not changing:
    if strcmp(leftbctype,'fixedicethick') && strcmp(leftbcparam,'file')==0
        icethickchangerate_l=0;
    end
    if strcmp(rightbctype,'fixedicethick') && strcmp(rightbcparam,'file')==0
        icethickchangerate_r=0;
    end
    
    % Assign sea level:
    if strcmp(rightbctype,'front') && strcmp(rightbcparam,'file')==0
        sealevel=rightbcparam;
    elseif strcmp(rightbctype,'front')
        sealevel=interp1(Time_input,RightBCParam_input,time,interpstyle,RightBCParam_input(end));
    else
        sealevel=min(BedElev_input)-1;
    end
    
    % Define input ice surface:
    SurfElev_input=max(BedElev_input+Icethick_input,sealevel+(1-rho_i/rho_sw)*Icethick_input);
    SurfElev_input(isnan(SurfElev_input))=0;
    
    % Create first grid:
    firstgrid=1;
    GridSetUp_2D_v3;
    firstgrid=0;
    
    % Compute time-variable climate anomalies:
    if strcmp(climatetimedependence,'elev')
        % Interpolate elevation anomaly:
        elevanomaly=interp1(Time_input,ElevAnomaly_input,time,interpstyle,ElevAnomaly_input(end));
        % Check if spatial dependence also goes by elevation:
        if strcmp(climatespatialdependence,'z')
            % Direct anomalies are zero:
            accumanomaly=0;
            annualmeltanomaly=0;
            surftempanomaly=0;
        else
            % Special  behavior to be compatible with x-dependence:
            % Constant accum and surftemp anomalies:
            accumanomaly=elevanomaly*elevanomaccumlapserate;
            surftempanomaly=elevanomaly*elevanomtemplapserate;
            % Annual melt anomaly is spatially variable:
            AnnualMeltAnomaly_u=zeros(1,xsize);
            % Linear profile below elevation anomaly:
            AnnualMeltAnomaly_u(SurfElev_c<=elevanomaly)=elevanomablatelapserate*(elevanomaly-SurfElev_c(SurfElev_c<=elevanomaly));
        end
        % Accum multiplier is one:
        accummultiplier=1;
    elseif strcmp(climatetimedependence,'direct')
        % Interpolate direct anomalies:
        surftempanomaly=interp1(Time_input,SurfTempAnomaly_input,time,interpstyle,SurfTempAnomaly_input(end));
        accumanomaly=interp1(Time_input,AccumAnomaly_input,time,interpstyle,AccumAnomaly_input(end));
        annualmeltanomaly=interp1(Time_input,AnnualMeltAnomaly_input,time,interpstyle,AnnualMeltAnomaly_input(end));
        % Elevation anomaly is zero:
        elevanomaly=0;
        % Accum multiplier is one:
        accummultiplier=1;
    elseif strcmp(climatetimedependence,'mix')
        % Interpolate elevation anomaly:
        elevanomaly=interp1(Time_input,ElevAnomaly_input,time,interpstyle,ElevAnomaly_input(end));
        % Interpolate accumulation multiplier:
        accummultiplier=interp1(Time_input,AccumMultiplier_input,time,interpstyle,AccumMultiplier_input(end));
        % Check if spatial dependence also goes by elevation:
        if strcmp(climatespatialdependence,'z')
            % Direct anomalies are zero:
            accumanomaly=0;
            annualmeltanomaly=0;
            surftempanomaly=0;
        else
            % Special  behavior to be compatible with x-dependence:
            % Constant surftemp anomaly:
            surftempanomaly=elevanomaly*elevanomtemplapserate;
            % Accum anomaly is zero:
            accumanomaly=0;
            % Annual melt anomaly is spatially variable:
            AnnualMeltAnomaly_u=zeros(1,xsize);
            % Linear profile below elevation anomaly:
            AnnualMeltAnomaly_u(SurfElev_c<=elevanomaly)=elevanomablatelapserate*(elevanomaly-SurfElev_c(SurfElev_c<=elevanomaly));
        end
    else
        % All anomalies are zero:
        accumanomaly=0;
        annualmeltanomaly=0;
        surftempanomaly=0;
        elevanomaly=0;
        % Accum multiplier is one:
        accummultiplier=1;
    end
    
    % Adjust surface temperature anomaly for seasonal cycle:
    surftempanomaly=surftempanomaly-(seasonaltempcycle/(2*pi*dt_yr))*(sin(2*pi*(time_yr+dt_yr))-sin(2*pi*time_yr)); % averages a sinusoid across the timestep
    
    % Compute climate variables:
    if strcmp(climatespatialdependence,'x')
        % Interpolate w.r.t. x:
        SurfTemp_u=min(interp1(X_input,SurfTemp_input,X_c,interpstyle,SurfTemp_input(end))+surftempanomaly,tmelt);
        Accum_u=accummultiplier*(interp1(X_input,Accum_input,X_c,interpstyle,Accum_input(end))+accumanomaly);
        if strcmp(climatetimedependence,'elev') || strcmp(climatetimedependence,'mix')
            AnnualMelt_u=interp1(X_input,AnnualMelt_input,X_c,interpstyle,AnnualMelt_input(end))+AnnualMeltAnomaly_u;
            annualmelt_r=interp1(X_input,AnnualMelt_input,domainwidth,interpstyle,AnnualMelt_input(end))+elevanomablatelapserate*max(0,elevanomaly-.5*surfelev_r);
        else
            AnnualMelt_u=interp1(X_input,AnnualMelt_input,X_c,interpstyle,AnnualMelt_input(end))+annualmeltanomaly;
            annualmelt_r=interp1(X_input,AnnualMelt_input,domainwidth,interpstyle,AnnualMelt_input(end))+annualmeltanomaly;
        end
    else
        % Interpolate w.r.t. z:
        SurfTemp_u=interp1(Z_input,SurfTemp_input,SurfElev_c-elevanomaly,interpstyle)+surftempanomaly; % positive elevanomaly means warmer conditions.
        if strcmp(climatetimedependence,'mix')
            Accum_u=accummultiplier*(interp1(Z_input,Accum_input,SurfElev_c,interpstyle)+accumanomaly); % no vertical shift for accum in mixed version
        else
            Accum_u=accummultiplier*(interp1(Z_input,Accum_input,SurfElev_c-elevanomaly,interpstyle)+accumanomaly);
        end
        AnnualMelt_u=interp1(Z_input,AnnualMelt_input,SurfElev_c-elevanomaly,interpstyle)+annualmeltanomaly;
        annualmelt_r=interp1(Z_input,AnnualMelt_input,.5*(surfelev_r+sealevel),interpstyle)+annualmeltanomaly; % neglects curvature in melt/elevation profile
        % Use lapse rate for out-of-range low values:
        SurfTemp_u(SurfElev_c-elevanomaly<Z_input(1))=SurfTemp_input(1)+surftempanomaly+surftempextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1));   % assumes Z_input is monotonically increasing
        if strcmp(climatetimedependence,'mix')
            Accum_u(SurfElev_c<Z_input(1))=accummultiplier*(Accum_input(1)+accumanomaly+accumextrapslope*(SurfElev_c(SurfElev_c<Z_input(1))-Z_input(1)));   % no vertical shift for accum in mixed version
        else
            Accum_u(SurfElev_c-elevanomaly<Z_input(1))=accummultiplier*(Accum_input(1)+accumanomaly+accumextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1)));   % assumes Z_input is monotonically increasing
        end
        AnnualMelt_u(SurfElev_c-elevanomaly<Z_input(1))=AnnualMelt_input(1)+annualmeltanomaly+annualmeltextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1));   % assumes Z_input is monotonically increasing
        if .5*(surfelev_r+sealevel)<Z_input(1)
            annualmelt_r=AnnualMelt_input(1)+annualmeltanomaly+annualmeltextrapslope*(.5*(surfelev_r+sealevel)-elevanomaly-Z_input(1));   % neglects curvature in melt/elevation profile, assumes Z_input is monotonically increasing
        end
        % Truncate temperature at the melting point:
        SurfTemp_u=min(SurfTemp_u,tmelt);
    end
    
    % Convert mean annual melt into peak summer melt:
    PeakMelt_u=AnnualMelt_u*pi/(2*meltseasonduration);
    
    % Compute fractional year:
    yearfracstart=time_yr-floor(time_yr);
    yearfracend=time_yr+dt_yr-floor(time_yr); % values greater than 1 indicate timestep goes into the next year
    
    % Compute timestep melt based on mean annual melt:
    % Check timestep length relative to annual cycle and melt season length:
    if dt_yr>=1 % timestep longer than a year
        % Instantaneous melt rate equals annual melt rate:
        MeltRate_u=AnnualMelt_u;
    elseif dt_yr<1 && dt_yr>=meltseasonduration % timestep longer than the melt season but shorter than a year
        if (yearfracend<=.5*(1-meltseasonduration) && yearfracstart<=.5*(1-meltseasonduration) && yearfracend-1<=.5*(1-meltseasonduration)) || (yearfracend>=.5*(1+meltseasonduration) && yearfracstart>=.5*(1+meltseasonduration))
            % Timestep is completely outside of the melt season
            % Instantaneous melt rate is zero:
            MeltRate_u=zeros(xsize,1);
        elseif yearfracstart<=.5*(1-meltseasonduration) && yearfracend>=.5*(1+meltseasonduration)
            % Timestep completely contains the melt season
            % Average the melt season out through the timestep:
            MeltRate_u=PeakMelt_u*(2/pi)*(meltseasonduration/dt_yr);
        elseif yearfracstart<.5*(1-meltseasonduration) && yearfracend>=.5*(1-meltseasonduration)
            % End of timestep straddles the beginning of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
        elseif yearfracstart<=.5*(1+meltseasonduration) && yearfracend>.5*(1+meltseasonduration)
            % Begining of timestep straddles the end of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)+1);
            % Check if end of timestep straddles the beginning of the next
            % melt season:
            if yearfracend-1>=.5*(1-meltseasonduration)
                % Add the average of an additional partial sinusoid:
                MeltRate_u=MeltRate_u+PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(yearfracend-1-.5*(1-meltseasonduration))/meltseasonduration));
            end
        end
    else % timestep shorter than the melt season
        if (yearfracend<=.5*(1-meltseasonduration) && yearfracstart<=.5*(1-meltseasonduration)) || (yearfracend>=.5*(1+meltseasonduration) && yearfracstart>=.5*(1+meltseasonduration))
            % Timestep is completely outside of the melt season
            % Instantaneous melt rate is zero:
            MeltRate_u=zeros(xsize,1);
        elseif yearfracstart>=.5*(1-meltseasonduration) && yearfracend<=.5*(1+meltseasonduration)
            % Timestep is completely inside of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
        elseif yearfracstart<.5*(1-meltseasonduration) && yearfracend>=.5*(1-meltseasonduration)
            % End of timestep straddles the beginning of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
        elseif yearfracstart<=.5*(1+meltseasonduration) && yearfracend>.5*(1+meltseasonduration)
            % Begining of timestep straddles the end of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)+1);
        end
    end
    
    % Set side temperature BC:
    temp_ld=Temp_d(1);
    temp_rd=Temp_d(end);
    
    % Compute left side BC shape functions:
    % Interpolate temperature to grid edge corners:
    Temp_lud=zeros(zsize+1,1);
    Temp_lud(end)=SurfTemp_u(1);
    Temp_lud(2:end-1)=.5*(Temp_l(1:end-1)+Temp_l(2:end));
    Temp_lud(1)=temp_ld;
    % Evaulate rheological constant on grid edge corners:
    A_lud=zeros(zsize+1,1);
    A_lud(Temp_lud>=t0)=a0*exp(-(q_big/r)*((1./Temp_lud(Temp_lud>=t0))-(1/t0)));
    A_lud(Temp_lud<t0)=a0*exp(-(q_small/r)*((1./Temp_lud(Temp_lud<t0))-(1/t0)));
    % Get normalized strain rate on grid corners:
    StrainRate_lud=A_lud.*((1-Zhat_ud).^n);
    % Integrate for normalized horizontal velocity:
    U_l=cumsum(DZhat_ud(1:end-1).*StrainRate_lud(1:end-1));
    % Define deformation shape function for horizontal velocity:
    ShapeFunctionU_l=U_l/sum(U_l.*DZhat_c);
    % Compute horizontal shape function:
    ShapeFunctionU_l=leftbcslidefrac*ones(zsize,1)+(1-leftbcslidefrac)*ShapeFunctionU_l;
    % Compute vertical shape function:
    ShapeFunctionW_lud=[0;cumsum(ShapeFunctionU_l.*DZhat_c)];
    
    % Compute right side BC shape functions:
    % Interpolate temperature to grid edge corners:
    Temp_rud=zeros(zsize+1,1);
    Temp_rud(end)=SurfTemp_u(end);
    Temp_rud(2:end-1)=.5*(Temp_r(1:end-1)+Temp_r(2:end));
    Temp_rud(1)=temp_rd;
    % Evaulate rheological constant on grid edge corners:
    A_rud=zeros(zsize+1,1);
    A_rud(Temp_rud>=t0)=a0*exp(-(q_big/r)*((1./Temp_rud(Temp_rud>=t0))-(1/t0)));
    A_rud(Temp_rud<t0)=a0*exp(-(q_small/r)*((1./Temp_rud(Temp_rud<t0))-(1/t0)));
    % Get normalized strain rate on grid corners:
    StrainRate_rud=A_rud.*((1-Zhat_ud).^n);
    % Integrate for normalized horizontal velocity:
    U_r=cumsum(DZhat_ud(1:end-1).*StrainRate_rud(1:end-1));
    % Define deformation shape function for horizontal velocity:
    ShapeFunctionU_r=U_r/sum(U_r.*DZhat_c);
    % Compute horizontal shape function:
    ShapeFunctionU_r=rightbcslidefrac*ones(zsize,1)+(1-rightbcslidefrac)*ShapeFunctionU_r;
    % Compute vertical shape function:
    ShapeFunctionW_rud=[0;cumsum(ShapeFunctionU_r.*DZhat_c)];
    
    % Compute vertically averaged flow direction:
    FlowDir_lr=sign(sum(ModelRecord.U_lr.*repmat(DZhat_c,[1,xsize+1]),1));
    
    % Compute first guess of side drag coefficient:
    if dosidedrag
        amin=a0*exp(-(q_small/r)*((1./min(SurfTemp_u))-(1/t0)));
        amax=a0*exp(-(q_big/r)*((1./tmelt)-(1/t0)));
        if strcmp(velocitytype,'SSA')
            DragCoefficient_lr=((exp(mean(log([amin,amax])))*Width_lr).^(-1/n)).*(U_lr(1,:).^((1-n)/n));
        else
            DragCoefficient_lr=((exp(mean(log([amin,amax])))*repmat(Width_lr,[zsize,1])).^(-1/n)).*(U_lr.^((1-n)/n));
        end
    else
        if strcmp(velocitytype,'SSA')
            DragCoefficient_lr=zeros(1,xsize+1);
        else
            DragCoefficient_lr=zeros(zsize,xsize+1);
        end
    end
    
    % Initial time variables:
    spinuptimedone=1;
    nextrecordtime=time+recordinterval;
    if branchrun==0
        recordcounter=numrecords(end)+1;
    else
        recordcounter=floor(resumetime_yr/recordinterval_yr+.5)+1;
    end
    
    % Compute new numrecords and numdigits:
    numrecords=[numrecords;numrecords(end)+1];
    numdigits=[numdigits;floor(log10(floor(runtime/recordinterval)+1))+1];
    
    % Pre-allocate extra space in the model time series:
    totalnumrecords=ceil(runtime_yr/recordinterval_yr);
    if numtotalrecords>length(ModelTimeSeries.Time)
        ModelTimeSeries.Time=[ModelTimeSeries.Time;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.X_gl=[ModelTimeSeries.X_gl;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.Domainwidth=[ModelTimeSeries.Domainwidth;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.Volume=[ModelTimeSeries.Volume;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.VAF=[ModelTimeSeries.VAF;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.Accumulation=[ModelTimeSeries.Accumulation;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.Ablation=[ModelTimeSeries.Ablation;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.OceanMelt=[ModelTimeSeries.OceanMelt;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.Calving=[ModelTimeSeries.Calving;zeros(totalnumrecords-recordcounter+1,1)];
        ModelTimeSeries.SideInflux=[ModelTimeSeries.SideInflux;zeros(totalnumrecords-recordcounter+1,1)];
    end
    
    % Initialize new record file for a branched run:
    if branchrun==1
        save(outputfile,'*Parameters','InitialConditions','*_input*','Zhat_*','DZhat_*','numdigits','inputfile','resumeoldrun','resumetime_yr')
    else
        save(outputfile,'resumeoldrun','resumetime_yr','-append')
    end
    
    
else
    error('Parameter "resumeoldrun" must be equal to 0 or 1')
end



% Create (or reset) model record structure:
% Basic model record:
ModelRecord=struct('time_yr',0 ... % elapsed time record
    ,'Icethick_c',zeros(1,xsize),'IceBottom_c',zeros(1,xsize),'domainwidth',0,'x_gl',0 ... % ice thickness, domain size records
    ,'U_lr',zeros(zsize,xsize+1),'U_lrd',zeros(1,xsize+1),'u_r',0,'calvingrate_r',0,'subaerialmeltrate_r',0,'oceanmeltrate_r',0,'W_ud',zeros(zsize+1,xsize),'Drag_lrd',zeros(1,xsize+1),'Viscosity_c',zeros(strcmp(velocitytype,'SSA')+zsize*(strcmp(velocitytype,'SSA')==0),xsize)... % velocity related records
    ,'Accum_u',zeros(1,xsize),'MeltRate_u',zeros(1,xsize)... % surface climate records
    ,'numstokesiterations',0,'numunstablecolumns',0,'normalizedvolumeerror',0);  % numerical records
% Add thermal and hydrology records:
if dotemp
    ModelRecord.Temp_c=zeros(zsize,xsize);
    ModelRecord.Temp_d=zeros(1,xsize);
    ModelRecord.MeltRate_c=zeros(zsize,xsize);
    ModelRecord.StrainHeat_c=zeros(zsize,xsize);
    ModelRecord.SlideHeat_d=zeros(1,xsize);
    ModelRecord.MeltRate_d=zeros(1,xsize);
    ModelRecord.HydroHeat_d=zeros(1,xsize);
    ModelRecord.WaterFlux_lrd=zeros(1,xsize+1);
    ModelRecord.hasbasins=0;
    ModelRecord.SurfTemp_u=zeros(1,xsize);
    ModelRecord.normalizedenergyerror=0;
end
% Add tracer records:
if dotracers
    ModelRecord.ID_t=zeros(numtracers,1);
    ModelRecord.X_t=zeros(numtracers,1);
    ModelRecord.Zhat_t=zeros(numtracers,1);
    ModelRecord.Age_t_yr=zeros(numtracers,1);
    ModelRecord.InDomain_t=false(numtracers,1);
    ModelRecord.IsOriginal_t=false(numtracers,1);
    ModelRecord.IsAcc_t=false(numtracers,1);
    ModelRecord.IsLooseEnd_t=false(numtracers,1);
    ModelRecord.Connectivity_t=zeros(numconnectors,2);
    ModelRecord.DeadConnector_t=false(numconnectors,1);
end
% Add plume records:
if doplume
    ModelRecord.Thick_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
    ModelRecord.U_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
    ModelRecord.Temp_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
    ModelRecord.Salinity_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
    ModelRecord.MeltRate_c_plume=zeros(1,xsize_plume+zsize_plume);
    if strcmp(plumesolvertype,'custom')
        ModelRecord.meanplumeiterations=0;
    end
    ModelRecord.plumeseparationx=0;
    ModelRecord.plumeseparationz=0;
end
% Add dynamic pressure:
if strcmp(velocitytype,'fullstokes')
    ModelRecord.PressureDynamic_c=zeros(zsize,xsize);
end
% Add sill thickness:
if dosill && dosillerosion
    ModelRecord.SillThick_input=zeros(1,inputxsize);
end

%% Run Model:

% Initial debugging display: (insert plotting commands here)
% figure(1)
% hold off
% plot(0,0,'.k')
% hold on
% xlabel('Time (yr)')
% ylabel('Error (unitless)')
% plot(X_c/1000,W_ud(end,:)*secondsperyear,'k')
% hold on
% plot(x_gl*[1,1]/1000,[-10,20],'--k')
% xlim([0,domainwidth/1000])
% ylim([-10,20])
% imagesc(X_c/1000,Zhat_c,Temp_c-MeltPoint_c)
% hold on
% plot([x_gl,x_gl]/1000,[0,1],'w')
% set(gca,'ydir','normal')
% caxis([min(SurfTemp_u-tmelt),0])
% xlim([0,500])
% ylim([0,1])
% xlabel('Distance (km)')
% ylabel('Elevation (unitless)')
% title(['Temperature, T=0 yr'])
% plot(0,x_gl/1000,'dk','MarkerFaceColor','k')
% xlim([0,1000])
% xlabel('Time (yr)')
% ylabel('Position (km)')
% title('Grounding Line Position, T=0 yr')
% hold on
% spot_disp=4.5e5;
% GLind_c_disp=NaN*zeros(1000,1);
% NumIterations_disp=NaN*zeros(1000,1);
% T_yr_disp=zeros(1000,1);
% Ubar_lr=sum(U_lr.*repmat(DZhat_c,[1,xsize+1]),1);
% GLind_c_disp(1)=find(IsGrounded_c,1,'last');
% plot(T_yr_disp,Icethick_spot_disp,'k')
% xlim([0,23])
% ylim([400,1000])
% plot([0,X_c,domainwidth]/1000,[bedelev_l,BedElev_c,bedelev_r],'k')
% hold on
% plot(X_input/1000,BedElev_input,'k')
% plot([0,X_c(X_c<=x_gl),x_gl,X_c(X_c>x_gl),domainwidth]/1000,[surfelev_l,SurfElev_c(X_c<=x_gl),surfelev_gl,SurfElev_c(X_c>x_gl),surfelev_r],'k')
% plot([0,X_c(X_c<=x_gl),x_gl,X_c(X_c>x_gl),domainwidth]/1000,[icebottom_l,IceBottom_c(X_c<=x_gl),bedelev_gl,IceBottom_c(X_c>x_gl),icebottom_r],'k')
% plot(domainwidth*[1,1]/1000,[surfelev_r,icebottom_r],'k')
% plot(x_gl*[1,1]/1000,[surfelev_gl,bedelev_gl],'k')
% plot([0,maxdomainwidth]/1000,[0,0],'--k')
% title('Model Geometry, T=0 yr')
% xlabel('Distance (km)')
% ylabel('Elevation (m)')
% imagesc(X_lr/domainwidth,Zhat_c,U_lr*secondsperyear)
% hold on
% plot([1,1]*x_gl/domainwidth,[0,1],'w')
% set(gca,'ydir','normal')
% caxis([0,4e3])
% subplot(2,1,1)
% plot(T_yr_disp,GLind_c_disp,'k')
% xlim([0,23])
% ylim([60,80])
% title('Grounding Line Index, T=0 yr')
% xlabel('Time (yr)')
% ylabel('Index')
% subplot(2,1,2)
% plot(T_yr_disp,NumIterations_disp,'k')
% xlim([0,23])
% ylim([0,30])
% xlabel('Time (yr)')
% ylabel('Number')
% title('Number of stokes iterations, T=0 yr')
% drawnow
% pause(.03)
% moviename='/home/mike/Documents/Research/FullTempModel/FTM_2D/FlatBedResTest_test1_v5_viscmisfit_v1.avi';
% Movie=VideoWriter(moviename);
% open(Movie);

% % Plume calibration code: (runs with FlowlinePlumeCalibrate_v1.m)
% NormalizedX=(X_c(lastgroundedind+1:end)-x_gl)/(domainwidth-x_gl);
% NormalizedY=(icethick_gl-Icethick_c(lastgroundedind+1:end))/(icethick_gl-icethick_r);
% A=[ones(sum(NormalizedX<=.2),1),NormalizedX(NormalizedX<=.2)'];
% equation=(A'*A)\A'*NormalizedY(NormalizedX<=.2)';
% initialslope=equation(2);

% Initial command-line display:
if verbose
    disp('...')
    if resumeoldrun==0
        disp(['record interval=',num2str(recordinterval_yr),' yr'])
        disp(['spinup period=',num2str(spinuptime_yr),' yr'])
        disp(['dt=',num2str(dt_yr),' yr'])
    else
        disp(['Resuming old model run at t=',num2str(time_yr),' yr'])
    end
    disp('...')
    disp('Running Model')
end

% Initial scalar variables:
done=0;
timestep=1;
displaycounter=1;
recordedtracers=0;
hasbasins=0;
blewup=0;

% % Pre-allocate variables for calving calibration:
% calibrationtime_yr=1;
% calibrationinterval_yr=1;
% Meltrate=zeros(ceil(calibrationtime_yr/dt_yr),1);
% %CalvingStress=zeros(ceil(calibrationinterval_yr/dt_yr),1);
% FrontThick=zeros(ceil(calibrationtime_yr/dt_yr),1);
% FrontU=zeros(ceil(calibrationtime_yr/dt_yr),1);

% Loop through timesteps:
while done==0
    
    % Set ice thickness on sides, if fixed from file:
    if strcmp(leftbcparam,'file') && strcmp(leftbctype,'fixedicethick')
        Icethick_c(1)=interp1(Time_input,LeftBCParam_input,time,interpstyle,LeftBCParam_input(end));
    end
    if strcmp(rightbcparam,'file') && strcmp(rightbctype,'fixedicethick')
        Icethick_c(end)=interp1(Time_input,RightBCParam_input,time,interpstyle,RightBCParam_input(end));
    end
    
    % Check for time-variable side BC, interpolate in time:
    % Left side ice BC:
    if strcmp(leftbcparam,'file') && strcmp(leftbctype,'flux')
        influx_l=interp1(Time_input,LeftBCParam_input,time,interpstyle,LeftBCParam_input(end))/secondsperyear;
    elseif strcmp(leftbcparam,'file')
        icethickchangerate_l=(interp1(Time_input,LeftBCParam_input,time+.5*dt,interpstyle,LeftBCParam_input(end))-interp1(Time_input,LeftBCParam_input,time-.5*dt,interpstyle,LeftBCParam_input(end)))/dt;
    end
    % Right side ice BC:
    if strcmp(rightbcparam,'file') && strcmp(rightbctype,'flux')
        influx_r=interp1(Time_input,RightBCParam_input,time,interpstyle,RightBCParam_input(end))/secondsperyear;
    elseif strcmp(rightbcparam,'file') && strcmp(rightbctype,'fixedicethick')
        icethickchangerate_r=(interp1(Time_input,RightBCParam_input,time+.5*dt,interpstyle,RightBCParam_input(end))-interp1(Time_input,RightBCParam_input,time-.5*dt,interpstyle,RightBCParam_input(end)))/dt;
    elseif strcmp(rightbcparam,'file') && strcmp(rightbctype,'front')
        sealevel=interp1(Time_input,RightBCParam_input,time,interpstyle,RightBCParam_input(end));
    end
    % Left side water BC:
    if strcmp(waterinflux_l_yr,'file')
        waterinflux_l=interp1(Time_input,WaterInflux_l_yr_input,time,interpstyle,WaterInflux_l_yr_input(end))/secondsperyear;
    end
    
    % Set up new grid:
    GridSetUp_2D_v3;
    
    % Compute time-variable climate anomalies:
    if strcmp(climatetimedependence,'elev')
        % Interpolate elevation anomaly:
        elevanomaly=interp1(Time_input,ElevAnomaly_input,time,interpstyle,ElevAnomaly_input(end));
        % Check if spatial dependence also goes by elevation:
        if strcmp(climatespatialdependence,'z')
            % Direct anomalies are zero:
            accumanomaly=0;
            annualmeltanomaly=0;
            surftempanomaly=0;
        else
            % Special  behavior to be compatible with x-dependence:
            % Constant accum and surftemp anomalies:
            accumanomaly=elevanomaly*elevanomaccumlapserate;
            surftempanomaly=elevanomaly*elevanomtemplapserate;
            % Annual melt anomaly is spatially variable:
            AnnualMeltAnomaly_u=zeros(1,xsize);
            % Linear profile below elevation anomaly:
            AnnualMeltAnomaly_u(SurfElev_c<=elevanomaly)=elevanomablatelapserate*(elevanomaly-SurfElev_c(SurfElev_c<=elevanomaly));
        end
        % Accum multiplier is one:
        accummultiplier=1;
    elseif strcmp(climatetimedependence,'direct')
        % Interpolate direct anomalies:
        surftempanomaly=interp1(Time_input,SurfTempAnomaly_input,time,interpstyle,SurfTempAnomaly_input(end));
        accumanomaly=interp1(Time_input,AccumAnomaly_input,time,interpstyle,AccumAnomaly_input(end));
        annualmeltanomaly=interp1(Time_input,AnnualMeltAnomaly_input,time,interpstyle,AnnualMeltAnomaly_input(end));
        % Elevation anomaly is zero:
        elevanomaly=0;
        % Accum multiplier is one:
        accummultiplier=1;
    elseif strcmp(climatetimedependence,'mix')
        % Interpolate elevation anomaly:
        elevanomaly=interp1(Time_input,ElevAnomaly_input,time,interpstyle,ElevAnomaly_input(end));
        % Interpolate accumulation multiplier:
        accummultiplier=interp1(Time_input,AccumMultiplier_input,time,interpstyle,AccumMultiplier_input(end));
        % Check if spatial dependence also goes by elevation:
        if strcmp(climatespatialdependence,'z')
            % Direct anomalies are zero:
            accumanomaly=0;
            annualmeltanomaly=0;
            surftempanomaly=0;
        else
            % Special  behavior to be compatible with x-dependence:
            % Constant surftemp anomaly:
            surftempanomaly=elevanomaly*elevanomtemplapserate;
            % Accum anomaly is zero:
            accumanomaly=0;
            % Annual melt anomaly is spatially variable:
            AnnualMeltAnomaly_u=zeros(1,xsize);
            % Linear profile below elevation anomaly:
            AnnualMeltAnomaly_u(SurfElev_c<=elevanomaly)=elevanomablatelapserate*(elevanomaly-SurfElev_c(SurfElev_c<=elevanomaly));
        end
    else
        % All anomalies are zero:
        accumanomaly=0;
        annualmeltanomaly=0;
        surftempanomaly=0;
        elevanomaly=0;
        % Accum multiplier is one:
        accummultiplier=1;
    end
    
    % Adjust surface temperature anomaly for seasonal cycle:
    surftempanomaly=surftempanomaly-(seasonaltempcycle/(2*pi*dt_yr))*(sin(2*pi*(time_yr+dt_yr))-sin(2*pi*time_yr)); % averages a sinusoid across the timestep
    
    % Compute climate variables:
    if strcmp(climatespatialdependence,'x')
        % Interpolate w.r.t. x:
        SurfTemp_u=min(interp1(X_input,SurfTemp_input,X_c,interpstyle,SurfTemp_input(end))+surftempanomaly,tmelt);
        Accum_u=accummultiplier*(interp1(X_input,Accum_input,X_c,interpstyle,Accum_input(end))+accumanomaly);
        if strcmp(climatetimedependence,'elev') || strcmp(climatetimedependence,'mix')
            AnnualMelt_u=interp1(X_input,AnnualMelt_input,X_c,interpstyle,AnnualMelt_input(end))+AnnualMeltAnomaly_u;
            annualmelt_r=interp1(X_input,AnnualMelt_input,domainwidth,interpstyle,AnnualMelt_input(end))+elevanomablatelapserate*max(0,elevanomaly-.5*surfelev_r);
        else
            AnnualMelt_u=interp1(X_input,AnnualMelt_input,X_c,interpstyle,AnnualMelt_input(end))+annualmeltanomaly;
            annualmelt_r=interp1(X_input,AnnualMelt_input,domainwidth,interpstyle,AnnualMelt_input(end))+annualmeltanomaly;
        end
    else
        % Interpolate w.r.t. z:
        SurfTemp_u=interp1(Z_input,SurfTemp_input,SurfElev_c-elevanomaly,interpstyle)+surftempanomaly; % positive elevanomaly means warmer conditions.
        if strcmp(climatetimedependence,'mix')
            Accum_u=accummultiplier*(interp1(Z_input,Accum_input,SurfElev_c,interpstyle)+accumanomaly); % no vertical shift for accum in mixed version
        else
            Accum_u=accummultiplier*(interp1(Z_input,Accum_input,SurfElev_c-elevanomaly,interpstyle)+accumanomaly);
        end
        AnnualMelt_u=interp1(Z_input,AnnualMelt_input,SurfElev_c-elevanomaly,interpstyle)+annualmeltanomaly;
        annualmelt_r=interp1(Z_input,AnnualMelt_input,.5*(surfelev_r+sealevel),interpstyle)+annualmeltanomaly; % neglects curvature in melt/elevation profile
        % Use lapse rate for out-of-range low values:
        SurfTemp_u(SurfElev_c-elevanomaly<Z_input(1))=SurfTemp_input(1)+surftempanomaly+surftempextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1));   % assumes Z_input is monotonically increasing
        if strcmp(climatetimedependence,'mix')
            Accum_u(SurfElev_c<Z_input(1))=accummultiplier*(Accum_input(1)+accumanomaly+accumextrapslope*(SurfElev_c(SurfElev_c<Z_input(1))-Z_input(1)));   % no vertical shift for accum in mixed version
        else
            Accum_u(SurfElev_c-elevanomaly<Z_input(1))=accummultiplier*(Accum_input(1)+accumanomaly+accumextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1)));   % assumes Z_input is monotonically increasing
        end
        AnnualMelt_u(SurfElev_c-elevanomaly<Z_input(1))=AnnualMelt_input(1)+annualmeltanomaly+annualmeltextrapslope*(SurfElev_c(SurfElev_c-elevanomaly<Z_input(1))-elevanomaly-Z_input(1));   % assumes Z_input is monotonically increasing
        if .5*(surfelev_r+sealevel)<Z_input(1)
            annualmelt_r=AnnualMelt_input(1)+annualmeltanomaly+annualmeltextrapslope*(.5*(surfelev_r+sealevel)-elevanomaly-Z_input(1));   % neglects curvature in melt/elevation profile, assumes Z_input is monotonically increasing
        end
        % Truncate temperature at the melting point:
        SurfTemp_u=min(SurfTemp_u,tmelt);
    end
    
    % Convert mean annual melt into peak summer melt:
    PeakMelt_u=AnnualMelt_u*pi/(2*meltseasonduration);
    peakmelt_r=annualmelt_r*pi/(2*meltseasonduration);
    
    % Compute fractional year:
    yearfracstart=time_yr-floor(time_yr);
    yearfracend=time_yr+dt_yr-floor(time_yr); % values greater than 1 indicate timestep goes into the next year
    
    % Compute timestep melt based on mean annual melt:
    % Check timestep length relative to annual cycle and melt season length:
    if dt_yr>=1 % timestep longer than a year
        % Instantaneous melt rate equals annual melt rate:
        MeltRate_u=AnnualMelt_u;
        subaerialmeltrate_r=annualmelt_r;
    elseif dt_yr<1 && dt_yr>=meltseasonduration % timestep longer than the melt season but shorter than a year
        if (yearfracend<=.5*(1-meltseasonduration) && yearfracstart<=.5*(1-meltseasonduration) && yearfracend-1<=.5*(1-meltseasonduration)) || (yearfracend>=.5*(1+meltseasonduration) && yearfracstart>=.5*(1+meltseasonduration))
            % Timestep is completely outside of the melt season
            % Instantaneous melt rate is zero:
            MeltRate_u=zeros(1,xsize);
            subaerialmeltrate_r=0;
        elseif yearfracstart<=.5*(1-meltseasonduration) && yearfracend>=.5*(1+meltseasonduration)
            % Timestep completely contains the melt season
            % Average the melt season out through the timestep:
            MeltRate_u=PeakMelt_u*(2/pi)*(meltseasonduration/dt_yr);
            subaerialmeltrate_r=peakmelt_r*(2/pi)*(meltseasonduration/dt_yr);
        elseif yearfracstart<.5*(1-meltseasonduration) && yearfracend>=.5*(1-meltseasonduration)
            % End of timestep straddles the beginning of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
            subaerialmeltrate_r=peakmelt_r*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
        elseif yearfracstart<=.5*(1+meltseasonduration) && yearfracend>.5*(1+meltseasonduration)
            % Begining of timestep straddles the end of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)+1);
            subaerialmeltrate_r=peakmelt_r*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)+1);
        end
        % Check if end of timestep straddles the beginning of the next
        % melt season:
        if yearfracend-1>=.5*(1-meltseasonduration)
            % Add the next year's melt season to the existing melt:
            MeltRate_u=MeltRate_u+PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(min(yearfracend-1,.5*(1+meltseasonduration))-.5*(1-meltseasonduration))/meltseasonduration));
            subaerialmeltrate_r=subaerialmeltrate_r+peakmelt_r*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(min(yearfracend-1,.5*(1+meltseasonduration))-.5*(1-meltseasonduration))/meltseasonduration));
        end
    else % timestep shorter than the melt season
        if (yearfracend<=.5*(1-meltseasonduration) && yearfracstart<=.5*(1-meltseasonduration)) || (yearfracend>=.5*(1+meltseasonduration) && yearfracstart>=.5*(1+meltseasonduration))
            % Timestep is completely outside of the melt season
            % Instantaneous melt rate is zero:
            MeltRate_u=zeros(1,xsize);
            subaerialmeltrate_r=0;
        elseif yearfracstart>=.5*(1-meltseasonduration) && yearfracend<=.5*(1+meltseasonduration)
            % Timestep is completely inside of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
            subaerialmeltrate_r=peakmelt_r*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
        elseif yearfracstart<.5*(1-meltseasonduration) && yearfracend>=.5*(1-meltseasonduration)
            % End of timestep straddles the beginning of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
            subaerialmeltrate_r=peakmelt_r*(meltseasonduration/(pi*dt_yr))*(1-cos(pi*(yearfracend-.5*(1-meltseasonduration))/meltseasonduration));
        elseif yearfracstart<=.5*(1+meltseasonduration) && yearfracend>.5*(1+meltseasonduration)
            % Begining of timestep straddles the end of the melt season
            % Average a partial sinusoid within this timestep:
            MeltRate_u=PeakMelt_u*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)+1);
            subaerialmeltrate_r=peakmelt_r*(meltseasonduration/(pi*dt_yr))*(cos(pi*(yearfracstart-.5*(1-meltseasonduration))/meltseasonduration)+1);
        end
    end
    
    % Set moving front temperature BC:
    if strcmp(rightbctype,'front') && dotemp
        if bedelev_r>sealevel
            temp_rd=min(refsurftemp+templapserate*(bedelev_r-refelev_temp),MeltPoint_d(end));
            Temp_r=min(refsurftemp+templapserate*(bedelev_r+Zhat_c*icethick_r-refelev_temp),MeltPoint_r);
        else
            temp_rd=MeltPoint_d(end);
            if sum(GridElev_lr(:,end)>sealevel)>0
                Temp_r(GridElev_lr(:,end)>sealevel)=min(refsurftemp+templapserate*(GridElev_lr(GridElev_lr(:,end)>sealevel,end)-refelev_temp),MeltPoint_r(GridElev_lr(:,end)>sealevel));
            end
            Temp_r(GridElev_lr(:,end)<=sealevel)=MeltPoint_r(GridElev_lr(:,end)<=sealevel);
        end
    end
    
    % Run grounded hydrology model:
    if dotemp
        BasalHydrology_2D_quasistatic_v2;
    else
        % Constant basal melting:
        MeltRate_d(:)=constbasalmelt;
        % Every grid cell at the melting point:
        IsTied_d=true(1,xsize);
        % Basal temperature at the melting point:
        Temp_d=MeltPoint_d;
        % Integrate basal and surface melt:
        WaterFlux_lrd(1)=waterinflux_l;
        for d2=1:xsize
            if IsGrounded_c(d2)
                WaterFlux_lrd(d2+1)=(Width_lr(d2)*WaterFlux_lrd(d2)+dx*Width_c(d2)*(MeltRate_d(d2)+MeltRate_u(d2)))/Width_lr(d2+1);
            else
                WaterFlux_lrd(d2+1)=0;
            end
        end
    end
    
    % Run plume model:
    if doplume
        if strcmp(plumesolvertype,'custom')
             PlumeModel_v2;
        else
             PlumeModel_v3;
        end
    else
        oceanmeltrate_r=0;
    end
    
    % Combine ocean melt and subaerial melt at the front:
    meltrate_r=subaerialmeltrate_r*min(1,(surfelev_r-sealevel)/icethick_r)+oceanmeltrate_r*max(0,(sealevel-icebottom_r)/icethick_r);
    
    % Interpolate basal temperature and melt point to grid edges:
    Temp_lrd=[Temp_d(1),.5*(Temp_d(1:end-1)+Temp_d(2:end)),Temp_d(end)];
    MeltPoint_lrd=[MeltPoint_d(1),.5*(MeltPoint_d(1:end-1)+MeltPoint_d(2:end)),MeltPoint_d(end)];
    
    % Compute sliding constant:
    if slidingtempscale==0
        TempEffect_lrd(2:end-1)=Temp_lrd(2:end-1)==MeltPoint_lrd(2:end-1);
    else
        TempEffect_lrd(2:end-1)=exp((Temp_lrd(2:end-1)-MeltPoint_lrd(2:end-1))/slidingtempscale);
    end
    SlidingConstant_lrd(2:end-1)=max(slidingstabilizer,SlidingVelocityScale_lrd(2:end-1).*TempEffect_lrd(2:end-1))./(SlidingStressScale_lrd(2:end-1).^M_lrd(2:end-1));
    
    % Calculate velocity:
    if strcmp(velocitytype,'SIA')
        ShallowIceApproximation_2D_v1
    elseif strcmp(velocitytype,'SSA')
        ShallowShelfApproximation_2D_v1
    elseif strcmp(velocitytype,'horzforce')
        HorizontalForceBalance_2D_v5;
    else % full stokes
        FullStokes_2D_v5a;
    end
    error('end here')
    
    %     % Calibrate plume model:
    %     % Target thinning rate: (taken from Pritchard et al., [2012],
    %     % supplementary table 1, PIG=-5.59,Thwaites=-6.14):
    %     targetthickeningrate=-6.14/secondsperyear;
    %     % Create vector of parameters:
    %     StantonNumbers=exp(linspace(log(1e-5),log(3e-3),100)');
    %     MeanW=zeros(100,1);
    %     % Loop through:
    %     for ii=1:100
    %         % Run plume model:
    %         stantonnumber=StantonNumbers(ii);
    %         PlumeModel_v2;
    %         % Compute mean thickening rate:
    %         ThickeningRate_c=W_ud(end,:)-W_ud(1,:)-MeltRate_d+Accum_u-MeltRate_u+SideInflux_c./Width_c;
    %         MeanW(ii)=mean(ThickeningRate_c(IsGrounded_c==0));
    %     end
    %     % Solve for best parameter:
    %     stopind=find(diff(MeanW)<=0,1,'last');
    %     bests=interp1(MeanW(1:stopind),StantonNumbers(1:stopind),targetthickeningrate,interpstyle)
    %     % Make a figure showing the calibration:
    %     figure(1)
    %     hold off
    %     semilogx(StantonNumbers,MeanW*secondsperyear,'k')
    %     hold on
    %     plot([1e-5,1e-2],[0,0],'--k')
    %     plot([1e-5,1e-2],targetthickeningrate*[1,1]*secondsperyear,'k')
    %     plot(bests,targetthickeningrate*secondsperyear,'Marker','o','MarkerFaceColor','k','Color','k')
    %     text(bests,targetthickeningrate*secondsperyear-3,['\Gamma_{TS}=',num2str(1e-5*round(1e5*bests))],'HorizontalAlignment','right','VerticalAlignment','top','Color','k')
    %     text(.1,targetthickeningrate*secondsperyear-3,'Pritchard et al., [2012] thinning rate','HorizontalAlignment','left','VerticalAlignment','top','Color','k')
    %     xlabel('Stanton Number (unitless)')
    %     ylabel('Mean Shelf Thickening Rate (m/yr)')
    %     title(['Plume Model Calibration in ',inputfile{1}(47:end-4)],'interpreter','none')
    %     figname=['/home/mjw/Documents/FjordSillPaper/Figures/',inputfile{1}(47:end-4),'_plumecalibration_v4.png'];
    %     set(gcf,'PaperSize',[8,6])
    %     set(gcf,'PaperPosition',[0,0,8,6])
    %     print('-dpng',figname,'-r300')
    %     % Quit here:
    %     error('end here')

    % Compute calving rate:
    if strcmp(rightbctype,'front')
        if strcmp(calvingtype,'fixed')
            calvingrate_r=sum(U_lr(:,end).*DZhat_c)-meltrate_r;
        elseif strcmp(calvingtype,'thick')
            calvingrate_r=refcalvingrate*refcalvingicethick/icethick_r;
        elseif strcmp(calvingtype,'uthick')
            calvingrate_r=sum(U_lr(:,end).*DZhat_c)*(refcalvingicethick/icethick_r);
        elseif strcmp(calvingtype,'meltthick')
            calvingrate_r=(refcalvingrate*refcalvingicethick/refcalvingmeltrate)*(max(MeltRate_c_plume(xsize_plume+1:end))/icethick_r);
        elseif strcmp(calvingtype,'vonmises')
            % Compute principle strain rates:
            strainratetensor_r=[sum(DZhat_c.*StrainRate_xx_c(:,end)),sum(DZhat_c.*StrainRate_xy_c(:,end));sum(DZhat_c.*StrainRate_xy_c(:,end)),sum(DZhat_c.*StrainRate_yy_c(:,end))];
            principlestrainrates_r=eig(strainratetensor_r);
            % Compute von mises stress:
            vonmisesstress_r=sqrt(3)*sum((A_c(:,end).^(-1/n)).*DZhat_c)*(sqrt(mean(max(0,principlestrainrates_r).^2))^(1/n));
            % Modify von mises stress for subaerial cliff face:
            %calvingstress_r=sqrt(vonmisesstress_r^2+(rho_i*g*min(surfelev_r-sealevel,icethick_r))^2);
            calvingstress_r=vonmisesstress_r;
            % Compute calving rate:
            calvingrate_r=sum(U_lr(:,end).*DZhat_c)*(calvingstress_r/icetensilestrength);
        elseif strcmp(calvingtype,'meltmultiplier')
            calvingrate_r=meltmultiplier*max(MeltRate_c_plume(xsize_plume+1:end));
        else % position-dependent from file
            calvingrate_r=interp1(X_input,CalvingRate_input,domainwidth,interpstyle,CalvingRate_input(end));
        end
        % Add calving for overhanging ice cliff:
        if strcmp(calvingtype,'fixed')==0
            calvingrate_r=calvingrate_r+(max(MeltRate_c_plume(xsize_plume+1:end))-meltrate_r);
        end
        % Reduce calving rate linearly for height above flotation:
        % (this step is in here to produce a thin land-based terminus)
        if strcmp(calvingtype,'fixed')==0
            calvingrate_r=calvingrate_r*max(0,1-(surfelev_r-(sealevel+(1-rho_i/rho_sw)*icethick_r))/((rho_i/rho_sw)*icethick_r));
        end
        % Ensure calving removes more mass than ablation and basal melt in final cell:
        % (this behavior deactivated now that the model can jump back when
        % thickness is too small)
        % calvingrate_r=sqrt(calvingrate_r^2+max(0,(MeltRate_d(end)+MeltRate_u(end))*dx/Icethick_c(end))^2);
        % Ensure that the ice front does not advance beyond the input
        % region:
        if domainwidth+dt*(sum(U_lr(:,end).*DZhat_c)-calvingrate_r-meltrate_r)>maxdomainwidth
            % Flag that the calving front is fixed:
            forcefixed=true(1);
            % Compute calving rate required to keep the front fixed:
            calvingrate_r=(maxdomainwidth-domainwidth)/dt+sum(U_lr(:,end).*DZhat_c)-meltrate_r;
        else
            % Flag that the calving front is not fixed:
            forcefixed=false(1);
        end
    else
        calvingrate_r=NaN;
    end
    
    %     % Calving model calibration:
    %     % Compute von mises stress:
    %     %strainratetensor_r=[StrainRate_xx_c(end),StrainRate_xy_c(end);StrainRate_xy_c(end),StrainRate_yy_c(end)];
    %     %principalstrainrates_r=eig(strainratetensor_r);
    %     %vonmisesstress_r=sqrt(3)*sum((A_c(:,end).^(-1/n)).*DZhat_c)*(sqrt(mean(max(0,principalstrainrates_r).^2))^(1/n));
    %     %calvingstress_r=sqrt(vonmisesstress_r^2+(rho_i*g*min(surfelev_r-sealevel,icethick_r))^2);
    %     %calvingstress_r=vonmisesstress_r;
    %     % Record relevant front variables:
    %     %Meltrate(timestep)=meltrate_r;
    %     Meltrate(timestep)=max(MeltRate_c_plume(xsize_plume+1:end));
    %     FrontThick(timestep)=icethick_r;
    %     %CalvingStress(timestep)=calvingstress_r;
    %     FrontU(timestep)=U_lr(end);
    %     % Check if the specified time has elapsed:
    %     if timestep==ceil(calibrationtime_yr/dt_yr)
    %         m0=exp(mean(log(Meltrate)))*secondsperyear
    %         error('end here')
    %         % Pre-allocate:
    %         %Sigma0=exp(linspace(log(1e3),log(1e7),100));
    %         %H0=exp(linspace(log(10),log(1e3),100));
    %         %MeltMultiplier=linspace(1,100,100);
    %         M0=exp(linspace(log(1),log(1000),100))/secondsperyear;
    %         MeanU=zeros(1,100);
    %         % Loop through potential parameter values:
    %         for ii=1:100
    %             % Compute calving rate and mean over last 5 years:
    %             %ThisCalvingRate=FrontU.*(CalvingStress/Sigma0(ii));
    %             %ThisCalvingRate=FrontU.*((FrontThick/H0(ii)).^-1);
    %             %ThisCalvingRate=Meltrate*MeltMultiplier(ii);
    %             ThisCalvingRate=FrontU.*(Meltrate/M0(ii)).*(calvingparam1./FrontThick);
    %             MeanU(ii)=mean(FrontU(end-ceil(calibrationinterval_yr/dt_yr)+1:end)-ThisCalvingRate(end-ceil(calibrationinterval_yr/dt_yr)+1:end)-Meltrate(end-ceil(calibrationinterval_yr/dt_yr)+1:end));
    %         end
    %         % Compute best-fit yield stress:
    %         %bestsigma0=interp1(MeanU,Sigma0,0,interpstyle);
    %         %besth0=interp1(MeanU,H0,0,interpstyle);
    %         %bestmultiplier=interp1(MeanU,MeltMultiplier,0,interpstyle)
    %         bestm0=interp1(MeanU,M0,0,interpstyle);
    %         %BestCalvingRate=FrontU.*(CalvingStress/bestsigma0);
    %         %BestCalvingRate=FrontU.*((FrontThick/besth0).^-1);
    %         %BestCalvingRate=Meltrate*bestmultiplier;
    %         BestCalvingRate=FrontU.*(Meltrate/bestm0).*(calvingparam1./FrontThick);
    %         NetU=FrontU-BestCalvingRate-Meltrate;
    %         % Make a figure showing the calibration:
    %         figure(1)
    %         hold off
    %         plot(linspace(0,calibrationtime_yr,length(FrontU)),FrontU*secondsperyear,'b')
    %         hold on
    %         plot(linspace(0,calibrationtime_yr,length(FrontU)),Meltrate*secondsperyear,'r')
    %         plot(linspace(0,calibrationtime_yr,length(FrontU)),BestCalvingRate*secondsperyear,'g')
    %         plot(linspace(0,calibrationtime_yr,length(FrontU)),NetU*secondsperyear,'k')
    %         plot([0,calibrationtime_yr],[0,0],'--k')
    %         xlabel('Time (yr)')
    %         ylabel('Velocity (m/yr)')
    %         %legend('Ice Velocity','Frontal Melt Rate',['Best-Fit Calving Rate, \sigma_0=',num2str(round(bestsigma0/1000)),' kPa'],'Net Frontal Velocity','location','Best')
    %         %legend('Ice Velocity','Frontal Melt Rate',['Best-Fit Calving Rate, H_0=',num2str(round(besth0)),' m'],'Net Frontal Velocity','location','Best')
    %         %legend('Ice Velocity','Frontal Melt Rate',['Best-Fit Calving Rate, K=',num2str(round(bestmultiplier))],'Net Frontal Velocity','location','Best')
    %         legend('Ice Velocity','Frontal Melt Rate',['Best-Fit Calving Rate, m_{0}=',num2str(round(bestm0*secondsperyear)),' m/yr, H_{0}=',num2str(round(calvingparam1)),' m'],'Net Frontal Velocity','location','Best')
    %         title(['Calibration of Melt and Thickness Dependent Calving Law in ',inputfile{1}(47:end-4)],'interpreter','none')
    %         figname=['/home/mjw/Documents/FjordSillPaper/Figures/',inputfile{1}(47:end-4),'_calvingcalibration_v8.png'];
    %         set(gcf,'PaperSize',[8,6])
    %         set(gcf,'PaperPosition',[0,0,8,6])
    %         print('-dpng',figname,'-r300')
    %         % Quit here:
    %         error('end here')
    %     end
    
    % Compute horizontal grid velocity in moving front case:
    if strcmp(rightbctype,'front') && strcmp(calvingtype,'fixed')==0 && forcefixed==0
        % Compute front velocity:
        u_r=sum(U_lr(:,end).*DZhat_c)-calvingrate_r-meltrate_r;
        % Calculate grid velocity:
        Ugrid_lr=repmat(u_r*X_lr/domainwidth,[zsize,1]);
    else
        % Grid is not moving horizontally:
        u_r=0;
        Ugrid_lr=zeros(zsize,xsize+1);
    end
    
    % Compute an ice thickness correction for horizontal grid motion:
    % This correction is upwind differencing for the "advection" of the
    % moving grid.  Without this correction, the model experiences
    % erroneous volume gain during front advance and erroneous volume loss
    % during front retreat.  Thus, the uncorrected model is prone to a
    % positive feedback of runaway loss or runaway growth.  With the
    % correction, residual mass conservation errors are reduced by several
    % orders of magnitude and reversed in sign: the residual errors are
    % positive (erroneous growth) during front retreat and negative
    % (erroneous loss) during front advance.  Thus, the residual errors are
    % negative feedbacks and do not tend to destabilize the model.  I
    % performed a rudimentary resolution test and found that the residual
    % errors scale roughly linearly with horizontal grid cell size.
    if u_r>0
        % Relative ice motion is backward:
        IcethickCorrectionRate_c=(u_r*X_c/domainwidth).*[Icethick_c(2:end)-Icethick_c(1:end-1),.5*(icethick_r-Icethick_c(end))]/dx;
    elseif u_r<0
        % Relative ice motion is forward:
        IcethickCorrectionRate_c=(u_r*X_c/domainwidth).*[.5*(Icethick_c(1)-icethick_l),Icethick_c(2:end)-Icethick_c(1:end-1)]/dx;
    else
        IcethickCorrectionRate_c=zeros(1,xsize);
    end
    
    % Compute a similar grid motion correction for the ice bottom:
    if dodynamicbottom
        if u_r>0
            % Relative ice motion is backward:
            IceBottomCorrectionRate_c=(u_r*X_c/domainwidth).*[IceBottom_c(2:end)-IceBottom_c(1:end-1),.5*(icebottom_r-IceBottom_c(end))]/dx;
        elseif u_r<0
            % Relative ice motion is forward:
            IceBottomCorrectionRate_c=(u_r*X_c/domainwidth).*[.5*(IceBottom_c(1)-icebottom_l),IceBottom_c(2:end)-IceBottom_c(1:end-1)]/dx;
        else
            IceBottomCorrectionRate_c=zeros(1,xsize);
        end
    end
    
    % Calculate grid vertical velocity:
    Wgrid_ud=repmat(W_ud(1,:)+MeltRate_d,[zsize+1,1])+repmat(W_ud(end,:)+Accum_u-MeltRate_u+IcethickCorrectionRate_c,[zsize+1,1]).*repmat(Zhat_ud,[1,xsize]);
    
    % Calculate diffusive temperature terms:
    if dotemp
        DiffusiveTerm2D_v3;
    end
    
    % Calculate advection terms:
    if dotemp
        AdvectiveTerm2D_v4;
    end
    
    % Compute strain heating:
    if dostrainheat
        % Englacial strain heating:
        if strcmp(velocitytype,'SSA')
            %StrainHeat_c=repmat(Viscosity_c.*(2*StrainRate_xy_c.^2+StrainRate_xx_c.^2+StrainRate_yy_c.^2+StrainRate_zz_c.^2),[zsize,1]);
            StrainHeat_c=repmat(Viscosity_c.*StrainRate_xx_c.^2,[zsize,1]);
        else
            StrainHeat_c=Viscosity_c.*(2*StrainRate_xz_c.^2+2*StrainRate_xy_c.^2+StrainRate_xx_c.^2+StrainRate_yy_c.^2+StrainRate_zz_c.^2);
        end
        % Basal shear heating:
        SlideHeat_lrd=Drag_lrd.*U_lrd;
        SlideHeat_lrd(1)=SlideHeat_lrd(2);
        SlideHeat_lrd(end)=SlideHeat_lrd(end-1);
        SlideHeat_d=.5*(SlideHeat_lrd(1:end-1)+SlideHeat_lrd(2:end)); % this doesn't feed back into basal temperature until the next timestep
    else
        StrainHeat_c=zeros(zsize,xsize);
        SlideHeat_d=zeros(1,xsize);
    end
    
    % Assign upwind boundary temperatures: (needed for energy conservation
    % check)
    if dotemp
        % Left boundary:
        if sum(U_lr(:,1).*DZhat_c)>=0
            Temp_l_upwind=Temp_l;
        else
            Temp_l_upwind=Temp_c(:,1);
        end
        % Right boundary:
        if sum((U_lr(:,end)-Ugrid_lr(:,end)).*DZhat_c)>=0
            Temp_r_upwind=Temp_r;
        else
            Temp_r_upwind=Temp_c(:,end);
        end
        % Top boundary:
        Temp_u_upwind=zeros(1,xsize);
        Temp_u_upwind(W_ud(end,:)-Wgrid_ud(end,:)>=0)=Temp_c(end,W_ud(end,:)-Wgrid_ud(end,:)>=0);
        Temp_u_upwind(W_ud(end,:)-Wgrid_ud(end,:)<0)=SurfTemp_u(W_ud(end,:)-Wgrid_ud(end,:)<0);
        % Bottom boundary:
        Temp_d_upwind=zeros(1,xsize);
        Temp_d_upwind(W_ud(1,:)-Wgrid_ud(1,:)>=0)=Temp_d(W_ud(1,:)-Wgrid_ud(1,:)>=0);
        Temp_d_upwind(W_ud(1,:)-Wgrid_ud(1,:)<0)=Temp_c(1,W_ud(1,:)-Wgrid_ud(1,:)<0);
        % Compute total energy inflows: (for energy conservation check later)
        totalenergyinflow=dt*(sum(rho_i*specheat_i*(Temp_l_upwind-tmelt).*U_lr(:,1).*DZ_lr_upwind(:,1)*Width_lr(1)+2*cond_i*(Temp_l-Temp_c(:,1)).*DZ_lr_center(:,1)*Width_lr(1)/dx)... % left side boundary
            +sum(-rho_i*specheat_i*(Temp_r_upwind-tmelt).*(U_lr(:,end)-Ugrid_lr(:,end)).*DZ_lr_upwind(:,end)*Width_lr(end)+2*cond_i*(Temp_r-Temp_c(:,end)).*DZ_lr_center(:,end)*Width_lr(end)/dx)... % right side boundary
            +sum(dx*Width_c.*(rho_i*specheat_i*(Temp_d_upwind-tmelt).*(W_ud(1,:)-Wgrid_ud(1,:))+cond_i*(Temp_d-Temp_c(1,:))./DZ_ud(1,:)))... % bottom boundary
            +sum(dx*Width_c.*(-rho_i*specheat_i*(Temp_u_upwind-tmelt).*(W_ud(end,:)-Wgrid_ud(end,:))+cond_i*(SurfTemp_u-Temp_c(end,:))./DZ_ud(end,:)))); % top boundary
        % Compute total energy in the model domain before changes:
        totalenergy=rho_i*specheat_i*sum(sum((Temp_c-tmelt).*DZ_c,1).*Width_c*dx,2);
    end
    
    % Compute total volume in the model domain before changes:
    totalvolume=sum(Icethick_c.*Width_c*dx);
    
    % Advance temperature to end of timestep:
    if dotemp
        Temp_c=Temp_c+dt*(AdvectiveTerm_c+DiffusiveTerm_c+StrainHeat_c)/(rho_i*specheat_i);
    else
        Temp_c(:)=consttemp;
    end
    
    % Detect unstable ice thickness and kill run:
    if max(abs(del2(Icethick_c,dx)))>maxicethickcurvature
        %error('Ice thickness is exploding.')
        disp('Ice Thickness is Exploding')
        disp(['Model time=',num2str(time_yr),' yr'])
        blewup=1;
        break
    end
    
    % Use a steady-state solution below the diffusion stability limit:
    if dotemp
        % Identify columns with grid cells below the stability limit:
        BelowStabilityLimit_c=Icethick_c<=zsize*maxdensify*2*sqrt(dt*cond_i/(rho_i*specheat_i));
        numunstablecolumns=sum(BelowStabilityLimit_c);
        if numunstablecolumns~=0
            % Loop through grid cells:
            for d2=2:xsize
                % Check that this grid column is unstable:
                if BelowStabilityLimit_c(d2)==0
                    continue
                end
                % Solve for steady state in this column:
                [Temp_c(:,d2),Temp_d(d2),~]=TempProfile1D_column_v4a(Temp_c(:,d2-1),StrainHeat_c(:,d2),U_lr(:,d2:d2+1),W_ud(:,d2),...
                    IsTied_d(d2),MeltRate_d(d2),DZ_c(:,d2),DZ_lr_upwind(:,d2:d2+1),dx,Gflux_c(d2)+SlideHeat_d(d2),SurfTemp_u(d2),MeltPoint_d(d2),ThermalParameters);
            end
        end
    else
        % Set number of unstable columns to zero:
        numunstablecolumns=0;
    end
    
    % Check whether or not to include englacial melting:
    if doenglacialmelting==1 && dotemp
        % Turn sensible heat above the melting point into latent heat:
        % (this is an enthalpy formulation)
        MeltRate_c=(specheat_i/latentheat)*max(0,Temp_c-MeltPoint_c)/dt; % units: 1/s
        % Correct vertical velocity for englacial melting:
        W_ud=W_ud+cumsum([zeros(1,xsize);-MeltRate_c.*DZ_c],1);
    else
        MeltRate_c=zeros(zsize,xsize);
    end
    
    % Truncate temperature at melting point:
    Temp_c=min(Temp_c,MeltPoint_c);
    
    % Advance domain geometry to end of timestep:
    Icethick_c=max(Icethick_c+dt*(W_ud(end,:)-W_ud(1,:)-MeltRate_d+Accum_u-MeltRate_u+IcethickCorrectionRate_c+SideInflux_c./Width_c),minicethick);
    domainwidth=domainwidth+dt*u_r;
    
    % Advance ice bottom to end of timestep:
    if dodynamicbottom
        IceBottom_c(IsGrounded_c==0)=IceBottom_c(IsGrounded_c==0)+dt*(W_ud(1,IsGrounded_c==0)-MeltRate_d(IsGrounded_c==0)+IceBottomCorrectionRate_c(IsGrounded_c==0));
    end
    
    % Do tracer operations:
    if dotracers
        % Interpolate velocity to tracer locations:
        U_t(InDomain_t)=interp2(X_lr,[0;Zhat_c;1],[U_lrd;U_lr;U_lr(end,:)],X_t(InDomain_t),Zhat_t(InDomain_t),interpstyle_t);
        What_t(InDomain_t)=interp2([0,X_c,domainwidth],Zhat_ud,[W_ud(:,1)-Wgrid_ud(:,1),W_ud-Wgrid_ud,W_ud(:,end)-Wgrid_ud(:,end)]./repmat([icethick_l,Icethick_c,icethick_r],[zsize+1,1]),X_t(InDomain_t),Zhat_t(InDomain_t),interpstyle_t);
        % Advance tracer positions to end of timestep:
        X_t(InDomain_t)=X_t(InDomain_t)+dt*U_t(InDomain_t);
        Zhat_t(InDomain_t)=Zhat_t(InDomain_t)+dt*What_t(InDomain_t);
        % Flag tracers that are still in the domain:
        InDomain_t=X_t>=0&X_t<=domainwidth&Zhat_t>=0&Zhat_t<=1;
        % Advance left edge tracer spawn points:
        What_t_l=interp1(Zhat_ud,-ShapeFunctionW_lud*(Accum_u(1)-MeltRate_u(1))/icethick_l,Zhat_t_l,interpstyle_t);
        Zhat_t_l=Zhat_t_l+dt*What_t_l;
        Age_t_l_yr=Age_t_l_yr+dt_yr;
        % Advance tracer age:
        Age_t_yr=Age_t_yr+dt_yr;
        % Spawn new tracers:
        SpawnTracers_v3;
    end
    
    % Advance time to end of timestep:
    time=time+dt;
    time_yr=time/secondsperyear;
    
    % Remember old dx and width:
    olddx=dx;
    OldWidth_c=Width_c;
    OldWidth_lr=Width_lr;
    
    % Check if the calving front needs to jump back:
    % If ice thickness hits the minimum thickness, the calving front jumps
    % back to that location to prevent unrealistic "pinch-out" of the model
    % of the ice shelf.
    if sum(Icethick_c==minicethick)>0
        % Flag jumpback:
        jumpback=1;
        % Remember old grid dimensions:
        olddomainwidth=domainwidth;
        OldX_c=X_c;
        % Find new domain width:
        domainwidth=X_c(find(Icethick_c==minicethick,1,'first')-1);
        % Flag tracers that are no longer in the domain:
        if dotracers
            InDomain_t(X_t>domainwidth)=0;
        end
    else
        % Flag jumpback:
        jumpback=0;
    end
    
    % Compute sill thickness profile:
    if dosill && time<=sillstarttime+sillconstructiontime
        % Forced sill growth:
        SillThick_input=sillamp*max(0,min(1,(time-sillstarttime)/(sillconstructiontime+1)))*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2)); % note stabilizer
    elseif dosill && dosillerosion && time>sillstarttime+sillconstructiontime
        % Sill erosion:
        % interpolate grounded velocity onto input grid:
        U_input=interp1(X_lr,U_lrd.*GroundedFraction_lr,X_input,interpstyle,0);
        % compute erosion rate:
        ErosionRate_input=abs(sillerosionparam*U_input);
        % erode the sill:
        SillThick_input=max(0,SillThick_input-dt*ErosionRate_input);
    end
    
    % Regrid for moving front:
    if strcmp(rightbctype,'front')
        % Regrid:
        dx=domainwidth/xsize;
        X_lr=linspace(0,domainwidth,xsize+1);
        X_c=linspace(.5*dx,domainwidth-.5*dx,xsize);
        % Interpolate bed:
        BedElev_c=interp1(X_input,BedElev_input+SillThick_input,X_c,interpstyle,BedElev_input(end));
        bedelev_l=interp1(X_input,BedElev_input+SillThick_input,0,interpstyle,BedElev_input(1)+SillThick_input(1));
        bedelev_r=interp1(X_input,BedElev_input+SillThick_input,domainwidth,interpstyle,BedElev_input(end)+SillThick_input(end));
        BedElev_lr=[bedelev_l,.5*(BedElev_c(1:end-1)+BedElev_c(2:end)),bedelev_r];
        % Interpolate ice geometry and temperature:
        if jumpback
            Icethick_c=interp1([0,OldX_c,olddomainwidth],[icethick_l,Icethick_c,icethick_r],X_c,interpstyle);
            IceBottom_c=interp1([0,OldX_c,olddomainwidth],[icebottom_l,IceBottom_c,icebottom_r],X_c,interpstyle);
            if dotemp
                Temp_c=interp1([0,OldX_c,olddomainwidth]',[Temp_l,Temp_c,Temp_r]',X_c',interpstyle)';
            end
        end
        % Interpolate width:
        Width_lr=interp1(X_input,Width_input,X_lr,interpstyle,Width_input(end));
        Width_c=.5*(Width_lr(1:end-1)+Width_lr(2:end));
        % Interpolate velocity scale:
        if strcmp(slidingvelocityscale_yr,'file')
            SlidingVelocityScale_lrd=interp1(X_input,SlidingVelocityScale_yr_input,X_lr,interpstyle,SlidingVelocityScale_yr_input(end))/secondsperyear;
        end
        % Interpolate slip exponent:
        if strcmp(m,'file')
            M_lrd=interp1(X_input,M_input,X_lr,interpstyle,M_input(end));
        end
        % Interpolate stress scale:
        if strcmp(slidingstressscale,'file')
            SlidingStressScale_lrd=interp1(X_input,SlidingStressScale_input,X_lr,interpstyle,SlidingStressScale_input(end));
        else
            SlidingStressScale_lrd=slidingstressscale*ones(1,xsize+1);
        end
        % Adjust stress scale for sill:
        if dosill && dosilldrag
            SillThick_lrd=interp1(X_input,SillThick_input,X_lr,interpstyle,0);
            SlidingStressScale_lrd=(1-exp(-SillThick_lrd/silldragthickscale))*sillstressscale+...
                exp(-SillThick_lrd/silldragthickscale).*SlidingStressScale_lrd;
        end
        % Interpolate geothermal flux:
        if strcmp(gflux,'file')
            Gflux_c=interp1(X_input,Gflux_input,X_c,interpstyle,Gflux_input(end));
        end
        % Interpolate side influx:
        if usesidemassinput==1
            % Interpolate flux:
            SideInflux_c=interp1(X_input,SideInflux_input,X_c,interpstyle,SideInflux_input(end));
            % Interpolate original surface:
            S0_c=interp1(X_input(1:lasticeind),SurfElev_input(1:lasticeind),X_c,interpstyle,SurfElev_input(lasticeind));
            % Adjust flux for surface changes:
            SideInflux_c=SideInflux_c.*max(0,1+(S0_c-SurfElev_c)./(S0_c+1)); % note stabilizer
        end
    elseif dosill
        % Interpolate bed: (sill thickness may have changed)
        BedElev_c=interp1(X_input,BedElev_input+SillThick_input,X_c,interpstyle,BedElev_input(end));
        bedelev_l=interp1(X_input,BedElev_input+SillThick_input,0,interpstyle,BedElev_input(1)+SillThick_input(1));
        bedelev_r=interp1(X_input,BedElev_input+SillThick_input,domainwidth,interpstyle,BedElev_input(end)+SillThick_input(end));
        BedElev_lr=[bedelev_l,.5*(BedElev_c(1:end-1)+BedElev_c(2:end)),bedelev_r];
        % Interpolate stress scale and adjust for sill:
        if dosill && dosilldrag
            % Interpolate stress scale:
            if strcmp(slidingstressscale,'file')
                SlidingStressScale_lrd=interp1(X_input,SlidingStressScale_input,X_lr,interpstyle,SlidingStressScale_input(end));
            else
                SlidingStressScale_lrd=slidingstressscale*ones(1,xsize+1);
            end
            % Adjust stress scale for sill:
            SillThick_lrd=interp1(X_input,SillThick_input,X_lr,interpstyle,0);
            SlidingStressScale_lrd=(1-exp(-SillThick_lrd/silldragthickscale))*sillstressscale+...
                exp(-SillThick_lrd/silldragthickscale).*(SlidingStressScale_lrd-sillstressscale);
        end
    end
    
    % Check volume (mass) conservation:
    % Compute total volume in the model domain after changes:
    newtotalvolume=sum(Icethick_c.*Width_c*dx); % uses new dx and width
    % Compute total volume inflows and sources/sinks:
    totalvolumeinflow=dt*(sum(U_lr(:,1).*DZ_lr_upwind(:,1))*OldWidth_lr(1)-sum((U_lr(:,end)-Ugrid_lr(:,end)).*DZ_lr_upwind(:,end))*OldWidth_lr(end)... % side boundaries (uses old width)
        +sum((Accum_u-MeltRate_u-MeltRate_d).*OldWidth_c*olddx)); % top/bottom boundaries (uses old dx and width)
    totalvolumesources=dt*sum(sum(-MeltRate_c.*DZ_c,1).*OldWidth_c*olddx,2); % old dx and width
    % Compute volume error:
    normalizedvolumeerror=(newtotalvolume-(totalvolume+totalvolumeinflow+totalvolumesources))/totalvolume;
    %disp(['Normalized mass misfit=',num2str(normalizedvolumeerror)])
    
    % Check energy conservation:
    if dotemp
        % Compute total energy in the model domain after changes:
        newtotalenergy=rho_i*specheat_i*sum(sum((Temp_c-tmelt).*DZ_c,1).*Width_c*dx,2); % uses new dx and width
        % Compute energy sources/sinks:
        totalenergysources=dt*(sum(sum(StrainHeat_c.*DZ_c,1).*OldWidth_c*dx,2)-sum(sum(rho_i*latentheat*MeltRate_c.*DZ_c,1).*OldWidth_c*dx,2));
        % Compute energy error:
        normalizedenergyerror=(newtotalenergy-(totalenergy+totalenergyinflow+totalenergysources))/totalenergy;
        %disp(['Normalized energy misfit=',num2str(normalizedenergyerror)])
    end
    
    % Record tracer locations:
    if time>=nextrecordtime-.5*recordinterval && recordedtracers==0 && dotracers
        ModelRecord.ID_t=ID_t;
        ModelRecord.X_t=X_t;
        ModelRecord.Zhat_t=Zhat_t;
        ModelRecord.Age_t_yr=Age_t_yr;
        ModelRecord.InDomain_t=InDomain_t;
        ModelRecord.IsOriginal_t=IsOriginal_t;
        ModelRecord.IsAcc_t=IsAcc_t;
        ModelRecord.IsLooseEnd_t=IsLooseEnd_t;
        ModelRecord.Connectivity_t=Connectivity_t;
        ModelRecord.DeadConnector_t=DeadConnector_t;
        recordedtracers=1;
    end
    
    % Compute quantities needed for time series output:
    % compute model volume: (m^3)
    totalvolume=sum(dx*Width_c.*Icethick_c);
    % compute volume above flotation: (m^3)
    FlotationThickness_c=max(0,-(rho_sw/rho_i)*BedElev_c);
    vaf=sum(dx*Width_c.*max(0,Icethick_c-FlotationThickness_c));
    % compute volume fluxes: (m^3/yr)
    accuminflux_yr=sum(dx*Width_c.*Accum_u)*secondsperyear;
    ablationoutflux_yr=(sum(dx*Width_c.*MeltRate_u)+Width_lr(end)*(surfelev_r-max(icebottom_r,0))*subaerialmeltrate_r)*secondsperyear;
    if hasshelf
        oceanmeltoutflux_yr=(sum(DS_c_plume(1:xsize_plume).*Width_c_plume(1:xsize_plume).*MeltRate_c_plume(1:xsize_plume))+...
            Width_lr(end)*(-icebottom_r)*oceanmeltrate_r)*secondsperyear;
    else
        oceanmeltoutflux_yr=Width_lr(end)*(-icebottom_r)*oceanmeltrate_r*secondsperyear;
    end
    calvingflux_yr=Width_lr(end)*icethick_r*calvingrate_r*secondsperyear;
    sideinflux_yr=sum(dx*SideInflux_c)*secondsperyear;
    
    % Compute overlap between timestep and record interval:
    fracinrecord=1-max(0,(time-nextrecordtime)/dt)-max(0,((nextrecordtime-recordinterval)-(time-dt))/dt);
    
    % Record this timestep:
    % Thermal records:
    if dotemp
        ModelRecord.Temp_c=ModelRecord.Temp_c+Temp_c*dt*fracinrecord/recordinterval;
        ModelRecord.Temp_d=ModelRecord.Temp_d+Temp_d*dt*fracinrecord/recordinterval;
        ModelRecord.MeltRate_c=ModelRecord.MeltRate_c+MeltRate_c*dt*fracinrecord/recordinterval;
        ModelRecord.StrainHeat_c=ModelRecord.StrainHeat_c+StrainHeat_c*dt*fracinrecord/recordinterval;
        ModelRecord.SlideHeat_d=ModelRecord.SlideHeat_d+SlideHeat_d*dt*fracinrecord/recordinterval;
    end
    % Domain geometry records:
    ModelRecord.Icethick_c=ModelRecord.Icethick_c+Icethick_c*dt*fracinrecord/recordinterval;
    ModelRecord.IceBottom_c=ModelRecord.IceBottom_c+IceBottom_c*dt*fracinrecord/recordinterval;
    ModelRecord.domainwidth=ModelRecord.domainwidth+domainwidth*dt*fracinrecord/recordinterval;
    ModelRecord.x_gl=ModelRecord.x_gl+x_gl*dt*fracinrecord/recordinterval;
    % Basal hydrology records:
    if dotemp
        ModelRecord.MeltRate_d=ModelRecord.MeltRate_d+MeltRate_d*dt*fracinrecord/recordinterval;
        ModelRecord.HydroHeat_d=ModelRecord.HydroHeat_d+HydroHeat_d*dt*fracinrecord/recordinterval;
        ModelRecord.WaterFlux_lrd=ModelRecord.WaterFlux_lrd+WaterFlux_lrd*dt*fracinrecord/recordinterval;
    end
    % Plume model records:
    if doplume
        ModelRecord.Thick_lr_plume=ModelRecord.Thick_lr_plume+Thick_lr_plume*dt*fracinrecord/recordinterval;
        ModelRecord.U_lr_plume=ModelRecord.U_lr_plume+U_lr_plume*dt*fracinrecord/recordinterval;
        ModelRecord.Temp_lr_plume=ModelRecord.Temp_lr_plume+Temp_lr_plume*dt*fracinrecord/recordinterval;
        ModelRecord.Salinity_lr_plume=ModelRecord.Salinity_lr_plume+Salinity_lr_plume*dt*fracinrecord/recordinterval;
        ModelRecord.MeltRate_c_plume=ModelRecord.MeltRate_c_plume+MeltRate_c_plume*dt*fracinrecord/recordinterval;
        ModelRecord.plumeseparationx=ModelRecord.plumeseparationx+plumeseparationx*dt*fracinrecord/recordinterval;
        ModelRecord.plumeseparationz=ModelRecord.plumeseparationz+plumeseparationz*dt*fracinrecord/recordinterval;
        if strcmp(plumesolvertype,'custom')
            ModelRecord.meanplumeiterations=ModelRecord.meanplumeiterations+meanplumeiterations*dt*fracinrecord/recordinterval;
        end
    end
    % Sill records:
    if dosill && dosillerosion
        ModelRecord.SillThick_input=ModelRecord.SillThick_input+SillThick_input*dt*fracinrecord/recordinterval;
    end
    % Velocity related records:
    ModelRecord.U_lr=ModelRecord.U_lr+U_lr*dt*fracinrecord/recordinterval;
    ModelRecord.U_lrd=ModelRecord.U_lrd+U_lrd*dt*fracinrecord/recordinterval;
    ModelRecord.u_r=ModelRecord.u_r+u_r*dt*fracinrecord/recordinterval;
    ModelRecord.calvingrate_r=ModelRecord.calvingrate_r+calvingrate_r*dt*fracinrecord/recordinterval;
    ModelRecord.subaerialmeltrate_r=ModelRecord.subaerialmeltrate_r+subaerialmeltrate_r*dt*fracinrecord/recordinterval;
    ModelRecord.oceanmeltrate_r=ModelRecord.oceanmeltrate_r+oceanmeltrate_r*dt*fracinrecord/recordinterval;
    ModelRecord.W_ud=ModelRecord.W_ud+W_ud*dt*fracinrecord/recordinterval;
    ModelRecord.Drag_lrd=ModelRecord.Drag_lrd+Drag_lrd*dt*fracinrecord/recordinterval;
    ModelRecord.Viscosity_c=ModelRecord.Viscosity_c+Viscosity_c*dt*fracinrecord/recordinterval;
    if strcmp(velocitytype,'fullstokes')
        ModelRecord.PressureDynamic_c=ModelRecord.PressureDynamic_c+PressureDynamic_c*dt*fracinrecord/recordinterval;
    end
    % Surface climate records:
    ModelRecord.Accum_u=ModelRecord.Accum_u+Accum_u*dt*fracinrecord/recordinterval;
    ModelRecord.MeltRate_u=ModelRecord.MeltRate_u+MeltRate_u*dt*fracinrecord/recordinterval;
    if dotemp
        ModelRecord.SurfTemp_u=ModelRecord.SurfTemp_u+SurfTemp_u*dt*fracinrecord/recordinterval;
    end
    % Numerical records:
    ModelRecord.numstokesiterations=ModelRecord.numstokesiterations+(triedreset*maxiterations+iteration_stokes)*dt*fracinrecord/recordinterval;
    ModelRecord.numunstablecolumns=ModelRecord.numunstablecolumns+numunstablecolumns*dt*fracinrecord/recordinterval;
    ModelRecord.normalizedvolumeerror=ModelRecord.normalizedvolumeerror+normalizedvolumeerror*dt*fracinrecord/recordinterval;
    if dotemp
        ModelRecord.normalizedenergyerror=ModelRecord.normalizedenergyerror+normalizedenergyerror*dt*fracinrecord/recordinterval;
    end
    % Time series:
    ModelTimeSeries.Volume(recordcounter)=ModelTimeSeries.Volume(recordcounter)+totalvolume*dt*fracinrecord/recordinterval;
    ModelTimeSeries.VAF(recordcounter)=ModelTimeSeries.VAF(recordcounter)+vaf*dt*fracinrecord/recordinterval;
    ModelTimeSeries.Accumulation(recordcounter)=ModelTimeSeries.Accumulation(recordcounter)+accuminflux_yr*dt*fracinrecord/recordinterval;
    ModelTimeSeries.Ablation(recordcounter)=ModelTimeSeries.Ablation(recordcounter)+ablationoutflux_yr*dt*fracinrecord/recordinterval;
    ModelTimeSeries.OceanMelt(recordcounter)=ModelTimeSeries.OceanMelt(recordcounter)+oceanmeltoutflux_yr*dt*fracinrecord/recordinterval;
    ModelTimeSeries.Calving(recordcounter)=ModelTimeSeries.Calving(recordcounter)+calvingflux_yr*dt*fracinrecord/recordinterval;
    ModelTimeSeries.SideInflux(recordcounter)=ModelTimeSeries.SideInflux(recordcounter)+sideinflux_yr*dt*fracinrecord/recordinterval;
    
    % Save record:
    if time>=nextrecordtime
        
        % Record elapsed time:
        ModelRecord.time_yr=(nextrecordtime-.5*recordinterval)/secondsperyear;  % time is midpoint of record
        
        % Record presence of closed basins:
        ModelRecord.hasbasins=hasbasins;
        
        % Create specific model record structure:
        prefix='0'*ones(1,numdigits(end)-floor(log10(recordcounter))-1);
        idnumber=num2str(recordcounter);
        eval(['ModelRecord_',prefix,idnumber,'=ModelRecord;'])
        
        % Add specific model record structure to output file:
        numrecords(end)=recordcounter;
        save(outputfile,['ModelRecord_',prefix,idnumber],'numrecords','numdigits','-append')
        clear(['ModelRecord_',prefix,idnumber])
        
        % Take elements of time series that are duplicates of model record:
        ModelTimeSeries.Time(recordcounter)=ModelRecord.time_yr;
        ModelTimeSeries.X_gl(recordcounter)=ModelRecord.x_gl;
        ModelTimeSeries.Domainwidth(recordcounter)=ModelRecord.domainwidth;
        
        % Reset record counting variables:
        nextrecordtime=nextrecordtime+recordinterval;
        hasbasins=0;
        recordcounter=recordcounter+1;
        recordedtracers=0;
        
        % Compute overlap between major timestep and (next) record interval:
        fracinrecord=1-max(0,(time-nextrecordtime)/dt)-max(0,((nextrecordtime-recordinterval)-(time-dt))/dt);
        
        % Reset records:  (record contribution of this timestep to next record interval):
        % Thermal records:
        if dotemp
            ModelRecord.Temp_c=Temp_c*dt*fracinrecord/recordinterval;
            ModelRecord.Temp_d=Temp_d*dt*fracinrecord/recordinterval;
            ModelRecord.MeltRate_c=MeltRate_c*dt*fracinrecord/recordinterval;
            ModelRecord.StrainHeat_c=StrainHeat_c*dt*fracinrecord/recordinterval;
            ModelRecord.SlideHeat_d=SlideHeat_d*dt*fracinrecord/recordinterval;
        end
        % Domain geometry records:
        ModelRecord.Icethick_c=Icethick_c*dt*fracinrecord/recordinterval;
        ModelRecord.IceBottom_c=IceBottom_c*dt*fracinrecord/recordinterval;
        ModelRecord.domainwidth=domainwidth*dt*fracinrecord/recordinterval;
        ModelRecord.x_gl=x_gl*dt*fracinrecord/recordinterval;
        if dotemp
            % Basal hydrology records:
            ModelRecord.MeltRate_d=MeltRate_d*dt*fracinrecord/recordinterval;
            ModelRecord.HydroHeat_d=HydroHeat_d*dt*fracinrecord/recordinterval;
            ModelRecord.WaterFlux_lrd=WaterFlux_lrd*dt*fracinrecord/recordinterval;
        end
        % Plume model records:
        if doplume
            ModelRecord.Thick_lr_plume=Thick_lr_plume*dt*fracinrecord/recordinterval;
            ModelRecord.U_lr_plume=U_lr_plume*dt*fracinrecord/recordinterval;
            ModelRecord.Temp_lr_plume=Temp_lr_plume*dt*fracinrecord/recordinterval;
            ModelRecord.Salinity_lr_plume=Salinity_lr_plume*dt*fracinrecord/recordinterval;
            if strcmp(plumesolvertype,'custom')
                ModelRecord.meanplumeiterations=meanplumeiterations*dt*fracinrecord/recordinterval;
            end
            ModelRecord.MeltRate_c_plume=MeltRate_c_plume*dt*fracinrecord/recordinterval;
            ModelRecord.plumeseparationx=plumeseparationx*dt*fracinrecord/recordinterval;
            ModelRecord.plumeseparationz=plumeseparationz*dt*fracinrecord/recordinterval;
        end
        % Sill records:
        if dosill && dosillerosion
            ModelRecord.SillThick_input=SillThick_input*dt*fracinrecord/recordinterval;
        end
        % Velocity related records:
        ModelRecord.U_lr=U_lr*dt*fracinrecord/recordinterval;
        ModelRecord.U_lrd=U_lrd*dt*fracinrecord/recordinterval;
        ModelRecord.u_r=u_r*dt*fracinrecord/recordinterval;
        ModelRecord.calvingrate_r=calvingrate_r*dt*fracinrecord/recordinterval;
        ModelRecord.subaerialmeltrate_r=subaerialmeltrate_r*dt*fracinrecord/recordinterval;
        ModelRecord.oceanmeltrate_r=oceanmeltrate_r*dt*fracinrecord/recordinterval;
        ModelRecord.W_ud=W_ud*dt*fracinrecord/recordinterval;
        ModelRecord.Drag_lrd=Drag_lrd*dt*fracinrecord/recordinterval;
        ModelRecord.Viscosity_c=Viscosity_c*dt*fracinrecord/recordinterval;
        if strcmp(velocitytype,'fullstokes')
            ModelRecord.PressureDynamic_c=PressureDynamic_c*dt*fracinrecord/recordinterval;
        end
        % Surface climate records:
        ModelRecord.Accum_u=Accum_u*dt*fracinrecord/recordinterval;
        ModelRecord.MeltRate_u=MeltRate_u*dt*fracinrecord/recordinterval;
        if dotemp
            ModelRecord.SurfTemp_u=SurfTemp_u*dt*fracinrecord/recordinterval;
        end
        % Numerical records:
        ModelRecord.numstokesiterations=(triedreset*maxiterations+iteration_stokes)*dt*fracinrecord/recordinterval;
        ModelRecord.numunstablecolumns=numunstablecolumns*dt*fracinrecord/recordinterval;
        ModelRecord.normalizedvolumeerror=normalizedvolumeerror*dt*fracinrecord/recordinterval;
        if dotemp
            ModelRecord.normalizedenergyerror=normalizedenergyerror*dt*fracinrecord/recordinterval;
        end
    end
    
    % Display progress on command line:
    if displaycounter>=displayinterval
        if verbose
            disp('...')
            disp(strcat('Model Time=',num2str(time_yr),' yr'))
            disp(strcat('Real Time=',num2str(toc(t2)),' s'))
            if spinuptimedone==0
                disp(['dt=',num2str(dt_yr),'yr'])
            end
        end
        displaycounter=1;
    else
        displaycounter=displaycounter+1;
    end
    
    % Lengthen timesteps during the spinup period:
    if spinuptimedone==0 && time<spinuptime
        % Adjust timesteps:
        dt=dt_basic*((1+1/rapidfactor)/2)+.5*(dt_basic-dt_basic/rapidfactor)*(-cos(pi*time/spinuptime));
        dt_yr=dt/secondsperyear;
    elseif spinuptimedone==0 % finish the spniup period:
        % Adjust timesteps:
        dt=dt_basic;
        dt_yr=dt/secondsperyear;
        spinuptimedone=1;
        displaycounter=1;
        % Display on command line:
        if verbose
            disp('...')
            disp(strcat('Time=',num2str(time_yr),'yr'))
            disp('Spinup Period Complete:')
            disp(['dt=',num2str(dt_yr),'yr'])
        end
    end
    
    % Debugging Display Section:  (uncomment plotting commands as needed)
    % if rem(time_yr,.25)<dt_yr*.999
        %     figure(1)
        %     if timestep==1
        %         hold off
        %         plot(X_c/1000,(Icethick_c-InitialIcethick_c)/time_yr,'r')
        %         hold on
        %         xlim([x_gl,domainwidth]/1000)
        %         xlabel('Distance (km)')
        %         ylabel('Average Thickness Change Rate (m/yr)')
        %     elseif rem(timestep,10)==0
        %         plot(X_c/1000,(Icethick_c-InitialIcethick_c)/time_yr,'k')
        %         title(['t=',num2str(time_yr),' yr'])
        %     end
        %     if exist('h3','var')
        %         set(h3,'XData',[get(h3,'XData'),time_yr])
        %         set(h3,'YData',[get(h3,'YData'),normalizedvolumeerror])
        %     else
        %         h3=plot(time_yr,normalizedvolumeerror,'.r');
        %     end
        %     title(['Normalized Volume Error, T=',num2str(time_yr),' yr'])
        %     if timestep>=2
        %         delete(h1)
        %         delete(hgl)
        %     end
        %     h1=plot(X_c/1000,W_ud(end,:)*secondsperyear,'r');
        %     hgl=plot(x_gl*[1,1]/1000,[-1,1],'--r');
        %     title(['Vertical Velocity, T=',num2str(time_yr),' yr'])
        %     xlim([0,domainwidth/1000])
        %     plot(time_yr,x_gl/1000,'dk','MarkerFaceColor','k')
        %     title(['Grounding Line Position, T=',num2str(time_yr),' yr'])
        %     T_yr_disp(timestep+1)=time_yr;
        %     Icethick_spot_disp(timestep+1)=interp1(X_c,Icethick_c,spot_disp);
        %     plot(T_yr_disp,Icethick_spot_disp,'k')
        %     xlabel('Time (yr)')
        %     ylabel('Thickness (m)')
        %     xlim([0,23])
        %     ylim([400,1000])
        %     hold off
        %     imagesc(X_c/1000,Zhat_c,Temp_c-MeltPoint_c)
        %     hold on
        %     plot([x_gl,x_gl]/1000,[0,1],'w')
        %     set(gca,'ydir','normal')
        %     caxis([min(SurfTemp_u-tmelt),0])
        %     xlim([0,500])
        %     ylim([0,1])
        %     xlabel('Distance (km)')
        %     ylabel('Elevation (unitless)')
        %     title(['Temperature, T=',num2str(time_yr),' yr'])
        %         if exist('h1','var')
        %             delete([h1,h2,h3,h4])
        %         end
        %         if exist('h5','var')
        %             delete(h5)
        %         end
        %         h1=plot([0,X_c(X_c<=x_gl),x_gl,X_c(X_c>x_gl),domainwidth]/1000,[surfelev_l,SurfElev_c(X_c<=x_gl),surfelev_gl,SurfElev_c(X_c>x_gl),surfelev_r],'r');
        %         h2=plot([0,X_c(X_c<=x_gl),x_gl,X_c(X_c>x_gl),domainwidth]/1000,[icebottom_l,IceBottom_c(X_c<=x_gl),bedelev_gl,IceBottom_c(X_c>x_gl),icebottom_r],'r');
        %         h3=plot(domainwidth*[1,1]/1000,[surfelev_r,icebottom_r],'r');
        %         h4=plot(x_gl*[1,1]/1000,[surfelev_gl,bedelev_gl],'r');
        %         if dosill
        %             h5=plot(X_input(SillThick_input>silldragthickscale)/1000,BedElev_input(SillThick_input>silldragthickscale)+SillThick_input(SillThick_input>silldragthickscale),'b');
        %         end
        %         title(['Model Geometry, T=',num2str(time_yr),' yr'])
        %     hold off
        %     imagesc(X_lr/domainwidth,Zhat_c,U_lr*secondsperyear)
        %     hold on
        %     plot([1,1]*x_gl/domainwidth,[0,1],'w')
        %     set(gca,'ydir','normal')
        %     caxis([0,4e3])
        %     Ubar_lr=sum(U_lr.*repmat(DZhat_c,[1,xsize+1]),1);
        %     GLind_c_disp(timestep+1)=interp1(X_lr,Ubar_lr,spot_disp);
        %     plot(T_yr_disp,GLind_c_disp*secondsperyear,'k')
        %     xlim([0,23])
        %     ylim([0,4000])
        %     xlabel('Time (yr)')
        %     ylabel('Velocity, (m/yr)')
        %     GLind_c_disp(timestep+1)=find(IsGrounded_c,1,'last');
        %     NumIterations_disp(timestep+1)=iteration_stokes;
        %     subplot(2,1,1)
        %     plot(T_yr_disp,GLind_c_disp,'k')
        %     xlim([0,23])
        %     ylim([60,80])
        %     title(['Grounding Line Index, T=',num2str(time_yr),' yr'])
        %     xlabel('Time (yr)')
        %     ylabel('Index')
        %     subplot(2,1,2)
        %     plot(T_yr_disp,NumIterations_disp,'k')
        %     xlim([0,23])
        %     ylim([0,30])
        %     xlabel('Time (yr)')
        %     ylabel('Number')
        %     title(['Number of stokes iterations, T=',num2str(time_yr),' yr'])
        % drawnow
        % pause(.03)
        %     frame=getframe(gcf);
        %     writeVideo(Movie,frame);
    % end
    
    %     % Plume calibration code: (runs with FlowlinePlumeCalibrate_v1.m)
    %     if time_yr>=plumecalibrationtime_yr
    %         lastind=find(X_c>=x_gl&Icethick_c==max(Icethick_c(X_c>=x_gl)),1,'first');
    %         if lastind~=lastgroundedind+1
    %             refx=X_c(lastind);
    %             reficethick=Icethick_c(lastind);
    %         else
    %             refx=x_gl;
    %             reficethick=icethick_gl;
    %         end
    %         NormalizedX=(X_c(lastind+1:end)-refx)/(domainwidth-refx);
    %         NormalizedY=(reficethick-Icethick_c(lastind+1:end))/(reficethick-icethick_r);
    %         A=[ones(sum(NormalizedX<=.2),1),NormalizedX(NormalizedX<=.2)'];
    %         equation=(A'*A)\A'*NormalizedY(NormalizedX<=.2)';
    %         finalslope=equation(2);
    %         return
    %     end
    
    % Finish Model Run:
    realtime=clock;
    thishour=realtime(4)+realtime(5)/60+realtime(6)/3600;
    if time>=runtime || (thishour>shutoffhour && thishour<restarthour)
        % Set "done" variable:
        done=1;
    else
        timestep=timestep+1;
    end
    
end

%% Finish Up:

% Record the final model conditions:
% Basic final state:
FinalConditions=struct('time_yr',time_yr ... % elapsed time record
    ,'Icethick_c',Icethick_c,'IceBottom_c',IceBottom_c,'domainwidth',domainwidth,'x_gl',x_gl ... % ice thickness, domain size records
    ,'U_lr',U_lr,'U_lrd',U_lrd,'u_r',u_r,'calvingrate_r',calvingrate_r,'subaerialmeltrate_r',subaerialmeltrate_r,'oceanmeltrate_r',oceanmeltrate_r,'W_ud',W_ud,'Drag_lrd',Drag_lrd,'Viscosity_c',Viscosity_c... % velocity related records
    ,'Accum_u',Accum_u,'MeltRate_u',MeltRate_u... % surface climate records
    ,'numstokesiterations',triedreset*maxiterations+iteration_stokes,'numunstablecolumns',numunstablecolumns,'normalizedvolumeerror',normalizedvolumeerror);  % numerical records
% Add thermal and hydrology records:
if dotemp
    FinalConditions.Temp_c=Temp_c;
    FinalConditions.Temp_d=Temp_d;
    FinalConditions.MeltRate_c=MeltRate_c;
    FinalConditions.StrainHeat_c=StrainHeat_c;
    FinalConditions.SlideHeat_d=SlideHeat_d;
    FinalConditions.MeltRate_d=MeltRate_d;
    FinalConditions.HydroHeat_d=HydroHeat_d;
    FinalConditions.WaterFlux_lrd=WaterFlux_lrd;
    FinalConditions.hasbasins=hasbasins;
    FinalConditions.SurfTemp_u=SurfTemp_u;
    FinalConditions.normalizedenergyerror=normalizedenergyerror;
end
% Add tracer records:
if dotracers
    FinalConditions.ID_t=ID_t;
    FinalConditions.X_t=X_t;
    FinalConditions.Zhat_t=Zhat_t;
    FinalConditions.Age_t_yr=Age_t_yr;
    FinalConditions.InDomain_t=InDomain_t;
    FinalConditions.IsOriginal_t=IsOriginal_t;
    FinalConditions.IsAcc_t=IsAcc_t;
    FinalConditions.IsLooseEnd_t=IsLooseEnd_t;
    FinalConditions.Connectivity_t=Connectivity_t;
    FinalConditions.DeadConnector_t=DeadConnector_t;
end
% Add plume records:
if doplume
    FinalConditions.Thick_lr_plume=Thick_lr_plume;
    FinalConditions.U_lr_plume=U_lr_plume;
    FinalConditions.Temp_lr_plume=Temp_lr_plume;
    FinalConditions.Salinity_lr_plume=Salinity_lr_plume;
    FinalConditions.MeltRate_c_plume=MeltRate_c_plume;
    if strcmp(plumesolvertype,'custom')
        FinalConditions.meanplumeiterations=meanplumeiterations;
    end
    FinalConditions.plumeseparationx=plumeseparationx;
    FinalConditions.plumeseparationz=plumeseparationz;
end
% Sill records:
if dosill && dosillerosion
    FinalConditions.SillThick_input=SillThick_input;
end
% Add dynamic pressure:
if strcmp(velocitytype,'fullstokes')
    FinalConditions.PressureDynamic_c=PressureDynamic_c;
end

% Compute total computer time:
computertime=computertime+toc(t2);

% Save the final state:
numrecords(end)=recordcounter-1;
save(outputfile,'FinalConditions','ModelTimeSeries','blewup','numrecords','computertime','-append')

% Final display for this model run:
if verbose
    disp(strcat('This Run Computation Time=',num2str(computertime),' s'))
end

% catch errorid
%     
%     % Save final things:
%     if exist(outputfile,'file') && exist('recordcounter','var')
%         computertime=computertime+toc(t2);
%         numrecords(end)=recordcounter-1;
%         save(outputfile,'ModelTimeSeries','numrecords','computertime','errorid','-append')
%     end
%     
%     % Display the error message and express my anger:
%     disp('!!!')
%     disp('Hit an Error Message.')
%     disp(errorid.identifier)
%     disp(errorid.message)
%     disp('...')
%     disp('Error Found In:')
%     for ii=1:length(errorid.stack)
%         disp(['Script=',errorid.stack(ii).name,', Line=',num2str(errorid.stack(ii).line)])
%     end
%     if exist('time_yr','var')
%         disp('...')
%         disp(['Model time=',num2str(time_yr),' yr'])
%     end
%     disp('...')
%     disp('Model Run Failed.')
%     disp('DAMMIT!')
%     disp(strcat('This Run Computation Time=',num2str(toc(t2)),' s'))
%     disp('...')
%     %quit
%     
% end

%end

% Final display:
disp('Finished This Run')
%disp(strcat('This Run Real Time=',num2str(toc(t1)),' s'))
%quit