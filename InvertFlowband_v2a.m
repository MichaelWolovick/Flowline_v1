% InvertFlowband_v1

% Mike Wolovick, 4/5/2016

% This script takes the output of ProduceFlowband_v2 and inverts for basal
% sliding and side drag parameters.  The script also puts the flowband
% variables into the format that the model is expecting to see, flips the
% distance coordinate (so it is now pointing in the direction of flow).

% Velocity inversion: 
% Solve for constant side drag coefficient, constant ice rheology
% coefficient (for long stresses), and variable sliding coefficient.  The
% inversion assumes a linear sliding law for simplicity.  The constraints
% are the velocity and strain rate observations, a smoothness constraint on
% the sliding coefficient, and expected values of the longitudinal rheology
% coefficient.  I do not constrain the side drag coefficient to allow for
% shear marign weakening.

% The inversion strategy is a genetic algorithm (evolution!).  

% The inversion calls ShallowShelfApproximation_2D_v1 in it's forward
% model.

% The inputs from the flowband are assumed to lie at grid edges.

% v2a:  the evolutionary algorithm computes a single weighted mean model at
% the end of each generation to serve as the seed of the next generation,
% rather than pairing off models into couples and doing sexual
% reproduction.  The old way worked, but I hope the new way has a better
% signal-to-noise ratio and therefore converges faster.  


clear all
tic


%% Parameters:

% File namesand paths:
inputfolder='/home/mjw/Documents/FjordSillPaper/Flowbands/';
outputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/WeightingTest2/';
%outputfolder='/work/Michael.Wolovick/FjordSillPaper/ModelInput/WeightingTest/';
lastfile='/work/Michael.Wolovick/FjordSillPaper/ModelInput/WeightingTest/KangerFlowband_v2_inverted2.mat';
% inputfiles={'HelheimFlowband_v2.mat';...
%     'HelheimFlowband_v3.mat';...
%     'HelheimFlowband_v4.mat';...
%     'KangerFlowband_v1.mat';...
%     'KangerFlowband_v2.mat';...
%     'KangerFlowband_v3.mat';...
%     'JakobshavnFlowband_v2.mat';...
%     'JakobshavnFlowband_v3.mat';...
%     'JakobshavnFlowband_v4.mat';...
%     'PetermannFlowband_v2.mat';...
%     'PetermannFlowband_v3.mat';...
%     'PetermannFlowband_v4.mat';...
%     'PineIslandFlowband_v1.mat';...
%     'PineIslandFlowband_v2.mat';...
%     'PineIslandFlowband_v3.mat';...
%     'PineIslandFlowband_v4.mat';...
%     'PineIslandFlowband_v5.mat';...
%     'ThwaitesFlowband_v1.mat';...
%     'ThwaitesFlowband_v2.mat'};
inputfiles={'KangerFlowband_v2.mat'};

% Physical parameters:
rho_i=917;                        % kg/m^3
rho_sw=1028;                      % kg/m^3
g=9.8;                            % m/s^2
n=3;                              % unitless
secondsperyear=60*60*24*365.25;   % s/yr
tmelt=273.15;                     % K
defaulta_long=5e-25;              % Pa^-n*s^-1
defaulta_side=5e-25;              % Pa^-n*s^-1
deltaa_long=30;                   % unitless>>1

% Inversion parameters:
AllDataWeights=[.3:.1:1]';        % unitless (0,1)
inversionwavelength=1e4;          % m
strainratewavelength=1e5;         % m (for scale only)
populationsize=100;               % integer
numgenerations=200;               % integer
startmutationfactor=3;            % unitless>1
minmutationfactor=1.05;           % unitless>=1
efoldinggen=25;                   % unitless
annealingfactor=4;                % unitless
forceconstlonga=0;                % logical
costweightingexponent=2;          % unitless

% Forward model parameters:
dosidedrag=1;                     % 0 or 1
useinputsidea=1;                  % 0 or 1
useinputa=1;                      % 0 or 1
doimplicitsurface=0;              % 0 or 1
calvingtype='fixed';              % 'fixed', 'thick', or 'file'
leftbctype='flux';                % 'flux' 
rightbctype='front';              % 'front'
tolerance=1e-3;                   % unitless
damping_visc=.75;                 % unitless [0,1)
damping_drag=.5;                  % unitless [0,1)
miniterations=2;                  % integer
maxiterations=100;                % integer
solvertype='direct';              % 'direct', 'pcg', 'bicg', 'bicgstab', 'bicgstabl', 'cgs', 'gmres', 'lsqr', 'minres', 'qmr', 'symmlq', 'tfqmr'  (names of matlab iterative solvers)
viscosityguess=1e15;              % Pa*s
dragcoefficientguess=1e10;        % Pa/(m/s)
sidedragcoefficientguess=1e10;    % Pa/(m/s)
strainratestabilizer_yr=1e-6;     % 1/yr
slidingstabilizer_yr=1e-3;        % m/yr
singularwarning='on';             % 'on' or 'off'
nearsingularwarning='off';        % 'on' or 'off'
ForwardModelParameters=struct('n',n,'dosidedrag',dosidedrag,'useinputsidea',useinputsidea,'useinputa',useinputa,'doimplicitsurface',doimplicitsurface,'calvingtype',calvingtype,...
    'leftbctype',leftbctype,'rightbctype',rightbctype,'tolerance',tolerance,'damping_visc',damping_visc,'damping_drag',damping_drag,'miniterations',miniterations,'maxiterations',maxiterations,...
    'solvertype',solvertype,'viscosityguess',viscosityguess,'dragcoefficientguess',dragcoefficientguess,'sidedragcoefficientguess',sidedragcoefficientguess,...
    'strainratestabilizer',strainratestabilizer_yr/secondsperyear,'slidingstabilizer',slidingstabilizer_yr/secondsperyear,'singularwarning',singularwarning,'nearsingularwarning',nearsingularwarning);



%% Preparation:

% Pre-allocate costs:
% DataCost=zeros(size(AllDataWeights));
% PriorCost=zeros(size(AllDataWeights));

% Load last output file:
load(lastfile,'AllDataWeights','DataCost','PriorCost')

% Loop through input files:
for filenumber=1:length(AllDataWeights)
    
    if DataCost(filenumber)~=0
        continue
    end
    
% Set this data weight:
dataweight=AllDataWeights(filenumber);
    
% % Define output file name:
% outputfile=[outputfolder,inputfiles{filenumber}(1:end-4),'_inverted.mat'];
% 
% % Display file name:
% inputfiles{filenumber}

% Define output file name:
outputfile=[outputfolder,inputfiles{1}(1:end-4),'_inverted',num2str(filenumber),'.mat'];

% Display weighting parameter:
disp(['Data Weight=',num2str(dataweight)])

% try/catch pair:
try

% Control the ill-conditioned warnings:
warning(nearsingularwarning,'MATLAB:nearlySingularMatrix')
warning(singularwarning,'MATLAB:SingularMatrix')

% Control the negative data warning:
warning('off','MATLAB:Axes:NegativeDataInLogAxis')

% Ensure that population size is divisible by 4:
populationsize=4*ceil(populationsize/4);

% Load input file:
load([inputfolder,inputfiles{1}])

% Flip and rename variables: (some get _input suffix, some don't)
X_input=fliplr(Distance_flowband');
BedElev_input=fliplr(BedElev_flowband');
SurfElev=fliplr(SurfElev_flowband');
Icethick_input=fliplr(IceThick_flowband');
VelMag_yr_input=fliplr(VelMag_flowband');
SurfTemp=fliplr(SurfTemp_flowband')+tmelt;
Accum_yr=fliplr(Accum_flowband');
Gflux_input=fliplr(Gflux_flowband');
IceFrac=fliplr(IceFrac_flowband');
GroundedFrac=fliplr(GroundedFrac_flowband');
Influx_input=fliplr(Influx_flowband');

% Deal with width variables:
if exist('AreaConservingWidth_flowband','var')
    Width_ac_input=fliplr(AreaConservingWidth_flowband');
    Width_input=fliplr(FluxConservingWidth_flowband');
else
    Width_ac_input=fliplr(Width_flowband');
    Width_input=fliplr(Width_flowband');
end

% Recompute distance coordinate after flipping:
X_input=X_input(1)-X_input;

% Clear extraneous variables:
clear *_flowband

% Locate calving front:
lasticeind=find(isnan(Icethick_input)==0,1,'last')-1; % last grid center
domainwidth=X_input(lasticeind+1);

% Interpolate to grid centers:
X_c=.5*(X_input(1:end-1)+X_input(2:end));
BedElev_c=.5*(BedElev_input(1:end-1)+BedElev_input(2:end));
Icethick_c=.5*(Icethick_input(1:end-1)+Icethick_input(2:end));

% Compute hydraulic head:
HydroHead_c=BedElev_c+(rho_i/rho_sw)*Icethick_c;

% Locate grounding line:
if HydroHead_c(lasticeind)<0
    hasshelf=1;
    x_gl=interp1(HydroHead_c(isnan(HydroHead_c)==0),X_c(isnan(HydroHead_c)==0),0,'spline');
    lastgroundedind=find(X_c<=x_gl,1,'last');
else
    hasshelf=0;
    x_gl=X_input(lasticeind+1);
    lastgroundedind=lasticeind;
end

% Compute surface and bottom elevations:
SurfElev_c=max(BedElev_c+Icethick_c,(1-rho_i/rho_sw)*Icethick_c);
IceBottom_c=SurfElev_c-Icethick_c;

%%  Create the environment needed by ShallowShelfApproximation:

% Grid dimensions:
Environment.xsize=lasticeind;
Environment.zsize=1;
Environment.domainwidth=domainwidth;
Environment.dx=steplength;
Environment.lastgroundedind=lastgroundedind;

% Ice geometry:
Environment.Icethick_c=Icethick_c(1:lasticeind);
Environment.SurfElev_c=SurfElev_c(1:lasticeind);
Environment.IceBottom_c=IceBottom_c(1:lasticeind);
Environment.Icethick_lr_center=Icethick_input(1:lasticeind+1);
Environment.Icethick_lr_upwind=Icethick_input(1:lasticeind+1);
Environment.Width_c=.5*(Width_input(1:lasticeind)+Width_input(2:lasticeind+1));
Environment.Width_lr=Width_input(1:lasticeind+1);

% DZ variables:
Environment.DZhat_c=1;
Environment.DZ_c=Icethick_c(1:lasticeind);
Environment.DZ_lr_center=Icethick_input(1:lasticeind+1);
Environment.DZ_lr_upwind=Icethick_input(1:lasticeind+1);

% Boundary conditions:
Environment.icethick_l=Icethick_input(1);
Environment.icethick_r=Icethick_input(lasticeind+1);
Environment.NormalStress_r=-.5*rho_i*g*(1-rho_i/rho_sw)*Icethick_input(lasticeind+1);
Environment.U_l=VelMag_yr_input(1)/secondsperyear;

% Melt rates:
Environment.MeltRate_d=zeros(1,lasticeind);
Environment.MeltRate_c=zeros(1,lasticeind);

% Driving stress gradient:
Environment.DrivingStressGrad_lr=-rho_i*g*[SurfElev_c(2)-SurfElev_c(1),SurfElev_c(2:lasticeind)-SurfElev_c(1:lasticeind-1),0]/steplength;

% Grounding info:
if hasshelf
    Environment.IsGrounded_c=[true(1,lastgroundedind),false(1,lasticeind-lastgroundedind)];
    Environment.IsGrounded_lr=[ones(1,lastgroundedind),.5,zeros(1,lasticeind-lastgroundedind)];
    Environment.GroundedFraction_lr=[ones(1,lastgroundedind),(x_gl-X_c(lastgroundedind))/steplength,zeros(1,lasticeind-lastgroundedind)];
else
    Environment.IsGrounded_c=true(1,lasticeind);
    Environment.IsGrounded_lr=ones(1,lasticeind+1);
    Environment.GroundedFraction_lr=ones(1,lasticeind+1);
end

% Guess of side drag coefficient:
DragCoefficient_lr=((defaulta_side*Width_input(1:lasticeind+1)).^(-1/n)).*((VelMag_yr_input(1:lasticeind+1)/secondsperyear).^((1-n)/n));

% Target velocity:
Environment.U_lr_data=VelMag_yr_input(1:lasticeind+1)/secondsperyear;

% Drag coefficient guess:
DragCoefficientGuess_lrd=abs(Environment.DrivingStressGrad_lr.*Environment.Icethick_lr_center./Environment.U_lr_data);

% Characteristic scales:
Environment.viscositylogscale=sqrt(log(((defaulta_long^(-1/n))*(((VelMag_yr_input(lasticeind)-VelMag_yr_input(1))/(secondsperyear*domainwidth))^((1-n)/n)))/((defaulta_long^-1)*(mean(abs(Environment.DrivingStressGrad_lr.*Environment.Icethick_lr_center))^(1-n))))^2+1);
Environment.dragcoefflogscale=log(max(DragCoefficientGuess_lrd(2:lastgroundedind))/min(DragCoefficientGuess_lrd(2:lastgroundedind))); 
Environment.sidedragcoefflogscale=sqrt(log(max(max(DragCoefficient_lr(:,2:end-1)))/min(min(DragCoefficient_lr(:,2:end-1))))^2+1);
Environment.drivingstressscale=mean(Environment.DrivingStressGrad_lr(2:lastgroundedind).*Environment.Icethick_lr_center(2:lastgroundedind));

% Compute a variable strain rate scale:
Environment.StrainRate_xx_c_data=(Environment.U_lr_data(1,2:end)-Environment.U_lr_data(1,1:end-1))/steplength;
Environment.StrainRateScale_c=exp(intuitive_lowpass(log(abs(Environment.StrainRate_xx_c_data)),strainratewavelength/steplength));

% A_long parameters:
Environment.defaulta_long=defaulta_long;
Environment.deltaa_long=deltaa_long;

%% Solve the Inverse Problem:

% Note: mutation amplitude for A is multiplied by n. (not any more)

% Create the initial input vector:
InitialVector=[log(defaulta_long);log(defaulta_side);log(DragCoefficientGuess_lrd(2:lastgroundedind+1)')];

% Create an initial population:
PopulationVectors=zeros(lastgroundedind+2,populationsize);
PopulationVectors(1:2,:)=repmat(InitialVector(1:2),[1,populationsize])+log(startmutationfactor)*randn(2,populationsize);
for ii=1:populationsize
    PopulationVectors(3:end,ii)=InitialVector(3:end)+log(startmutationfactor)*(randn(1)+makerednoise(lastgroundedind,inversionwavelength/steplength)');
end

% Pre-allocate things:
ViscosityGuess_c=viscosityguess*ones(1,lasticeind); 
DragCoefficientGuess_lr=sidedragcoefficientguess*ones(1,lasticeind+1); 
OutputViscosity_c=viscosityguess*ones(lasticeind,populationsize); % x in d1, pop in d2
OutputDragCoefficient_lr=sidedragcoefficientguess*ones(lasticeind+1,populationsize); % x in d1, pop in d2
Cost=zeros(1,populationsize);

% Create figure 1:
figure(1)

% Loop through generations:
for generation=1:numgenerations
    % Force a constant longitudinal rheology:
    if forceconstlonga
        PopulationVectors(1,:)=log(defaulta_long);
    end
    % Evaluate cost function for each population member:
    for ii=1:populationsize
        [Cost(ii),ThisViscosityGuess_c,ThisDragCoefficientGuess_lr,~]=InvertFlowbandCostFunction_v1(PopulationVectors(:,ii),Environment,ForwardModelParameters,...
            ViscosityGuess_c,DragCoefficientGuess_lr,dataweight,inversionwavelength);
        OutputViscosity_c(:,ii)=ThisViscosityGuess_c';
        OutputDragCoefficient_lr(:,ii)=ThisDragCoefficientGuess_lr';
    end
    % Produce a single consensus result: (cost-weighted average)
    totalweight=sum(Cost.^(-costweightingexponent));
    FinalVector=(1/totalweight)*sum(PopulationVectors./repmat(Cost.^costweightingexponent,[lastgroundedind+2,1]),2);
    ViscosityGuess_c=(1/totalweight)*sum(OutputViscosity_c./repmat(Cost.^costweightingexponent,[lasticeind,1]),2)';
    DragCoefficientGuess_lr=(1/totalweight)*sum(OutputDragCoefficient_lr./repmat(Cost.^costweightingexponent,[lasticeind+1,1]),2)';
    % Quite here on the final generation:
    if generation==numgenerations
        break
    end
    % Create a new population using the consensus result plus random noise:
    PopulationVectors(1:2,:)=repmat(FinalVector(1:2,:),[1,populationsize])+max(log(minmutationfactor),log(startmutationfactor)*exp(-(generation-efoldinggen*annealingfactor*floor(generation/(efoldinggen*annealingfactor)))/efoldinggen))*randn(2,populationsize);
    for ii=1:populationsize
        PopulationVectors(3:end,ii)=FinalVector(3:end)+max(log(minmutationfactor),log(startmutationfactor)*exp(-(generation-efoldinggen*annealingfactor*floor(generation/(efoldinggen*annealingfactor)))/efoldinggen))*(randn(1)+makerednoise(lastgroundedind,inversionwavelength/steplength)');
    end
    % Sort by cost:
    [~,SortedInd]=sort(Cost,'ascend');
    % Find inverted index:
    InvertedInd=linspace(1,populationsize,populationsize);
    InvertedInd=InvertedInd(SortedInd);
    % Quite here on the final generation:
    if generation==numgenerations
        break
    end
    % Plot cost function:
    %figure(1)
    plot(generation*ones(1,populationsize),Cost,'.k')
    hold on
    ylim([0,1])
    xlim([0,10*ceil(generation/10)])
    drawnow
end

% Evaluate consensus misfit and velocity field:
[finalcost,~,~,U_lr]=InvertFlowbandCostFunction_v1(FinalVector,Environment,ForwardModelParameters,...
    ViscosityGuess_c,DragCoefficientGuess_lr,dataweight,inversionwavelength);

% Finish plotting cost function:
he=plot(generation*ones(1,populationsize),Cost,'.k');
hm=plot([0,generation],finalcost*[1,1],'--k');
ylim([0,1])
xlim([0,10*ceil(generation/10)])
legend([he(1);hm],'Ensemble Members','Final Ensemble Mean','location','NorthEast')
xlabel('Generation Number')
ylabel('Cost Function')
title('Convergence of Evolutionary Algorithm')
drawnow

% Save figure 1:
set(gcf,'PaperSize',[7,6])
set(gcf,'PaperPosition',[0,0,7,6])
figname=[outputfile(1:end-4),'_costfunction.png'];
print('-dpng',figname,'-r300')

%% Organize and Save Results of Inversion:

% Compute drag and strain rate:
unpack(Environment)
DragCoefficient_lrd=zeros(1,xsize+1);
DragCoefficient_lrd(2:lastgroundedind+1)=exp(FinalVector(3:end)');
Drag_lrd=U_lr.*DragCoefficient_lrd;
DrivingStress_lr=DrivingStressGrad_lr.*Icethick_lr_center;
StrainRate_xx_c=(U_lr(1,2:end)-U_lr(1,1:end-1))/dx;

% Convert drag to characteristic stress scale:
SlidingStressScale_input=zeros(size(X_input));
SlidingStressScale_input(2:lastgroundedind)=Drag_lrd(2:lastgroundedind);
SlidingStressScale_input(1)=Drag_lrd(2);

% Extrapolate stress scale beyond the grounding line:
glstress=Drag_lrd(lastgroundedind);
glstressgrad=(Drag_lrd(lastgroundedind)-Drag_lrd(lastgroundedind-1))/dx;
farfieldstress=glstress+inversionwavelength*glstressgrad;
SlidingStressScale_input(lastgroundedind+1:end)=farfieldstress+(glstress-farfieldstress)*exp(-(X_input(lastgroundedind+1:end)-X_input(lastgroundedind))/inversionwavelength);

% Convert velocity to characteristic scale:
SlidingVelocityScale_yr_input=zeros(size(X_input));
SlidingVelocityScale_yr_input(1:lastgroundedind)=U_lr(1:lastgroundedind)*secondsperyear;

% Extrapolate stress scale beyond the grounding line:
glu=U_lr(lastgroundedind);
glugrad=(U_lr(lastgroundedind)-U_lr(lastgroundedind-1))/dx;
farfieldu=glu+inversionwavelength*glugrad;
SlidingVelocityScale_yr_input(lastgroundedind+1:end)=secondsperyear*(farfieldu+(glu-farfieldu)*exp(-(X_input(lastgroundedind+1:end)-X_input(lastgroundedind))/inversionwavelength));

%% Make Profile Figure:

% Compute ensemble velocity and drag:
U_lr_population=zeros(xsize+1,populationsize);
Drag_lrd_population=zeros(xsize+1,populationsize);
SideDrag_lr_population=zeros(xsize+1,populationsize);
for ii=1:populationsize
    % Compute velocity:
    [~,~,~,ThisU]=InvertFlowbandCostFunction_v1(PopulationVectors(:,ii),Environment,ForwardModelParameters,...
        ViscosityGuess_c,DragCoefficientGuess_lr,dataweight,inversionwavelength);
    U_lr_population(:,ii)=ThisU';
    % Compute drag:
    Drag_lrd_population(2:lastgroundedind,ii)=exp(PopulationVectors(3:end-1,ii)).*ThisU(2:lastgroundedind)';
end
% Compute strain rate:
StrainRate_xx_c_population=(U_lr_population(2:end,:)-U_lr_population(1:end-1,:))/dx;
% Compute side drag:
SideDrag_lr=(U_lr./(exp(FinalVector(2))*Width_lr)).^(1/n);
SideDrag_lr_population=(U_lr_population./(repmat(exp(PopulationVectors(2,:)),[xsize+1,1]).*repmat(Width_lr',[1,populationsize]))).^(1/n);

%% Make figure 2:
figure(2)
% Velocity:
subplot(2,2,1)
h11=semilogy(X_input(1:lasticeind+1)/1000,U_lr_population*secondsperyear,'Color',[.5,.5,.5]);
hold on
h12=semilogy(X_input(1:lasticeind+1)/1000,U_lr*secondsperyear,'k');
h13=semilogy(X_input(1:lasticeind+1)/1000,U_lr_data*secondsperyear,'r');
legend([h11(1);h12;h13],'Ensemble Members','Ensemble Mean','Observations','location','SouthEast')
xlim([0,domainwidth/1000])
ylim([1,1e4])
xlabel('Distance (km)')
ylabel('Velocity (m/yr)')
title('Ice Velocity')
% Strain rate:
subplot(2,2,2)
h21=semilogy(.5*(X_input(1:lasticeind)+X_input(2:lasticeind+1))/1000,abs(StrainRate_xx_c_population)*secondsperyear*1000,'Color',[.5,.5,.5]);
hold on
h22=semilogy(.5*(X_input(1:lasticeind)+X_input(2:lasticeind+1))/1000,abs(StrainRate_xx_c)*secondsperyear*1000,'k');
h23=semilogy(.5*(X_input(1:lasticeind)+X_input(2:lasticeind+1))/1000,abs(StrainRate_xx_c_data)*secondsperyear*1000,'r');
legend([h21(1);h22;h23],'Ensemble Members','Ensemble Mean','Observations','location','NorthWest')
xlim([0,domainwidth/1000])
ylim([1e-3,1e3])
xlabel('Distance (km)')
ylabel('Strain Rate (1/ka)')
title('Longitudinal Strain Rate')
% Stress plot:
subplot(2,2,3)
h31=plot(X_input(2:lastgroundedind)/1000,Drag_lrd_population(2:lastgroundedind,:)/1000,'Color',[.5,.5,.5]);
hold on
h32=plot(X_input(1:lastgroundedind)/1000,Drag_lrd(1:lastgroundedind)/1000,'k');
plot(X_input(1:lasticeind+1)/1000,(2*SideDrag_lr_population.*repmat(Icethick_lr_center'./Width_lr',[1,populationsize]))/1000,'Color',[.75,.75,1])
h33=plot(X_input(1:lasticeind+1)/1000,(2*SideDrag_lr.*Icethick_lr_center./Width_lr)/1000,'b');
h34=plot(X_input(2:lasticeind+1)/1000,DrivingStress_lr(2:end)/1000,'r');
legend([h32;h33;h34],'Basal Drag','Side Drag (*2H/W)','Driving Stress','location','NorthWest')
xlim([0,domainwidth/1000])
ylim([0,1e5*ceil(max([max(DrivingStress_lr),max(max(Drag_lrd_population)),max(max(2*SideDrag_lr_population.*repmat(Icethick_lr_center'./Width_lr',[1,populationsize])))])/1e5)]/1000)
xlabel('Distance (km)')
ylabel('Stress (kPa)')
title('Force Balance')
% Drag coefficient plot:
subplot(2,2,4)
h41=semilogy(X_input(2:lastgroundedind)/1000,exp(PopulationVectors(3:end-1,:)),'Color',[.5,.5,.5]);
hold on
h42=semilogy(X_input(2:lastgroundedind)/1000,exp(FinalVector(3:end-1)),'k');
h43=semilogy(X_input(2:lastgroundedind)/1000,DragCoefficientGuess_lrd(2:lastgroundedind),'r');
legend([h41(1);h42;h43],'Ensemble Members','Ensemble Mean','\tau_{d}/u_{obs}','location','SouthWest')
xlim([0,domainwidth/1000])
xlabel('Distance (km)')
ylabel('Drag Coefficient (Pa/(m/s))')
title('Drag Coefficient')

%% Save figure 2:
set(gcf,'PaperSize',[7,6])
set(gcf,'PaperPosition',[0,0,7,6])
figname=[outputfile(1:end-4),'_profiles.png'];
print('-dpng',figname,'-r300')

% Display some things on the command line:
disp(['a_long=',num2str(exp(FinalVector(1)))])
disp(['a_side=',num2str(exp(FinalVector(2)))])
velocitycost=max(sqrt(mean((log(U_lr./U_lr_data)/log(max(U_lr_data)/min(U_lr_data))).^2)),...
    sqrt(mean(((U_lr-U_lr_data)./U_lr_data).^2))) % assumes velocity is always positive
strainratecost=sqrt(mean(((StrainRate_xx_c-StrainRate_xx_c_data)./StrainRateScale_c).^2))
priorcost=(1/3)*(sqrt(mean((del2(Drag_lrd(2:lastgroundedind),dx)*(inversionwavelength^2)/drivingstressscale).^2))...
    +sqrt(mean((del2(FinalVector(3:end-1),dx)*(inversionwavelength^2)/dragcoefflogscale).^2))...
    +((FinalVector(1)-log(defaulta_long))/log(deltaa_long))^2)

% Record costs:
DataCost(filenumber)=.5*(velocitycost+strainratecost);
PriorCost(filenumber)=priorcost;

% Save output:
save(outputfile)

% Close figures:
close all

catch errorid
    
    disp('Hit an Error Message.')
    disp(errorid.identifier)
    disp(errorid.message)
    disp('...')
    disp('Error Found In:')
    for ii=1:length(errorid.stack)
        disp(['Script=',errorid.stack(ii).name,', Line=',num2str(errorid.stack(ii).line)])
    end
    continue
    
end

end

% Plot data cost against prior cost:
figure(3)
plot(PriorCost,DataCost,'Color','k','MarkerFaceColor','k','Marker','o')
for ii=1:length(DataCost)
    text(PriorCost(ii),DataCost(ii),num2str(AllDataWeights(ii)),'HorizontalAlignment','left','VerticalAlignment','bottom','Color','k')
end
xlim([0,1])
ylim([0,1])
xlabel('Prior Cost (unitless)')
ylabel('Data Cost (unitless)')
title({['Inversion Tradeoff for ',inputfiles{1}(1:end-4)];'(numbers indicate data weighting)'},'Interpreter','none')
set(gcf,'PaperSize',[7,6])
set(gcf,'PaperPosition',[0,0,7,6])
figname=[outputfolder,inputfiles{1}(1:end-4),'_inversion_tradeoff.png'];
print('-dpng',figname,'-r300')

% Final display:
disp('Done!')
toc