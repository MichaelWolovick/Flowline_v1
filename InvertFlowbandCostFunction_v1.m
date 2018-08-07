function [cost,Viscosity_c,DragCoefficient_lr,U_lr]=InvertFlowbandCostFunction_v1(InputVector,Environment,ForwardModelParameters,Viscosity_c,DragCoefficient_lr,dataweight,xwavelength)

% This is the cost function for the sliding inversion in InvertFlowband_v2.

% Note that the values within InputVector are assumed to be on a log scale.

% The function also outputs viscosity and side drag coefficient to make the
% iteration go faster.

% Unpack structures:
unpack(Environment)
unpack(ForwardModelParameters)

% Pre-allocate variables:
DragCoefficient_lrd=zeros(1,xsize+1);
U_lr=zeros(1,xsize+1);
Drag_lrd=zeros(1,xsize+1);

% Parse the input:
a_input=exp(InputVector(1));
a_side_input=exp(InputVector(2));
DragCoefficient_lrd(2:lastgroundedind+1)=exp(InputVector(3:end)');

% Adjust drag coefficient for grounded fraction:
DragCoefficient_lrd=DragCoefficient_lrd.*GroundedFraction_lr;

% Unnecessary sliding constant variable:
SlidingConstant_lrd=1./DragCoefficient_lrd;

% Set variables:
timestep=1;
FlowDir_lr=true(1,xsize+1);
IsInternal_lrd=[false(1),true(1,xsize-1),false(1)];
NewlyGrounded_lr=false(1,xsize+1);
m=1;

% Compute velocity field:
ShallowShelfApproximation_2D_v1

% Evaluate prior cost:
dragcost=sqrt(mean((del2(Drag_lrd(2:lastgroundedind),dx)*(xwavelength^2)/drivingstressscale).^2));
dragcoeffcost=sqrt(mean((del2(InputVector(3:end-1),dx)*(xwavelength^2)/dragcoefflogscale).^2));
acost=((log(a_input)-log(defaulta))/log(deltaa))^2+((log(min(a_side_input,defaulta))-log(defaulta))/log(deltaa))^2;
priorcost=(1/3)*(dragcost+dragcoeffcost+acost);

% Evaluate data cost:
velocitycost=max(sqrt(mean((log(U_lr./U_lr_data)/log(max(U_lr_data)/min(U_lr_data))).^2)),...
    sqrt(mean(((U_lr-U_lr_data)./U_lr_data).^2))); % assumes velocity is always positive
strainratecost=sqrt(mean(((StrainRate_xx_c-StrainRate_xx_c_data)./StrainRateScale_c).^2));
datacost=.5*(velocitycost+strainratecost);

% Compute total cost:
cost=dataweight*datacost+(1-dataweight)*priorcost;
