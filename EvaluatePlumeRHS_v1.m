function Derivs=EvaluatePlumeRHS_v1(s,Values,OtherInputs)

% EvaluatePlumeRHS_v1

% Mike Wolovick, 6/8/2016

% This function returns the right-hand-side for the plume model when I'm
% using Matlab's built-in ODE solvers.

% This function can be modified to only produce the RHS or to solve the
% "mass matrix" problem and output individual derivatives.

% Unpack other inputs:
unpack(OtherInputs)
unpack(PlumeParameters)
unpack(ThermalParameters)
% Unpack Values:
thick=Values(1);
u=Values(2);
temp=Values(3)+tmelt;
sal=Values(4);
% Interpolate forcings:
ambienttemp=TempInterpolant(s);
ambientsal=SalInterpolant(s);
pressure=PressureInterpolant(s);
sintheta=SinThetaInterpolant(s);
if dotemp
    condflux=CondFluxInterpolant(s);
end
% Compute freezing temperature of plume water:
meltpoint=tmelt-meltingpointslope*pressure-salmeltcoeff*sal;
% Compute melt rate:
if dotemp
    meltrate=(rho_sw*specheat_sw*u*stantonnumber*(temp-meltpoint)-condflux)/(rho_sw*latentheat);
else
    meltrate=specheat_sw*u*stantonnumber*(temp-meltpoint)/(latentheat+specheat_i*(meltpoint-consttemp));
end
% Compute entrainment rate:
entrainment=e0*abs(u*sintheta);
% Compute normalized density contrast: (negative for buoyant plume)
deltarho=PlumeParameters.beta*(sal-ambientsal)-PlumeParameters.alpha*(temp-ambienttemp);
% Compute mass-gradient:
massgradient=entrainment+meltrate;
% Compute momentum gradient:
momentumgradient=-thick*deltarho*g*sintheta-plumedragcoeff*(u^2);
% Compute heat gradient:
heatgradient=ambienttemp*entrainment+meltpoint*meltrate+stantonnumber*u*(meltpoint-temp);
% Compute salt gradient:
saltgradient=ambientsal*entrainment;
% Assemble the right-hand side:
RHS=[massgradient;momentumgradient;heatgradient;saltgradient];
% Set up system of equations to solve for individual derivatives:
M=[u,thick,0,0;u^2,2*u*thick,0,0;temp*u,temp*thick,u*thick,0;sal*u,sal*thick,0,u*thick];
% Solve for derivatives:
Derivs=M\RHS;
end


