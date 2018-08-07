function [deltarho,isterminal,direction]=TerminatePlumeModel_v1(s,Values,OtherInputs)

% TerminatePlumeModel_v1

% Mike Wolovick, 6/8/2016

% This function determines whether the plume model has reached neutral
% buoyancy and detached from the ice face when I'm using Matlab's built-in
% ODE solvers to do the plume model.  All it does is compute the density
% contrast; Matlab considers the zero-crossing to be an "event" that
% terminates the plume.

% Unpack other inputs:
unpack(OtherInputs)
unpack(PlumeParameters)
unpack(ThermalParameters)
% Unpack Values:
temp=Values(3)+tmelt;
sal=Values(4);
% Interpolate ambient properties:
ambienttemp=TempInterpolant(s);
ambientsal=SalInterpolant(s);
% Compute normalized density contrast: (negative for buoyant plume)
deltarho=PlumeParameters.beta*(sal-ambientsal)-PlumeParameters.alpha*(temp-ambienttemp);
% This is terminal:
isterminal=1;
direction=0;
end