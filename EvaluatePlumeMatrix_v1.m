function M=EvaluatePlumeMatrix_v1(s,Values)

% EvaluatePlumeMatrix_v1

% Mike Wolovick, 6/8/2016

% This function returns the "mass matrix" for the plume model when I'm
% using Matlab's built-in ODE solvers.  The "mass matrix" converts the
% right-hand side into individual derivatives for D,U,T,S

% Unpack Values:
thick=Values(1);
u=Values(2);
temp=Values(3);
sal=Values(4);
% Set up system of equations to solve for individual derivatives:
M=[u,thick,0,0;u^2,2*u*thick,0,0;temp*u,temp*thick,u*thick,0;sal*u,sal*thick,0,u*thick];
end