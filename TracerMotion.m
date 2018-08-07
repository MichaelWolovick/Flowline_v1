function TracerVelocity=TracerMotion(t,TracerPosition,X_c,Y_c,Zhat_c,DZhat_c,U_c,V_c,What_c,SlidingFraction_c)

% TracerMotion

% Mike Wolovick, 12/3/2011

% This function evaluates the RHS for tracer advection by linearly
% interpolating a velocity field onto the positions of the tracers.
% Horizontal position and velocity have units of m and m/s, but vertical
% position is unitless and vertical velocity has units of 1/s.  

% This code expects velocities to be defined on grid centers.  It uses
% no-gradient velocity BC for the upper surface and it determines velocity
% along the base from a given sliding fraction.

% What_c has two expanded fields in the third dimension for velocity at the
% surface and the bed.

% Because this function expects to run inside of one of matlab's ODE
% integrators, it expects to see this indexing scheme:
% TracerPosition=[TracerX;TracerY;TracerZhat];
% TracerVelocity=[TracerU;TracerV;TracerWhat];

% In other words, the TracerPosition and TracerVelocity vectors must have a
% length divisible by 3.

% Check input:
if rem(length(TracerPosition),3)~=0
    error('Input "TracerPosition" must have a length divisible by 3')
end
numtracers=length(TracerPosition)/3;

% Unpackage input:
TracerX=TracerPosition(1:numtracers);
TracerY=TracerPosition(numtracers+1:2*numtracers);
TracerZhat=TracerPosition(2*numtracers+1:3*numtracers);

% Find grid sizes:
[ysize,xsize,zsize]=size(U_c);
xlims=[X_c(1,1),X_c(1,end)];
ylims=[Y_c(1,1),Y_c(end,1)];

% Find nans:
NanInd=isnan(TracerX)|isnan(TracerY)|isnan(TracerZhat)|TracerX>xlims(2)|TracerX<xlims(1)|TracerY>ylims(2)|TracerY<ylims(1)|TracerZhat>1|TracerZhat<0;

% Pre-allocate velocity:
TracerU=zeros(numtracers,1);
TracerV=zeros(numtracers,1);
TracerWhat=zeros(numtracers,1);
TracerU(NanInd==1)=NaN;
TracerV(NanInd==1)=NaN;
TracerWhat(NanInd==1)=NaN;

% Interpolate velocity:
TracerU(NanInd==0)=interp3(repmat(X_c,[1,1,zsize+2]),repmat(Y_c,[1,1,zsize+2]),repmat(cat(3,0,Zhat_c,1),[ysize,xsize,1]),cat(3,SlidingFraction_c.*sum(U_c.*repmat(DZhat_c,[ysize,xsize,1]),3),U_c,U_c(:,:,end)),TracerX(NanInd==0),TracerY(NanInd==0),TracerZhat(NanInd==0),'linear');
TracerV(NanInd==0)=interp3(repmat(X_c,[1,1,zsize+2]),repmat(Y_c,[1,1,zsize+2]),repmat(cat(3,0,Zhat_c,1),[ysize,xsize,1]),cat(3,SlidingFraction_c.*sum(V_c.*repmat(DZhat_c,[ysize,xsize,1]),3),V_c,V_c(:,:,end)),TracerX(NanInd==0),TracerY(NanInd==0),TracerZhat(NanInd==0),'linear');
TracerWhat(NanInd==0)=interp3(repmat(X_c,[1,1,zsize+2]),repmat(Y_c,[1,1,zsize+2]),repmat(cat(3,0,Zhat_c,1),[ysize,xsize,1]),What_c,TracerX(NanInd==0),TracerY(NanInd==0),TracerZhat(NanInd==0),'linear');

% Package velocity:
TracerVelocity=[TracerU;TracerV;TracerWhat];


