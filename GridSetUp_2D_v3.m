% GridSetUp_2D_v2

% Mike Wolovick, 10/21/2011

% This script builds the skewed grid (parallel to dimensionless elevation)
% for FTM_2D.  The grid is a linear spacing between the surface and the
% bed.  This means that dz is not a constant value, and that grid edge
% incidence angles need to be calculated (although cos(theta)>.99 in most
% cases).

% This script also interpolates variables to grid cell edges, using either
% upwind interpolation or linear interpolation as appropriate.

% GridElev_ud(1,:) represents the top of the water layer.

% Changes for v2:  the script now handles variable side BC and sediments.
% It also has the option for 'loosefixedicethick' right side BC

% Changes for v3:  Runs with Flowline_v1 now.  No more sediment grid.
% Deals with grounding lines and such.


%% Deal with Grounding Lines and Ice Surface:

% Check if we're doing dynamic or hydrostatic bottom:
if dodynamicbottom==0 || firstgrid
    
    % Compute and interpolate hydraulic head: (no gradient BC) (ocean density)
    HydroHead_c=BedElev_c+(rho_i/rho_sw)*Icethick_c;
    HydroHead_lr=[HydroHead_c(1),.5*(HydroHead_c(1:end-1)+HydroHead_c(2:end)),HydroHead_c(end)];
    
    % Record old grounded status at grid edges:
    if firstgrid==0
        IsGrounded_lr_last=IsGrounded_lr;
    end
    
    % Identify grounded cells:
    IsGrounded_c=HydroHead_c>sealevel;
    IsGrounded_lr=[IsGrounded_c(1),.5*(IsGrounded_c(1:end-1)+IsGrounded_c(2:end)),IsGrounded_c(end)];
    
    % Count grounded cells:
    numgroundedinds=sum(IsGrounded_c);
    
    % Identify newly grounded cells: (important for drag coefficient)
    if firstgrid==0
        NewlyGrounded_lr=IsGrounded_lr~=0&IsGrounded_lr_last==0;
    end
    
    % Compute ice bottom:
    IceBottom_c(IsGrounded_c)=BedElev_c(IsGrounded_c);
    if firstgrid || dodynamicbottom==0
        IceBottom_c(IsGrounded_c==0)=sealevel-(rho_i/rho_sw)*Icethick_c(IsGrounded_c==0);
    end
    
    % Compute grounded fraction: (staggered for velocity solver grid)
    % Pre-allocate:
    GroundedFraction_lr=zeros(1,xsize+1);
    % Do fully grounded cells:
    GroundedFraction_lr(IsGrounded_lr==1)=1;
    % Identify partially grounded cells:
    partiallygroundedinds=find(IsGrounded_lr==.5);
    % Loop through partially grounded cells:
    if isempty(partiallygroundedinds)==0
        for thisind=1:length(partiallygroundedinds) % must loop forward!
            % Check grounding line interpolation style:
            if strcmp(glinterp,'cubic') && partiallygroundedinds(thisind)<xsize
                % Evaluate 3 key hydraulic heat differences:
                hydroheaddiffupstream=HydroHead_c(partiallygroundedinds(thisind)-1)-HydroHead_c(partiallygroundedinds(thisind)-2);
                hydroheaddiffcenter=HydroHead_c(partiallygroundedinds(thisind))-HydroHead_c(partiallygroundedinds(thisind)-1);
                hydroheaddiffdownstream=HydroHead_c(partiallygroundedinds(thisind)+1)-HydroHead_c(partiallygroundedinds(thisind));
                % Compute cubic polynomial coefficients:
                thesecoeffs=[hydroheaddiffupstream-2*hydroheaddiffcenter+hydroheaddiffdownstream,...
                    -2*hydroheaddiffupstream+3*hydroheaddiffcenter-hydroheaddiffdownstream,...
                    hydroheaddiffupstream,...
                    HydroHead_c(partiallygroundedinds(thisind)-1)];
                % Compute roots of polynomial:
                theseroots=roots(thesecoeffs);
                % Check roots:
                if sum(theseroots>=0&theseroots<=1&imag(theseroots)==0)==1
                    % Use cubic function for grounding line:
                    x_gl=X_c(partiallygroundedinds(thisind)-1)+dx*theseroots(theseroots>=0&theseroots<=1&imag(theseroots)==0);
                else
                    % Use linear interpolation instead:
                    x_gl=X_c(partiallygroundedinds(thisind)-1)+dx*(-HydroHead_c(partiallygroundedinds(thisind)-1)/hydroheaddiffcenter);
                end
            else
                % Linear interpolation:
                hydroheaddiffcenter=HydroHead_c(partiallygroundedinds(thisind))-HydroHead_c(partiallygroundedinds(thisind)-1);
                x_gl=X_c(partiallygroundedinds(thisind)-1)+dx*(-HydroHead_c(partiallygroundedinds(thisind)-1)/hydroheaddiffcenter);
            end
            % Evaluate grounded fraction:
            if IsGrounded_c(partiallygroundedinds(thisind)-1)==1 && IsGrounded_c(partiallygroundedinds(thisind))==0
                GroundedFraction_lr(partiallygroundedinds(thisind))=(x_gl-X_c(partiallygroundedinds(thisind)-1))/dx;
            else
                GroundedFraction_lr(partiallygroundedinds(thisind))=(X_c(partiallygroundedinds(thisind))-x_gl)/dx;
            end
        end
        % Determine last grounded index: (assumes last GL is an
        % ungrounding)
        lastgroundedind=partiallygroundedinds(end)-1;
    else
        % Grounding line is terminus:
        x_gl=domainwidth;
        lastgroundedind=xsize;
    end
    
else
    % Ensure that grounded ice bottom sits on the bed:
    IceBottom_c(IsGrounded_c)=BedElev_c(IsGrounded_c);
end

% Compute surface elevation:
SurfElev_c=IceBottom_c+Icethick_c;

% Compute water thickness:
WaterThick_c=IceBottom_c-BedElev_c;

%% Address Side BC:

% Linear extrapolation of ice thickness to the sides:
icethick_l=1.5*Icethick_c(1)-.5*Icethick_c(2);
icethick_r=1.5*Icethick_c(end)-.5*Icethick_c(end-1);

% Apply minimum ice thickness to BC:
icethick_l=max(icethick_l,minicethick);
icethick_r=max(icethick_r,minicethick);

% Compute ice bottom and surface elevation on side boundaries:
if IsGrounded_lr(1)==1
    icebottom_l=bedelev_l;
else
    icebottom_l=sealevel-(rho_i/rho_sw)*icethick_l;
end
surfelev_l=icebottom_l+icethick_l;
if IsGrounded_lr(end)==1
    icebottom_r=bedelev_r;
else
    icebottom_r=sealevel-(rho_i/rho_sw)*icethick_r;
end
surfelev_r=icebottom_r+icethick_r;

% Compute side strain rates and vertical velocities:
if strcmp(leftbctype,'fixedicethick') && firstgrid==0
    W_ud(:,1)=-(MeltRate_d(1)+ShapeFunctionW_lud*(Accum_c(1)+icethickchangerate_l-MeltRate_d(1)));
end
if strcmp(rightbctype,'fixedicethick') && firstgrid==0
    W_ud(:,end)=-(MeltRate_d(end)+ShapeFunctionW_rud*(Accum_c(end)+icethickchangerate_r-MeltRate_d(end)));
end

% Compute side velocities:
if strcmp(leftbctype,'flux') && firstgrid==0
    U_l=ShapeFunctionU_l*influx_l/icethick_l;
end
if strcmp(rightbctype,'flux') && firstgrid==0
    U_r=-ShapeFunctionU_r*influx_r/icethick_r;
end

%% Compute a Fictitious Monotonic Ice Bottom for the Floating Shelf:

% Set the monotonic bottom to the real ice bottom:
MonotonicIceBottom_c=IceBottom_c;

% Check if we're running the plume model and have a shelf:
if doplume && sum(IsGrounded_c==0)>1
    
    % Detect downward sloping ice bottom in the shelf:
    DownwardSlopingShelf_lr=IsGrounded_lr~=1&([IceBottom_c,icebottom_r]-[icebottom_l,IceBottom_c])<0;
    
    % Check if there are any cells that need to be fixed:
    if sum(DownwardSlopingShelf_lr)>0
        
        % Identify cells that need to be fixed:
        needfixinginds=find(DownwardSlopingShelf_lr);
        
        % Pre-allocate a record of whether they have been fixed:
        hasbeenfixed=false(size(needfixinginds));
        
        % Loop through from the calving front to the grounding line:
        for thisind=length(needfixinginds):-1:1
            
            % Check if this grid cell has been fixed:
            if hasbeenfixed(thisind)
                continue
            end
            
            % Identify real index of this grid edge:
            d2=needfixinginds(thisind);
            
            % Identify downstream ice bottom:
            if d2<xsize+1
                downstreamicebottomx=X_c(d2);
                downstreamicebottomz=IceBottom_c(d2);
            else
                downstreamicebottomx=domainwidth;
                downstreamicebottomz=icebottom_r;
            end
            
            % Identify next ice bottom below that point:
            nextind=find(IceBottom_c(1:d2-1)<downstreamicebottomz-minbotgrad*(downstreamicebottomx-X_c(1:d2-1)),1,'last');
            if isempty(nextind)
                nextind=lastgroundedind;
                nexticebottom=downstreamicebottomz-minbotgrad*(downstreamicebottomx-X_c(nextind));
            else
                nexticebottom=IceBottom_c(nextind);
            end
            
            % Create a linear profile to the next ice bottom:
            if d2<xsize+1
                MonotonicIceBottom_c(nextind:d2-1)=linspace(nexticebottom,downstreamicebottomz-(downstreamicebottomz-nexticebottom)/(d2-nextind+1),d2-nextind);
            else
                MonotonicIceBottom_c(nextind:d2-1)=linspace(nexticebottom,downstreamicebottomz-.5*(downstreamicebottomz-nexticebottom)/(d2-nextind+1),d2-nextind);
            end
            
            % Flag cells that have been fixed:
            hasbeenfixed(needfixinginds>=nextind+1)=1;
            
        end
        
    end
    
end

%% Build Grid:

% Interpolate geometry onto the grounding line: (hydrostatic assumption)
bedelev_gl=interp1(X_input,BedElev_input+SillThick_input,x_gl,interpstyle,BedElev_lr(end));
if lastgroundedind<xsize
    icethick_gl=(sealevel-bedelev_gl)*(rho_sw/rho_i);
else
    icethick_gl=icethick_r;
end
surfelev_gl=bedelev_gl+icethick_gl;

% Interpolate ice bottom to grid edges:
IceBottom_lr=[icebottom_l,.5*(IceBottom_c(1:end-1)+IceBottom_c(2:end)),icebottom_r];

% Interpolate ice thickness to grid edges, centered style:
Icethick_lr_center=[icethick_l,(Icethick_c(1:end-1)+Icethick_c(2:end))/2,icethick_r];

% Compute grid elevation matrices:
GridElev_c=repmat(IceBottom_c,[zsize,1])+repmat(Zhat_c,[1,xsize]).*repmat(Icethick_c,[zsize,1]);
GridElev_lr=repmat(IceBottom_lr,[zsize,1])+repmat(Zhat_c,[1,xsize+1]).*repmat(Icethick_lr_center,[zsize,1]);
GridElev_ud=repmat(IceBottom_c,[zsize+1,1])+repmat(Zhat_ud,[1,xsize]).*repmat(Icethick_c,[zsize+1,1]);
GridElev_lrud=repmat(IceBottom_lr,[zsize+1,1])+repmat(Zhat_ud,[1,xsize+1]).*repmat(Icethick_lr_center,[zsize+1,1]);

% Define DZ at centers:
DZ_c=repmat(DZhat_c,[1,xsize]).*repmat(Icethick_c,[zsize,1]);

% First guess flow direction is all forward:
if firstgrid
    FlowDir_lr=ones(1,xsize+1);
end

% Interpolate ice thickness to grid edges, upwind style:
Icethick_lr_upwind=zeros(1,xsize+1);
CandidateIcethick=[icethick_l,Icethick_c];
Icethick_lr_upwind(FlowDir_lr>=0)=CandidateIcethick(FlowDir_lr>=0);
CandidateIcethick=[Icethick_c,icethick_r];
Icethick_lr_upwind(FlowDir_lr<0)=CandidateIcethick(FlowDir_lr<0);
clear CandidateIcethick

% Define DZ at grid edges:
DZ_lr_center=repmat(DZhat_c,[1,xsize+1]).*repmat(Icethick_lr_center,[zsize,1]);
DZ_lr_upwind=repmat(DZhat_c,[1,xsize+1]).*repmat(Icethick_lr_upwind,[zsize,1]);
DZ_ud=repmat(DZhat_ud,[1,xsize]).*repmat(Icethick_c,[zsize+1,1]); % first and last rows are half-cells 
DZ_lrud=repmat(DZhat_ud,[1,xsize+1]).*repmat([icethick_l,.5*(Icethick_c(1:end-1)+Icethick_c(2:end)),icethick_r],[zsize+1,1]); % first and last rows are half-cells

% Compute static pressure:
PressureStatic_c=flipud(cumsum(rho_i*g*flipud(DZ_ud(2:end,:)),1));
PressureStatic_d=PressureStatic_c(1,:)+rho_i*g*DZ_ud(1,:);
PressureStatic_l=flipud(cumsum(rho_i*g*flipud(DZ_lrud(2:end,1)),1));
PressureStatic_r=flipud(cumsum(rho_i*g*flipud(DZ_lrud(2:end,end)),1));

% Compute grid slope angle: (positive for upward slope)
SinTheta_lr=sin(atan([2*(GridElev_c(:,1)-GridElev_lr(:,1))/dx,(GridElev_c(:,2:end)-GridElev_c(:,1:end-1))/dx,2*(GridElev_lr(:,end)-GridElev_c(:,end))/dx]));

% Compute driving stress gradient: (points downhill)
DrivingStressGrad_lr=[zeros(zsize,1),-rho_i*g*SinTheta_lr(:,2:end-1)-(PressureStatic_c(:,2:end)-PressureStatic_c(:,1:end-1))/dx,zeros(zsize,1)];

% Compute side stress:  
if strcmp(rightbctype,'front') && firstgrid==0
    if strcmp(velocitytype,'SSA')
        % Constant normal stress:
        NormalStress_r=-.5*rho_i*g*(1-rho_i/rho_sw)*icethick_r*ones(zsize,1);
    else
        % Vertically variable normal stress:
        PressureWater_r=max(0,sealevel-GridElev_lr(:,end))*rho_sw*g;
        NormalStress_r=PressureWater_r-PressureStatic_r;
    end
end

% Define melting point: (using ice pressure, not water pressure)
MeltPoint_c=tmelt-meltingpointslope*PressureStatic_c;
MeltPoint_d=tmelt-meltingpointslope*PressureStatic_d;
MeltPoint_l=tmelt-meltingpointslope*PressureStatic_l;
MeltPoint_r=tmelt-meltingpointslope*PressureStatic_r;

%%

% Here are some commands to visualize the cubic grounding line interpolant:
% Testx=linspace(0,1,1000);
% TestHydroHead=thesecoeffs(1)*Testx.^3+thesecoeffs(2)*Testx.^2+thesecoeffs(3)*Testx+thesecoeffs(4);
% plot(X_c/1000,HydroHead_c,'k')
% hold on
% plot((X_c(lastgroundedind)+dx*Testx)/1000,TestHydroHead,'r')
% plot([0,domainwidth]/1000,[0,0],'--k')
% xlim([X_c(lastgroundedind)-dx,X_c(lastgroundedind)+2*dx]/1000)

