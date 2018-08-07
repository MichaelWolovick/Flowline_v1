% MakeSchoofIC_2D_v1

% Mike Wolovick, 3/24/2016

% This script builds the bed profile used by Schoof (2007).  It computes
% the steady-state grounding line (you can choose which one to use, when
% multiple solutions are possible), then computes the steady-state surface
% profile.  It uses the simpler form of the asymptotic approximation (B) to
% choose the grounding line.  You even have the option of starting with the
% unstable profile that has a grounding line in the overdeepening.

% The floating shelf profile (not included in Schoof's equations) uses an
% arbitrary power-law profile.

% The grounded ice surface profile begins at the last grounded grid cell
% (just before the grounding line), and it has a thickness slightly more
% than the flotation thickness there.  The integration inland of the
% grounding line starts with a subgrid interpolated grounding line just
% beyond this last grounded cell.

% The upstream integration is implicit.

% Ice rheology is isothermal, specified through a temperature rather than a
% value of A.  The material parameters and sliding parameters need not be
% the same as in Schoof's paper (the same methodology can still compute a
% grounding line position and surface profile).  The sliding law is
% specified the same way I specify it in my model, then converted to
% Schoof's parameters when implementing the equations.

% Schoof's sliding law:  tau=C*u^m
% My sliding law:        u=u_0*(tau/tau_0)^m

% Conversion:  m_schoof=1/m_me 
%              C=tau_0*u_0^-1/m_me
%              tau_0=arbitrary

% The script also makes a time-dependent accumulation rate.  The accum
% moves in steps up and down.  

% The accumulation steps have been manually modified to have double the
% duration when the grounding line crosses the overdeepening.
              

clear all
tic

%% Parameters:

% File name:
outputfile='/home/mjw/Documents/FjordSillPaper/ModelInput/SchoofTestInput_v3.mat';

% Bed profile parameters:
bedconst1=729;                    % m
bedconst2=-2184.8;                % m
bedconst3=1031.72;                % m
bedconst4=-151.72;                % m
bedx0=7.5e5;                      % m

% Other domain-building parameters:
domainwidth=1.75e6;               % m
icethick_r=300;                   % m
shelfthickpower=10;               % unitless
numpoints=1.5e4;                  % integer
numiterations=5;                  % integer
glinterpbuffer=10;                % integer

% Which grounding line are we using:
gltouse=3;                        % 1, 2, or 3 

% Physical parameters:
rho_i=917;                        % kg/m^3
rho_w=1028;                       % kg/m^3
g=9.8;                            % m/s^2
n=3;                              % unitless
a0=4.9e-25;                       % Pa^-n*s^-1
t0=263;                           % K
q_big=1.39e5;                     % J/mol
q_small=6e4;                      % J/mol
T=248;                            % K
m=3;                              % unitless
tau0=1e5;                         % Pa 
u0_yr=300;                        % m/yr
accum_yr=.2;                      % m/yr

% Other parameters:
r=8.314;                          % J/(mol*K) (ideal gas constant)
secondsperyear=60*60*24*365.25;   % s/yr

% Forcing parameters: 
amin_yr=.2;                       % m/yr
amax_yr=1;                        % m/yr
runtime_yr=1.1e5;                 % yr
numtimesamples=1e4;               % integer
% astep_yr=.2;                      % m/yr 
% steplength_yr=1e4;                % yr
% veryshorttime_yr=1e-6;            % yr

%% Work:

% Convert from years to seconds:
u0=u0_yr/secondsperyear;
accum=accum_yr/secondsperyear;

% Convert between the way I specify a sliding law and the way Schoof does:
m_schoof=1/m;
C_schoof=tau0*(u0^(-1/m));

% Make x-coordinate:
X_input=linspace(0,domainwidth,numpoints);
dx=X_input(2)-X_input(1);

% Compute bed profile:
BedElev_input=bedconst1+bedconst2*(X_input/bedx0).^2+bedconst3*(X_input/bedx0).^4+bedconst4*(X_input/bedx0).^6;

% COmpute bed gradient:
BedGradient=gradient(BedElev_input,dx);

% Compute balance flux:
BalanceFlux=accum*X_input;

% Compute rheological constant:
if T>=t0
    A=a0*exp(-(q_big/r)*((1/T)-(1/t0)));
else
    A=a0*exp(-(q_small/r)*((1/T)-(1/t0)));
end

% Compute grounding line thickness:
GroundingLineThick=-(rho_w/rho_i)*BedElev_input;
GroundingLineThick(GroundingLineThick<0)=0;

% Compute grounding line flux:
GroundingLineFlux=(((A*((rho_i*g)^(n+1))*((1-rho_i/rho_w)^n))/((4^n)*C_schoof))^(1/(m_schoof+1)))*(GroundingLineThick.^((m_schoof+n+3)/(m_schoof+1)));

% Compute grounding line position:
% locate flux intersections:
glinds=find(diff(sign(BalanceFlux-GroundingLineFlux)));
% discard first intersection (at first grid cell)
glinds=glinds(2:end);
% Check how many intersections we have, and which one we're using:
if length(glinds)>1
    glind=glinds(gltouse);
else
    glind=glinds;
end
% locate grounding line:
x_gl=interp1(BalanceFlux(glind-glinterpbuffer:glind+glinterpbuffer)-GroundingLineFlux(glind-glinterpbuffer:glind+glinterpbuffer),...
    X_input(glind-glinterpbuffer:glind+glinterpbuffer),0,'spline');

% Find ice thickness at grounding line:
bedelev_gl=interp1(X_input,BedElev_input,x_gl,'spline');
icethick_gl=-(rho_w/rho_i)*bedelev_gl;

% Pre-allocate ice thickness:
Icethick_input=zeros(1,numpoints);

% Make shelf ice thickness:
Icethick_input(glind+1:end)=icethick_r+(icethick_gl-icethick_r)*(((domainwidth-X_input(glind+1:end))/(domainwidth-x_gl)).^shelfthickpower);

% Assign first guess of last grounded cell thickness:
Icethick_input(glind)=GroundingLineThick(glind);

% Integrate for grounded ice thickness:
for ii=glind:-1:1
    % Check if we're in the first grid cell:
    if ii==glind
        % Iterate for this ice thickness:
        for jj=1:numiterations
            % Compute ice thickness gradient:
            icethickgrad=-BedGradient(ii)-(C_schoof/(rho_i*g))*(BalanceFlux(ii)^m_schoof)/(Icethick_input(ii)^(m_schoof+1));
            % Advance ice thickness:
            Icethick_input(ii)=icethick_gl-(x_gl-X_input(ii))*icethickgrad;
        end
    else
        % Iterate for this ice thickness:
        for jj=1:numiterations
            % Compute ice thickness gradient:
            icethickgrad=-BedGradient(ii)-(C_schoof/(rho_i*g))*(BalanceFlux(ii)^m_schoof)/(Icethick_input(ii)^(m_schoof+1));
            % Advance ice thickness:
            Icethick_input(ii)=Icethick_input(ii+1)-dx*icethickgrad;
        end
    end
    % Assign first guess of next ice thickness (linear extrapolation):
    if ii>1
        Icethick_input(ii-1)=Icethick_input(ii)-dx*icethickgrad;
    end
end

% Compute surface elevation:
SurfElev=BedElev_input+Icethick_input;
SurfElev(SurfElev<(1-rho_i/rho_w)*Icethick_input)=(1-rho_i/rho_w)*Icethick_input(SurfElev<(1-rho_i/rho_w)*Icethick_input);

% % Make time-dependent forcing:
% % Compute number of steps and number of time points:
% numsteps=2*floor((amax_yr-amin_yr)/astep_yr);
% nt=2*numsteps;
% % Create time vector:
% Time_yr_input=reshape([steplength_yr*linspace(0,numsteps+1,numsteps+2)-veryshorttime_yr;steplength_yr*linspace(0,numsteps+1,numsteps+2)+veryshorttime_yr],[],1);
% % Create accumulation vector:
% Accum_yr_input=zeros(size(Time_yr_input));
% Accum_yr_input(1:nt/2+1)=amin_yr+astep_yr*floor(Time_yr_input(1:nt/2+1)/steplength_yr);
% Accum_yr_input(nt/2+2:end)=amin_yr+astep_yr*floor((Time_yr_input(end)-veryshorttime_yr-Time_yr_input(nt/2+2:end))/steplength_yr);
% Accum_yr_input(1)=amin_yr;
% Accum_yr_input(end)=amin_yr;
% 
% % Lengthen the duration of the steps that cross the overdeepening:
% Time_yr_input(11:end)=Time_yr_input(11:end)+steplength_yr;
% Time_yr_input(19:end)=Time_yr_input(19:end)+steplength_yr;

% Make time-dependent forcing:
Time_yr_input=linspace(0,runtime_yr,numtimesamples)';
Accum_yr_input=amin_yr+(amax_yr-amin_yr)*(.5-.5*cos(2*pi*Time_yr_input/runtime_yr));

% Save output:
save(outputfile)

% Final display:
disp('Done!')
toc

% % Plot geometry:
% figure(1)
% %hold off
% plot(X_input/1000,SurfElev,'k')
% hold on
% plot(X_input/1000,SurfElev-Icethick_input,'k')
% plot(X_input/1000,BedElev_input,'k')
% xlim([0,domainwidth/1000])
% ylim([-1000,4000])
% 
% % Plot fluxes:
% figure(2)
% %hold off
% plot(X_input/1000,BalanceFlux*secondsperyear,'k')
% hold on
% plot(X_input/1000,GroundingLineFlux*secondsperyear,'k')
% xlim([0,domainwidth/1000])
% ylim([0,1]*3e6)