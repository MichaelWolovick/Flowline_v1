% ProcessResolutionTests_v1

% Mike Wolovick, 3/29/2016

% This script loads the model output from the Schoof resolution tests,
% processes their output, and saves the results in a single summary file.

% Time series:
% Time [numrecordsx1]
% Grounding line position [numrecordsxnumruns]
% Grounding line velocity [numrecordsxnumruns]
% Volume above flotation [numrecordsxnumruns]
% VAF change rate [numrecordsxnumruns]

% Steady state datasets:
% Accumulation rate [numstepsx1]
% Analytic grounding line position [numstepsx1]
% Model grounding line position [numstepsxnumruns]
% Volume above flotation [numstepsxnumruns]

% Dynamic datasets:
% Advance rate [1xnumruns]
% VAF growth rate [1xnumruns]
% Retreat rate [1xnumruns]
% VAF loss rate [1xnumruns]

% The steady state datasets are the mean within a set time period before
% the forcing changes.

% The dynamic datasets are the average rates as the grounding line 
% traverses the overdeepening.  Only the rate through the steepest part of
% the overdeepening is considered.  

% The script assumes that all the runs have the same parameters, bed
% topography, and forcing.  The script also assumes that the forcing went
% low->high->low, and that all the runs have the same recordinterval.  The
% script assumes that sealevel=0;

clear all
close all
tic

%% Parameters:

% File names and conventions:
inputfileformat='/net/mjw/FjordSillPaper/ModelOutput/SchoofTest_v2_run';
outputfile='/home/mjw/Documents/FjordSillPaper/ModelOutput/SchoofTest_v2_processedruns.mat';

% Number of runs:
numruns=9;                        % integer
minrun=2;                         % integer

% For compiling one additional run after doing a partial compilation:
thisrun=2;                        % integer (can be commented)

% Arbitrary analysis parameters:
ssavgtime=500;                    % yr
bedslopecutoff=.5;                % unitless [0,1)

% Does the forcing have discrete steps or not?
dosteps=1;                        % logical

% Grounding line interpolation bufer:
glinterpbuffer=10;                % integer

% Timescale for low- and high-pass filtering the rates:
timescale_yr=5e2;                 % yr

%% Compute Steady State Solutions:

% Check if we're only compiling one run:
if exist('thisrun','var')==0
    
    % Communicate:
    disp('Doing preliminary stuff')
    
    % Load one set of inputs:
    load([inputfileformat,num2str(numruns),'.mat'],'*_input','*Parameters','numdigits')
    
    % Unpack needed parameters:
    unpack(BottomParameters)
    unpack(RheologyParameters)
    unpack(OtherParameters)
    unpack(TimingParameters)
    unpack(ThermalParameters)
    rho_sw=PlumeParameters.rho_sw;
    
    % Convert between the way I specify a sliding law and the way Schoof does:
    m_schoof=1/m;
    C_schoof=slidingstressscale*((slidingvelocityscale_yr/secondsperyear)^(-1/m));
    
    % Compute rheological constant:
    if consttemp>=t0
        A=a0*exp(-(q_big/r)*((1/consttemp)-(1/t0)));
    else
        A=a0*exp(-(q_small/r)*((1/consttemp)-(1/t0)));
    end
    
    % Compute grounding line thickness:
    GroundingLineThick=-(rho_sw/rho_i)*BedElev_input;
    GroundingLineThick(GroundingLineThick<0)=0;
    
    % Compute grounding line flux:
    GroundingLineFlux=(((A*((rho_i*g)^(n+1))*((1-rho_i/rho_sw)^n))/((4^n)*C_schoof))^(1/(m_schoof+1)))*(GroundingLineThick.^((m_schoof+n+3)/(m_schoof+1)));
    
    % Does the forcing have discrete steps?
    if dosteps
        
        % Identify accumulation steps:
        % Identify jumps in accum forcing:
        jumpinds=find(abs(diff(Accum_yr_input))>verysmallnumber);
        numjumps=length(jumpinds);
        jumpinds=[1;jumpinds];
        numsteps=numjumps+1;
        % Identify times of the jumps:
        Time_yr_jumps=[.5*(Time_yr_input(jumpinds)+Time_yr_input(jumpinds+1));runtime_yr];
        
        % Pre-allocate:
        Accum_yr=zeros(numsteps,1);
        SteadyGroundingLine_analytic=zeros(numsteps,1);
        
        % Loop through accumulation steps:
        for ii=1:numsteps
            
            % Identify accumulation rate:
            Accum_yr(ii)=Accum_yr_input(jumpinds(ii)+1);
            
            % Compute balance flux:
            BalanceFlux=(Accum_yr(ii)/secondsperyear)*X_input;
            
            % Locate flux intersections:
            glinds=find(diff(sign(BalanceFlux-GroundingLineFlux)));
            
            % Discard first intersection (at first grid cell)
            glinds=glinds(2:end);
            
            % Check how many intersections we have, and which one we're using:
            if length(glinds)>1 && ii>numsteps/2
                glind=glinds(3);
            else
                glind=glinds(1);
            end
            
            % Locate grounding line:
            SteadyGroundingLine_analytic(ii)=interp1(BalanceFlux(glind-glinterpbuffer:glind+glinterpbuffer)-GroundingLineFlux(glind-glinterpbuffer:glind+glinterpbuffer),...
                X_input(glind-glinterpbuffer:glind+glinterpbuffer),0,'spline');
            
        end
        
    end
    
else
    
    % Communicate:
    disp('Loading previous output file')
    
    % Load existing output:
    thisrun1=thisrun;
    minrun1=minrun;
    load(outputfile)
    thisrun=thisrun1;
    minrun=minrun1;
    clear thisrun1 minrun1
    
end

%% Load and Analyze Model Output:

% Check if we're only compiling one run:
if exist('thisrun','var')==0
    
    % Compute number of records:
    numrecords=floor(runtime_yr/recordinterval_yr);
    
    % Pre-allocate time series:
    Time_yr=.5*recordinterval_yr+recordinterval_yr*linspace(0,numrecords-1,numrecords)';
    GroundingLine=zeros(numrecords,numruns);
    GroundingLineRate=zeros(numrecords,numruns);
    GroundingLineRate_lowpass=zeros(numrecords,numruns);
    GroundingLineRate_highpass=zeros(numrecords,numruns);
    VAF=zeros(numrecords,numruns);
    VAFRate=zeros(numrecords,numruns);
    VAFRate_lowpass=zeros(numrecords,numruns);
    VAFRate_highpass=zeros(numrecords,numruns);
    
    % Pre-allocate model numeric attributes:
    DX=zeros(1,numruns);
    NX=zeros(1,numruns);
    DT=zeros(1,numruns);
    NT=zeros(1,numruns);
    ComputerTime=zeros(1,numruns);
    
    % Set startrun and stoprun:
    startrun=minrun;
    stoprun=numruns;
    
else
    
    % Set startrun and stoprun:
    startrun=thisrun;
    stoprun=thisrun;
    
end

% Communicate:
disp('Building model time series')

% Loop through model runs:
for jj=startrun:stoprun
    
    % Communicate:
    disp(['Model run=',num2str(jj),'/',num2str(numruns)])
    
    % Loop through model records:
    for ii=1:numrecords
        
        % Communicate:
        if rem(ii,200)==0
            disp(['record=',num2str(ii),'/',num2str(numrecords)])
        end
        
        % Construct name of this model record:
        prefix='0'*ones(1,numdigits(end)-floor(log10(ii))-1);
        idnumber=num2str(ii);
        ThisName=['ModelRecord_',prefix,idnumber];
        
        % Load this model record:
        if ii==1
            load([inputfileformat,num2str(jj),'.mat'],ThisName,'TimingParameters','GridParameters','computertime')
        else
            load([inputfileformat,num2str(jj),'.mat'],ThisName)
        end
        
        % Assign model attributes and bed topography:
        if ii==1
            % Assign domainwidth:
            eval(['domainwidth=',ThisName,'.domainwidth;'])
            % Assign numeric attributes:
            NX(jj)=GridParameters.xsize;
            DX(jj)=domainwidth/NX(jj);
            DT(jj)=TimingParameters.dt_yr;
            NT(jj)=ceil(runtime_yr/DT(jj));
            ComputerTime(jj)=computertime;
            % Produce horizontal grid:
            X_lr=linspace(0,domainwidth,NX(jj)+1);
            X_c=linspace(.5*DX(jj),domainwidth-.5*DX(jj),NX(jj));
            % Interpolate bed:
            BedElev_c=interp1(X_input,BedElev_input,X_c,'linear',BedElev_input(end));
        end
        
        % Assign values from this model record:
        eval(['GroundingLine(ii,jj)=',ThisName,'.x_gl;'])
        eval(['Icethick_c=',ThisName,'.Icethick_c;'])
        
        % Compute height above flotation:
        HeightAboveFlotation=Icethick_c-max(0,min(Icethick_c,-(rho_sw/rho_i)*BedElev_c));
        
        % Compute volume above flotation:
        VAF(ii,jj)=sum(DX(jj)*HeightAboveFlotation);
        
        % Clear this model record:
        eval(['clear ',ThisName])
        
    end
    
    % Compute rates:
    GroundingLineRate(:,jj)=gradient(GroundingLine(:,jj),recordinterval_yr);
    VAFRate(:,jj)=gradient(VAF(:,jj),recordinterval_yr);
    
    % Low-pass and high-pass filter the rates:
    GroundingLineRate_lowpass(:,jj)=intuitive_lowpass(GroundingLineRate(:,jj),timescale_yr/recordinterval_yr);
    VAFRate_lowpass(:,jj)=intuitive_lowpass(VAFRate(:,jj),timescale_yr/recordinterval_yr);
    GroundingLineRate_highpass(:,jj)=GroundingLineRate(:,jj)-GroundingLineRate_lowpass(:,jj);
    VAFRate_highpass(:,jj)=VAFRate(:,jj)-VAFRate_lowpass(:,jj);
    
end

% Communicate:
disp('Analyzing model time series')

% Check if we're only compiling one run:
if exist('thisrun','var')==0 && dosteps
    % Pre-allocate steady-state properties:
    SteadyGroundingLine_model=zeros(numsteps,numruns);
    SteadyVAF=zeros(numsteps,numruns);
end

% Compute steady-state properties:
if dosteps
    for ii=1:numsteps
        
        % Identify records within the steady-state averaging window:
        InWindow=Time_yr>=Time_yr_jumps(ii+1)-ssavgtime&Time_yr<=Time_yr_jumps(ii+1);
        
        % Loop through model runs:
        for jj=startrun:stoprun
            % Compute averages:
            SteadyGroundingLine_model(ii,jj)=mean(GroundingLine(InWindow,jj));
            SteadyVAF(ii,jj)=mean(VAF(InWindow,jj));
        end
        
    end
end

% Check if we're only compiling one run:
if exist('thisrun','var')==0
    
    % Pre-allocate dynamic properties:
    AdvanceRate=zeros(1,numruns);
    GrowthRate=zeros(1,numruns);
    RetreatRate=zeros(1,numruns);
    LossRate=zeros(1,numruns);
    
    % Compute bed gradient:
    BedGradient=gradient(BedElev_input);
    
    % Identify steep section of overdeepening:
    ind1=find(BedGradient>bedslopecutoff*max(BedGradient),1,'first');
    ind2=find(BedGradient>bedslopecutoff*max(BedGradient),1,'last');
    x1=X_input(ind1);
    x2=X_input(ind2);
    
end

% Loop through model runs to compute dynamic properties:
for jj=startrun:stoprun
    
    % Identify times when the grounding line was in the steep section of
    % the overdeepening:
    InOverdeepening=GroundingLine(:,jj)>=x1&GroundingLine(:,jj)<=x2;
    
    % Compute dynamic properties:
    AdvanceRate(jj)=mean(GroundingLineRate(InOverdeepening&GroundingLineRate(:,jj)>0,jj));
    GrowthRate(jj)=mean(VAFRate(InOverdeepening&VAFRate(:,jj)>0,jj));
    RetreatRate(jj)=mean(GroundingLineRate(InOverdeepening&GroundingLineRate(:,jj)<0,jj));
    LossRate(jj)=mean(VAFRate(InOverdeepening&VAFRate(:,jj)<0,jj));
    
end

% Save output:
save(outputfile)

% Final display:
disp('Done!')
toc

