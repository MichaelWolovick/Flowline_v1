% TracerTester_2D_v1

% Mike Wolovick, 6/11/2012

% This script tests the tracer tracking routines for FTM_2D.

% It tests them inside a model domain that is a parallel-sided slab where
% horizontal velocity is constant with depth and strain rate is in
% equilibrium with the surface accumulation.

clear all
tic

%% Parameters:

% Output file:
outputfile='/home/mike/Documents/Research/FullTempModel/FTM_2D/TracerTest_v2.mat';

% Test Parameters:
domainwidth=1e5;                  % m
icethick=2.5e3;                   % m
u_l_yr=5;                         % m/yr
accum_yr=.05;                     % m/yr
runtime_yr=1e5;                   % yr
amp=1;                            % unitless
stdev=1.5e4;                      % m
n=3;                              % unitless
TestParameters=struct('domainwidth',domainwidth,'icethick',icethick,'u_l_yr',u_l_yr,'accum_yr',accum_yr,'runtime_yr',runtime_yr,'amp',amp,'stdev',stdev,'n',n);

% Numeric Parameters:
xsize=100;                        % integer
zsize=20;                         % integer
zpower=2;                         % unitless
dt_major_yr=25;                   % yr
numiterations=3;                  % integer 
recordinterval_yr=50;             % yr
NumericParameters=struct('xsize',xsize,'zsize',zsize,'zpower',zpower,'dt_major_yr',dt_major_yr,'numiterations',numiterations,'recordinterval_yr',recordinterval_yr);

% Tracer Parameters:
efzheight=.2;                     % unitless (0,1)
tracerzsize=30;                   % integer
tracerxsize=100;                  % integer
interpstyle_t='linear';           % valid interp2 method
storagebuffer=2;                  % integer
TracerParameters=struct('efzheight',efzheight,'tracerzsize',tracerzsize,'tracerxsize',tracerxsize,'spawninterval_u_yr',[],'spawninterval_l_yr',[],'interpstyle_t',interpstyle_t,'storagebuffer',storagebuffer);

% Other Parameters:
secondsperyear=60*60*24*365.25;   % s/yr
displayinterval=250;              % integer
OtherParameters=struct('secondsperyear',secondsperyear,'displayinterval',displayinterval);

%% Preparation:

% Define dependent parameters:
dt_major=dt_major_yr*secondsperyear;
recordinterval=recordinterval_yr*secondsperyear;
accum=accum_yr/secondsperyear;
runtime=runtime_yr*secondsperyear;
u_l=u_l_yr/secondsperyear;

% Compute melt rate profile:
dx=domainwidth/xsize;
X_c=domainwidth*linspace(.5/xsize,1-.5/xsize,xsize);
X_lr=domainwidth*linspace(0,1,xsize+1);
MeltRate_c=-accum*amp*exp(-.5*(((X_c-.5*domainwidth)/stdev).^2));

% Define dimensionless elevation:
Zhat_ud=linspace(0,1,zsize+1)'.^zpower;
Zhat_c=.5*(Zhat_ud(1:end-1)+Zhat_ud(2:end));
DZhat_c=Zhat_ud(2:end)-Zhat_ud(1:end-1);
DZhat_ud=[Zhat_c(1);Zhat_c(2:end)-Zhat_c(1:end-1);1-Zhat_c(end)];  

% Define shape function:
ShapeFunctionU_c=cumsum(DZhat_ud(1:end-1).*((1-Zhat_ud(1:end-1)).^n));
ShapeFunctionU_c=ShapeFunctionU_c/sum(ShapeFunctionU_c.*DZhat_c);
ShapeFunctionW_ud=[0;cumsum(ShapeFunctionU_c.*DZhat_c)];

% Define horizontal velocity field:
Ubar_lr=u_l+[0,cumsum((dx/icethick)*(accum-MeltRate_c))];
U_lr=repmat(ShapeFunctionU_c,[1,xsize+1]).*repmat(Ubar_lr,[zsize,1]);

% Define vertical velocity field:
W_ud=-accum*repmat(ShapeFunctionW_ud,[1,xsize])+repmat(-MeltRate_c,[zsize+1,1]).*repmat(1-ShapeFunctionW_ud,[1,xsize]);
What_ud=W_ud/icethick;

% Create left edge tracer spawn points:
Zhat_t_l=exp(linspace(0,log(efzheight),tracerzsize)');
Age_t_l_yr=(icethick/accum_yr)*log(1./Zhat_t_l); 
dzhat_u_l=Zhat_t_l(1)-Zhat_t_l(2);

% Compute spawn intervals:
spawninterval_u=dzhat_u_l*icethick/accum;
spawninterval_l=(domainwidth/tracerxsize)/u_l;
spawninterval_l=spawninterval_u/ceil(spawninterval_u/spawninterval_l);
TracerParameters.spawninterval_u_yr=spawninterval_u/secondsperyear;
TracerParameters.spawninterval_l_yr=spawninterval_l/secondsperyear;
nextspawntime_u=spawninterval_u;
nextspawntime_l=spawninterval_l;

% Create initial tracer distribution:
numlivetracers=tracerxsize*tracerzsize;
numliveconnectors=(tracerxsize-1)*tracerzsize;
numtracers=storagebuffer*numlivetracers;
numconnectors=storagebuffer*numliveconnectors;
Inds_t=linspace(1,numtracers,numtracers)';
X_t_start=[reshape(repmat(linspace(0,domainwidth,tracerxsize),[tracerzsize,1]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
Zhat_t_start=[reshape(repmat(Zhat_t_l,[1,tracerxsize]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
Age_t_yr=[reshape(repmat(Age_t_l_yr,[1,tracerxsize]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
IsOriginal_t=[true(numlivetracers,1);false((storagebuffer-1)*numlivetracers,1)];
IsAcc_t=false(numtracers,1);
InDomain_t=[true(numlivetracers,1);false((storagebuffer-1)*numlivetracers,1)];
IsLooseEnd_t=X_t_start==0;

% Compute tracer connectivity and sort by age:
Connectivity_t=[[linspace(1,numliveconnectors,numliveconnectors)',linspace(tracerzsize+1,numlivetracers,numliveconnectors)'];NaN*ones((storagebuffer-1)*numliveconnectors,2)];
ConnectorAge_t_yr=[Age_t_yr(Connectivity_t(1:numliveconnectors,1));NaN*ones((storagebuffer-1)*numliveconnectors,1)];
[ConnectorAge_t_yr,ConnectorAgeOrder_t]=sort(ConnectorAge_t_yr);
Connectivity_t=Connectivity_t(ConnectorAgeOrder_t,:);
DeadConnector_t=[false(numliveconnectors,1);true((storagebuffer-1)*numliveconnectors,1)];
clear ConnectivityAgeOrder_t

% Sort tracers by age:  (R_t=reindexing, Rinv_t=inverse reindexing)
[~,R_t]=sort(Age_t_yr);
X_t_start=X_t_start(R_t);
Zhat_t_start=Zhat_t_start(R_t);
Age_t_yr=Age_t_yr(R_t);
IsOriginal_t=IsOriginal_t(R_t);
IsAcc_t=IsAcc_t(R_t);
InDomain_t=InDomain_t(R_t);
IsLooseEnd_t=IsLooseEnd_t(R_t);
[~,Rinv_t]=sort(R_t);
Connectivity_t(DeadConnector_t==0,:)=Rinv_t(Connectivity_t(DeadConnector_t==0,:));

% Pre-allocate memory:
U_t=zeros(numtracers,1);
What_t=zeros(numtracers,1);

% Key to reindexing connectivity after tracer lists are reindexed by R_t:
% [~,Rinv_t]=sort(R_t);
% Connectivity_t=Rinv_t(Connectivity_t);

% Rules for adding tracers:
% 1.  New tracers are added to the top of the list, old tracers are removed
%     from the bottom of the list.  
% 2.  Connectivity=Connectivty+numnewtracers

% Initialize record file:
numrecords=floor(runtime/recordinterval)+1;
numdigits=floor(log10(numrecords));
prefix='0'*ones(1,numdigits-1);
idnumber='0';
eval(['time_yr_',prefix,idnumber,'=0;'])
eval(['X_t_',prefix,idnumber,'=X_t_start;'])
eval(['Zhat_t_',prefix,idnumber,'=Zhat_t_start;'])
eval(['Age_t_yr_',prefix,idnumber,'=Age_t_yr;'])
eval(['IsOriginal_t_',prefix,idnumber,'=IsOriginal_t;'])
eval(['IsAcc_t_',prefix,idnumber,'=IsAcc_t;'])
eval(['InDomain_t_',prefix,idnumber,'=InDomain_t;'])
eval(['IsLooseEnd_t_',prefix,idnumber,'=IsLooseEnd_t;'])
eval(['Connectivity_t_',prefix,idnumber,'=Connectivity_t;'])
eval(['DeadConnector_t_',prefix,idnumber,'=DeadConnector_t;'])
wildcard=['*',prefix,idnumber];
save(outputfile,'*Parameters','X_c','X_lr','Zhat_c','Zhat_ud','U_lr','W_ud','DZhat*','numrecords',wildcard)
clear(wildcard)

%% Run Model:

% Initial variables:
disp('...')
disp(['dt_major=',num2str(dt_major_yr),'yr'])
disp(['record interval=',num2str(recordinterval_yr),'yr'])
disp('...')
disp('Running Model')
done=0;
timestep=1;
time=0;
time_yr=0;
nextrecordtime=recordinterval;
recordcounter=1;
displaycounter=1;

% Loop through timesteps:
while done==0
    
    % Assign begining values as first guesses of the middle and end (ie, euler):
    X_t_mid=X_t_start;
    Zhat_t_mid=Zhat_t_start;
    
    % Perform the iterative method:
    for iteration_major=1:numiterations
        
        % Interpolate velocity to tracer locations:
        U_t(InDomain_t)=interp2(X_lr,[0;Zhat_c;1],[zeros(1,xsize+1);U_lr;U_lr(end,:)],X_t_mid(InDomain_t),Zhat_t_mid(InDomain_t),interpstyle_t);
        What_t(InDomain_t)=interp2([0,X_c,domainwidth],Zhat_ud,[What_ud(:,1),What_ud,What_ud(:,end)],X_t_mid(InDomain_t),Zhat_t_mid(InDomain_t),interpstyle_t);
        
        % Advance tracer positions to middle of timestep:
        X_t_mid(InDomain_t)=X_t_start(InDomain_t)+.5*dt_major*U_t(InDomain_t);
        Zhat_t_mid(InDomain_t)=Zhat_t_start(InDomain_t)+.5*dt_major*What_t(InDomain_t);
        
        % Flag tracers which are still in the domain:
        InDomain_t=X_t_mid>=0&X_t_mid<=domainwidth&Zhat_t_mid>=0&Zhat_t_mid<=1;
        
    end
    
    % Advance nonphysics variables to end of major timestep:
    X_t_end=X_t_start+dt_major*U_t;
    Zhat_t_end=Zhat_t_start+dt_major*What_t;
    
    % Advance left edge tracer spawn points:
    What_t_l=interp1(Zhat_ud,What_ud(:,1),Zhat_t_l,interpstyle_t);
    Zhat_t_l=Zhat_t_l+dt_major*What_t_l;
    Age_t_l_yr=Age_t_l_yr+dt_major_yr;
    
    % Advance end variables to start variables:
    X_t_start=X_t_end;
    Zhat_t_start=Zhat_t_end;
    
    % Advance time:
    time=time+dt_major;
    time_yr=time/secondsperyear;
    Age_t_yr=Age_t_yr+dt_major_yr;
    
    % Spawn new tracers:
    if time>=nextspawntime_l
        
        % Erase tracers which have left the domain:
        numindomain=sum(InDomain_t);
        X_t_start(InDomain_t==0)=NaN;
        Zhat_t_start(InDomain_t==0)=NaN;
        Age_t_yr(InDomain_t==0)=NaN;
        
        % Consolidate young in-domain tracers at the top of the list:
        [~,R_t]=sort(Age_t_yr);
        X_t_start=X_t_start(R_t);
        Zhat_t_start=Zhat_t_start(R_t);
        Age_t_yr=Age_t_yr(R_t);
        InDomain_t=InDomain_t(R_t);
        IsOriginal_t=IsOriginal_t(R_t);
        IsAcc_t=IsAcc_t(R_t);
        IsLooseEnd_t=IsLooseEnd_t(R_t);
        [~,Rinv_t]=sort(R_t);
        Connectivity_t(DeadConnector_t==0,:)=Rinv_t(Connectivity_t(DeadConnector_t==0,:));
        
        % Check if new tracers are on left edge only:
        if time>=nextspawntime_u % all edges
            
            % Kill an old spawn point, add a new spawn point:
            Zhat_t_l=[1;Zhat_t_l(1:end-1)];
            Age_t_l_yr=[0;Age_t_l_yr(1:end-1)];
            IsLooseEnd_t(Age_t_yr==max(Age_t_yr(IsLooseEnd_t)))=0;
            
            % Compute locations of new tracers:  [surf;leftBC;acc]
            X_t_acc=linspace(0,domainwidth,tracerxsize)';
            IsLive_t_acc=interp1([0,X_c,domainwidth],[MeltRate_c(1),MeltRate_c,MeltRate_c(end)],X_t_acc,interpstyle_t)<0;
            X_t_new=[X_t_acc;zeros(tracerzsize-1,1);X_t_acc(IsLive_t_acc)];
            Zhat_t_new=[ones(tracerxsize,1);Zhat_t_l(2:end);zeros(sum(IsLive_t_acc),1)];
            Age_t_new_yr=[zeros(tracerxsize,1);Age_t_l_yr(2:end);zeros(sum(IsLive_t_acc),1)];
            numnewtracers=length(X_t_new);
            numnewconnectors=tracerxsize-1+tracerzsize-1+sum(IsLive_t_acc)-1;
            
            % Compute logicals of new tracers:
            InDomain_t_new=true(numnewtracers,1);
            IsOriginal_t_new=false(numnewtracers,1);
            IsAcc_t_new=[false(tracerxsize+tracerzsize-1,1);true(sum(IsLive_t_acc),1)];
            IsLooseEnd_t_new=[true(1,1);false(tracerxsize-1,1);true(tracerzsize-1,1);false(sum(IsLive_t_acc),1)];
            
            % Compute connectivity of new tracers:
            Connectivity_t_new=[linspace(1,tracerxsize-1,tracerxsize-1)',linspace(2,tracerxsize,tracerxsize-1)';...
                linspace(tracerxsize+1,tracerxsize+tracerzsize-1,tracerzsize-1)',Inds_t(IsLooseEnd_t)+numnewtracers;
                linspace(tracerxsize+tracerzsize,numnewtracers-1,sum(IsLive_t_acc)-1)',linspace(tracerxsize+tracerzsize+1,numnewtracers,sum(IsLive_t_acc)-1)'];
            
            % Update spawn timer:
            nextspawntime_u=nextspawntime_u+spawninterval_u;
            nextspawntime_l=nextspawntime_l+spawninterval_l;
            
        else % left edge only
            
            % Compute locations of new tracers:
            X_t_new=zeros(tracerzsize,1);
            Zhat_t_new=Zhat_t_l;
            Age_t_new_yr=Age_t_l_yr;
            numnewtracers=length(X_t_new);
            numnewconnectors=tracerzsize;
            
            % Compute logicals of new tracers:
            InDomain_t_new=true(numnewtracers,1);
            IsOriginal_t_new=false(numnewtracers,1);
            IsAcc_t_new=false(numnewtracers,1);
            IsLooseEnd_t_new=true(numnewtracers,1);
            
            % Compute connectivity of new tracers:
            Connectivity_t_new=[linspace(1,tracerzsize,tracerzsize)',Inds_t(IsLooseEnd_t)+numnewtracers];
            
            % Update spawn timer:
            nextspawntime_l=nextspawntime_l+spawninterval_l;
            
        end
        
        % Add new tracers to the list:
        X_t_start=[X_t_new;X_t_start(1:end-numnewtracers)];
        Zhat_t_start=[Zhat_t_new;Zhat_t_start(1:end-numnewtracers)];
        Age_t_yr=[Age_t_new_yr;Age_t_yr(1:end-numnewtracers)];
        InDomain_t=[InDomain_t_new;InDomain_t(1:end-numnewtracers)];
        IsOriginal_t=[IsOriginal_t_new;IsOriginal_t(1:end-numnewtracers)];
        IsAcc_t=[IsAcc_t_new;IsAcc_t(1:end-numnewtracers)];
        IsLooseEnd_t=[IsLooseEnd_t_new;false(numtracers-numnewtracers,1)];
        
        % Flag dead connectors:
        DeadConnector_t(DeadConnector_t==0)=isnan(Age_t_yr(Connectivity_t(DeadConnector_t==0,1)))|isnan(Age_t_yr(Connectivity_t(DeadConnector_t==0,2)))|max(Connectivity_t(DeadConnector_t==0,:),[],2)+numnewtracers>numtracers;
        ConnectorAge_t_yr(DeadConnector_t)=NaN;
        numdeadconnectors=sum(DeadConnector_t);
        
        % Sort connectors by age:
        [ConnectorAge_t_yr,ConnectorAgeOrder_t]=sort(ConnectorAge_t_yr);
        Connectivity_t=Connectivity_t(ConnectorAgeOrder_t,:);
        DeadConnector_t=DeadConnector_t(ConnectorAgeOrder_t);
        clear ConnectivityAgeOrder_t
        
        % Update connectivity:
        Connectivity_t=[Connectivity_t_new;Connectivity_t(1:end-numnewconnectors,:)+numnewtracers];
        ConnectorAge_t_yr=[zeros(numnewconnectors,1);ConnectorAge_t_yr(1:end-numnewconnectors)];
        DeadConnector_t=[max(Connectivity_t_new,[],2)>numtracers;DeadConnector_t(1:end-numnewconnectors)];
        
        % Check connectivity:
        if max(max(Connectivity_t(DeadConnector_t==0,:)))>numtracers
            error('Live connector linked to dead tracer.')
        end
        if max(max(abs(Connectivity_t-round(Connectivity_t))))>1e-6
            error('Connector index not an integer.')
        end
        
        % Compute number of live tracers and connectors:
        numlivetracers=min([numindomain+numnewtracers,numtracers]);
        numliveconnectors=sum(DeadConnector_t==0);
        
        % Clear extraneous variables:
        clear *t_new*
        
    end
    
    % Display elapsed time on command line:
    if rem(displaycounter,displayinterval)==1
        disp('...')
        disp(strcat('Time=',num2str(time_yr),'yr'))
    end
    displaycounter=displaycounter+1;
    
    % Record tracer locations:
    if time>=nextrecordtime
        prefix='0'*ones(1,numdigits-floor(log10(recordcounter))-1);
        idnumber=num2str(recordcounter);
        eval(['time_yr_',prefix,idnumber,'=time_yr;'])
        eval(['X_t_',prefix,idnumber,'=X_t_start;'])
        eval(['Zhat_t_',prefix,idnumber,'=Zhat_t_start;'])
        eval(['Age_t_yr_',prefix,idnumber,'=Age_t_yr;'])
        eval(['IsOriginal_t_',prefix,idnumber,'=IsOriginal_t;'])
        eval(['IsAcc_t_',prefix,idnumber,'=IsAcc_t;'])
        eval(['InDomain_t_',prefix,idnumber,'=InDomain_t;'])
        eval(['IsLooseEnd_t_',prefix,idnumber,'=IsLooseEnd_t;'])
        eval(['Connectivity_t_',prefix,idnumber,'=Connectivity_t;'])
        eval(['DeadConnector_t_',prefix,idnumber,'=DeadConnector_t;'])
        wildcard=['*',prefix,idnumber];
        save(outputfile,wildcard,'-append')
        clear(wildcard)
        recordcounter=recordcounter+1;
        nextrecordtime=nextrecordtime+recordinterval;
    end
    
    % Break from loop:
    if time>=runtime
        done=1;
    end
    
end

% Final Display:
disp('Done!')
toc