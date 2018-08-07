% SpawnTracers_v1

% Mike Wolovick, 6/17/2012

% This script contains the commands for spawning new tracers in FTM_2D_v2.

% Key to reindexing connectivity after tracer lists are reindexed by R_t:
% [~,Rinv_t]=sort(R_t);
% Connectivity_t=Rinv_t(Connectivity_t);

% Rules for adding tracers:
% 1.  New tracers are added to the top of the list, old tracers are removed
%     from the bottom of the list.  
% 2.  Connectivity=Connectivity+numnewtracers

% Check if this is the IC:
if exist('Zhat_t_l','var')==0  % Spawn initial tracers
    
    % Create left edge tracer spawn points:
    Zhat_t_l=exp(linspace(0,log(efzheight),tracerzsize)');
    Age_t_l_yr=(icethick_l/(Accum_c(1)*secondsperyear))*log(1./Zhat_t_l);
    dzhat_u_l=Zhat_t_l(1)-Zhat_t_l(2);
    
    % Compute spawn intervals:
    spawninterval_u=dzhat_u_l*icethick_l/Accum_c(1);
    spawninterval_l=(domainwidth_start/tracerxsize)/sum(U_l.*DZhat_c);
    spawninterval_l=spawninterval_u/ceil(spawninterval_u/spawninterval_l);
    TracerParameters.spawninterval_u_yr=spawninterval_u/secondsperyear;
    TracerParameters.spawninterval_l_yr=spawninterval_l/secondsperyear;
    nextspawntime_u=spawninterval_u;
    nextspawntime_l=spawninterval_l;
    
    % Create initial tracer distribution:
    Inds_t=linspace(1,numtracers,numtracers)';
    X_t_start=[reshape(repmat(linspace(0,domainwidth_start,tracerxsize),[tracerzsize+1,1]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
    Zhat_t_start=[reshape(repmat([Zhat_t_l;0],[1,tracerxsize]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
    Age_t_yr=[reshape(repmat([Age_t_l_yr;0],[1,tracerxsize]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
    IsOriginal_t=[true(numlivetracers,1);false((storagebuffer-1)*numlivetracers,1)];
    IsAcc_t=false(numtracers,1);
    IsAcc_t(Zhat_t_start==0)=1;
    InDomain_t=[true(numlivetracers,1);false((storagebuffer-1)*numlivetracers,1)];
    IsLooseEnd_t=X_t_start==0&Zhat_t_start~=0;
    
    % Compute tracer connectivity and sort by age:
    Connectivity_t=[[linspace(1,numliveconnectors,numliveconnectors)',linspace(tracerzsize+2,numlivetracers,numliveconnectors)'];NaN*ones((storagebuffer-1)*numliveconnectors,2)]; 
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
    
elseif time>=nextspawntime_l     % Spawn new tracers
    
    % Erase tracers that have left the domain:
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
        X_t_acc=linspace(0,domainwidth_mid,tracerxsize)';
        X_t_new=[X_t_acc;zeros(tracerzsize-1,1);X_t_acc];
        Zhat_t_new=[ones(tracerxsize,1);Zhat_t_l(2:end);zeros(tracerxsize,1)];
        Age_t_new_yr=[zeros(tracerxsize,1);Age_t_l_yr(2:end);zeros(tracerxsize,1)];
        numnewtracers=length(X_t_new);
        numnewconnectors=tracerxsize-1+tracerzsize-1+tracerxsize-1;
        
        % Compute logicals of new tracers:
        InDomain_t_new=true(numnewtracers,1);
        IsOriginal_t_new=false(numnewtracers,1);
        IsAcc_t_new=[false(tracerxsize+tracerzsize-1,1);true(tracerxsize,1)];
        IsLooseEnd_t_new=[true(1,1);false(tracerxsize-1,1);true(tracerzsize-1,1);false(tracerxsize,1)];
        
        % Compute connectivity of new tracers:
        Connectivity_t_new=[linspace(1,tracerxsize-1,tracerxsize-1)',linspace(2,tracerxsize,tracerxsize-1)';...
            linspace(tracerxsize+1,tracerxsize+tracerzsize-1,tracerzsize-1)',Inds_t(IsLooseEnd_t)+numnewtracers;
            linspace(tracerxsize+tracerzsize,numnewtracers-1,tracerxsize-1)',linspace(tracerxsize+tracerzsize+1,numnewtracers,tracerxsize-1)'];
        
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