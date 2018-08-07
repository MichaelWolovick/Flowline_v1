% SpawnTracers_v1

% Mike Wolovick, 6/17/2012

% This script contains the commands for spawning new tracers in FTM_2D_v2.

% v2:  the script is now able to handle creating tracers at the upper
% surface without simultaneously creating them at the left edge.  ie, it is
% now safe to run the script with zero influx.  IC tracer distribution
% includes accreted ice at an elevation of "icacczhat".  This layer of
% basal ice is also a loose end that is advected in from the left edge.
% Tracers intering as part of this layer always have zero age.  In
% addition, unique tracer ID numbers are stored along with the position and
% connectivity information.

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
    
    % Compute next spawn times:
    nextspawntime_u=spawninterval_u;
    nextspawntime_d=spawninterval_d;
    nextspawntime_l=spawninterval_l;
    
    % Create initial tracer distribution:
    % Index within the list:
    Inds_t=linspace(1,numtracers,numtracers)';  
    % Unique ID number:
    ID_t=[linspace(1,numlivetracers,numlivetracers)';NaN*ones((storagebuffer-1)*numlivetracers,1)]; 
    nextid=numlivetracers+1;
    % X position:
    X_t=[reshape(repmat(linspace(0,domainwidth_start,tracerxsize),[tracerzsize+1,1]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
    % Zhat position:
    Zhat_t=[reshape(repmat([Zhat_t_l;icacczhat],[1,tracerxsize]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
    % Age:
    Age_t_yr=[reshape(repmat([Age_t_l_yr;0],[1,tracerxsize]),[],1);NaN*ones((storagebuffer-1)*numlivetracers,1)];
    % Logicals:
    IsOriginal_t=[true(numlivetracers,1);false((storagebuffer-1)*numlivetracers,1)];
    IsAcc_t=[reshape([false(tracerzsize,tracerxsize);true(1,tracerxsize)],[],1);false((storagebuffer-1)*numlivetracers,1)];
    IsAcc_t(Zhat_t==0)=1;
    InDomain_t=[true(numlivetracers,1);false((storagebuffer-1)*numlivetracers,1)];
    IsLooseEnd_t=X_t==0;
    
    % Compute tracer connectivity and sort by age:
    Connectivity_t=[[linspace(1,numliveconnectors,numliveconnectors)',linspace(tracerzsize+2,numlivetracers,numliveconnectors)'];NaN*ones((storagebuffer-1)*numliveconnectors,2)]; 
    ConnectorAge_t_yr=[Age_t_yr(Connectivity_t(1:numliveconnectors,1));NaN*ones((storagebuffer-1)*numliveconnectors,1)];
    [ConnectorAge_t_yr,ConnectorAgeOrder_t]=sort(ConnectorAge_t_yr);
    Connectivity_t=Connectivity_t(ConnectorAgeOrder_t,:);
    DeadConnector_t=[false(numliveconnectors,1);true((storagebuffer-1)*numliveconnectors,1)];
    clear ConnectivityAgeOrder_t
    
    % Sort tracers by age:  (R_t=reindexing, Rinv_t=inverse reindexing)
    [~,R_t]=sort(Age_t_yr);
    X_t=X_t(R_t);
    Zhat_t=Zhat_t(R_t);
    Age_t_yr=Age_t_yr(R_t);
    IsOriginal_t=IsOriginal_t(R_t);
    IsAcc_t=IsAcc_t(R_t);
    InDomain_t=InDomain_t(R_t);
    IsLooseEnd_t=IsLooseEnd_t(R_t);
    [~,Rinv_t]=sort(R_t);
    Connectivity_t(DeadConnector_t==0,:)=Rinv_t(Connectivity_t(DeadConnector_t==0,:));
    
elseif time>=nextspawntime_l || time>=nextspawntime_u || time>=nextspawntime_d   % Spawn any new tracers
    
    % Erase tracers that have left the domain:
    numindomain=sum(InDomain_t);
    ID_t(InDomain_t==0)=NaN;
    X_t(InDomain_t==0)=NaN;
    Zhat_t(InDomain_t==0)=NaN;
    Age_t_yr(InDomain_t==0)=NaN;
    
    % Consolidate young in-domain tracers at the top of the list:
    [~,R_t]=sort(Age_t_yr);
    ID_t=ID_t(R_t);
    X_t=X_t(R_t);
    Zhat_t=Zhat_t(R_t);
    Age_t_yr=Age_t_yr(R_t);
    InDomain_t=InDomain_t(R_t);
    IsOriginal_t=IsOriginal_t(R_t);
    IsAcc_t=IsAcc_t(R_t);
    IsLooseEnd_t=IsLooseEnd_t(R_t);
    [~,Rinv_t]=sort(R_t);
    Connectivity_t(DeadConnector_t==0,:)=Rinv_t(Connectivity_t(DeadConnector_t==0,:));
    
    % Check where to create new tracers:
    if time>=nextspawntime_u && time>=nextspawntime_l   % Upper and left surfaces
        
        % Kill an old spawn point, add a new spawn point:
        Zhat_t_l=[1;Zhat_t_l(1:end-1)];
        Age_t_l_yr=[0;Age_t_l_yr(1:end-1)];
        IsLooseEnd_t(Age_t_yr==max(Age_t_yr(IsLooseEnd_t)))=0;
        
        % Compute locations of new tracers:  [surf;left]
        X_t_acc=linspace(0,domainwidth_mid,tracerxsize)';
        X_t_new=[X_t_acc;zeros(tracerzsize,1)];
        Zhat_t_new=[ones(tracerxsize,1);Zhat_t_l(2:end);icacczhat];
        Age_t_new_yr=[zeros(tracerxsize,1);Age_t_l_yr(2:end);0];
        numnewtracers=length(X_t_new);
        numnewconnectors=tracerxsize-1+tracerzsize;
        
        % Compute unique ID numbers of new tracers:
        ID_t_new=linspace(nextid,nextid+numnewtracers-1,numnewtracers)';
        nextid=nextid+numnewtracers;
        
        % Compute logicals of new tracers:
        InDomain_t_new=true(numnewtracers,1);
        IsOriginal_t_new=false(numnewtracers,1);
        IsAcc_t_new=[false(tracerxsize+tracerzsize-1,1);true(1)];
        IsLooseEnd_t_new=[true(1);false(tracerxsize-1,1);true(tracerzsize,1)];
        
        % Compute connectivity of new tracers:
        Temp1=Inds_t(IsLooseEnd_t);
        [~,Temp2]=sort(Zhat_t(Temp1),'descend');
        Temp3=Temp1(Temp2);
        Connectivity_t_new=[linspace(1,tracerxsize-1,tracerxsize-1)',linspace(2,tracerxsize,tracerxsize-1)';...
            linspace(tracerxsize+1,tracerxsize+tracerzsize,tracerzsize)',Temp3+numnewtracers];
        clear Temp1 Temp2 Temp3
        
        % Tie off old loose ends:
        IsLooseEnd_t(:)=0;
        
        % Update spawn timer:
        nextspawntime_u=time+spawninterval_u;
        nextspawntime_l=time+spawninterval_l;
        
    elseif time>=nextspawntime_u    % Upper surface only
        
        % Kill an old spawn point, add a new spawn point:
        Zhat_t_l=[1;Zhat_t_l(1:end-1)];
        Age_t_l_yr=[0;Age_t_l_yr(1:end-1)];
        IsLooseEnd_t(Age_t_yr==max(Age_t_yr(IsLooseEnd_t)))=0;
        
        % Compute locations of new tracers:  
        X_t_acc=linspace(0,domainwidth_mid,tracerxsize)';
        X_t_new=X_t_acc;
        Zhat_t_new=ones(tracerxsize,1);
        Age_t_new_yr=zeros(tracerxsize,1);
        numnewtracers=length(X_t_new);
        numnewconnectors=tracerxsize-1;
        
        % Compute unique ID numbers of new tracers:
        ID_t_new=linspace(nextid,nextid+numnewtracers-1,numnewtracers)';
        nextid=nextid+numnewtracers;
        
        % Compute logicals of new tracers:
        InDomain_t_new=true(numnewtracers,1);
        IsOriginal_t_new=false(numnewtracers,1);
        IsAcc_t_new=false(numnewtracers,1);
        IsLooseEnd_t_new=[true(1,1);false(tracerxsize-1,1)];
        
        % Compute connectivity of new tracers:
        Connectivity_t_new=[linspace(1,tracerxsize-1,tracerxsize-1)',linspace(2,tracerxsize,tracerxsize-1)'];
        
        % Update spawn timer:
        nextspawntime_u=time+spawninterval_u;
        
    elseif time>=nextspawntime_d % lower surface only
        
        % Compute locations of new tracers:  
        X_t_acc=linspace(0,domainwidth_mid,tracerxsize)';
        X_t_new=X_t_acc;
        Zhat_t_new=zeros(tracerxsize,1);
        Age_t_new_yr=zeros(tracerxsize,1);
        numnewtracers=length(X_t_new);
        numnewconnectors=tracerxsize-1;
        
        % Compute unique ID numbers of new tracers:
        ID_t_new=linspace(nextid,nextid+numnewtracers-1,numnewtracers)';
        nextid=nextid+numnewtracers;
        
        % Compute logicals of new tracers:
        InDomain_t_new=true(numnewtracers,1);
        IsOriginal_t_new=false(numnewtracers,1);
        IsAcc_t_new=true(tracerxsize,1);
        IsLooseEnd_t_new=false(tracerxsize,1);
        
        % Compute connectivity of new tracers:
        Connectivity_t_new=[linspace(1,tracerxsize-1,tracerxsize-1)',linspace(2,tracerxsize,tracerxsize-1)'];
        
        % Update spawn timer:
        nextspawntime_d=time+spawninterval_d;
        
    else % left edge only
        
        % Compute locations of new tracers:
        X_t_new=zeros(tracerzsize+1,1);
        Zhat_t_new=[Zhat_t_l;icacczhat];
        Age_t_new_yr=[Age_t_l_yr;0];
        numnewtracers=length(X_t_new);
        numnewconnectors=tracerzsize+1;
        
        % Compute unique ID numbers of new tracers:
        ID_t_new=linspace(nextid,nextid+numnewtracers-1,numnewtracers)';
        nextid=nextid+numnewtracers;
        
        % Compute logicals of new tracers:
        InDomain_t_new=true(numnewtracers,1);
        IsOriginal_t_new=false(numnewtracers,1);
        IsAcc_t_new=[false(tracerzsize,1);true(1)];
        IsLooseEnd_t_new=true(numnewtracers,1);
        
        % Compute connectivity of new tracers:
        Temp1=Inds_t(IsLooseEnd_t);
        [~,Temp2]=sort(Zhat_t(Temp1),'descend');
        Temp3=Temp1(Temp2);
        Connectivity_t_new=[linspace(1,tracerzsize+1,tracerzsize+1)',Temp3+numnewtracers];
        clear Temp1 Temp2 Temp3
        
        % Tie off old loose ends:
        IsLooseEnd_t(:)=0;
        
        % Update spawn timer:
        nextspawntime_l=time+spawninterval_l;
        
    end
    
    % Add new tracers to the list:
    ID_t=[ID_t_new;ID_t(1:end-numnewtracers)];
    X_t=[X_t_new;X_t(1:end-numnewtracers)];
    Zhat_t=[Zhat_t_new;Zhat_t(1:end-numnewtracers)];
    Age_t_yr=[Age_t_new_yr;Age_t_yr(1:end-numnewtracers)];
    InDomain_t=[InDomain_t_new;InDomain_t(1:end-numnewtracers)];
    IsOriginal_t=[IsOriginal_t_new;IsOriginal_t(1:end-numnewtracers)];
    IsAcc_t=[IsAcc_t_new;IsAcc_t(1:end-numnewtracers)];
    IsLooseEnd_t=[IsLooseEnd_t_new;IsLooseEnd_t(1:end-numnewtracers)];
    
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
    
    if length(Zhat_t)~=length(X_t) || length(IsLooseEnd_t)~=numtracers
        error('wtf?')
    end
    
    % Clear extraneous variables:
    clear *t_new*
    
end
