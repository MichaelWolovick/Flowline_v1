% Flowline_GetMultipleGL_v1

% Mike Wolovick, 3/13/2017

% This script computes a time series of multiple grounding line positions
% based on model output of Flowline_v1.  This script exists because the
% built-in time series records of Flowline_v1 only record the most seaward
% grounding line position, but the model can have multiple grounding lines.

% The script tries to identify grounding lines that are continuous from one
% output timestep to the next.  The goal is for each grounding line to be
% listed as a continuous time series.  The algorithm works by looking at
% the grounding lines in the old timestep first, and then identifying the
% closest new grounding line.  If a new point is the closest one to
% multiple old points, priority is given to the shortest distance and the
% other old grounding lines are assumed to terminate within the intervening
% timestep.

% This script operates on one model run at a time.

% The runs are saved in a structure, so that the different time series
% don't have to be the same length.

% Note: the script assumes that sealevel=0.  It assumes that the plume
% model is on.

clear all
t1=tic;

%% Parameters:

% Files and folders:
inputfiles={'/work/Michael.Wolovick/FjordSillPaper/ModelOutput/ThwaitesC_v2_s4_m03_ch_v3a.mat';...
    '/work/Michael.Wolovick/FjordSillPaper/ModelOutput/ThwaitesC_v2_w_m03_ch_v3.mat';...
    '/work/Michael.Wolovick/FjordSillPaper/ModelOutput/ThwaitesC_v2_c_m03_ch_v3.mat'};  % cell array of strings
outputfolder='/net/mjw/FjordSillPaper/ModelOutput/';
outputsuffix='multiplegls'; % do not include underscore or .mat

% Number of grounding lines pre-allocated in this script:
numallocatedgls=20;           % integer

% Max grounding line jump between adjacent records:
maxglchange=7.5e3;            % m


%% Work:

% Create units note:
NOTE_units='Time is in yr, position is in m';

% Compute number of input files:
numinputfiles=length(inputfiles);

% loop through input files:
for thisfile=1:numinputfiles
    
    % Load common elements from this model run:
    load(inputfiles{thisfile},'numrecords','numdigits','*Parameters','*_input')
    
    % Pre-allocate output:
    Time=zeros(numrecords(end),1);
    X_gls=NaN*zeros(numrecords(end),numallocatedgls);
    Thesexgls=zeros(1,numallocatedgls);
    
    % Assume that sealevel=0:
    sealevel=0;
    
    % Loop through model records:
    for ii=1:numrecords(end)
        
        % Communicate:
        if rem(ii-1+1,100)==1
            disp(['record=',num2str(ii),'/',num2str(numrecords(end)-1+1)])
        end
        
        % Construct name of this model record:
        prefix='0'*ones(1,numdigits(end)-floor(log10(ii))-1);
        idnumber=num2str(ii);
        ThisName=['ModelRecord_',prefix,idnumber];
        
        % Load this model record:
        load(inputfiles{thisfile},ThisName)
        
        % Unpack model record:
        eval(['unpack(',ThisName,')'])
        
        % Define sill profile:
        if SillParameters.dosill
            % Compute sill thickness profile:
            if SillParameters.sillstarttime_yr==0 && SillParameters.sillconstructiontime_yr==0 && time_yr==0
                SillThick_input=SillParameters.sillamp*exp(-.5*(((X_input-SillParameters.sillcentx)/(.5*SillParameters.sillwidth)).^2));
            else
                SillThick_input=SillParameters.sillamp*max(0,min(1,(time_yr-SillParameters.sillstarttime_yr)/(SillParameters.sillconstructiontime_yr+1)))*exp(-.5*(((X_input-SillParameters.sillcentx)/(.5*SillParameters.sillwidth)).^2)); % note stabilizer
            end
        else
            SillThick_input=zeros(size(X_input));
        end
        
        % Assign time:
        Time(ii)=time_yr;
        
        % Compute dx:
        dx=domainwidth/GridParameters.xsize;
        
        % Produce horizontal grid:
        X_lr=linspace(0,domainwidth,GridParameters.xsize+1);
        X_c=.5*(X_lr(1:end-1)+X_lr(2:end));
        
        % Interpolate width and bed elevation:
        BedElev_c=interp1(X_input,BedElev_input+SillThick_input,X_c,ModelParameters.interpstyle,BedElev_input(end));
        bedelev_r=interp1(X_input,BedElev_input+SillThick_input,domainwidth,ModelParameters.interpstyle,BedElev_input(end));
        
        % Extrapolate ice thickness to terminus:
        icethick_r=1.5*Icethick_c(end)-.5*Icethick_c(end-1);
        
        % Compute hydraulic head:
        HydroHead_c=BedElev_c+(ThermalParameters.rho_i/PlumeParameters.rho_sw)*Icethick_c;
        hydrohead_r=bedelev_r+(ThermalParameters.rho_i/PlumeParameters.rho_sw)*icethick_r;
        
        % Concatenate terminus hydrohead onto grid center vector:
        HydroHead_cr=[HydroHead_c,hydrohead_r];
        X_cr=[X_c,domainwidth];
        
        % Locate sea level crossings:
        sealevelcrossinginds=find(diff(sign(HydroHead_cr-sealevel))); % index of last grid cell before the crossing
        
        % Compute number of grounding lines right now:
        thisnumgls=length(sealevelcrossinginds);
        
        % Locate all grounding lines this timestep:
        % Check number of grounding lines:
        if thisnumgls>=1
            % Loop through grounding lines:
            for thisgl=1:thisnumgls
                % Locate this grounding line precisely:
                Thesexgls(thisgl)=interp1(HydroHead_cr(max(sealevelcrossinginds(thisgl)-4,1):min(sealevelcrossinginds(thisgl)+4,GridParameters.xsize+1)),...
                    X_cr(max(sealevelcrossinginds(thisgl)-4,1):min(sealevelcrossinginds(thisgl)+4,GridParameters.xsize+1)),sealevel,'pchip',domainwidth);
            end
        else
            % No grounding line means that the terminus is the grounding line:
            Thesexgls(1)=domainwidth;
            thisnumgls=1;
        end
        % Set remaining values to NaN:
        Thesexgls(thisnumgls+1:end)=NaN;
        
        % Figure out relationship to previous timestep's grounding lines:
        if ii>1
            % Pre-allocate:
            ActiveInds=false(1,numallocatedgls); % active members of running list
            UnusedInds=false(1,numallocatedgls); % new gls that have not been used for an old gl
            UnusedInds(1:thisnumgls)=1;
            NewlyUsedInds=false(1,numallocatedgls); % new gls that have just been assigned to an old gl
            IndexingKey=NaN*zeros(1,numallocatedgls); % index into temporary list. in order of running list
            % Iterate until no changes have been made:
            madechanges=1;
            counter=1;
            while madechanges==1
                % Reset madechanges variable:
                madechanges=0;
                % Loop through old grounding lines:
                for oldgl=1:numallocatedgls
                    % Check whether to continue:
                    if isnan(X_gls(ii-1,oldgl)) || ActiveInds(oldgl) || sum(UnusedInds)==0
                        continue
                    end
                    % Compute distances to new gls:
                    Distances=abs(X_gls(ii-1,oldgl)-Thesexgls);
                    % Find min unused distance:
                    mindist=min(Distances(UnusedInds));
                    % Check against maximum distance:
                    if mindist>maxglchange
                        continue
                    end
                    % Find the closest unused new point to the old point:
                    IndexingKey(oldgl)=find(UnusedInds&Distances==mindist,1,'first');
                    % Mark this new point as used:
                    NewlyUsedInds(IndexingKey(oldgl))=1;
                    % Flag this old point as active:
                    ActiveInds(oldgl)=1;
                    % Changes have been made:
                    madechanges=1;
                end
                % Record newly used new gls:
                UnusedInds(NewlyUsedInds)=0;
                % Loop through new grounding lines:
                for newgl=1:thisnumgls
                    % Check that each new gl is used exactly once:
                    if sum(IndexingKey==newgl)>1 % used multiple times
                        % Identify competing indices:
                        CompetingInds=find(IndexingKey==newgl);
                        % Compute distances to competing points:
                        Distances=abs(X_gls(ii-1,CompetingInds)-Thesexgls(newgl));
                        % Identify the shortest distance:
                        winner=find(Distances==min(Distances),1,'first');
                        % Unflag all but the shortest distance:
                        Losers=1:length(CompetingInds)~=winner;
                        ActiveInds(CompetingInds(Losers))=0;
                        IndexingKey(CompetingInds(Losers))=NaN;
                        % Changes have been made:
                        madechanges=1;
                    elseif sum(IndexingKey==newgl&isnan(X_gls(ii-1,:)))==1 % used but unattached
                        % Identify competing indices:
                        CompetingInds=find(ActiveInds==0&isnan(X_gls(ii-1,:))==0);
                        % Compute distance to competing points:
                        Distances=abs(X_gls(ii-1,CompetingInds)-Thesexgls(newgl));
                        mindist=min(Distances);
                        % Check that there is a nearby old point to attach to:
                        if isempty(Distances) || mindist>maxglchange
                            continue
                        end
                        % Identify the shortest distance:
                        winner=find(Distances==mindist,1,'first');
                        % Remove the new gl from its unattached position:
                        ActiveInds(IndexingKey==newgl)=0;
                        IndexingKey(IndexingKey==newgl)=NaN;
                        % Attach the new gl to the winner:
                        ActiveInds(CompetingInds(winner))=1;
                        IndexingKey(CompetingInds(winner))=newgl;
                        % Mark this new gl as being used:
                        UnusedInds(newgl)=0;
                        % Changes have been made:
                        madechanges=1;
                    elseif sum(IndexingKey==newgl)==0 % completely unused
                        % Assign lowest available index: (from running list)
                        for dummyind=1:numallocatedgls
                            if isnan(X_gls(ii-1,dummyind)) && ActiveInds(dummyind)==0
                                ActiveInds(dummyind)=1;
                                IndexingKey(dummyind)=newgl;
                                break
                            end
                        end
                        % Changes have been made:
                        madechanges=1;
                    end
                end
                % Force break from loop:
                if counter>100
                    error('Unable to converge on a grounding line indexing solution.')
                else
                    counter=counter+1;
                end
            end
            % Record grounding lines using indexing keys:
            X_gls(ii,ActiveInds)=Thesexgls(IndexingKey(isnan(IndexingKey)==0));
        else
            % Record in order:
            X_gls(ii,:)=sort(Thesexgls);
        end
        
        % Clear this model record:
        eval(['clear ',ThisName])
        
    end
    
    % Save output:
    slashnums=strfind(inputfiles{thisfile},'/');
    outputfile=[outputfolder,inputfiles{thisfile}(slashnums(end)+1:end-4),'_',outputsuffix,'.mat'];
    save(outputfile,'NOTE_units','Time','X_gls')
    
end

% Final display:
disp('Done!')
toc(t1)