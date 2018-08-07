% FjordSillSnapshotsFigure_v1

% Mike Wolovick, 9/14/2016

% This figure makes a series of model snapshots for the artificial sill
% paper.

% v2: modified to make snapshots of all the model runs.  Snapshots shown
% every 100 years (11 total snapshots), on a 4x3 subplot array.  Axes ticks
% are set to 50 km if domainwidth<300km, to 100 km otherwise.  The script
% assumes that output records are at least a year long.

% Further modified to always do the same number of snapshots (11, counting
% IC) regardless of how long the  model run lasted.


clearvars
close all
tic

%% Parameters:

% File names and paths:
inputfolder='/work/Michael.Wolovick/FjordSillPaper/ModelOutput/';
outputfolder='/net/mjw/FjordSillPaper/Figures/';
figsuffix='_Snapshots.png';

% Color ribbon thickness:
ribbonthick=100;                % m (vertical thickness)

% Colormaps and limits:
oceantempcmap='parula';         % valid colormap string or 'fadedjet'
oceanmeltcmap='jet';            % valid colormap string or 'fadedjet'
velcmap='hsv';                  % valid colormap string or 'fadedjet'
oceantemplims='auto';           % [1x2] deg C or 'auto'
meltlims=[1,1e3];               % [1x2] m/yr
vellims=[1,1e4];                % [1x2] m/yr
meltscale='log';                % 'linear' or 'log'
velscale='log';                 % 'linear' or 'log'
includecalving=0;               % 0 or 1 (include calving rate in the melt ribbon)
showmonotonicbottom=0;          % 0 or 1 (put the melt color ribbon on the monotonic ice bottom)

% Line widths:
bedlinewidth=.5;                % points
surflinewidth=.5;               % points
gllinewidth=.5;                 % points
frontlinewidth=.5;              % points
elalinewidth=.5;                % points
otherlinewidth=.5;              % points

% Sky, bedrock, ice color, sill color, and ocean settings:
skycolor=[.75,.75,1];           % [R,G,B]
skyzsize=10;                    % integer
bedcolor=[.63,.32,.18];         % [R,G,B]
icecolor=[.45,.62,.78];         % [R,G,B]
sillcolor=[.5,.5,.5];           % [R,G,B]
oceanxsize=200;                 % integer
oceanzsize=50;                  % integer

% Other Settings:
ticklength=.02;                 % unitless
fontsize=12;                    % points

% Margins:
verttextbuffer=.07;             % unitless
titlebuffer=.05;                % unitless
horztextbuffer=.07;             % unitless
horznotextbuffer=.03;           % unitless

% Output size and resolution:
pagesize=[12,8];                % [1x2] inches
resolution=300;                 % dpi


% File naming convention:
% Flowband names:
% flowbandnames={'HelheimA_v2';...
%     'HelheimB_v2';...
%     'HelheimC_v2';...
%     'KangerA_v2';...
%     'KangerB_v2';...
%     'KangerC_v2';...
%     'JakobshavnA_v2';...
%     'JakobshavnB_v2';...
%     'JakobshavnC_v2';...
%     'PetermannA_v2';...
%     'PetermannB_v2';...
%     'PetermannC_v2';...
%     'PineIslandA_v2';...
%     'PineIslandB_v2';...
%     'PineIslandC_v2';...
%     'ThwaitesA_v2';...
%     'ThwaitesB_v2';...
%     'ThwaitesC_v2'};
flowbandnames={'ThwaitesA_v2';...
    'ThwaitesC_v2'};
% Scenarios:
scenariokeys={'s6'};
% Sliding rules:
mkeys={'m01';'m03';'m10'};
% Calving rules:
ckeys={'ch';'chu';'chm'};
% Version:
versionkey='v3a';

%% Build input file set:

% Compute number of runs:
numruns=size(flowbandnames,1)*size(scenariokeys,1)*size(mkeys,1)*size(ckeys,1);

% Pre-allocate input files:
inputfiles=cell(numruns,1);

% Set the counters:
flowbandcounter=1;
scenariocounter=1;
mcounter=1;
thisflowband=1;
thisscenario=1;
thism=1;
thisc=1;

% Loop through input files:
for ii=1:numruns
    
    % Construct the file name:
    inputfiles{ii}=[flowbandnames{thisflowband},'_',scenariokeys{thisscenario},'_',mkeys{thism},'_',ckeys{thisc},'_',versionkey,'.mat'];
    
    % Advance counters:
    flowbandcounter=flowbandcounter+1;
    scenariocounter=scenariocounter+1;
    mcounter=mcounter+1;
    
    % Advance this calving rule:
    thisc=thisc+1;
    if thisc>size(ckeys,1)
        thisc=1;
    end
    
    % Advance this sliding rule:
    if mcounter>size(ckeys,1)
        thism=thism+1;
        mcounter=1;
        if thism>size(mkeys,1)
            thism=1;
        end
    end
    
    % Advance this flowband:
    if flowbandcounter>size(ckeys,1)*size(mkeys,1)
        thisflowband=thisflowband+1;
        flowbandcounter=1;
        if thisflowband>size(flowbandnames,1)
            thisflowband=1;
        end
    end
    
    % Advance this scenario:
    if scenariocounter>size(flowbandnames,1)*size(mkeys,1)*size(ckeys,1)
        thisscenario=thisscenario+1;
        scenariocounter=1;
    end
    
end

%% Preparation:

% Hardcoded number of snapshots:
numsnapshots=11;

% Hardcoded 4x3 plot arrangement:
numvertplots=3;
numhorzplots=4;

% Compute subplot positions:
availablehorzspace=1-horztextbuffer-numhorzplots*horznotextbuffer;
availablevertspace=1-verttextbuffer-numvertplots*titlebuffer;
Boxes=cell(numhorzplots,numvertplots);
% Loop through boxes: (listed in transpose order)
for ii=1:numvertplots
    for jj=1:numhorzplots
        Boxes{jj,ii}=[horztextbuffer+(jj-1)*(horznotextbuffer+availablehorzspace/numhorzplots),...
            verttextbuffer+(numvertplots-ii)*(titlebuffer+availablevertspace/numvertplots),...
            availablehorzspace/numhorzplots,availablevertspace/numvertplots];
    end
end

% Make a figure:
figure(1)

% Make faded jet colormap:
FadedJet=colormap('jet');
FadedJet(3*round(size(FadedJet,1)/8)+1:round(size(FadedJet,1)/2),1)=linspace(0,1,round(size(FadedJet,1)/8))';
FadedJet(3*round(size(FadedJet,1)/8)+1:round(size(FadedJet,1)/2),3)=1;
FadedJet(round(size(FadedJet,1)/2)+1:5*round(size(FadedJet,1)/8),1)=1;
FadedJet(round(size(FadedJet,1)/2)+1:5*round(size(FadedJet,1)/8),3)=linspace(1,0,round(size(FadedJet,1)/8))';

% Assign colormaps:
if strcmp(oceantempcmap,'fadedjet')
    OceanTempcmap=FadedJet;
else
    OceanTempcmap=colormap(oceantempcmap);
end
if strcmp(oceanmeltcmap,'fadedjet')
    OceanMeltcmap=FadedJet;
else
    OceanMeltcmap=colormap(oceanmeltcmap);
end
if strcmp(velcmap,'fadedjet')
    Velcmap=FadedJet;
else
    Velcmap=colormap(velcmap);
end

% Convert color limits to log scale:
if strcmp(meltscale,'log')
    meltlims=log(meltlims);
elseif strcmp(meltscale,'linear')==0
    error('Parameter "meltscale" must be a string equal to "linear" or "log".')
end
if strcmp(velscale,'log')
    vellims=log(vellims);
elseif strcmp(velscale,'linear')==0
    error('Parameter "velscale" must be a string equal to "linear" or "log".')
end

%% Make Figures:

% Loop through model runs:
startrun=1;                          
plothandles=zeros(numsnapshots,1);
for runnum=startrun:numruns   
    
    % Delete all old subplots:
    for ii=1:numsnapshots
        try
            delete(plothandles(ii));
        end
    end
    
    % Load and unpack model output:
    load([inputfolder,inputfiles{runnum}],'*Parameters','*_input','Zhat_*','DZhat_*','numrecords','numdigits')
    unpack(ModelParameters)
    unpack(GridParameters)
    unpack(SideParameters)
    unpack(PlumeParameters)
    unpack(TimingParameters)
    unpack(TracerParameters)
    unpack(ThermalParameters)
    unpack(BottomParameters)
    unpack(OtherParameters)
    if exist('SillParameters','var')
        unpack(SillParameters)
        if exist('sillblockingfraction','var')==0
            sillblockingfraction=1;
        end
    else
        dosill=0;
    end
    if exist('interpmethod','var')==0
        interpmethod='linear';
    end
    
    % Get total runtime for this run:
    thisruntime=(numrecords-1)*recordinterval_yr;
    
    % Compute timing of snapshots:
    snapshottimes=linspace(0,thisruntime,11);
    
    % Kludge to show ocean melt at the first snapshot:
    snapshottimes(1)=.49;
    
    % Define ocean temperature limits:
    if strcmp(oceantemplims,'auto')
        if oceantimedependence
            theseoceantemplims=[min(min(Temperature_input)),max(max(Temperature_input))]-tmelt;
        else
            theseoceantemplims=[min(Temperature_input(:,1)),max(Temperature_input(:,1))]-tmelt;
        end
    else
        theseoceantemplims=oceantemplims;
    end
    
    % Define maximum domain width:
    maxdomainwidth=X_input(end);
    
    % Define input xsize:
    inputxsize=length(X_input);
    
    % Set x-limits:
    xlims=[0,maxdomainwidth]/1000;
    
    % Loop through snapshots:
    for snapshot=1:numsnapshots
        
        % Get this model record:
        if snapshottimes(snapshot)==0
            % This record is the IC:
            thisrecord=0;
            % Load IC:
            load([inputfolder,inputfiles{runnum}],'InitialConditions')
            ModelRecord=InitialConditions;
            clear InitialConditions
            % Assign things that aren't recorded in the IC:
            % time is zero:
            ModelRecord.time_yr=0;
            % zero surface melt:
            ModelRecord.MeltRate_u=ModelRecord.AnnualMelt_u;
            % interpolate basal drag:
            if strcmp(slidingstressscale,'file')
                X_lr=linspace(0,ModelRecord.domainwidth,xsize+1);
                ModelRecord.Drag_lrd=interp1(X_input,SlidingStressScale_input,X_lr,interpmethod);
            else
                ModelRecord.Drag_lrd=slidingstressscale*ones(1,xsize+1);
            end
            % interpolate velocity:
            ModelRecord.U_lr=interp1(X_input,VelMag_yr_input/secondsperyear,X_lr,interpmethod);
            % set plume melt rate to bottom of color limits:
            ModelRecord.MeltRate_c_plume=meltlims(1)*ones(1,xsize_plume+zsize_plume);
            % set calving rate to zero:
            ModelRecord.calvingrate_r=0;
        elseif snapshottimes(snapshot)==runtime_yr
            % This record is the final state:
            load([inputfolder,inputfiles{runnum}],'FinalConditions')
            if exist('FinalConditions','var')
                ModelRecord=FinalConditions;
                clear FinalConditions
            else
                % Break from snapshot loop:
                break
            end
        else
            % Determine this record number:
            thisrecord=round(snapshottimes(snapshot)/recordinterval_yr+.5);  % assumes that record time is the middle of the record
            % Load this model record:
            prefix='0'*ones(1,numdigits(end)-floor(log10(thisrecord))-1);
            idnumber=num2str(thisrecord);
            load([inputfolder,inputfiles{runnum}],['ModelRecord_',prefix,idnumber])
            % Check if model record exists:
            if exist(['ModelRecord_',prefix,idnumber],'var')
                % Rename model record:
                eval(['ModelRecord=ModelRecord_',prefix,idnumber,';'])
                clear(['ModelRecord_',prefix,idnumber])
            else
                % Break from snapshot loop:
                break
            end
        end
        
        % Assign sea level:
        if strcmp(rightbctype,'front') && strcmp(rightbcparam,'file')==0
            sealevel=rightbcparam;
        elseif strcmp(rightbctype,'front')
            sealevel=interp1(Time_yr_input,RightBCParam_input,snapshottimes(snapshot),interpmethod,RightBCParam_input(1));
        else
            sealevel=min(BedElev_input)-1;
        end
        
        % Define sill profile:
        if dosill && isfield(ModelRecord,'SillThick_input')
            SillThick_input=ModelRecord.SillThick_input;
        elseif dosill
            % Compute sill thickness profile:
            if sillstarttime_yr==0 && sillconstructiontime_yr==0 && snapshottimes(snapshot)==0
                SillThick_input=sillamp*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2));
            else
                SillThick_input=sillamp*max(0,min(1,(snapshottimes(snapshot)-sillstarttime_yr)/(sillconstructiontime_yr+1)))*exp(-.5*(((X_input-sillcentx)/(.5*sillwidth)).^2)); % note stabilizer
            end
        else
            SillThick_input=zeros(size(X_input));
        end
        
        % Define initial domain:
        dx=ModelRecord.domainwidth/xsize;
        X_lr=linspace(0,ModelRecord.domainwidth,xsize+1);
        X_c=linspace(.5*dx,ModelRecord.domainwidth-.5*dx,xsize);
        BedElev_lr=interp1(X_input,BedElev_input+SillThick_input,X_lr,interpmethod);
        BedElev_c=.5*(BedElev_lr(1:end-1)+BedElev_lr(2:end));
        
        % Define vertical grid:
        Zhat_ud=linspace(0,1,zsize+1)'/maxdensify+(1-1/maxdensify)*linspace(0,1,zsize+1)'.^densifypower;
        Zhat_c=.5*(Zhat_ud(1:end-1)+Zhat_ud(2:end));
        DZhat_c=Zhat_ud(2:end)-Zhat_ud(1:end-1);
        DZhat_ud=[Zhat_c(1);Zhat_c(2:end)-Zhat_c(1:end-1);1-Zhat_c(end)];  % first and last are half-cells
        
        % Calculate surface elevation:
        SurfElev_c=ModelRecord.IceBottom_c+ModelRecord.Icethick_c;
        SurfElev_lr=[1.5*SurfElev_c(1)-.5*SurfElev_c(2),.5*(SurfElev_c(1:end-1)+SurfElev_c(2:end)),1.5*SurfElev_c(end)-.5*SurfElev_c(end-1)]; % linear extrapolation
        
        % Compute ticks and limits:
        if snapshot==1
            % x-ticks:
            if maxdomainwidth<3e5
                xtick=50;
            else
                xtick=100;
            end
            % z-lims:
            elevlims=[500*floor(min(BedElev_input)/500),500*ceil(max(SurfElev_c)/500)];
            % z-tick:
            elevtick=500;
        end
        
        % Compute ELA:
        [MonotonicSurface,MonotonicSurfInd]=sort(SurfElev_c);
        ela=interp1(ModelRecord.Accum_u(MonotonicSurfInd)-ModelRecord.MeltRate_u(MonotonicSurfInd),MonotonicSurface,0,'linear');
        
        % Interpolate ice thickness to grid edges:
        Icethick_lr=[1.5*ModelRecord.Icethick_c(1)-.5*ModelRecord.Icethick_c(2),.5*(ModelRecord.Icethick_c(1:end-1)+ModelRecord.Icethick_c(2:end)),1.5*ModelRecord.Icethick_c(end)-.5*ModelRecord.Icethick_c(end-1)];
        
        % Interpolate ice bottom to grid edges:
        IceBottom_lr=[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(2),.5*(ModelRecord.IceBottom_c(1:end-1)+ModelRecord.IceBottom_c(2:end)),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)];
        
        % Interpolate geometry onto the grounding line: (hydrostatic assumption)
        lastgroundedind=find(X_c<=ModelRecord.x_gl,1,'last');
        bedelev_gl=interp1(X_lr,BedElev_lr,ModelRecord.x_gl,interpmethod,BedElev_lr(end));
        if lastgroundedind<xsize
            icethick_gl=(sealevel-bedelev_gl)*(rho_sw/rho_i);
            hasshelf=1;
        else
            icethick_gl=Icethick_lr(end);
            hasshelf=0;
        end
        surfelev_gl=bedelev_gl+icethick_gl;
        
        % Compute vertical exaggeration:
        length_h=pagesize(1)*Boxes{snapshot}(3);
        length_v=pagesize(2)*Boxes{snapshot}(4);
        density_h=(xlims(2)-xlims(1))/length_h;
        density_v=.001*(elevlims(2)-elevlims(1))/length_v;
        exag=density_h/density_v;
        
        % Compute surface and bed gradients and angles:
        SurfGrad_lr=gradient(SurfElev_lr,dx);
        BedGrad_c=[.5*(BedElev_c(1)-BedElev_lr(1)),(BedElev_lr(2:end)-BedElev_lr(1:end-1)),.5*(BedElev_lr(end)-BedElev_c(end))]/dx; % also includes domain edges
        SurfAngle_lr=atan(exag*SurfGrad_lr); % angle on exaggerated plot
        BedAngle_c=atan(exag*BedGrad_c); % angle on exaggerated plot
        
        % Compute monotonic floating ice bottom:
        % Set the monotonic bottom to the real ice bottom:
        MonotonicIceBottom_c=ModelRecord.IceBottom_c;
        % Check if we have a shelf and are using the monotonic bottom:
        if hasshelf && showmonotonicbottom
            % Detect downward sloping ice bottom in the shelf:
            DownwardSlopingShelf_lr=X_lr>ModelRecord.x_gl&([ModelRecord.IceBottom_c,1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)]-[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(2),ModelRecord.IceBottom_c])<0;
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
                        downstreamicebottomz=ModelRecord.IceBottom_c(d2);
                    else
                        downstreamicebottomz=1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1);
                    end
                    % Identify next ice bottom below that point:
                    nextind=find(ModelRecord.IceBottom_c(1:d2-1)<downstreamicebottomz,1,'last');
                    if isempty(nextind)
                        continue
                    end
                    nexticebottom=ModelRecord.IceBottom_c(nextind);
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
        
        
        
        % Create or modify color patches:
        if snapshot==1
            
            %             % ICE PATCH:
            %             IcePatch=struct('Vertices',[[X_lr,fliplr(X_lr)]'/1000,[IceBottom_lr,fliplr(SurfElev_lr)]'],...
            %                 'Faces',linspace(1,2*(xsize+1),2*(xsize+1)),'FaceVertexCData',icecolor);
            % Create patch face indexing:
            Faces=zeros((xsize+1)*zsize,4);
            Faces(1:zsize,:)=[linspace(1,zsize,zsize)',linspace(2,zsize+1,zsize)',linspace(zsize+3,2*zsize+2,zsize)',linspace(zsize+2,2*zsize+1,zsize)'];
            for ii=2:xsize+1
                Faces((ii-1)*zsize+1:ii*zsize,:)=repmat(Faces((ii-1)*zsize,:)+1,[zsize,1])+repmat(linspace(1,zsize,zsize)',[1,4]);
            end
            % Create patch structure:
            GridElev_ud=[ModelRecord.IceBottom_c;SurfElev_c];
            IcePatch=struct('Vertices',[reshape(repmat([0,X_c,ModelRecord.domainwidth]/1000,[zsize+1,1]),[],1),reshape([IceBottom_lr(1)+Zhat_ud*Icethick_lr(1),GridElev_ud,IceBottom_lr(end)+Zhat_ud*Icethick_lr(end)],[],1)],...
                'Faces',Faces,...
                'FaceVertexCData',zeros((xsize+1)*zsize,3));
            % Interpolate RGB data for ice patch structure:
            if strcmp(velscale,'log')
                IcePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,1),max(0,min(1,(log(abs(ModelRecord.U_lr'*secondsperyear))-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
                IcePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,2),max(0,min(1,(log(abs(ModelRecord.U_lr'*secondsperyear))-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
                IcePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,3),max(0,min(1,(log(abs(ModelRecord.U_lr'*secondsperyear))-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
            else
                IcePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,1),max(0,min(1,(ModelRecord.U_lr'*secondsperyear-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
                IcePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,2),max(0,min(1,(ModelRecord.U_lr'*secondsperyear-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
                IcePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,3),max(0,min(1,(ModelRecord.U_lr'*secondsperyear-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
            end
            
            % OCEAN PATCH:
            % Create the ocean grid:
            % Vertical grid:
            Zhat_ud_ocean=linspace(0,1,oceanzsize+1)';
            Zhat_c_ocean=.5*(Zhat_ud_ocean(1:end-1)+Zhat_ud_ocean(2:end));
            % Check whether bed topography at the grounding line is above sea level:
            ind1=find(BedElev_input+SillThick_input<sealevel&X_input>ModelRecord.x_gl,1,'first');
            ind2=find(X_input>ModelRecord.x_gl,1,'first');
            if ind1==ind2
                % The ocean starts at the grounding line:
                oceanstartx=ModelRecord.x_gl;
            else
                % The oceans starts where the bed goes below sea level:
                oceanstartx=interp1(BedElev_input(ind1-1:ind1)+SillThick_input(ind1-1:ind1),X_input(ind1-1:ind1),sealevel,'linear');
            end
            % Create ocean horizontal grid:
            X_lr_ocean=linspace(oceanstartx,maxdomainwidth,oceanxsize+1);
            % Create ocean top and bottom:
            BedElev_lr_ocean=interp1(X_input,BedElev_input+SillThick_input,X_lr_ocean);
            % Interpolate ocean properties in time:
            if oceantimedependence
                ThisTemperature_input=interp1(Time_yr_input,Temperature_input',snapshottimes(snapshot),'linear');
                if snapshottimes(snapshot)>Time_yr_input(end)
                    ThisTemperature_input=Temperature_input(:,end);
                end
            else
                ThisTemperature_input=Temperature_input(:,1);
            end
            % Pre-allocate:
            Temp_lrud_ocean=zeros(oceanzsize+1,oceanxsize+1);
            % Loop through grid columns to assign ocean temperature:
            ThisTemperature_input1=ThisTemperature_input;
            for ii=1:oceanxsize+1
                % Find sill depth:
                silldepth=sealevel-max(BedElev_lr_ocean(min(ii+1,oceanxsize+1):end));
                abovesillind=find(Depth_input<=silldepth,1,'last');
                % Interpolate water properties:
                if abovesillind>1
                    % Compute partial blocking:
                    ThisTemperature_input1(Depth_input>silldepth)=...
                        sillblockingfraction*ThisTemperature_input(abovesillind)+(1-sillblockingfraction)*ThisTemperature_input(Depth_input>silldepth);
                    % Interpolate:
                    Temp_lrud_ocean(:,ii)=interp1(Depth_input,ThisTemperature_input1-tmelt,sealevel-(BedElev_lr_ocean(ii)+Zhat_ud_ocean*(sealevel-BedElev_lr_ocean(ii))),'linear',ThisTemperature_input1(end)-tmelt);
                else
                    Temp_lrud_ocean(:,ii)=tmelt;
                end
            end
            % Create patch face indexing for ocean patch:
            Faces=zeros(oceanxsize*oceanzsize,4);
            Faces(1:oceanzsize,:)=[linspace(1,oceanzsize,oceanzsize)',linspace(2,oceanzsize+1,oceanzsize)',linspace(oceanzsize+3,2*oceanzsize+2,oceanzsize)',linspace(oceanzsize+2,2*oceanzsize+1,oceanzsize)'];
            for ii=2:oceanxsize
                Faces((ii-1)*oceanzsize+1:ii*oceanzsize,:)=repmat(Faces((ii-1)*oceanzsize,:)+1,[oceanzsize,1])+repmat(linspace(1,oceanzsize,oceanzsize)',[1,4]);
            end
            % Create ocean patch structure:
            OceanPatch=struct('Vertices',[reshape(repmat(X_lr_ocean/1000,[oceanzsize+1,1]),[],1),reshape(repmat(BedElev_lr_ocean,[oceanzsize+1,1])+repmat(Zhat_ud_ocean,[1,oceanxsize+1]).*repmat(sealevel-BedElev_lr_ocean,[oceanzsize+1,1]),[],1),zeros((oceanxsize+1)*(oceanzsize+1),1)],...
                'Faces',Faces,'FaceVertexCData',zeros((oceanxsize+1)*(oceanzsize+1),3));
            % Interpolate RGB data for ocean patch structure:
            OceanPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,1),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
            OceanPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,2),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
            OceanPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,3),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
            
            % SKY PATCH:
            % Interpolate surface elevation onto input grid:
            SurfElev_input=interp1(X_lr,SurfElev_lr,X_input,'linear',sealevel);
            % Create patch face indexing for sky patch:
            Faces=zeros((inputxsize-1)*skyzsize,4);
            Faces(1:skyzsize,:)=[linspace(1,skyzsize,skyzsize)',linspace(2,skyzsize+1,skyzsize)',linspace(skyzsize+3,2*skyzsize+2,skyzsize)',linspace(skyzsize+2,2*skyzsize+1,skyzsize)'];
            for ii=2:inputxsize-1
                Faces((ii-1)*skyzsize+1:ii*skyzsize,:)=repmat(Faces((ii-1)*skyzsize,:)+1,[skyzsize,1])+repmat(linspace(1,skyzsize,skyzsize)',[1,4]);
            end
            % Create sky patch structure:
            %SkyPatch=struct('Vertices',zeros(skyzsize+3,2),'Faces',Faces,'FaceVertexCData',repmat(skycolor,[(skyzsize+1)*inputxsize,1])+repmat([1,1,1]-skycolor,[(skyzsize+1)*inputxsize,1]).*repmat(reshape(repmat(linspace(0,1,skyzsize+1)',[1,inputxsize]),[],1),[1,3]));
            SkyPatch=struct('Vertices',zeros((skyzsize+1)*inputxsize,2),'Faces',Faces,'FaceVertexCData',zeros((skyzsize+1)*inputxsize,3));
            % Assign values to sky patch:
            SkyPatch.Vertices=[reshape(repmat(X_input/1000,[skyzsize+1,1]),[],1),reshape(repmat(SurfElev_input,[skyzsize+1,1])+repmat(linspace(0,1,skyzsize+1)',[1,inputxsize]).*repmat(elevlims(2)-SurfElev_input,[skyzsize+1,1]),[],1)];
            SkyPatch.FaceVertexCData=repmat(skycolor,[(skyzsize+1)*inputxsize,1])+repmat([1,1,1]-skycolor,[(skyzsize+1)*inputxsize,1]).*repmat((SkyPatch.Vertices(:,2)-sealevel)/(elevlims(2)-sealevel),[1,3]);
            
            % BED PATCH:
            % Create bed patch:
            BedPatch=struct('Vertices',[[X_input,maxdomainwidth,0]'/1000,[BedElev_input,elevlims(1),elevlims(1)]'],...
                'Faces',linspace(1,length(X_input)+2,length(X_input)+2),'FaceVertexCData',bedcolor);
            
            % SILL PATCH:
            if dosill
                SillPatch=struct('Vertices',[[X_input,fliplr(X_input)]'/1000,[BedElev_input,fliplr(BedElev_input+SillThick_input)]'],...
                    'Faces',linspace(1,2*length(X_input),2*length(X_input)),'FaceVertexCData',sillcolor);
            end
            
            % OCEAN MELT PATCH:
            % Compute plume model grid size:
            dx_plume=(ModelRecord.domainwidth-ModelRecord.x_gl)/xsize_plume;
            % Create plume model horizontal grid:
            X_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
            Z_lr_plume=zeros(1,xsize_plume+zsize_plume+1);
            X_lr_plume(1:xsize_plume+1)=linspace(ModelRecord.x_gl,ModelRecord.domainwidth,xsize_plume+1);
            X_lr_plume(xsize_plume+2:xsize_plume+zsize_plume+1)=ModelRecord.domainwidth;
            X_c_plume=.5*(X_lr_plume(1:end-1)+X_lr_plume(2:end));
            % Compute plume model elevation:
            if hasshelf
                % Interpolate elevation under the floating shelf:
                Z_lr_plume(1:xsize_plume+1)=interp1([ModelRecord.x_gl,X_c(X_c>ModelRecord.x_gl),ModelRecord.domainwidth],...
                    [bedelev_gl,MonotonicIceBottom_c(X_c>ModelRecord.x_gl),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],...
                    X_lr_plume(1:xsize_plume+1),'linear');
            else
                Z_lr_plume(1:xsize_plume+1)=1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1);
            end
            % Compute elevation along the vertical face:
            Z_lr_plume(xsize_plume+1:xsize_plume+zsize_plume+1)=linspace(1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1),sealevel,zsize_plume+1);
            % Interpolate to grid centers:
            Z_c_plume=.5*(Z_lr_plume(1:end-1)+Z_lr_plume(2:end));
            % Compute plume grid gradient and angle:
            PlumeGrad_lr_plume=gradient(Z_lr_plume,dx_plume); % also includes domain edges
            PlumeAngle_lr_plume=atan(exag*PlumeGrad_lr_plume); % angle on exaggerated plot
            % Force vertical face to be vertical:
            PlumeAngle_lr_plume(xsize_plume+2:end)=pi/2;
            % Fix the angle just before the vertical face:
            PlumeAngle_lr_plume(xsize_plume+1)=PlumeAngle_lr_plume(xsize_plume);
            % Compute plume ribbon coordinates: (one extra point to accomodate right
            % angle turn at ice front)
            PlumeRibbonBotZ_lr=[Z_lr_plume-ribbonthick*cos(PlumeAngle_lr_plume),Z_lr_plume(xsize_plume+1)];
            PlumeRibbonBotX_lr=[X_lr_plume+exag*ribbonthick*sin(PlumeAngle_lr_plume),X_lr_plume(xsize_plume+1)+exag*ribbonthick];
            % Reduce under-shelf ribbon to zero size if there is no shelf:
            if hasshelf==0
                PlumeRibbonBotZ_lr(1:xsize_plume+1)=Z_lr_plume(1:xsize_plume+1);
                PlumeRibbonBotX_lr(1:xsize_plume+1)=X_lr_plume(1:xsize_plume+1);
            end
            % Create plume ribbon patch structure:
            PlumePatch=struct('Vertices',[[X_lr_plume,PlumeRibbonBotX_lr]'/1000,...
                [Z_lr_plume,PlumeRibbonBotZ_lr]'],...
                'Faces',[linspace(1,xsize_plume+zsize_plume,xsize_plume+zsize_plume)',...
                linspace(2,xsize_plume+zsize_plume+1,xsize_plume+zsize_plume)',...
                linspace(xsize_plume+zsize_plume+3,2*(xsize_plume+zsize_plume)+2,xsize_plume+zsize_plume)',...
                [linspace(xsize_plume+zsize_plume+2,2*xsize_plume+zsize_plume+1,xsize_plume),2*(xsize_plume+zsize_plume)+3,linspace(2*xsize_plume+zsize_plume+3,2*(xsize_plume+zsize_plume)+1,zsize_plume-1)]'],...% inserted extra point at right angle
                'FaceVertexCData',zeros(xsize_plume+zsize_plume,3));
            % Create combined melt and calving vector:
            MeltPlusCalving_c_plume=ModelRecord.MeltRate_c_plume;
            if includecalving==1
                MeltPlusCalving_c_plume(xsize_plume+1:end)=MeltPlusCalving_c_plume(xsize_plume+1:end)+ModelRecord.calvingrate_r;
            end
            % Interpolate plume patch RGB data:
            if strcmp(meltscale,'log')
                PlumePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,1),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltlims(1))))-meltlims(1))/(meltlims(2)-meltlims(1)))),'linear')));
                PlumePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,2),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltlims(1))))-meltlims(1))/(meltlims(2)-meltlims(1)))),'linear')));
                PlumePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,3),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltlims(1))))-meltlims(1))/(meltlims(2)-meltlims(1)))),'linear')));
            else
                PlumePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,1),(MeltPlusCalving_c_plume'*secondsperyear-meltlims(1))/(meltlims(2)-meltlims(1)),'linear')));
                PlumePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,2),(MeltPlusCalving_c_plume'*secondsperyear-meltlims(1))/(meltlims(2)-meltlims(1)),'linear')));
                PlumePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,3),(MeltPlusCalving_c_plume'*secondsperyear-meltlims(1))/(meltlims(2)-meltlims(1)),'linear')));
            end
            
        else
            
            % ICE PATCH:
            % IcePatch.Vertices=[[X_lr,fliplr(X_lr)]'/1000,[IceBottom_lr,fliplr(SurfElev_lr)]'];
            % Modify vertices:
            GridElev_ud=[ModelRecord.IceBottom_c;SurfElev_c];
            IcePatch.Vertices=[reshape(repmat([0,X_c,ModelRecord.domainwidth]/1000,[zsize+1,1]),[],1),reshape([IceBottom_lr(1)+Zhat_ud*Icethick_lr(1),GridElev_ud,IceBottom_lr(end)+Zhat_ud*Icethick_lr(end)],[],1)];
            % Interpolate RGB data for ice patch structure:
            if strcmp(velscale,'log')
                IcePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,1),max(0,min(1,(log(abs(ModelRecord.U_lr'*secondsperyear))-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
                IcePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,2),max(0,min(1,(log(abs(ModelRecord.U_lr'*secondsperyear))-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
                IcePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,3),max(0,min(1,(log(abs(ModelRecord.U_lr'*secondsperyear))-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
            else
                IcePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,1),max(0,min(1,(ModelRecord.U_lr'*secondsperyear-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
                IcePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,2),max(0,min(1,(ModelRecord.U_lr'*secondsperyear-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
                IcePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(Velcmap,1))',Velcmap(:,3),max(0,min(1,(ModelRecord.U_lr'*secondsperyear-vellims(1))/(vellims(2)-vellims(1)))),'linear')));
            end
            
            % OCEAN PATCH:
            % Check whether bed topography at the grounding line is above sea level:
            ind1=find(BedElev_input+SillThick_input<sealevel&X_input>ModelRecord.x_gl,1,'first');
            ind2=find(X_input>ModelRecord.x_gl,1,'first');
            if ind1==ind2
                % The ocean starts at the grounding line:
                oceanstartx=ModelRecord.x_gl;
            else
                % The oceans starts where the bed goes below sea level:
                oceanstartx=interp1(BedElev_input(ind1-1:ind1)+SillThick_input(ind1-1:ind1),X_input(ind1-1:ind1),sealevel,'linear');
            end
            % Create ocean horizontal grid:
            X_lr_ocean=linspace(oceanstartx,maxdomainwidth,oceanxsize+1);
            % Create ocean top and bottom:
            BedElev_lr_ocean=interp1(X_input,BedElev_input+SillThick_input,X_lr_ocean);
            % Interpolate ocean properties in time:
            if oceantimedependence
                ThisTemperature_input=interp1(Time_yr_input,Temperature_input',ModelRecord.time_yr,'linear');
                if ModelRecord.time_yr>Time_yr_input(end)
                    ThisTemperature_input=Temperature_input(:,end);
                end
            else
                ThisTemperature_input=Temperature_input(:,1);
            end
            % Loop through grid columns to assign ocean temperature:
            ThisTemperature_input1=ThisTemperature_input;
            for ii=1:oceanxsize+1
                % Find sill depth:
                silldepth=sealevel-max(BedElev_lr_ocean(min(ii+1,oceanxsize+1):end));
                abovesillind=find(Depth_input<=silldepth,1,'last');
                % Interpolate water properties:
                if abovesillind>1
                    % Compute partial blocking:
                    ThisTemperature_input1(Depth_input>silldepth)=...
                        sillblockingfraction*ThisTemperature_input(abovesillind)+(1-sillblockingfraction)*ThisTemperature_input(Depth_input>silldepth);
                    % Interpolate:
                    Temp_lrud_ocean(:,ii)=interp1(Depth_input,ThisTemperature_input1-tmelt,sealevel-(BedElev_lr_ocean(ii)+Zhat_ud_ocean*(sealevel-BedElev_lr_ocean(ii))),'linear',ThisTemperature_input1(end)-tmelt);
                else
                    Temp_lrud_ocean(:,ii)=tmelt;
                end
            end
            % Update ocean patch structure:
            OceanPatch.Vertices=[reshape(repmat(X_lr_ocean/1000,[oceanzsize+1,1]),[],1),reshape(repmat(BedElev_lr_ocean,[oceanzsize+1,1])+repmat(Zhat_ud_ocean,[1,oceanxsize+1]).*repmat(sealevel-BedElev_lr_ocean,[oceanzsize+1,1]),[],1),zeros((oceanxsize+1)*(oceanzsize+1),1)];
            % Interpolate RGB data for ocean patch structure:
            OceanPatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,1),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
            OceanPatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,2),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
            OceanPatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanTempcmap,1))',OceanTempcmap(:,3),max(0,min(1,(reshape(Temp_lrud_ocean,[],1)-theseoceantemplims(1))/(theseoceantemplims(2)-theseoceantemplims(1)))),'linear')));
            
            % SKY PATCH:
            % Interpolate surface elevation onto input grid:
            SurfElev_input=interp1(X_lr,SurfElev_lr,X_input,'linear',sealevel);
            % Update sky patch:
            SkyPatch.Vertices=[reshape(repmat(X_input/1000,[skyzsize+1,1]),[],1),reshape(repmat(SurfElev_input,[skyzsize+1,1])+repmat(linspace(0,1,skyzsize+1)',[1,inputxsize]).*repmat(elevlims(2)-SurfElev_input,[skyzsize+1,1]),[],1)];
            SkyPatch.FaceVertexCData=repmat(skycolor,[(skyzsize+1)*inputxsize,1])+repmat([1,1,1]-skycolor,[(skyzsize+1)*inputxsize,1]).*repmat((SkyPatch.Vertices(:,2)-sealevel)/(elevlims(2)-sealevel),[1,3]);
            
            % SILL PATCH:
            if dosill
                SillPatch.Vertices(:,2)=[BedElev_input,fliplr(BedElev_input+SillThick_input)]';
            end
            
            % OCEAN MELT PATCH:
            % Compute plume model grid size:
            dx_plume=(ModelRecord.domainwidth-ModelRecord.x_gl)/xsize_plume;
            % Create plume model horizontal grid:
            X_lr_plume(1:xsize_plume+1)=linspace(ModelRecord.x_gl,ModelRecord.domainwidth,xsize_plume+1);
            X_lr_plume(xsize_plume+2:xsize_plume+zsize_plume+1)=ModelRecord.domainwidth;
            X_c_plume=.5*(X_lr_plume(1:end-1)+X_lr_plume(2:end));
            % Compute plume model elevation:
            if hasshelf
                % Interpolate elevation under the floating shelf:
                Z_lr_plume(1:xsize_plume+1)=interp1([ModelRecord.x_gl,X_c(X_c>ModelRecord.x_gl),ModelRecord.domainwidth],...
                    [bedelev_gl,MonotonicIceBottom_c(X_c>ModelRecord.x_gl),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],...
                    X_lr_plume(1:xsize_plume+1),'linear');
            else
                Z_lr_plume(1:xsize_plume+1)=1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1);
            end
            % Compute elevation along the vertical face:
            Z_lr_plume(xsize_plume+1:xsize_plume+zsize_plume+1)=linspace(1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1),sealevel,zsize_plume+1);
            % Interpolate to grid centers:
            Z_c_plume=.5*(Z_lr_plume(1:end-1)+Z_lr_plume(2:end));
            % Compute plume grid gradient and angle:
            PlumeGrad_lr_plume=gradient(Z_lr_plume,dx_plume); % also includes domain edges
            PlumeAngle_lr_plume=atan(exag*PlumeGrad_lr_plume); % angle on exaggerated plot
            % Force vertical face to be vertical:
            PlumeAngle_lr_plume(xsize_plume+2:end)=pi/2;
            % Fix the angle just before the vertical face:
            PlumeAngle_lr_plume(xsize_plume+1)=PlumeAngle_lr_plume(xsize_plume);
            % Compute plume ribbon coordinates: (one extra point to accomodate right
            % angle turn at ice front)
            PlumeRibbonBotZ_lr=[Z_lr_plume-ribbonthick*cos(PlumeAngle_lr_plume),Z_lr_plume(xsize_plume+1)];
            PlumeRibbonBotX_lr=[X_lr_plume+exag*ribbonthick*sin(PlumeAngle_lr_plume),X_lr_plume(xsize_plume+1)+exag*ribbonthick];
            % Reduce under-shelf ribbon to zero size if there is no shelf:
            if hasshelf==0
                PlumeRibbonBotZ_lr(1:xsize_plume+1)=Z_lr_plume(1:xsize_plume+1);
                PlumeRibbonBotX_lr(1:xsize_plume+1)=X_lr_plume(1:xsize_plume+1);
            end
            % Upddate plume patch structure:
            PlumePatch.Vertices=[[X_lr_plume,PlumeRibbonBotX_lr]'/1000,[Z_lr_plume,PlumeRibbonBotZ_lr]'];
            % Create combined melt and calving vector:
            MeltPlusCalving_c_plume=ModelRecord.MeltRate_c_plume;
            if includecalving==1
                MeltPlusCalving_c_plume(xsize_plume+1:end)=MeltPlusCalving_c_plume(xsize_plume+1:end)+ModelRecord.calvingrate_r;
            end
            % Interpolate plume patch RGB data:
            if strcmp(meltscale,'log')
                PlumePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,1),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltlims(1))))-meltlims(1))/(meltlims(2)-meltlims(1)))),'linear')));
                PlumePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,2),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltlims(1))))-meltlims(1))/(meltlims(2)-meltlims(1)))),'linear')));
                PlumePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,3),max(0,min(1,(log(max(MeltPlusCalving_c_plume'*secondsperyear,exp(meltlims(1))))-meltlims(1))/(meltlims(2)-meltlims(1)))),'linear')));
            else
                PlumePatch.FaceVertexCData(:,1)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,1),(MeltPlusCalving_c_plume'*secondsperyear-meltlims(1))/(meltlims(2)-meltlims(1)),'linear')));
                PlumePatch.FaceVertexCData(:,2)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,2),(MeltPlusCalving_c_plume'*secondsperyear-meltlims(1))/(meltlims(2)-meltlims(1)),'linear')));
                PlumePatch.FaceVertexCData(:,3)=max(0,min(1,interp1(linspace(0,1,size(OceanMeltcmap,1))',OceanMeltcmap(:,3),(MeltPlusCalving_c_plume'*secondsperyear-meltlims(1))/(meltlims(2)-meltlims(1)),'linear')));
            end
            
        end
        
        % Create subplot:
        plothandles(snapshot)=subplot('Position',Boxes{snapshot});
        
        % Turn off hold:
        hold off
        
        % Plot ocean back patch:
        patch('Vertices',[[0;maxdomainwidth;maxdomainwidth;0],[sealevel;sealevel;elevlims(1);elevlims(1)]],...
            'Faces',linspace(1,4,4),'FaceColor',OceanTempcmap(1,:),'EdgeColor','none')
        hold on
        
        % Plot ocean patch:
        patch(OceanPatch,'FaceColor','interp','EdgeColor','none')
        
        % Plot bed patch:
        patch(BedPatch,'FaceColor','flat','EdgeColor','none')
        
        % Plot sill patch and line:
        if dosill
            patch(SillPatch,'FaceColor','flat','EdgeColor','none')
            plot(X_input/1000,BedElev_input+SillThick_input,'Color','k','LineWidth',bedlinewidth)
        end
        
        % Plot sky:
        patch(SkyPatch,'FaceColor','interp','EdgeColor','none')
        
        % Plot plume:
        patch(PlumePatch,'FaceColor','flat','EdgeColor','none')
        plot(PlumeRibbonBotX_lr(1:xsize_plume+1)/1000,PlumeRibbonBotZ_lr(1:xsize_plume+1),'Color','k','LineWidth',otherlinewidth)
        plot(PlumeRibbonBotX_lr(xsize_plume+2:end)/1000,PlumeRibbonBotZ_lr(xsize_plume+2:end),'Color','k','LineWidth',otherlinewidth)
        
        % Plot ice patch:
        patch(IcePatch,'FaceColor','flat','EdgeColor','none')
        
        % Plot bed line:
        plot(X_input/1000,BedElev_input,'Color','k','LineWidth',bedlinewidth)
        
        % Plot ice surface:
        plot([0,X_c,ModelRecord.domainwidth]/1000,[SurfElev_lr(1),SurfElev_c,SurfElev_lr(end)],'Color','k','LineWidth',surflinewidth)
        
        % Plot ice bottom:
        plot([0,X_c,ModelRecord.domainwidth]/1000,[1.5*ModelRecord.IceBottom_c(1)-.5*ModelRecord.IceBottom_c(1),ModelRecord.IceBottom_c,1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],'Color','k','LineWidth',bedlinewidth)
        
        % Plot ice front:
        plot(ModelRecord.domainwidth*[1,1]/1000,[SurfElev_lr(end),1.5*ModelRecord.IceBottom_c(end)-.5*ModelRecord.IceBottom_c(end-1)],'Color','k','LineWidth',frontlinewidth)
        
        % Plot grounding line:
        plot(ModelRecord.x_gl*[1,1]/1000,[bedelev_gl,surfelev_gl],'Color','k','LineWidth',gllinewidth)
        
        % Plot ELA:
        elainterceptx=X_c(find(SurfElev_c>=ela,1,'last'));
        plot([elainterceptx/1000,xlims(2)],ela*[1,1],'Color','k','LineStyle','--','LineWidth',elalinewidth)
        
        % Organize axes:
        set(gca,'XLim',xlims)
        set(gca,'XTick',xtick*[ceil(xlims(1)/xtick):1:floor(xlims(2)/xtick)])
        set(gca,'YLim',elevlims)
        set(gca,'YTick',elevtick*[ceil(elevlims(1)/elevtick):1:floor(elevlims(2)/elevtick)])
        if rem(snapshot,4)~=1
            set(gca,'YTickLabel',[])
        else
            ylabel('Elevation (m)','FontSize',fontsize)
        end
        if snapshot<=8
            set(gca,'XTickLabel',[])
        else
            xlabel('Distance (km)','FontSize',fontsize)
        end
        set(gca,'XMinorTick','off')
        set(gca,'YMinorTick','off')
        set(gca,'TickDir','out')
        set(gca,'TickLength',[ticklength,1])
        set(gca,'LineWidth',otherlinewidth)
        set(gca,'Box','on')
        set(gca,'FontSize',fontsize)
        
        % Make title:
        if snapshot==1
            title({inputfiles{runnum}(1:end-4);['Time= ',num2str(round(snapshottimes(snapshot))),' year']},'FontSize',fontsize,'interpreter','none')
        else
            title(['Time= ',num2str(round(snapshottimes(snapshot))),' year'],'FontSize',fontsize)
        end
        
    end
    
    % Create colorbars in remaining space:
    % Velocity colorbar:
    startx=Boxes{end}(1);
    width=Boxes{end}(3)/15;
    subplot('Position',[startx,Boxes{end}(2),width,Boxes{end}(4)])
    patch('Vertices',[[zeros(length(Velcmap),1);ones(length(Velcmap),1)],[linspace(vellims(1),vellims(2),length(Velcmap))';linspace(vellims(1),vellims(2),length(Velcmap))']],...
        'Faces',[linspace(1,length(Velcmap)-1,length(Velcmap)-1)',linspace(2,length(Velcmap),length(Velcmap)-1)',linspace(length(Velcmap)+2,2*length(Velcmap),length(Velcmap)-1)',linspace(length(Velcmap)+1,2*length(Velcmap)-1,length(Velcmap)-1)'],...
        'FaceVertexCData',[Velcmap;Velcmap],'FaceColor','interp','EdgeColor','none')
    xlim([0,1])
    ylim(vellims)
    set(gca,'XTick',[])
    if strcmp(velscale,'log')
        set(gca,'YTick',log([1,10,100,1000,10000]))
        set(gca,'YTickLabel',[1,10,100,1000,10000])
    end
    ylabel('Ice Velocity (m/yr)')
    set(gca,'YAxisLocation','right')
    % Melt colorbar:
    startx=Boxes{end}(1)+(1/3)*Boxes{end}(4);
    subplot('Position',[startx,Boxes{end}(2),width,Boxes{end}(4)])
    patch('Vertices',[[zeros(length(OceanMeltcmap),1);ones(length(OceanMeltcmap),1)],[linspace(meltlims(1),meltlims(2),length(OceanMeltcmap))';linspace(meltlims(1),meltlims(2),length(OceanMeltcmap))']],...
        'Faces',[linspace(1,length(OceanMeltcmap)-1,length(OceanMeltcmap)-1)',linspace(2,length(OceanMeltcmap),length(OceanMeltcmap)-1)',linspace(length(OceanMeltcmap)+2,2*length(OceanMeltcmap),length(OceanMeltcmap)-1)',linspace(length(OceanMeltcmap)+1,2*length(OceanMeltcmap)-1,length(OceanMeltcmap)-1)'],...
        'FaceVertexCData',[OceanMeltcmap;OceanMeltcmap],'FaceColor','interp','EdgeColor','none')
    xlim([0,1])
    ylim(meltlims)
    set(gca,'XTick',[])
    if strcmp(meltscale,'log')
        set(gca,'YTick',log([.001,.01,.1,1,10,100,1000,10000]))
        set(gca,'YTickLabel',[.001,.01,.1,1,10,100,1000,10000])
    end
    ylabel('Submarine Melt Rate (m/yr)')
    set(gca,'YAxisLocation','right')
    % Ocean temp colorbar:
    startx=Boxes{end}(1)+(2/3)*Boxes{end}(4);
    subplot('Position',[startx,Boxes{end}(2),width,Boxes{end}(4)])
    patch('Vertices',[[zeros(length(OceanTempcmap),1);ones(length(OceanTempcmap),1)],[linspace(theseoceantemplims(1),theseoceantemplims(2),length(OceanTempcmap))';linspace(theseoceantemplims(1),theseoceantemplims(2),length(OceanTempcmap))']],...
        'Faces',[linspace(1,length(OceanTempcmap)-1,length(OceanTempcmap)-1)',linspace(2,length(OceanTempcmap),length(OceanTempcmap)-1)',linspace(length(OceanTempcmap)+2,2*length(OceanTempcmap),length(OceanTempcmap)-1)',linspace(length(OceanTempcmap)+1,2*length(OceanTempcmap)-1,length(OceanTempcmap)-1)'],...
        'FaceVertexCData',[OceanTempcmap;OceanTempcmap],'FaceColor','interp','EdgeColor','none')
    xlim([0,1])
    ylim(theseoceantemplims)
    set(gca,'XTick',[])
    ylabel('Ocean Temperature (\circC)')
    set(gca,'YAxisLocation','right')
    
    % Construct figure name:
    figname=[outputfolder,inputfiles{runnum}(1:end-4),figsuffix];
    
    % Save figure:
    set(gcf,'PaperSize',pagesize)
    set(gcf,'PaperPosition',[0,0,pagesize])
    print('-dpng',figname,'-r300')
    
end


% Final Display:
disp('Done!')
toc