% GetFlowlineTimeSeries_v1

% Mike Wolovick, 10/25/2016

% This script loads the time series variables from multiple flowline model
% runs and combines them into a single file.

% The script is set up to list the files in the same order as they appear 
% in my notebook. (note: doesn't work for climate experiment)

% This script appends the timeseries to an existing output file if the
% output file already exists on disc

clearvars
t1=tic;

%% Parameters:

% Files and folders:
inputfolder='/work/Michael.Wolovick/FjordSillPaper/ModelOutput/';
outputfile='/net/mjw/FjordSillPaper/ModelOutput/ParamSweepRuns_v3_v3a.mat';

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


%% Work:

% Create units note:
NOTE_units='Time is in yr, position is in m, volume is in m^3, and fluxes are in m^3/yr';

% Pre-allocate output structures:
TimeSeries1=struct('filename',cell(numruns,1),'Time',cell(numruns,1),'X_gl',cell(numruns,1),...
    'Domainwidth',cell(numruns,1),'Volume',cell(numruns,1),'VAF',cell(numruns,1),...
    'Accumulation',cell(numruns,1),'Ablation',cell(numruns,1),'OceanMelt',cell(numruns,1),'Calving',cell(numruns,1),'SideInflux',cell(numruns,1));
Parameters1=struct('ModelParameters',cell(numruns,1),'SillParameters',cell(numruns,1),'TopParameters',cell(numruns,1),...
    'BottomParameters',cell(numruns,1),'SideParameters',cell(numruns,1),'CalvingParameters',cell(numruns,1),'PlumeParameters',cell(numruns,1),...
    'ThermalParameters',cell(numruns,1),'RheologyParameters',cell(numruns,1),'GridParameters',cell(numruns,1),'TimingParameters',cell(numruns,1),...
    'IterationParameters',cell(numruns,1),'VelocitySolverParameters',cell(numruns,1),'TracerParameters',cell(numruns,1),...
    'ScaleParameters',cell(numruns,1),'OtherParameters',cell(numruns,1));

% Loop through model runs:
for ii=1:numruns
    
    % Communicate:
    disp(['Model run=',num2str(ii),'/',num2str(numruns)])
    
    % Load this time series from this model run:
    load([inputfolder,inputfiles{ii}],'ModelTimeSeries','*Parameters')
    
    % Assign these time series:
    TimeSeries1(ii)=ModelTimeSeries;
    
    % Set filename properly:
    TimeSeries1(ii).filename=inputfiles{ii};
    
    % Assign parameters:
    Parameters1(ii).ModelParameters=ModelParameters;
    Parameters1(ii).SillParameters=SillParameters;
    Parameters1(ii).TopParameters=TopParameters;
    Parameters1(ii).BottomParameters=BottomParameters;
    Parameters1(ii).SideParameters=SideParameters;
    Parameters1(ii).CalvingParameters=CalvingParameters;
    Parameters1(ii).PlumeParameters=PlumeParameters;
    Parameters1(ii).ThermalParameters=ThermalParameters;
    Parameters1(ii).RheologyParameters=RheologyParameters;
    Parameters1(ii).GridParameters=GridParameters;
    Parameters1(ii).TimingParameters=TimingParameters;
    Parameters1(ii).IterationParameters=IterationParameters;
    Parameters1(ii).VelocitySolverParameters=VelocitySolverParameters;
    Parameters1(ii).TracerParameters=TracerParameters;
    Parameters1(ii).ScaleParameters=ScaleParameters;
    Parameters1(ii).OtherParameters=OtherParameters;
    
end

% Check if output file exists:
if exist(outputfile,'file')
    % Load existing variables:
    load(outputfile,'TimeSeries','Parameters')
    % Append new variables onto the old ones:
    TimeSeries=[TimeSeries;TimeSeries1];
    Parameters=[Parameters;Parameters1];
else
    % Rename variables:
    TimeSeries=TimeSeries1;
    Parameters=Parameters1;
end

% Save output:
save(outputfile,'NOTE_units','TimeSeries','Parameters')

% Final display:
disp('Done!')
toc(t1)

