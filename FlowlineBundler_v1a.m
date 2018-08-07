% FlowlineBundler_v1

% Mike Wolovick, 5/24/2016

% This script bundles multiple model runs of Flowline_v1.  It requires that
% Flowline_v1 be in function mode.  

% v1a:  modified for the climate experiment.  Scenario keys are in the
% input files now.

clear all
tic

%% Parameters:

% Input files:
inputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/ClimateExperimentInputs/';
inputfiles=[repmat({'ThwaitesA_v2'},[15,1]);repmat({'ThwaitesB_v2'},[15,1]);repmat({'ThwaitesC_v2'},[15,1])];

% Scenario keys:
scenariokeys=repmat({'c';'ac5';'ac10';'ac20';'ab3';'ab6';'ab9';'up1';'up2';'up3';'up4';'cdw.5';'cdw1';'cdw2';'upc'},[3,1]);

% Output suffix:
outputsuffix=repmat({'_mva_chm_clim1'},[45,1]);

% Climate dependencies:
climatespatialdependence=repmat('x',[45,1]);
climatetimedependence=[{'none'};repmat({'mix'},[14,1]);{'none'};repmat({'mix'},[14,1]);{'none'};repmat({'mix'},[14,1])];
oceantimedependence=[0;ones(14,1);0;ones(14,1);0;ones(14,1)];

% Upstream influx:
influx_l_yr=[2.14*ones(15,1);2.22*ones(15,1);4.04*ones(15,1)]*1e4;

% Calving:
calvingtype=repmat({'meltthick'},[45,1]);
calvingparam1=[300*ones(15,1);295*ones(15,1);299*ones(15,1)];
calvingparam2=[2082*ones(15,1);3355*ones(15,1);2825*ones(15,1)];
calvingparam3=[52*ones(15,1);54*ones(15,1);54*ones(15,1)];

% Plume parameter:
plumeparam=[3.25*ones(15,1);2.94*ones(15,1);4.31*ones(15,1)]*1e-4;

% Sill parameters:
dosill=zeros(45,1);
silltopz=-100*ones(45,1);
sillstressscale=1e4*ones(45,1);
sillerosionparam=1e-4*ones(45,1);

% Activate/deactivate side input:
usesidemassinput=zeros(45,1);

% Sliding exponent:
%m=10*ones(45,1);
m=repmat('file',[45,1]);

%% Run Model:

% Check sizes:
numruns=length(inputfiles);
if length(outputsuffix)~=numruns
    error('outputsuffix wrong length')
end
if length(scenariokeys)~=numruns
    error('scenariokeys wrong length')
end
if length(climatespatialdependence)~=numruns
    error('climatespatialdependence wrong length')
end
if length(climatetimedependence)~=numruns
    error('climatetimedependence wrong length')
end
if length(oceantimedependence)~=numruns
    error('oceantimedependence wrong length')
end
if length(influx_l_yr)~=numruns
    error('influx_l_yr wrong length')
end
if length(calvingtype)~=numruns
    error('calvingtype wrong length')
end
if length(calvingparam1)~=numruns
    error('calvingparam1 wrong length')
end
if length(calvingparam2)~=numruns
    error('calvingparam2 wrong length')
end
if length(calvingparam3)~=numruns
    error('calvingparam3 wrong length')
end
if length(plumeparam)~=numruns
    error('plume parameter wrong length')
end
if length(dosill)~=numruns
    error('dosill wrong length')
end
if length(silltopz)~=numruns
    error('silltopz wrong length')
end
if length(usesidemassinput)~=numruns
    error('usesidemassinput wrong length')
end
if length(m)~=numruns
    error('m wrong length')
end
if length(sillstressscale)~=numruns
    error('sillstressscale wrong length')
end
if length(sillerosionparam)~=numruns
    error('sillerosionparam wrong length')
end

% Loop through model runs and go!
for ii=1:numruns                               
    % Define this input file name:
    thisinputfile{1}=[inputfolder,inputfiles{ii},'_',scenariokeys{ii},'.mat'];
    
    % Check nomenclature:
    % Check scenario:
    if strcmp(scenariokeys{ii},'c') && (strcmp(climatetimedependence{ii},'none')==0 || oceantimedependence(ii)~=0 || dosill(ii)~=0)
        error('Scenario nomenclature disagreement')
    elseif strcmp(scenariokeys{ii},'c')==0 && (strcmp(climatetimedependence{ii},'mix')==0 || oceantimedependence(ii)~=1 || dosill(ii)~=0)
        error('Scenario nomenclature disagreement')
    end
    % Check exponent:
    if (str2double(outputsuffix{ii}(3:4))~=m(ii) && strcmp(m(ii,:),'file')==0) 
        error('Sliding exponent nomenclature disagreement')
    elseif (strcmp(m(ii,:),'file') && strcmp(outputsuffix{ii}(3:4),'va')==0)
        error('Sliding exponent nomenclature disagreement')
    end
    % Check calving law:
    if strcmp(outputsuffix{ii}(7:8),'h_') && strcmp(calvingtype{ii},'thick')==0
        error('Calving type nomenclature disagreement')
    elseif strcmp(outputsuffix{ii}(7:8),'hm') && strcmp(calvingtype{ii},'meltthick')==0
        error('Calving type nomenclature disagreement')
    elseif strcmp(outputsuffix{ii}(7:8),'hu') && strcmp(calvingtype{ii},'uthick')==0
        error('Calving type nomenclature disagreement')
    end
    
    % Run model:
    Flowline_v1(thisinputfile,outputsuffix{ii},climatespatialdependence(ii),climatetimedependence{ii},oceantimedependence(ii),...
        influx_l_yr(ii),calvingtype{ii},calvingparam1(ii),calvingparam2(ii),calvingparam3(ii),plumeparam(ii),...
        dosill(ii),silltopz(ii),usesidemassinput(ii),m(ii,:),sillstressscale(ii),sillerosionparam(ii));
end

% Final display:
disp('...')
disp('DONE WITH MODEL RUNS')
toc


