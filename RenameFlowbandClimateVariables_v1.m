% RenameFlowbandClimateVariables_v1

% Mike Wolovick, 5/12/2016

% This script renames the width-averaged surface temperature and
% accumulation variables into the "_input" format that the model is looking
% for.  It is assumed that annual melt is zero.  (This script is meant to
% be run on PIG and Thwaites only)

clear all
tic

%% Parameters:

% File names and paths:
inputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/SimpleInversions/';
inputfiles={'PineIslandFlowband_v1_simpleinversion.mat';...
    'PineIslandFlowband_v2_simpleinversion.mat';...
    'PineIslandFlowband_v3_simpleinversion.mat';...
    'PineIslandFlowband_v4_simpleinversion.mat';...
    'PineIslandFlowband_v5_simpleinversion.mat';...
    'ThwaitesFlowband_v1_simpleinversion.mat';...
    'ThwaitesFlowband_v2_simpleinversion.mat'};


%% Work:

% Loop through input files:
for thisfile=1:length(inputfiles)
    
    % Load input file:
    load([inputfolder,inputfiles{thisfile}],'Accum_yr','SurfTemp')
    
    % Rename variables:
    Accum_yr_input=Accum_yr;
    SurfTemp_input=SurfTemp;
    
    % Zero ablation:
    AnnualMelt_yr_input=zeros(size(Accum_yr));
    
    % Save results:
    save([inputfolder,inputfiles{thisfile}],'*_input','-append')
    
end


% Final display:
disp('Done!')
toc