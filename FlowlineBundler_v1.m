% FlowlineBundler_v1

% Mike Wolovick, 5/24/2016

% This script bundles multiple model runs of Flowline_v1.  It requires that
% Flowline_v1 be in function mode.  

clear all
tic

%% Parameters:

% Input files:
inputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/';
% inputfiles=repmat({'HelheimA_v2.mat';...
%     'HelheimB_v2.mat';...
%     'HelheimC_v2.mat';...
%     'KangerA_v2.mat';...
%     'KangerB_v2.mat';...
%     'KangerC_v2.mat';...
%     'JakobshavnA_v2.mat';...
%     'JakobshavnB_v2.mat';...
%     'JakobshavnC_v2.mat';...
%     'PetermannA_v2.mat';...
%     'PetermannB_v2.mat';...
%     'PetermannC_v2.mat';...
%     'PineIslandA_v2.mat';...
%     'PineIslandB_v2.mat';...
%     'PineIslandC_v2.mat';...
%     'ThwaitesA_v2.mat';...
%     'ThwaitesB_v2.mat';...
%     'ThwaitesC_v2.mat'},[3,1]);
inputfiles=repmat({'ThwaitesA_v2.mat';...
    'ThwaitesC_v2.mat'},[2,1]);

% Output suffix:
% outputsuffix=[repmat({'_c_m10_chm_v3'},[18,1]);repmat({'_w_m10_chm_v3'},[18,1]);repmat({'_s_m10_chm_v3'},[18,1])];
%outputsuffix=[repmat({'_s4_m01_chm_v3a'},[2,1]);repmat({'_s5_m01_chm_v3a'},[2,1]);repmat({'_s4_m03_chm_v3a'},[2,1]);repmat({'_s5_m03_chm_v3a'},[2,1])];
outputsuffix=[repmat({'_s4_m10_chm_v3a'},[2,1]);repmat({'_s5_m10_chm_v3a'},[2,1])];

% Climate dependencies:
% climatespatialdependence=repmat([repmat('z',[12,1]);repmat('x',[6,1])],[3,1]);
% climatetimedependence=[repmat({'none'},[18,1]);repmat({'mix'},[36,1])];
% oceantimedependence=[zeros(18,1);ones(36,1)];
climatespatialdependence=repmat('x',[4,1]);
climatetimedependence=repmat({'mix'},[4,1]);
oceantimedependence=ones(4,1);

% Upstream influx:
% influx_l_yr=repmat([1.33;4.79;1.65;2.03;2.72;2.48;1.05;20.0;1.08;1.74;3.55;1.98;3.19;3.32;6.40;2.14;2.22;4.04]*1e4,[3,1]);
influx_l_yr=repmat([2.14;4.04]*1e4,[2,1]);

% Calving:
% calvingtype=repmat({'meltthick'},[54,1]);
% calvingparam1=repmat([655;680;754;589;682;868;971;938;1043;77;77;92;449;453;456;300;295;299],[3,1]);
% calvingparam2=repmat([5523;5797;6274;5883;6228;7931;7813;7203;8469;1015;1000;1027;2852;3080;3128;2082;3355;2825],[3,1]);
% calvingparam3=repmat([316;290;347;197;198;247;598;429;658;8.3;8.3;8.5;70;73;76;52;54;54],[3,1]);
calvingtype=repmat({'meltthick'},[4,1]);
calvingparam1=repmat([300;299],[2,1]);
calvingparam2=repmat([2082;2825],[2,1]);
calvingparam3=repmat([52;54],[2,1]);

% Plume parameter:
% plumeparam=repmat([59*ones(9,1);7.35;4.68;8.00;28.5;25.3;50.1;32.5;29.4;43.1]*1e-5,[3,1]);
plumeparam=repmat([3.25;4.31]*1e-4,[2,1]);

% Sill parameters:
% dosill=[zeros(36,1);ones(18,1)];
% silltopz=repmat([-200;-200;-200;-100;-100;-100;-150;-150;-150;-100;-100;-100;-100;-100;-100;-100;-100;-100],[3,1]);
dosill=ones(4,1);
silltopz=-250*ones(4,1);
sillstressscale=1e4*ones(4,1);
sillerosionparam=1e-4*ones(4,1);
sillblockingfraction=repmat([.5;.5;0;0],[1,1]);

% Activate/deactivate side input:
% usesidemassinput=repmat([0;1;0],[18,1]);
usesidemassinput=zeros(4,1);

% Sliding exponent:
% m=10*ones(54,1);
m=10*ones(4,1);

%% Run Model:

% Check sizes:
numruns=length(inputfiles);
if length(outputsuffix)~=numruns
    error('outputsuffix wrong length')
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
if length(sillblockingfraction)~=numruns
    error('sillblockingfraction wrong length')
end

% Loop through model runs and go!
for ii=1:numruns                               
    % Define this input file name:
    thisinputfile{1}=[inputfolder,inputfiles{ii}];
    % Check nomenclature:
    % Check scenario:
    if strcmp(outputsuffix{ii}(2),'c') && (strcmp(climatetimedependence{ii},'none')==0 || oceantimedependence(ii)~=0 || dosill(ii)~=0)
        error('Scenario nomenclature disagreement')
    elseif strcmp(outputsuffix{ii}(2),'w') && (strcmp(climatetimedependence{ii},'mix')==0 || oceantimedependence(ii)~=1 || dosill(ii)~=0)
        error('Scenario nomenclature disagreement')
    elseif strcmp(outputsuffix{ii}(2),'s') && (strcmp(climatetimedependence{ii},'mix')==0 || oceantimedependence(ii)~=1 || dosill(ii)~=1)
        error('Scenario nomenclature disagreement')
    end
    % Check exponent:
    if str2double(outputsuffix{ii}(5:6))~=m(ii) && str2double(outputsuffix{ii}(6:7))~=m(ii)
        error('Sliding exponent nomenclature disagreement')
    end
    % Check calving law:
    if (strcmp(outputsuffix{ii}(9:10),'h_') || strcmp(outputsuffix{ii}(10:11),'h_')) && strcmp(calvingtype{ii},'thick')==0
        error('Calving type nomenclature disagreement')
    elseif (strcmp(outputsuffix{ii}(9),'m') || strcmp(outputsuffix{ii}(9:10),'_m')) && strcmp(calvingtype{ii},'meltmultiplier')==0
        error('Calving type nomenclature disagreement')
    elseif (strcmp(outputsuffix{ii}(9:10),'hm') || strcmp(outputsuffix{ii}(10:11),'hm')) && strcmp(calvingtype{ii},'meltthick')==0
        error('Calving type nomenclature disagreement')
    elseif (strcmp(outputsuffix{ii}(9:10),'hu') || strcmp(outputsuffix{ii}(10:11),'hu')) && strcmp(calvingtype{ii},'uthick')==0
        error('Calving type nomenclature disagreement')
    end
    % Check sill scenario:
    if strcmp(outputsuffix{ii}(2:3),'s1') && sillstressscale(ii)~=1e3
        error('Sill scenario nomenclature disagreement')
    elseif strcmp(outputsuffix{ii}(2:3),'s2') && sillstressscale(ii)~=1e4
        error('Sill scenario nomenclature disagreement')
    elseif strcmp(outputsuffix{ii}(2:3),'s3') && sillstressscale(ii)~=1e5
        error('Sill scenario nomenclature disagreement')
    end
    if strcmp(outputsuffix{ii}(2:3),'s1') && sillerosionparam(ii)~=1e-3
        error('Sill scenario nomenclature disagreement')
    elseif strcmp(outputsuffix{ii}(2:3),'s2') && sillerosionparam(ii)~=1e-4
        error('Sill scenario nomenclature disagreement')
    elseif strcmp(outputsuffix{ii}(2:3),'s3') && sillerosionparam(ii)~=1e-5
        error('Sill scenario nomenclature disagreement')
    end
    
    % Run model:
    Flowline_v1(thisinputfile,outputsuffix{ii},climatespatialdependence(ii),climatetimedependence{ii},oceantimedependence(ii),...
        influx_l_yr(ii),calvingtype{ii},calvingparam1(ii),calvingparam2(ii),calvingparam3(ii),plumeparam(ii),...
        dosill(ii),silltopz(ii),usesidemassinput(ii),m(ii),sillstressscale(ii),sillerosionparam(ii),sillblockingfraction(ii));
end

% Final display:
disp('...')
disp('DONE WITH MODEL RUNS')
toc


