% FlowlinePlumeCalibrater_v1

% Mike Wolovick, 8/11/2016

% This script calibrates the plume model for Flowline_v1.  It is based on 
% FlowlineBundler_v1.  It runs Flowline_v1 in function mode, and it assumes
% that you've uncommented the lines that compute "initialslope" and
% "finalslope".

% The calibration looks for the smallest change in the relative slope near
% the grounding line.  Relative slope is slope normalized by
% (H_groundingline-H_calvingfront)/(X_groundingline-X_calvingfront).  The
% idea is that I'm allowing the shelf thickness to vary over time, but I'm
% assuming that the present-day shelf shape is close to it's natural shape.

clear all
tic

%% Parameters:

% Calibration parameters:
stantonlims=[1e-5,1e-3];   
numsamples=10;

% Input files:
inputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/';
inputfiles={'PetermannA.mat';...
    'PetermannB.mat';...
    'PetermannC.mat';...
    'PineIslandA.mat';...
    'PineIslandB.mat';...
    'PineIslandC.mat';...
    'ThwaitesA.mat';...
    'ThwaitesB.mat';...
    'ThwaitesC.mat'};

% Output suffix:
outputsuffix=repmat({'_Test_v4'},[9,1]); % not used

% Climate dependencies:
climatespatialdependence=[repmat('z',[3,1]);repmat('x',[6,1])];
climatetimedependence=repmat({'none'},[9,1]);
oceantimedependence=zeros(9,1);

% Upstream influx:
influx_l_yr=[1.74;3.55;1.98;3.19;3.32;6.40;2.14;2.22;4.04]*1e4;

% Calving parameter:
calvingparam1=[71;75;94;449;453;456;300;295;299];

% Grid sizes:
xsize=[829;629;829;756;748;756;913;919;913];
xsize_plume=[132;132;120;111;111;108;78;81;78];
zsize_plume=[21;22;28;133;135;136;89;88;89];

% Timing:
dt_yr=.02*ones(9,1);
runtime_yr=300*ones(9,1); % not used

% Sill parameters:
% dosill=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
% silltopz=[-200;-200;-200;-100;-100;-100;-150;-150;-150;-100;-100;-100;-200;-200;-200;-200;-200;-200;-200];
dosill=zeros(9,1);
silltopz=[-100;-100;-100;-200;-200;-200;-200;-200;-200];

% Activate/deactivate side input:
usesidemassinput=repmat([0;1;0],[3,1]);

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
if length(calvingparam1)~=numruns
    error('calvingparam1 wrong length')
end
if length(xsize)~=numruns
    error('xsize wrong length')
end
if length(xsize_plume)~=numruns
    error('xsize_plume wrong length')
end
if length(zsize_plume)~=numruns
    error('zsize_plume wrong length')
end
if length(dt_yr)~=numruns
    error('dt_yr wrong length')
end
if length(runtime_yr)~=numruns
    error('runtime_yr wrong length')
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

% Loop through flowbands:
for ii=7:numruns   % NOTE 7!!!!!
    
    % Create grid for calibration:
    StantonNumbers=exp(linspace(log(stantonlims(1)),log(stantonlims(2)),numsamples)');
    InitialSlopes=zeros(numsamples,1);
    FinalSlopes=zeros(numsamples,1);
    
    % Get this input file:
    thisinputfile{1}=[inputfolder,inputfiles{ii}]
    
    % Loop through stanton numbers:
    for jj=1:numsamples
        % Communicate:
        disp(['Sample=',num2str(jj),'/',num2str(numsamples)])
        % Run model:
        try
            [InitialSlopes(jj),FinalSlopes(jj)]=Flowline_v1(thisinputfile,outputsuffix{ii},climatespatialdependence(ii),climatetimedependence{ii},oceantimedependence(ii),...
                influx_l_yr(ii),calvingparam1(ii),StantonNumbers(jj),xsize(ii),xsize_plume(ii),zsize_plume(ii),dt_yr(ii),runtime_yr(ii),dosill(ii),silltopz(ii),usesidemassinput(ii));
        catch errorid
            % Set slopes to NaN:
            InitialSlopes(jj)=NaN;
            FinalSlopes(jj)=NaN;
        end
    end
    
    % Check initial slopes:
    meaninitialslope=mean(InitialSlopes(isnan(InitialSlopes)==0));
    if max(abs(InitialSlopes(isnan(InitialSlopes)==0)-meaninitialslope))>1e-12
        error('Inconsistent Initial Slopes')
    end
    
    % Interpolate smoother vectors:
    StantonNumbers_highres=exp(linspace(log(stantonlims(1)),log(max(StantonNumbers(isnan(InitialSlopes)==0))),10*numsamples)');
    FinalSlopes_highres=interp1(log(StantonNumbers(isnan(InitialSlopes)==0)),FinalSlopes(isnan(InitialSlopes)==0),log(StantonNumbers_highres),'pchip');
    
    % Solve for best stanton number:
    % Find crossing index:
    ind1=find(diff(sign(FinalSlopes_highres-meaninitialslope)),1,'first');
    % Solve using pchip interpolation:
    bests=exp(interp1(FinalSlopes_highres(ind1-5:ind1+5)-meaninitialslope,log(StantonNumbers_highres(ind1-5:ind1+5)),0,'pchip'))
    
    % Make a figure showing the calibration:
    figure(1)
    hold off
    semilogx(StantonNumbers_highres,FinalSlopes_highres,'k')
    hold on
    plot(StantonNumbers,FinalSlopes,'Marker','s','MarkerFaceColor','k','Color','k','LineStyle','none')
    plot([1e-5,1e-3],meaninitialslope*[1,1],'r')
    plot(bests,meaninitialslope,'Marker','o','MarkerFaceColor','r','Color','r')
    text(bests,meaninitialslope-.05,['\Gamma_{TS}=',num2str(1e-6*round(1e6*bests))],'HorizontalAlignment','left','VerticalAlignment','top','Color','k')
    text(2e-5,meaninitialslope-.05,'Initial Slope','HorizontalAlignment','left','VerticalAlignment','top','Color','r')
    xlabel('Stanton Number (unitless)')
    ylabel('Slope Near Grounding Line After Ten Years, normalized by (H_{GL}-H_{C})/(X_{GL}-X_{C})')
    title(['Plume Model Calibration in ',thisinputfile{1}(47:end-4)],'interpreter','none')
    figname=['/home/mjw/Documents/FjordSillPaper/Figures/',thisinputfile{1}(47:end-4),'_plumecalibration_v5.png'];
    set(gcf,'PaperSize',[8,6])
    set(gcf,'PaperPosition',[0,0,8,6])
    print('-dpng',figname,'-r300')
    
end

% Final display:
disp('...')
disp('DONE WITH MODEL RUNS')
toc


