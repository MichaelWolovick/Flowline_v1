% FlowlineAnalyzer_v1

% Mike Wolovick, 6/15/2016

% This script loads the model output from Flowline_v1 and computes various
% time series.  The script saves the time series from multiple runs into a
% single file.

% Time series:
% Time (yr)
% Grounding line position (m)
% Grounding line elevation (m)
% Calving front position (m)
% Total volume (m^3)
% Volume above flotation (m^3)
% Accumulation (m^3/yr)
% Surface melt (m^3/yr)
% Ocean melt (m^3/yr)
% Calving flux (m^3/yr)
% Side influx (m^3/yr)

% The runs are saved in a structure, so that the different time series
% don't have to be the same length.

% Note: the script assuems that sealevel=0.  It assumes that the plume
% model is on.


clear all
t1=tic;

%% Parameters:

% Files and folders:
inputfolder='/net/mjw/FjordSillPaper/ModelOutput/';
inputfilekeys={'*ConstForcing_v4.mat';'*Warming_v4.mat';'*Sill_v4.mat'}; % cell array of strings, with wildcards (*) as appropriate
outputfile='/net/mjw/FjordSillPaper/ModelOutput/ConstForcing_Warming_Sill_v4_TimeSeries_v1.mat';

%% Work:

% Identify input files:
inputfiles=[];
for ii=1:length(inputfilekeys)
    thesefiles=dir([inputfolder,inputfilekeys{ii}]);
    inputfiles=[inputfiles;thesefiles];
end
clear thesefiles

% Determine number of runs to compile:
numruns=length(inputfiles);

% Create units note:
NOTE_units='Time is in yr, position is in m, volume is in m^3, and fluxes are in m^3/yr';

% Pre-allocate output structure:
ModelTimeSeries=struct('filename',cell(numruns,1),'Time',cell(numruns,1),'X_gl',cell(numruns,1),'Z_gl',cell(numruns,1),...
    'Domainwidth',cell(numruns,1),'Volume',cell(numruns,1),'VAF',cell(numruns,1),...
    'Accumulation',cell(numruns,1),'Ablation',cell(numruns,1),'OceanMelt',cell(numruns,1),'Calving',cell(numruns,1),'SideInflux',cell(numruns,1));

% % Load output file:
% load(outputfile,'ModelTimeSeries')

% Loop through model runs:
for jj=1:numruns
    
    %     % Skip appropriate runs:
    %     if jj<17 || (jj>19 && jj<48) || jj>50
    %         continue
    %     end
    
    % Communicate:
    disp(['Model run=',num2str(jj),'/',num2str(numruns)])
    
    % Load common elements from this model run:
    load([inputfolder,inputfiles(jj).name],'numrecords','numdigits','*Parameters','*_input')
    
    % Define input ice surface:
    SurfElev_input=max(BedElev_input+Icethick_input,(1-ThermalParameters.rho_i/PlumeParameters.rho_sw)*Icethick_input); % assumes sealevel=0
    SurfElev_input(isnan(SurfElev_input))=0;
    
    % Save file name:
    ModelTimeSeries(jj).filename=inputfiles(jj).name;
    
    % Pre-allocate time series within the output structure:
    ModelTimeSeries(jj).Time=zeros(numrecords(end),1);
    ModelTimeSeries(jj).X_gl=zeros(numrecords(end),1);
    ModelTimeSeries(jj).Z_gl=zeros(numrecords(end),1);
    ModelTimeSeries(jj).Domainwidth=zeros(numrecords(end),1);
    ModelTimeSeries(jj).Volume=zeros(numrecords(end),1);
    ModelTimeSeries(jj).VAF=zeros(numrecords(end),1);
    ModelTimeSeries(jj).Accumulation=zeros(numrecords(end),1);
    ModelTimeSeries(jj).Ablation=zeros(numrecords(end),1);
    ModelTimeSeries(jj).OceanMelt=zeros(numrecords(end),1);
    ModelTimeSeries(jj).Calving=zeros(numrecords(end),1);
    ModelTimeSeries(jj).SideInflux=zeros(numrecords(end),1);
    
    % Loop through model records:
    for ii=1:numrecords(end)
        
        % Communicate:
        if rem(ii,200)==0
            disp(['record=',num2str(ii),'/',num2str(numrecords(end))])
        end
        
        % Construct name of this model record:
        prefix='0'*ones(1,numdigits(end)-floor(log10(ii))-1);
        idnumber=num2str(ii);
        ThisName=['ModelRecord_',prefix,idnumber];
        
        % Load this model record:
        load([inputfolder,inputfiles(jj).name],ThisName)
        
        % Unpack model record:
        eval(['unpack(',ThisName,')'])
        
        % Compute dx:
        dx=domainwidth/GridParameters.xsize; 
        
        % Produce horizontal grid:
        X_lr=linspace(0,domainwidth,GridParameters.xsize+1);
        X_c=.5*(X_lr(1:end-1)+X_lr(2:end));
        
        % Interpolate width and bed elevation:
        BedElev_c=interp1(X_input,BedElev_input,X_c,'linear',BedElev_input(end));
        bedelev_r=interp1(X_input,BedElev_input,domainwidth,'linear',BedElev_input(end));
        Width_c=interp1(X_input,Width_input,X_c,'linear',Width_input(end));
        width_l=interp1(X_input,Width_input,0,'linear',Width_input(1));
        width_r=interp1(X_input,Width_input,domainwidth,'linear',Width_input(end));
        
        % Compute surface elevation:
        SurfElev_c=IceBottom_c+Icethick_c;
        
        % Compute flotation thickness:
        FlotationThickness_c=max(0,-(PlumeParameters.rho_sw/ThermalParameters.rho_i)*BedElev_c);
        
        % Compute boundary thicknesses and top/bottom elevations:
        icethick_l=max(1.5*Icethick_c(1)-.5*Icethick_c(2),GridParameters.minicethick);
        icethick_r=max(1.5*Icethick_c(end)-.5*Icethick_c(end-1),GridParameters.minicethick);
        icebottom_r=max(bedelev_r,-(ThermalParameters.rho_i/PlumeParameters.rho_sw)*icethick_r);
        surfelev_r=icebottom_r+icethick_r;
        
        % Determine if a floating shelf is present:
        if domainwidth-x_gl>2*dx
            % Mark that a shelf is present:
            hasshelf=1;
            % Compute plume grid:
            X_lr_plume=x_gl+(domainwidth-x_gl)*(linspace(0,1,GridParameters.xsize_plume+1)/GridParameters.plumemaxdensify+(1-1/GridParameters.plumemaxdensify)*linspace(0,1,GridParameters.xsize_plume+1).^GridParameters.plumedensifypower);
            X_c_plume=.5*(X_lr_plume(1:end-1)+X_lr_plume(2:end));
            % Compute dx of the plume grid:
            DX_c_plume=X_lr_plume(2:end)-X_lr_plume(1:end-1);
            % Interpolate width onto the plume grid:
            Width_c_plume=interp1(X_input,Width_input,X_c_plume,'linear',Width_input(end));
        else
            hasshelf=0;
        end
        
        % Interpolate side influx:
        if ModelParameters.usesidemassinput==1
            % Interpolate flux:
            SideInflux_yr_c=interp1(X_input,Influx_yr_input,X_c,'linear',Influx_yr_input(end));
            % Interpolate original surface:
            S0_c=interp1(X_input,SurfElev_input,X_c,'linear',SurfElev_input(end));
            % Adjust flux for surface changes:
            SideInflux_yr_c=SideInflux_yr_c.*max(0,1+(S0_c-SurfElev_c)./(S0_c+1)); % note stabilizer
        end
        
        % Assign time:
        ModelTimeSeries(jj).Time(ii)=time_yr;
        
        % Assign model dimensions:
        ModelTimeSeries(jj).X_gl(ii)=x_gl;
        ModelTimeSeries(jj).Z_gl(ii)=interp1(X_input,BedElev_input,x_gl,'linear',BedElev_input(end));
        ModelTimeSeries(jj).Domainwidth(ii)=domainwidth;
        
        % Compute model volume:
        ModelTimeSeries(jj).Volume(ii)=sum(dx*Width_c.*Icethick_c);
        ModelTimeSeries(jj).VAF(ii)=sum(dx*Width_c.*max(0,Icethick_c-FlotationThickness_c));
        
        % Compute fluxes:
        ModelTimeSeries(jj).Accumulation(ii)=sum(dx*Width_c.*Accum_u)*OtherParameters.secondsperyear;
        ModelTimeSeries(jj).Ablation(ii)=(sum(dx*Width_c.*MeltRate_u)+width_r*(surfelev_r-max(icebottom_r,0))*subaerialmeltrate_r)*OtherParameters.secondsperyear;
        if hasshelf
            ModelTimeSeries(jj).OceanMelt(ii)=(sum(DX_c_plume.*Width_c_plume.*MeltRate_c_plume(1:GridParameters.xsize_plume))+...
                width_r*(-icebottom_r)*oceanmeltrate_r)*OtherParameters.secondsperyear;
        else
            ModelTimeSeries(jj).OceanMelt(ii)=width_r*(-icebottom_r)*oceanmeltrate_r*OtherParameters.secondsperyear;
        end
        ModelTimeSeries(jj).Calving(ii)=width_r*icethick_r*calvingrate_r*OtherParameters.secondsperyear;
        ModelTimeSeries(jj).SideInflux(ii)=sum(dx*SideInflux_yr_c);
        
        % Clear this model record:
        eval(['clear ',ThisName])
        
    end
    
end

% Save output:
save(outputfile,'NOTE_units','ModelTimeSeries')

% Final display:
disp('Done!')
toc(t1)
return

%% Figures:

% This is an ad hoc section to make a bunch of figures showing the time
% series.

figfolder='/home/mjw/Documents/FjordSillPaper/Figures/';

% Loop through glaciers:
for flowband=1:12
    
    % Determine these files:
    constforcingfile=flowband;
    warmingfile=flowband+12;
    sillfile=flowband+24;
    
    % Determine glacier name:
    underscorenum=strfind(ModelTimeSeries(constforcingfile).filename,'_');
    glaciername=ModelTimeSeries(constforcingfile).filename(1:underscorenum(1)-1);
    
    % Call figure:
    figure(1)
    
    % Call first subplot:
    subplot(2,1,1)
    
    % Plot zero line:
    hold off
    plot([0,300],[0,0],'k')
    hold on
    
    % Plot grounding line:
    plot(ModelTimeSeries(constforcingfile).Time,(ModelTimeSeries(constforcingfile).X_gl-ModelTimeSeries(constforcingfile).X_gl(1))/1e3,'Color','b')
    plot(ModelTimeSeries(warmingfile).Time,(ModelTimeSeries(warmingfile).X_gl-ModelTimeSeries(warmingfile).X_gl(1))/1e3,'Color','r')
    plot(ModelTimeSeries(sillfile).Time,(ModelTimeSeries(sillfile).X_gl-ModelTimeSeries(sillfile).X_gl(1))/1e3,'Color','g')
    
    % Fix up axes and title:
    xlim([0,300])
    set(gca,'XTickLabel',[])
    ylim([min([ModelTimeSeries(constforcingfile).X_gl-ModelTimeSeries(constforcingfile).X_gl(1);ModelTimeSeries(warmingfile).X_gl-ModelTimeSeries(warmingfile).X_gl(1);ModelTimeSeries(sillfile).X_gl-ModelTimeSeries(sillfile).X_gl(1)]),...
        max([ModelTimeSeries(constforcingfile).X_gl-ModelTimeSeries(constforcingfile).X_gl(1);ModelTimeSeries(warmingfile).X_gl-ModelTimeSeries(warmingfile).X_gl(1);ModelTimeSeries(sillfile).X_gl-ModelTimeSeries(sillfile).X_gl(1)])]/1000)
    ylabel('Change in Grounding Line Position (km)')
    title(['v4 Runs for ',glaciername,', Grounding Line Position'],'interpreter','none')
    
    % Plot sill construction time:
    ylims=get(gca,'Ylim');
    patch('Vertices',[[100;100;110;110],[ylims(1);ylims(2);ylims(2);ylims(1)]],'Faces',1:4,'FaceVertexCData',[.5,.5,.5],'EdgeColor','none','FaceColor','flat')
    
    % Call second subplot:
    subplot(2,1,2)
    
    % Plot zero line:
    hold off
    plot([0,300],[0,0],'k')
    hold on
    
    % Plot volume above flotation:
    h1=plot(ModelTimeSeries(constforcingfile).Time,(ModelTimeSeries(constforcingfile).VAF-ModelTimeSeries(constforcingfile).VAF(1))/1e9,'Color','b');
    h2=plot(ModelTimeSeries(warmingfile).Time,(ModelTimeSeries(warmingfile).VAF-ModelTimeSeries(warmingfile).VAF(1))/1e9,'Color','r');
    h3=plot(ModelTimeSeries(sillfile).Time,(ModelTimeSeries(sillfile).VAF-ModelTimeSeries(sillfile).VAF(1))/1e9,'Color','g');
    
    % Fix up axes and title:
    xlim([0,300])
    xlabel('Time (yr)')
    ylim([min([ModelTimeSeries(constforcingfile).VAF-ModelTimeSeries(constforcingfile).VAF(1);ModelTimeSeries(warmingfile).VAF-ModelTimeSeries(warmingfile).VAF(1);ModelTimeSeries(sillfile).VAF-ModelTimeSeries(sillfile).VAF(1)]),...
        max([ModelTimeSeries(constforcingfile).VAF-ModelTimeSeries(constforcingfile).VAF(1);ModelTimeSeries(warmingfile).VAF-ModelTimeSeries(warmingfile).VAF(1);ModelTimeSeries(sillfile).VAF-ModelTimeSeries(sillfile).VAF(1)])]/1e9)
    ylabel('Change in Volume Above Flotation (km^3)')
    title('Volume Above Flotation')
    
    % Create legend:
    ylims=get(gca,'Ylim');
    if abs(ylims(1))>abs(ylims(2))
        legend([h1;h2;h3],'Constant Forcing','Warming','Warming+Sill','location','SouthWest')
    else
        legend([h1;h2;h3],'Constant Forcing','Warming','Warming+Sill','location','NorthWest')
    end
    
    % Plot sill construction time:
    patch('Vertices',[[100;100;110;110],[ylims(1);ylims(2);ylims(2);ylims(1)]],'Faces',1:4,'FaceVertexCData',[.5,.5,.5],'EdgeColor','none','FaceColor','flat')
    
    % Save figure:
    set(gcf,'PaperSize',[7,10])
    set(gcf,'PaperPosition',[0,0,7,10])
    figname=[figfolder,glaciername,'_ConstForcing_Warming_Sill_v4_vaf_gl.png'];
    print('-dpng',figname,'-r300')
    
end







