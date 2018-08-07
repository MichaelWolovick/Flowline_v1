% FlowbandInversionFigure_v1

% Mike Wolovick, 11/13/2017

% This script makes a figure showing the results of flowband inversion for
% the artificial sill paper.  It takes a lot of plotting commands from
% InvertFlowband_v2

clear all
close all
tic

%% Parameters

% File names:
inputfile='/home/wolovick/Dropbox/FjordSillPaper/ModelInput/ThwaitesC_v2.mat';
figname='/home/wolovick/Dropbox/FjordSillPaper/Figures/ThwaitesInversionFigure_v2.png';




%% Work

% Load input file:
thisfigname=figname;
load(inputfile)
figname=thisfigname;

% Make figure:
figure(1)

% Velocity:
subplot(2,2,1)
h11=semilogy(X_input(1:lasticeind+1)/1000,U_lr_population*secondsperyear,'Color',[.5,.5,.5]);
hold on
h12=semilogy(X_input(1:lasticeind+1)/1000,U_lr*secondsperyear,'k');
h13=semilogy(X_input(1:lasticeind+1)/1000,U_lr_data*secondsperyear,'r');
h14=semilogy(X_input/1000,SlidingVelocityScale_yr_input,'g');
plot([1,1]*X_input(lastgroundedind)/1000,[1,1e4],'--k')
legend([h11(1);h12;h13;h14],'Ensemble Members','Ensemble Mean','Observations','u_0 (model input)','location','SouthEast')
xlim([0,X_input(end)/1000])
ylim([1,1e4])
xlabel('Distance (km)')
ylabel('Velocity (m/yr)')
title('a) Ice Velocity, u')

% Strain rate:
subplot(2,2,2)
h21=semilogy(.5*(X_input(1:lasticeind)+X_input(2:lasticeind+1))/1000,abs(StrainRate_xx_c_population)*secondsperyear*1000,'Color',[.5,.5,.5]);
hold on
h22=semilogy(.5*(X_input(1:lasticeind)+X_input(2:lasticeind+1))/1000,abs(StrainRate_xx_c)*secondsperyear*1000,'k');
h23=semilogy(.5*(X_input(1:lasticeind)+X_input(2:lasticeind+1))/1000,abs(StrainRate_xx_c_data)*secondsperyear*1000,'r');
plot([1,1]*X_input(lastgroundedind)/1000,[1e-3,1e3],'--k')
legend([h21(1);h22;h23],'Ensemble Members','Ensemble Mean','Observations','location','NorthWest')
xlim([0,X_input(end)/1000])
ylim([1e-3,1e3])
xlabel('Distance (km)')
ylabel('Strain Rate (1/ka)')
title('b) Longitudinal Strain Rate, du/dx')

% Stress plot:
subplot(2,2,3)
h31=plot(X_input(2:lastgroundedind)/1000,Drag_lrd_population(2:lastgroundedind,:)/1000,'Color',[.5,.5,.5]);
hold on
h32=plot(X_input(1:lastgroundedind)/1000,Drag_lrd(1:lastgroundedind)/1000,'k');
plot(X_input(1:lasticeind+1)/1000,(2*SideDrag_lr_population.*repmat(Icethick_lr_center'./Width_lr',[1,populationsize]))/1000,'Color',[.75,.75,1])
h33=plot(X_input(1:lasticeind+1)/1000,(2*SideDrag_lr.*Icethick_lr_center./Width_lr)/1000,'b');
h34=plot(X_input(2:lasticeind+1)/1000,DrivingStress_lr(2:end)/1000,'r');
h35=plot(X_input/1000,SlidingStressScale_input/1000,'g');
plot([1,1]*X_input(lastgroundedind)/1000,[0,1e6],'--k')
legend([h32;h33;h34;h35],'Basal Drag','Side Drag (*2H/W)','Driving Stress','\tau_0 (model input)','location','NorthWest')
xlim([0,X_input(end)/1000])
ylim([0,1e5*ceil(max([max(DrivingStress_lr),max(max(Drag_lrd_population)),max(max(2*SideDrag_lr_population.*repmat(Icethick_lr_center'./Width_lr',[1,populationsize])))])/1e5)]/1000)
xlabel('Distance (km)')
ylabel('Stress (kPa)')
title('c) Force Balance')

% Drag coefficient plot:
subplot(2,2,4)
h41=semilogy(X_input(2:lastgroundedind)/1000,exp(PopulationVectors(3:end-1,:)),'Color',[.5,.5,.5]);
hold on
h42=semilogy(X_input(2:lastgroundedind)/1000,exp(FinalVector(3:end-1)),'k');
h43=semilogy(X_input(2:lastgroundedind)/1000,DragCoefficientGuess_lrd(2:lastgroundedind),'r');
ylims=get(gca,'Ylim');
plot([1,1]*X_input(lastgroundedind)/1000,ylims,'--k')
legend([h41(1);h42;h43],'Ensemble Members','Ensemble Mean','\tau_{d}/u_{obs}','location','NorthEast')
xlim([0,X_input(end)/1000])
xlabel('Distance (km)')
ylabel('Drag Coefficient (Pa/(m/s))')
title('d) Drag Coefficient, C')

% Save figure:
set(gcf,'PaperSize',[10,8])
set(gcf,'PaperPosition',[0,0,10,8])
print('-dpng',figname,'-r300')

% Final display:
disp('Done!')
toc

