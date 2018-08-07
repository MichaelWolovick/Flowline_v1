% FixUngroundedExtrapolation_v1

% Mike Wolovick

% This script fixes the extrapolation of the sliding parameters to
% ungrounded regions from the flowband inversions.  I need to ensure that
% the extrapolation is not biased by erroneous short-wavelength gradients
% at the grounding line, and also that basal drag is always positive.

% The script starts by imposing a minimum drag to the grounded ice, then
% low-pass filtering the sliding stress scale (at the inversion 
% wavelength).  After that the extrapolation for stress is repeated, but 
% with a minimum imposed on the far-field stress.  

% Sliding velocity scale is treated in a similar way, but the smoothing is
% imposed in a log scale and no minimum is used.

% NOTE: as a result of these modifications, the model will not reproduce
% the velocity field that was compared with data in the inversion process.
% The original results can be recovered by using the drag coefficient
% variable, or by using the drag variable to produce a sliding stress
% scale.  (I suspect the smoother results may produce a more stable model,
% especially when the sliding exponent is very high).

clear all
tic

%% Parameters:

% File names and paths:
ungroundedfixerinputfolder='/home/mjw/Documents/FjordSillPaper/ModelInput/';
ungroundedfixerinputfiles={'HelheimA.mat';...
    'HelheimB.mat';...
    'HelheimC.mat';...
    'KangerA.mat';...
    'KangerB.mat';...
    'KangerC.mat';...
    'JakobshavnA.mat';...
    'JakobshavnB.mat';...
    'JakobshavnC.mat';...
    'PetermannA.mat';...
    'PetermannB.mat';...
    'PetermannC.mat'};

% Minimum basal drag:
minbasaldrag=1e3;         % Pa


%% Work:

% Loop through input files:
for thisfile=1:length(ungroundedfixerinputfiles)
    
    % Load input file:
    load([ungroundedfixerinputfolder,ungroundedfixerinputfiles{thisfile}])
    
    % Define "a" parameters:
    a_input=exp(FinalVector(1));
    a_side_input=exp(FinalVector(2));
    
    % Convert drag to characteristic stress scale:
    SlidingStressScale_input(2:lastgroundedind)=Drag_lrd(2:lastgroundedind);
    SlidingStressScale_input(1)=Drag_lrd(2);
    
    % Impose minimum sliding stress scale:
    SlidingStressScale_input=max(SlidingStressScale_input,minbasaldrag);
    
    % Low-pass filter stress scale:
    SlidingStressScale_input(1:lastgroundedind)=intuitive_lowpass(SlidingStressScale_input(1:lastgroundedind),inversionwavelength/steplength);
    
    % Extrapolate stress scale beyond the grounding line:
    glstress=SlidingStressScale_input(lastgroundedind);
    glstressgrad=(SlidingStressScale_input(lastgroundedind)-SlidingStressScale_input(lastgroundedind-1))/dx;
    farfieldstress=max(glstress+inversionwavelength*glstressgrad,minbasaldrag);
    SlidingStressScale_input(lastgroundedind+1:end)=farfieldstress+(glstress-farfieldstress)*exp(-(X_input(lastgroundedind+1:end)-X_input(lastgroundedind))/inversionwavelength);
    
    % Convert velocity to characteristic scale:
    SlidingVelocityScale_yr_input(1:lastgroundedind)=U_lr(1:lastgroundedind)*secondsperyear;
    
    % Low-pass filter velocity scale:
    SlidingVelocityScale_yr_input(1:lastgroundedind)=exp(intuitive_lowpass(log(SlidingVelocityScale_yr_input(1:lastgroundedind)),inversionwavelength/steplength));
    
    % Extrapolate velocity scale beyond the grounding line:
    glu_yr=SlidingVelocityScale_yr_input(lastgroundedind);
    glugrad_yr=(SlidingVelocityScale_yr_input(lastgroundedind)-SlidingVelocityScale_yr_input(lastgroundedind-1))/dx;
    farfieldu_yr=glu_yr+inversionwavelength*glugrad_yr;
    SlidingVelocityScale_yr_input(lastgroundedind+1:end)=farfieldu_yr+(glu_yr-farfieldu_yr)*exp(-(X_input(lastgroundedind+1:end)-X_input(lastgroundedind))/inversionwavelength);
    
    % Save results:
    save([ungroundedfixerinputfolder,ungroundedfixerinputfiles{thisfile}],'SlidingStressScale_input','SlidingVelocityScale_yr_input','-append')
    
end


% Final display:
disp('Done!')
toc