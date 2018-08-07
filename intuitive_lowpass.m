function dataout=intuitive_lowpass(datain,minwavelength,varargin)
% Function by Mike Wolovick, updated on 10/18/2011

% This function lowpass filters an input in 1D by convolving it with a 
% gaussian weighting function.  

% The standard deviation of the gaussian is one-half the minwavelength. ie,
% full width at half max (FWHM) is roughly minwavelength.

% Minwavelength must be in terms of number of samples of datain. 

% A linear best-fit is removed before the data are filtered, then added
% back in at the end. 

% Edge effects are minimized by adding a buffer to either side of the input
% equal to the value at the endpoints.  This is done after the best-fit is
% removed.  IE, the boundary condition is an extrapolation of the endpoint
% values based on the best-fit slope of the entire input. 

% The function assumes that the input is evenly spaced.

% MOD on 10/18/2011:
% The function now offers an option to control how BC's are dealt with.
% The default is the old method described above, the alternate is an
% extrapolation of the endpoints based on the local slope.  A global
% bestfit is removed from the data in either case; this only affects the
% extrapolation method.

% Syntax:
% output=intuitive_lowpass(input,minwavelength);
% or
% output=intuitive_lowpass(input,minwavelength,'method');  'method' can be
% either 'global' (default) or 'local'

rowvec=0;
if isscalar(datain)==1
    error('Data Input must be a vector or a matrix.')
elseif size(minwavelength,1)~=1 || size(minwavelength,2)~=1
    error('Minimum Wavelength must be a scalar.')
elseif length(varargin)>1
    error('Function cannot have more than 3 inputs.')
elseif size(datain,1)==1
    rowvec=1;
    datain=datain';
end
if isempty(varargin)==0
    if strcmp(varargin,'local')==0 && strcmp(varargin,'global')==0
        error('method must be either "global" or "local".')
    elseif strcmp(varargin,'local')==1
        bcmethod='local';
    else
        bcmethod='global';
    end
else
    bcmethod='global';
end
L=length(datain);
integer_wavelength=max(1,round(minwavelength));
stdev=.5*minwavelength;
Lbig=L+6*integer_wavelength;
xs=linspace(1,Lbig,Lbig)';
xcent=Lbig/2;
A=[ones(L,1),linspace(1,L,L)'];
equation=(A'*A)\A'*datain;
bestline=equation(1)+equation(2)*linspace(1,L,L)';
datain=datain-bestline;
dataspread=zeros(Lbig,1);
dataspread(3*integer_wavelength+1:end-3*integer_wavelength)=datain;
if strcmp(bcmethod,'global')==1
    dataspread(1:3*integer_wavelength)=datain(1);
    dataspread(end-3*integer_wavelength+1:end)=datain(end);
else
    d_start=datain(1)-datain(2);
    dataspread(1:3*integer_wavelength)=datain(1)+d_start*linspace(3*integer_wavelength,1,3*integer_wavelength)';
    d_end=datain(end)-datain(end-1);
    dataspread(end-3*integer_wavelength+1:end)=datain(end)+d_end*linspace(1,3*integer_wavelength,3*integer_wavelength)';
end
filter=exp(-.5*((xs-xcent)/stdev).^2);
filter=filter/sum(filter);
filtertrans=fft(filter);
datatrans=fft(dataspread);
dataout=fftshift(ifft(filtertrans.*datatrans));
dataout=dataout(3*integer_wavelength:end-3*integer_wavelength-1);
dataout=dataout+bestline;
if rowvec==1
    dataout=dataout';
end
