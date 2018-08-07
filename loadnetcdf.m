function loadnetcdf(inputfile)

% loadnetcdf

% Mike Wolovick, 12/14/2011

% This function loads a netcdf file as if you were calling matlab's usual
% "load" command.  It is based on the more specialized script
% "getALBMAPv1.m" which I wrote to load the albmap dataset.

% Modification, 7/16/2013:
% The function automatically replaces "-" in variable names with "_".  This
% was a problem loading the Bamber DEM for Greenland, since "-" is not
% allowed in matlab variable names (matlab interprets it as a minus sign).
% In addition, the function replaces -9999 as the nodata value with NaN.

% Modification, 6/16/2016: The function puts the variables into the
% workspace of the function that called it, rather than the base
% workspace.  


% Open netcdf file and get number of variables:
ncid=netcdf.open(inputfile,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid]=netcdf.inq(ncid);

% Loop through variables::
for ii=1:nvars
    
    % Get information about the variable:
    [varname,xtype,dimids,natts]=netcdf.inqVar(ncid,ii-1);
    
    % Create permute vector:
    if length(dimids)<3
        permutevector=[2,1];
    else
        permutevector=[2,1,linspace(3,length(dimids),length(dimids)-2)];
    end
    
    % Load the variable:
    eval(['data=','double(permute(netcdf.getVar(ncid,',num2str(ii),'-1),[',num2str(permutevector),']));'])
    
    % Replace "-" in the variable name with "_"
    dashnums=strfind(varname,'-');
    if isempty(dashnums)==0
        varname(dashnums)='_';
    end
    
    % Replace -9999 as the nodata value with NaN:
    data(data==-9999)=NaN;
    
    % Put the variable in the base workspace:
    assignin('caller',varname,data)
    
end

% Close the netcdf file:
netcdf.close(ncid);
