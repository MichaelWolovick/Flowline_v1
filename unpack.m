function unpack(struct)

% unpack

% Mike Wolovick, 12/17/2011

% This function creates variables in the workspace out of the fields of a
% structure.  The structure must have only one element, but it can have
% many fields.



% Check input:
if isstruct(struct)==0
    error('Input must be a structure')
end
if length(struct)~=1
    error('Input must have a length of one')
end

% Unpack input:
fields=fieldnames(struct);
for ii=1:length(fields)
    eval(['data=struct.',cell2mat(fields(ii)),';'])
    assignin('caller',cell2mat(fields(ii)),data);
end



