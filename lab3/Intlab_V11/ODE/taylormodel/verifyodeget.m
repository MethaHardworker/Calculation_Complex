function res = verifyodeget(options, name, default)
%VERIFYODEGET   returns the value of the component "name" of the structure "options", 
%               analogous to MATLAB function odeget. 
%               res = verifyodeget(options,'name') % returns option.(name) if name is a valid field name. Otherwise [] is returned. 
%               res = verifyodeget(options,'name',default) % returns 'default' if the named property is not specified in options.
% 
%   res = verifyodeget(options, name, default)

% written  11/04/15     F. Buenger

if isfield(options,name)
    res = options.(name);
elseif nargin == 3 
    res = default;
else
    res = [];
end

end % function verifyodeget