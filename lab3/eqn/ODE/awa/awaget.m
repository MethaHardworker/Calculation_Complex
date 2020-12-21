function res = awaget(options, name, default)
%AWA_GET  returns the value of the component "name" of the structure "options", 
%         analogous to MATLAB function odeget. 
%         res = awaget(options,'name') % returns option.(name) if name is a valid fieldname. Otherwise [] is returned. 
%         res = awaget(options,'name',default) % returns 'default' if the named property is not specified in options.
% 
%   res = awa_get(options, name, default)

% written  07/22/17     F. Buenger
% modified 04/26/17     F. Buenger   lower case field names 'evalmeth','abstol','reltol' are also accepted.

switch name
    case 'evalmeth'
        name = 'EvalMeth';
    case 'abstol'
        name = 'AbsTol';
    case 'reltol'
        name = 'RelTol';
end

if isfield(options,name)
    res = options.(name);
elseif nargin == 3 
    res = default;
else
    res = [];
end

end % function awaget