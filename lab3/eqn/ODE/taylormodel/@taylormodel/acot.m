function r = acot(a)
%ACOT  Taylor model inverse cotangent acot(a)
%
%   r = acot(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@acot,'acot');
else
    r = intval('pi')/2 - atan(a); 
end
end % function acot
