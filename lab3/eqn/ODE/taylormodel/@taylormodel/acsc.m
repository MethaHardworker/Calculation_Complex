function r = acsc(a)
%ACSC  Taylor model inverse cosecans acsc(a)
%
%   r = acsc(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@acsc,'acsc');
else
    r = intval('pi')/2 - asec(a); 
end

end % function acsc 
  