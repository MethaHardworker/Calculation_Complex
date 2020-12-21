function r = acsc(a)
%ACSC  inverse cosecans of Taylor coefficients
%
%   r = acsc(a)

% written  08/02/17     F. Buenger

r = intval('pi')/2 - asec(a); 

end % function acsc 
  