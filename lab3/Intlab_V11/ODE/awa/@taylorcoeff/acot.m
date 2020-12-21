function r = acot(a)
%ACOT  arccotangent of Taylor coefficients 
%                    
%   r = acot(a)  
%
% Compare with Lohner's global function ARCCOT, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = 1+sqr(a);
r = stdfun2(a,b,@acot,'acot');

end % function acot