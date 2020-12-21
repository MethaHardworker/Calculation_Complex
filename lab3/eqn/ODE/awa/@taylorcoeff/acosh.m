function r = acosh(a)
%ACOSH  hyperbolic arccosine of Taylor coefficients 
%                    
%   r = cosh(a)  
%
% Compare with Lohner's global function ARCOSH, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = sqrt(sqr(a)-1); 
r = stdfun2(a,b,@acosh,'acosh');

end % function acosh