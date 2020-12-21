function r = acoth(a)
%ACOTH  hyperbolic arccotangent of Taylor coefficients 
%                    
%   r = acoth(a)  
%
% Compare with Lohner's global function ARCOTH, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = 1-sqr(a); 
r = stdfun2(a,b,@acoth,'acoth');

end % function acoth