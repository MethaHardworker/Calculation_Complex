function r = coth(a)
%COTH  hyperbolic cotangent of Taylor coefficients 
%                    
%   r = coth(a)  
%
% Compare with Lohner's global function COTH, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = sqr(sinh(a));
r = stdfun2(a,b,@coth,'coth');

end % function coth