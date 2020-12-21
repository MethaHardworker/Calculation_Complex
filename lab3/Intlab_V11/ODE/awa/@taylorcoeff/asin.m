function r = asin(a)
%ASIN  arcsine of Taylor coefficients 
%                    
%   r = asin(a)  
%
% Compare with Lohner's global function ARCSIN, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = sqrt(1-sqr(a));
r = stdfun2(a,b,@asin,'asin');

end % function asin