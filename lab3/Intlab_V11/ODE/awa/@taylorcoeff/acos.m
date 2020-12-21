function r = acos(a)
%ACOS  arccosine of Taylor coefficients 
%                    
%   r = acos(a)  
%
% Compare with Lohner's global function ARCCOS, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = sqrt(1-sqr(a));
r = stdfun2(a,b,@acos,'acos');

end % function acos