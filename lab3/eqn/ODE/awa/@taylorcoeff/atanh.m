function r = atanh(a)
%ATANH  hyperbolic arctangent of Taylor coefficients 
%                    
%   r = atanh(a)  
%
% Compare with Lohner's global function ARTANH, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = 1-sqr(a);
r = stdfun2(a,b,@atanh,'atanh');

end % function atanh