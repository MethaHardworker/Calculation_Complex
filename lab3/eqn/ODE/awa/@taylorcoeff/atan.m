function r = atan(a)
%ATAN  arctangent of Taylor coefficients 
%                    
%   r = atan(a)  
%
% Compare with Lohner's global function ARCTAN, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = 1+sqr(a);
r = stdfun2(a,b,@atan,'atan');

end % function atan