function r = tanh(a)
%TANH  tangens hyperbolicus of Taylor coefficients 
%                    
%   r = tanh(a)  
%
% Compare with Lohner's global function TANH, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = sqr(cosh(a));
r = stdfun2(a,b,@tanh,'tanh');

end % function tanh