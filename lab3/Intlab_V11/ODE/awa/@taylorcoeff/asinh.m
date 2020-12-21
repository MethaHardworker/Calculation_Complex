function r = asinh(a)
%ASINH  hyperbolic arcsine of Taylor coefficients 
%                    
%   r = asinh(a)  
%
% Compare with Lohner's global function ARSINH, file fknoten.p. 
 
% written  08/01/17     F. Buenger

b = sqrt(sqr(a)+1);
r = stdfun2(a,b,@asinh,'asinh');

end % function asinh