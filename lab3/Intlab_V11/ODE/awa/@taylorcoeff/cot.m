function r = cot(a)
%COT  cotangent of Taylor coefficients 
%                    
%   r = cot(a)  
%
% Compare with Lohner's global function COT, file fknoten.p. 

% written  08/01/17     F. Buenger

b = sqr(sin(a));
r = stdfun2(a,b,@cot,'cot');

end % function cot