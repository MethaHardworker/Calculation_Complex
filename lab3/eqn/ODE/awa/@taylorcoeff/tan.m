function r = tan(a)
%TAN  tangens of Taylor coefficients 
%                    
%   r = tan(a)  
%
% Compare with Lohner's global function TAN, file fknoten.p. 

% written  08/01/17     F. Buenger

b = sqr(cos(a));
r = stdfun2(a,b,@tan,'tan');

end % function tan