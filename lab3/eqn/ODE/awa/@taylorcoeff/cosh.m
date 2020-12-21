function r = cosh(a)
%COSH  hyperbolic cosine of Taylor coefficients 
%                    
%   r = cosh(a)  

% written  08/01/17     F. Buenger

r = stdfun1(a,@cosh,@sinh,'cosh','sinh');

end % function cosh