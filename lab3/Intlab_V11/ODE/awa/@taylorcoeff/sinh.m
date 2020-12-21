function r = sinh(a)
%SINH  sinus hyperbolicus of Taylor coefficients 
%                    
%   r = sinh(a)  

% written  08/01/17     F. Buenger

r = stdfun1(a,@sinh,@cosh,'sinh','cosh');

end % function sinh