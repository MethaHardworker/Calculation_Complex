function r = cos(a)
%COS  cosine of Taylor coefficients 
%                    
%   r = cos(a)  

% written  08/01/17     F. Buenger

r = stdfun1(a,@cos,@sin,'cos','sin');

end % function cos