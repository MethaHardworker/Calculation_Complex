function r = sin(a)
%SIN  sine of Taylor coefficients 
%                    
%   r = sin(a)  

% written  08/01/17     F. Buenger

r = stdfun1(a,@sin,@cos,'sin','cos');

end % function sin