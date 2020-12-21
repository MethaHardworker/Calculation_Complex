function r = tie(a,b)
%Tie     concatenation a o b of "left" and "right" Taylor models a and b
%
%   r = concat(a,b)  
% 
% It is assumed that the input arguments "a" and "b" come from a function call 
%
%   [T,a,b] = verifyode(odefun,tspan,y0,options)
% 
% with preconditioning switched on, see the documentation of the function verifyode.
% Then, the output argument r is the concatenation of the left Taylor models "a"
% and the right Taylor models "b". This means that 
%
%   r(i,j), i = 1,...,k := length(T)-1, j = 1,...,n := "dimension of the ODE" 
%
% is the Taylor model representing the enclosure of the j-th solution component 
% y_j of the ODE on the time interval [T(i),T(i+1)].
%
% If only a single solution component is of interest, say j = 1, then it suffices to call 
%  
%   r = tie(a(:,1),b)
% 
% Note carefully, that still b must be transferred and not(!) just b(:,1).

% written  07/16/18     F. Buenger

r = a; 
if ~isempty(b)
    for i = 1:size(r,1)
        r(i,:) = concatenate(a(i,:),b(i,:));
    end
end

end % function tie
