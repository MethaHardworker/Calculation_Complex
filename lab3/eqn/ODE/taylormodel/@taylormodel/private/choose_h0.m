function h0 = choose_h0(odefun,t0,y0,order,loc_err_tol)
%CHOOSE_H0  choose inital step size for solving an ODE
%
%  r = choose(odefun,t0,y0)
%
% The algorithm is taken from  
%   [HNW] Hairer,Norsett,Wanner, "Solving Ordinary Differential Equations I", 
%         p. 169.

% written  07/07/18     F. Buenger

%type = 2;    % Euclidean norm
type = inf; % max-norm

if type == 2
    n = length(y0.inf);
    c = 1/sqrt(n);
else
    c = 1;
end
c = c/loc_err_tol;

y0_ = iv2intval(y0);
d0 = c*sup(norm(y0,type));
f0 = odefun(t0,y0_);
d1 = c*sup(norm(f0,type));

if d0 < 1e-5 || d1 < 1e-5
    h0 = 1e-6;
else
    h0 = 0.01*(d0/d1);    
end

y1 = y0_ + h0*odefun(t0,y0_); % Euler step
f1 = odefun(t0+h0,y1);
d2 = c*sup(norm(f1-f0,type))/h0;      % estimate for second derivative of the solution.

if max(d1,d2) <= 1e-15
    h1 = max(1e-6,h0*1e-3);
else
    h1 = (0.01/max(d1,d2)).^(1/(order+1));
end

h0 = min(100*h0,h1);

end % function choose_h0

