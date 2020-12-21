function dydt = lorenz(t,y,i)
% LORENZ   ode-function for Lorenz attractor

% written  07/24/18   F. Buenger

sigma = 10;
rho = 28; 
if nargin == 2 || isempty(i) || i == 3
    if isfloat(y)
        beta = 8/3;
    else
        beta = iv2intval(iv_rdivide(8,3));
    end
end
if nargin == 2 || isempty(i)
    dydt = y;
    dydt(1) = sigma.*(y(2)-y(1));
    dydt(2) = (rho-y(3)).*y(1)-y(2);
    dydt(3) = y(1).*y(2)-beta.*y(3);
else
   switch i
       case 1
           dydt = sigma.*(y(2)-y(1));
       case 2
           dydt = (rho-y(3)).*y(1)-y(2);
       case 3
           dydt = y(1).*y(2)-beta.*y(3);
   end
end

end % function lorenz

