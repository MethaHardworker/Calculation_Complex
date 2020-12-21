function dydt = van_der_pol(t,y,i)
%VAN_DER_POL     ODE-function for Van der Pol's - differential equation

% written  07/26/18     F. Buenger
if nargin < 3 || isempty(i)
    dydt = [y(2);(1-sqr(y(1))).*y(2)-y(1)];
else
    switch i
        case 1
            dydt = y(2);
        case 2
           dydt = (1-sqr(y(1))).*y(2)-y(1); 
    end
end

end % function van_der_pol

