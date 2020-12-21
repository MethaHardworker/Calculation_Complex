function dydt = lotka_volterra(t,y,i)
%LOTKA_VOLTERRA   ode-function for Lotka-Volterra equation,  
%                 see [E], Chapter 5.1, p.146 et seq.  

% written  03/11/18   F. Buenger

if nargin == 2 || isempty(i)
    dydt = [2.*y(1).*(1-y(2));       % y1' = 2*y1*(1-y2)
            (y(1)-1).*y(2)];         % y2' = (y1-1)*y2
else
    switch i
        case 1                        % Compute only first component:
            dydt = 2.*y(1).*(1-y(2)); %   y1' = 2*y1*(1-y2)
        case 2                        % Compute only second component:
            dydt = (y(1)-1).*y(2);    %   y2' = (y1-1)*y2
    end
end

end % function lotka_volterra

