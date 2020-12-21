function r = tanh(a)
%ACOTH  Taylor model hyperbolic tangent tanh(a)
%
%   r = tanh(a)

% written  12/15/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@tanh,'tanh');
else
    r = sinh(a)./cosh(a);
    % r = 1 - 2./(exp(2.*a)+1); % alternative formula 
end

end % function tanh