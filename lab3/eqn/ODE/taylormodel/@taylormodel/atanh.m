function r = atanh(a)
%ATANH  Taylor model inverse hyperbolic tangent atanh(a)
%
%   r = atanh(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@atanh,'atanh');
else
    r = 1/2 .* log((1+a)./(1-a));
end

end % function atanh
