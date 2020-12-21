function r = acosh(a)
%ACOSH  Taylor model inverse hyperbolic cosine acosh(a)
%
%   r = acosh(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@acosh,'acosh');
else
    r = log(a+sqrt(sqr(a)-1)); % acosh(a) = log(a+sqrt(a^2-1))
end

end % function acosh