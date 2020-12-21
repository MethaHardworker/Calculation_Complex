function r = asinh(a)
%ASINH  Taylor model inverse hyperbolic sine asinh(a)
%
%   r = asinh(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@asinh,'asinh');
else
    r = log(a+sqrt(sqr(a)+1)); % asinh(a) = log(a+sqrt(a^2+1))
end

end % function asinh
