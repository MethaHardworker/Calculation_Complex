function r = coth(a)
%COTH  Taylor model hyperbolic cotangent coth(a)
%
%   r = coth(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@coth,'coth');
else
    r = cosh(a)./sinh(a);
end

end % function coth
