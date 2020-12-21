function r = acoth(a)
%ACOTH  Taylor model inverse hyperbolic cotangent acoth(a)
%
%   r = acoth(a

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@acoth,'acoth');
else
    r = 1/2 .* log((a+1)./(a-1));
end

end % function acoth

