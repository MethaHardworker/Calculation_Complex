function r = cot(a)
%COT  Taylor cotangent cot(a)
%
%   r = cot(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@cot,'cot');
else
    r = cos(a)./sin(a);
end

end % function cot