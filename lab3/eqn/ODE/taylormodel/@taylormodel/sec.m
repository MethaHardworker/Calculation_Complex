function r = sec(a)
%SEC  Taylor model secans sec(a)
%
%   r = sec(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@sec,'sec');
else
    r = inv(cos(a)); % sec(a) = 1/cos(a)
end

end % function sec