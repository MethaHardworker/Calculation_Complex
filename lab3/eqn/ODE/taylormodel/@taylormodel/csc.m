function r = csc(a)
%SEC  Taylor model cosecans csc(a)
%
% r = csc(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@csc,'csc');
else
    r = inv(sin(a)); % csc(a) = 1/sin(a)
end

end % function csc