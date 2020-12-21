function r = asec(a)
%ASEC  Taylor model inverse secans asec(a)
%
%   r = asec(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@asec,'asec');
else
    r = acos(1./a);
end

end % function asec
  