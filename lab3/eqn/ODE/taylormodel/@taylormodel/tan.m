function r = tan(a)
%TAN  Taylor model tangent tan(a)
%
%   r = tan(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@tan,'tan');
else
    r = sin(a)./cos(a);
end

end % function tan