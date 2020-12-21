function r = atan(a)
%ATAN  Taylor model inverse tangent atan(a)
%
%   r = atan(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

r = stdfun(a,@atan,'atan');

end % function atan