function r = acos(a)
%ACOS  Taylor model inverse cosine acos(a)
%
%   r = acos(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

    r = stdfun(a,@acos,'acos');

end % function acos