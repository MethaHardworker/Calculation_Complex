function r = asin(a)
%ASIN  Taylor model inverse sine asin(a)
%
%   r = asin(a)

% written  09/10/15     F. Buenger
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

r = stdfun(a,@asin,'asin');

end % function asin
