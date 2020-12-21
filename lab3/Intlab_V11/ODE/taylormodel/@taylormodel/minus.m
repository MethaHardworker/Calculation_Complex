function r = minus(a,b)
%MINUS  Taylor model minus a - b
%
%   r = minus(a,b)

% written  08/27/15     F. Buenger

if isiv(b)
    r = a + iv_uminus(b);
else
    r = a + (-b);
end

end % function minus
