function res =  iv2intval(a)
%IV2INTVAL  convert interval structure to intval
%
%   s = iv2intval(a)

% written  02/09/16     F. Buenger

res = intval(a.inf,a.sup,'infsup');

end % function iv2intval
