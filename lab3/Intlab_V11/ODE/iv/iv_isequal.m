function r =  iv_isequal(a,b)
%IV_ISEQUAL  isequal for interval structure
%
%   s = iv_isequal(a)

% written  05/04/16     F. Buenger

r = and(isequal(a.inf,b.inf),isequal(a.sup,b.sup));

end % function iv_isequal
