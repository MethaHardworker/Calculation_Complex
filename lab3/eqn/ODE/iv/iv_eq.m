function c =  iv_eq(a,b)
%IV_EQ  interval structure equality check a == b
%
%   c = iv_eq(a,b)

% written  02/09/16     F. Buenger

c = ( all(size(a.inf) == size(b.inf)) && all(a.inf(:) == b.inf(:)) && all(a.sup(:) == b.sup(:)) );

end % function iv_eq
