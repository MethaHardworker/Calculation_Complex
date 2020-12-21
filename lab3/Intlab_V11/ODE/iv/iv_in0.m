function c =  iv_in0(a,b)
%IV_IN0  interval structure strict inclusion check:
%        "a contained in the interior of b ?"
%
%   c = iv_in0(a,b)

% written  02/10/16     F. Buenger

if isfloat(a)
    c = ( (b.inf < a) && (a < b.sup) );
else
    c = ( (b.inf < a.inf) && (a.sup < b.sup) );
end

end % function iv_in0
