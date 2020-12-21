function c =  iv_in(a,b)
%IV_IN  interval structure inclusion check "a contained in b ?"
%
%   c = iv_in(a,b)

% written  02/10/16     F. Buenger

if isfloat(a)
    c = ( (b.inf <= a) & (a <= b.sup) );    
else
    c = ( (b.inf <= a.inf) & (a.sup <= b.sup) );
end

end % function iv_in
