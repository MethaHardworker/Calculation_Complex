function r = tcoeff2tc(a)
%TCOEFF2TC  converts generalized Taylor coefficient to tcoeff-structure
%
%    r = tcoeff2tc(a)

% written  12/22/17     F. Buenger

if isempty(a)
    r = [];
else
    k = length(a(1).inf);
    if k
        [m,n] = size(a);
        r.inf = reshape([a.inf].',m,n,k); 
        r.sup = reshape([a.sup].',m,n,k);
    else
        r = [];
    end
end

end % function tcoeff2tc


