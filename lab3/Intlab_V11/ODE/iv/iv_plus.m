function c = iv_plus(a,b)
%IV_PLUS  interval structure summation a + b
%
%   c = iv_plus(a,b)

% written  02/09/16     F. Buenger

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

if isfloat(a)
    x.inf = a;
    x.sup = a;
    a = x;
end
if isfloat(b)
    x.inf = b;
    x.sup = b;
    b = x;
end

c.inf = -( (-a.inf) - b.inf );
c.sup = a.sup + b.sup;

% s = [c.inf,c.sup];
% if any(isnan(s(:))) || any(isinf(s(:)))
%     error('NAN/INF occured.')
% end

if rndold ~= 1 
    setround(rndold)
end

end % function iv_plus
