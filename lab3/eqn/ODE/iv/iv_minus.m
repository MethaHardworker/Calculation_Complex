function c =  iv_minus(a,b)
%IV_MINUS  interval structure subtraction a - b
%
%   c = iv_minus(a,b)

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

c.inf = -( b.sup - a.inf ); % Recall that rounding is upwards.
c.sup = a.sup - b.inf;

if rndold ~= 1 
    setround(rndold)
end

end % function iv_minus
