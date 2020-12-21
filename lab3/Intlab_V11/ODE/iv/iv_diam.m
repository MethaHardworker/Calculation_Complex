function d = iv_diam(a)
%IV_DIAM  upper bound for diameter of interval structure "a" 
%
%   c = iv_diam(a)

% written  02/10/16     F. Buenger

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

d =  a.sup-a.inf;

if rndold ~= 1 
    setround(rndold)
end

end % function iv_diam
