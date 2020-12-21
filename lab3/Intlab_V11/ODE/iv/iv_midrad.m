function c = iv_midrad(a,r)
%IV_MIDRAD  create interval structure "c" for interval I:=[a-r,a+r] 
%           such that I is contained in "c".
%
%   c = iv_midrad(a,r) 

% written  02/11/16     F. Buenger

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

% Recall that rounding is upwards
c.inf = -((-a) + r);
c.sup = a + r;

if rndold ~= 1 
    setround(rndold)
end

end % function iv_midrad
