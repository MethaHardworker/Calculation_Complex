function c =  iv_sum(a,dim)
%IV_SUM  interval structure sum
%
%   c = iv_sum(a,dim)

% written  02/10/16     F. Buenger
% modified 02/10/16     F. Buenger    performance improvement

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

if nargin == 1
    % Recall that rounding is upwards.
    c.inf = -sum(-a.inf);
    c.sup = sum(a.sup);
else
    % Recall that rounding is upwards.
    c.inf = -sum(-a.inf,dim);
    c.sup = sum(a.sup,dim);
end

if rndold ~= 1
    setround(rndold)
end

end % function iv_sum
