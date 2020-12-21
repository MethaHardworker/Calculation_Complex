function r = iv_vecnorm(a,type)
%IV_VECNORM  norm of row or column vector a. 
%
%   r = iv_vecnorm(a,type) 
% 
% type  1    sum norm
%       2    Euclidean norm (default)
%       inf  max norm

% written  08/08/17     F. Buenger

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

if isfloat(a) && isreal(a)
    b.inf = a;
    b.sup = a;
    a = b;
end

if nargin < 2
    type = 2;
end

if isequal(type,'inf')
    type = inf; 
end

b = iv_abs(a);
r = b.sup(:); 

switch type
    case 1
        r = sum(r); % Recall that rounding is upwards.
    case 2
        r = sum(r.^2); % Recall that rounding is upwards.
        r = sup(sqrt(intval(r))); % r = sqrt(r)
    case inf
        r = max(r,[],'includenan');
    otherwise
        error('invalid type of norm')
end
if issparse(r)
    r = full(r);
end

if rndold ~= 1 
    setround(rndold)
end

end % function iv_vecnorm
