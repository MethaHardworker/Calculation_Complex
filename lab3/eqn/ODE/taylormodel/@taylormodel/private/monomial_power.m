function res =  monomial_power(x,M)
%MONOMIAL_POWER   nonnegative integer power a.^M of interval row vector x
%                 and nonnegative integer Matrix M. The length of x
%                 must agree with the number of columns of M.
%
%   res = monomial_power(x,M)

% written  05/02/16     F. Buenger

e = 1e-30;

if 1+e > 1 % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

m = max(max(M));
N = repmat((0:m)',1,size(M,2));

idx_pos = repmat((x.inf >= 0),m+1,1); % indices of nonnegative interval components

idx_neg = repmat((x.sup <= 0),m+1,1); % indices of nonpositive interval components
idx_neg_even = idx_neg;
idx_neg_even([1,2:2:m+1],:) = false;
idx_neg_odd = idx_neg;
idx_neg_odd(1:2:m+1,:) = false;

idx_0 = repmat(and(x.inf < 0,x.sup > 0),m+1,1); % indices of interval components containing zero
idx_0_even = idx_0;
idx_0_even([1,2:2:m+1],:) = false;
idx_0_odd = idx_0;
idx_0_odd(1:2:m+1,:) = false;

a_up = abs(x.inf).^N;
b_up = abs(x.sup).^N;

setround(-1)
a_down = abs(x.inf).^N;
b_down = abs(x.sup).^N;

setround(1)

r.inf = ones(size(N));
r.sup = r.inf;

% case [a,b]^k = [a^k,b^k], a,b>=0 , k>=0
r.inf(idx_pos) = a_down(idx_pos);
r.sup(idx_pos) = b_up(idx_pos);

% case [-a,-b]^(2k) = [b^(2k),a^(2k)] , a,b>=0 , k>=1
r.inf(idx_neg_even) = b_down(idx_neg_even);
r.sup(idx_neg_even) = a_up(idx_neg_even);

% case [-a,-b]^(2k+1) = [-a^(2k+1),-b^(2k+1)] , a,b>=0, k>=0
r.inf(idx_neg_odd) = -a_up(idx_neg_odd);
r.sup(idx_neg_odd) = -b_down(idx_neg_odd);

% case [-a,b]^(2k) = [0,max(a,b)^(2k)] , a,b>=0, k>=1
r.inf(idx_0_even) = 0;
r.sup(idx_0_even) = max(a_up(idx_0_even),b_up(idx_0_even));

% case [-a,b]^(2k+1) = [-a^(2k+1),b^(2k+1)] , a,b>=0, k>=0
r.inf(idx_0_odd) = -a_up(idx_0_odd);
r.sup(idx_0_odd) = b_up(idx_0_odd);

% exponent/index shift of M
I = M + 1 + repmat((m+1)*(0:size(M,2)-1),size(M,1),1);
res.inf = r.inf(I);
res.sup = r.sup(I);

if rndold ~= 1
    setround(rndold)
end

end % function monomial_power
