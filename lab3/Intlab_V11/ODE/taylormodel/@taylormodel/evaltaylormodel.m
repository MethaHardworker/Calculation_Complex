function r = evaltaylormodel(a,x)
% EVALTAYLORMODEL  computes verified enclosures of (image + error)         
%                  of the single Taylor model "a" on all subdomains 
%                  of "a" given in the rows of x. 
%
%   r = evaltaylormodel(a,x)
%                       
% Speaking roughly, the function computes res(i) := a(x(i,:)-a.center') + a.interval
% for all row indices i of x.  

% written  12/02/15     F. Buenger
% modified 02/10/16     F. Buenger  "intval"-components --> intval-like structures 
% modified 05/17/18     F. Buenger  redesign for performance increase 
 
e = 1e-30;

if 1+e > 1 % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

if numel(a) ~= 1
    error('Only single Taylor models can be evaluated.');
end

M = a.monomial;
m_M = size(M,1);
m = max(max(M));

xs = iv_minus(x,a.center'); % centering of the evaluation points xs := x - a.center'
if any(isnan(xs.inf(:))) || any(isnan(xs.sup(:))) % Later on we set xs^0 := 1. Thus we are taking care of NaN at this place.
    error('NaN occured!');
end
[m_x,n_x] = size(xs.inf);
abs_xs_inf = abs(xs.inf);
abs_xs_sup = abs(xs.sup);

idx_pos = repmat(xs.inf >= 0,1,1,m+1); % indices of nonnegative interval components

idx_neg = repmat(xs.sup <= 0,1,1,m+1); % indices of negative interval components
idx_neg_even = idx_neg;
idx_neg_even(:,:,2:2:m+1) = false;
idx_neg_odd = idx_neg;
idx_neg_odd(:,:,1:2:m+1) = false;

idx_0 = repmat(xs.inf < 0 & xs.sup > 0,1,1,m+1); % indices of interval components containing zero
idx_0_even = idx_0;
idx_0_even(:,:,2:2:m+1) = false;
idx_0_odd = idx_0;
idx_0_odd(:,:,1:2:m+1) = false;

a_up = ones(m_x,n_x,m+1); % some initialization
b_up = a_up;
a_down = a_up;
b_down = a_up;
r_inf = a_up;
r_sup = a_up;

u = repmat([abs_xs_inf,abs_xs_sup],1,1,m);

setround(-1)
v = cumprod(u,3);
a_down(:,:,2:m+1) = v(:,1:n_x,:);      % lower bound of |xs_inf|.^N;
b_down(:,:,2:m+1) = v(:,n_x+1:end,:);  % lower bound of |xs_sup|.^N;

setround(1)
v = cumprod(u,3);
a_up(:,:,2:m+1) = v(:,1:n_x,:);        % upper bound of |xs_inf|.^N;
b_up(:,:,2:m+1) = v(:,n_x+1:end,:);    % upper bound of |xs_sup|.^N;

% case 1: nonnegative interval [a,b]^k = [a^k,b^k], a,b>=0
r_inf(idx_pos) = a_down(idx_pos);
r_sup(idx_pos) = b_up(idx_pos);

% case 2: negative interval, even exponent [-a,-b]^(2k) = [b^(2k),a^(2k)] , a,b>=0
r_inf(idx_neg_even) = b_down(idx_neg_even);
r_sup(idx_neg_even) = a_up(idx_neg_even);

% case 3: negative interval odd exponent [-a,-b]^(2k+1) = [-a^(2k+1),-b^(2k+1)] , a,b>=0
r_inf(idx_neg_odd) = -a_up(idx_neg_odd);
r_sup(idx_neg_odd) = -b_down(idx_neg_odd);

% case 4: interval containing zero, even exponent  [-a,b]^(2k) = [0,max(a,b)^(2k)] , a,b>=0, k > 0
r_inf(idx_0_even) = 0; % 0^2k = 0, k > 0
r_inf(:,:,1) = 1; % xs^0 = 1, k = 0
r_sup(idx_0_even) = max(a_up(idx_0_even),b_up(idx_0_even),'includenan');

% case 5: interval containing zero, odd exponent: [-a,b]^(2k+1) = [-a^(2k+1),b^(2k+1)] , a,b>=0
r_inf(idx_0_odd) = -a_up(idx_0_odd);
r_sup(idx_0_odd) = b_up(idx_0_odd);

% exponent/index shift of M
%   Recall that the linear index of (i,j,k)-entry of an (m_x,n_x,m+1)-cube is (k-1)*(m_x*n_x) + (j-1)*n_x + i. 
%   The indices k represent exponents k-1 stored in M.

M2 = repmat(permute(M*(m_x*n_x),[3 2 1]),m_x,1); % Create m_x horizontal copies of M*(m_x*n_x). 
I = (1:m_x)' + (0:n_x-1)*m_x; 
M_ = M2 + I ;

r_inf = r_inf(M_);
r_sup = r_sup(M_);

r_inf = permute(r_inf,[3 2 1]);
r_sup = permute(r_sup,[3 2 1]);

c = repmat(a.coefficient,1,1,m_x); % c is duplicated m_x times.
r.inf = [r_inf,c];                  
r.sup = [r_sup,c];
r = iv_prod(r,2);                  % Build all products p(i,k) := c(i) * xs(k,:)^M(i,1) * ... * xs(k,n)^M(i,n), i = 1,...,len_c, k = 1,...,len_x, n := size(M,2). 

r.inf = reshape(r.inf,m_M,m_x);
r.sup = reshape(r.sup,m_M,m_x);
r = iv_transpose(iv_sum(r,1));     % Build sums s(k) := p(1,k) +...+ p(len_c,k), k = 1,...,len_x. 

r = iv_plus(r,a.interval);         % Add error interval.

if rndold ~= 1
    setround(rndold)
end

end % function evaltaylormodel