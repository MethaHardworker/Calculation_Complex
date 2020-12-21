 function r = sdotprod(c,a,Mt)
%SDOTPROD     direct, fast implementation of the dotproduct c'*(a.*(t.^Mt))
%
%   r = sdotprod(a)  
% 
%   c: double column vector 
%   a: Taylor model column vector of same length as "c". 
%      It is assumed that all Taylor models a(i) have the same type, order, 
%      domain and center. 
%      This is not checked again in order to increase performance.
%  Mt: vector of nonnegative integer exponents for time powers  t^k  

% written  05/19/17     F. Buenger

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards in verified mode
end

sparsity_tol = get_sparsity_tol;
n = length(c);
m = ones(n+1,1);
a_ = a(1);
d = a_.dim;
h = a_.domain.sup(d)-a_.domain.inf(d); % Upper bound for time interval length. Recall that rounding is upwards.

% At the end of the following for-loop, the error interval "err" shall contain an inclusion 
% of the sum of all errors c(i)*a(i).interval*t^k, where t is the time variable in [0,h] 
% with exponent k:= Mt_(i) is its exponent, i = 1,...,n.

err.inf = 0; 
err.sup = 0;
for i = 1:n
  a_ = a(i);  
  k = Mt(i);
  c_ = c(i);
  e_ = a_.interval;

  m(i+1) = m(i) + length(a_.coefficient); % Update the sum of coefficient lengths for later use.
  
  if k ~= 0  % Multiply a_ with t^k. If k = 0, nothing has to be done
    a_.monomial(:,d) = a_.monomial(:,d) + k; % direct multiplication a(i) := a(i)*t^k. Temporarily the new monomial order, which increased by k, 
                                             % is allowed to exceed the maximum order a_.order. Later, such higher order terms will be 
                                             % truncated and their image will be moved to the error interval.
    hk = h^k; % Upper bound for h^k. Recall that rounding is upwards.

    % Adapt error interval e_ := e_ * [0,h]^k    
    e_.sup =  max(0,e_.sup*hk);        % Recall that rounding is upwards.
    e_.inf =  min(0,-((-e_.inf)*hk));
    a_.interval = e_;
  
    % Adapt image im_ := im_ * [0,h]^k
    im_ = a_.image;
    im_.sup =  max(0,im_.sup*hk);      % Recall that rounding is upwards.
    im_.inf =  min(0,-((-im_.inf)*hk));
    a_.image = im_;
    a(i) = a_;
  end  
    
  % err = iv_plus(err,iv_times(c(i),a_.interval)); % err = err + c(i) * a_.interval  
  
  err.sup = err.sup +  max(c_.*e_.inf,c_.*e_.sup);  % Recall that rounding is upwards.  
  err.inf = - ( (-err.inf) - min(-((-c_).*e_.inf),-((-c_).*e_.sup)) );     
end

l = m(n+1)-1;         % total number of monomials of a
M = zeros(l,d);       % Initialize matrix for all monomials of a.
c_lower = zeros(l,1); % Initialize vector of lower bounds for all coefficients of c.*a.
c_upper = zeros(l,1); % Initialize vector of upper bounds for all coefficients of c.*a.
for i = 1:n
  a_ = a(i);
  M(m(i):m(i+1)-1,:) = a_.monomial; % M(m(i),:) becomes the first and M(m(i+1)-1,:) the last monomial of a(i).monomial.
  % Recall that rounding is switched to upwards. 
  a_coeff = a_.coefficient;
  c_lower(m(i):m(i+1)-1) = -( (-c(i)) * a_coeff ); 
  c_upper(m(i):m(i+1)-1) = c(i) * a_coeff;          
end

r = a(1); % Initialize result r.
%[U,iM,iU] = unique(M,'rows'); % Delete double entries and get corresponding indices: U(iU) = M, M(iM) = U.
[U,iU] = unique_rows(M); % Delete double entries and get corresponding indices: U(iU) = M.

% Recall that rounding is switched to upwards.
c_lower = -accumarray(iU,-c_lower); % Accumulate coefficient vector c, by summing up entries according to index iU,
                                    % that is summing up all coefficients having the same monomial, i.e.,
                                    % the same index in iU, and store that sum at this index in c_lower,
                                    % see MATLAB-documentation of the function accumarray.
c_upper = accumarray(iU,c_upper);   % Same accumulation as before for upper interval bound.
c_mid = 0.5 * (c_lower + c_upper);  % Estimate for interval midpoint. Floating-point operations are much faster than

% Find row indices i for which |c_mid(i)| >= sparsity_tol and degree(r.monomial(i)) <= r_.order. 
% The contribution of the terms of all other row indices will be moved to the error interval. See "rest" below.
row = and( (abs(c_mid) >= sparsity_tol) , (sum(U,2) <= r.order) ); 
if ~any(row)
    r.monomial = zeros(1,d);
    r.coefficient = 0;
else
    r.coefficient = c_mid(row) ; % Store the corresponding coefficients as those of the result r=a*b.
    r.monomial = U(row,:);       % Store also the corresponding monomials.
    
    % Recall that rounding is switched to upwards. 
    c_lower(row) = -(r.coefficient-c_lower(row)); % Compute centered lower bound for non-sparse coefficients.
    c_upper(row) = c_upper(row)-r.coefficient;    % Compute centered upper bound for non-sparse coefficients.
end
rest = r;
rest.monomial = U;
coeff.inf = c_lower;
coeff.sup = c_upper;
rest.coefficient = coeff;
rest.order = 2*rest.order; % maximum degree of higher order terms 
rest.image = image(rest);
r.interval = iv_plus(rest.image,err); % r.interval = rest.image + err  
r.image = image(r);

if rndold ~= 1
    setround(rndold)
end

end % function sdotprod

