function res =  iv_intpower(x,M)
%IV_INTPOWER  nonnegative integer power x.^M of float or interval matrix x
%             and nonnegative integer Matrix M. 
%             The implementation is adapted to the needs of evaluating 
%             multivariate polynomials p(x_1,...,x_n) with a small number n
%             of variables, say n <= 10, which appear in powers x_i^m_i 
%             with small exponents, say m_i <= 20. 
%             The number k := size(x,1) of points at at which p is evaluated
%             might be large. Also the number of monomials s := size(M,1)
%             might be large.
%
%   res = iv_intpower(x,M)

% written  05/02/16     F. Buenger
% modified 11/09/16     F. Buenger  x.^M for interval matrix x of same size as M.
% modified 08/01/17     F. Buenger  parameter 'includenan' in min-, max-calculations added
% modified 05/17/18     F. Buenger  float input x, size(x) and size(M) compatible.

e = 1e-30;

if 1+e > 1 % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

float_x = isfloat(x); 

if float_x
    x_inf = x;
    x_sup = x;
else
    x_inf = x.inf;
    x_sup = x.sup;
end

if length(size(x_inf)) ~= 2 || length(size(M)) ~= 2
    error('input arguments must be 2-dimensional.')
end

[m_x,n_x] = size(x_inf);
[m_M,n_M] = size(M);
m_r =  max(m_x,m_M);
n_r =  max(n_x,n_M);

if (m_x < m_r && m_x ~= 1) || (m_M < m_r && m_M ~= 1) || ...
   (n_x < n_r && n_x ~= 1) || (n_M < n_r && n_M ~= 1)

   error('incompatible size of input arguments.')
end
if n_x < n_M
    x_inf = repmat(x_inf,1,n_M);
    x_sup = repmat(x_sup,1,n_M);
    n_x = n_M;
end
if n_M < n_x
    M = repmat(M,1,n_x);
    n_M = n_x;
end
if m_M < m_x
    M = repmat(M,m_x,1);
    m_M = m_x;
end

if all(m_x == m_M && n_x == n_M)  
    
    %[r_inf,r_sup] = power1(x_inf,x_sup,M);   
    [r_inf,r_sup] = power1(x_inf(:),x_sup(:),M(:));   
    r_inf = reshape(r_inf,m_x,n_x);
    r_sup = reshape(r_sup,m_x,n_x);
    
else  % x is a row vector and M is a matrix of same row length.    
    
    abs_x_inf = abs(x_inf);
    abs_x_sup = abs(x_sup);
    
    m = max(max(M));
    
    idx_pos = repmat((x_inf >= 0),m+1,1); % indices of nonnegative interval components
    
    idx_neg = repmat((x_sup <= 0),m+1,1); % indices of negative interval components
    idx_neg_even = idx_neg;
    idx_neg_even(2:2:m+1,:) = false;
    idx_neg_odd = idx_neg;
    idx_neg_odd(1:2:m+1,:) = false;
    
    idx_0 = repmat(and(x_inf < 0,x_sup > 0),m+1,1); % indices of interval components containing zero
    idx_0_even = idx_0;
    idx_0_even(2:2:m+1,:) = false;
    idx_0_odd = idx_0;
    idx_0_odd(1:2:m+1,:) = false;
    
    w = ones(1,n_x);   
    if float_x
        u = repmat(abs_x_inf,m,1);
        a_up = [w;cumprod(u)];       % upper bound of |x|.^N;
        b_up = a_up;
        setround(-1)
        a_down = [w;cumprod(u)];     % lower bound of |x_inf|.^N;        
        b_down = a_down;
    else
        u = repmat([abs_x_inf,abs_x_sup],m,1);
        v = cumprod(u);
        a_up = [w;v(:,1:n_x)];       % upper bound of |x_inf|.^N;
        b_up = [w;v(:,n_x+1:end)];   % upper bound of |x_sup|.^N;
        setround(-1)
        v = cumprod(u);
        a_down = [w;v(:,1:n_x)];     % lower bound of |x_inf|.^N;
        b_down = [w;v(:,n_x+1:end)]; % lower bound of |x_sup|.^N;
    end    
    
    setround(1)
    
    r_inf = zeros(m+1,n_x); 
    r_sup = r_inf;
    
    % case [a,b]^k = [a^k,b^k], a,b>=0
    r_inf(idx_pos) = a_down(idx_pos);
    r_sup(idx_pos) = b_up(idx_pos);
    
    % case [-a,-b]^(2k) = [b^(2k),a^(2k)] , a,b>=0
    r_inf(idx_neg_even) = b_down(idx_neg_even);
    r_sup(idx_neg_even) = a_up(idx_neg_even);
    
    % case [-a,-b]^(2k+1) = [-a^(2k+1),-b^(2k+1)] , a,b>=0
    r_inf(idx_neg_odd) = -a_up(idx_neg_odd);
    r_sup(idx_neg_odd) = -b_down(idx_neg_odd);
    
    % case [-a,b]^(2k) = [0,max(a,b)^(2k)] , a,b>=0, k > 0    
    r_inf(idx_0_even) = 0; % 0^2k = 0, k > 0
    r_inf(1,:) = 1; % [-a,b]^0 = 1, k = 0   
    if float_x
        r_sup(idx_0_even) = a_up(idx_0_even);        
    else
        r_sup(idx_0_even) = max(a_up(idx_0_even),b_up(idx_0_even),'includenan');
    end
    
    % case [-a,b]^(2k+1) = [-a^(2k+1),b^(2k+1)] , a,b>=0
    r_inf(idx_0_odd) = -a_up(idx_0_odd);
    r_sup(idx_0_odd) = b_up(idx_0_odd);

    % exponent/index shift of M
    A = M + repmat((m+1)*(0:n_x-1),m_M,1) + 1;
    r_inf = r_inf(A);
    r_sup = r_sup(A);
    
end

res.inf = r_inf;
res.sup = r_sup;
    
if rndold ~= 1
    setround(rndold)
end

end % function iv_intpower

function [rinf,rsup] = power1(ainf,asup,b)
% interval a, integer b, size(a) = size(b)
% Rounding is already switched to upwards by caling function 'power'. 
  rinf = ones(size(ainf));
  rsup = rinf;
  b_0 = ( b == 0 );
  if any(b_0(:))
    rinf(b_0) = 1;
    rsup(b_0) = 1;
  end
  b_even = ( mod(b,2) == 0 ) & ~b_0;
  if any(b_even(:))
    a_0 = b_even & ( ainf < 0 ) & (asup > 0 );
    if any(a_0)                         % 0 in a  &  b even
      rinf(a_0) = 0;
      [lb,ub] = powerint(-ainf(a_0),asup(a_0),abs(b(a_0)),1,1);
      rsup(a_0) = max(lb,ub);
    end
    a_pos = b_even & (ainf >= 0 );
    if any(a_pos(:))                       % a >= 0  &  b even
      [rinf(a_pos) , rsup(a_pos)] = powerint(ainf(a_pos),asup(a_pos),abs(b(a_pos)),-1,1);
    end
    a_neg = b_even & (asup <= 0 );
    if any(a_neg(:))                       % a <= 0  &  b even
      [rinf(a_neg) , rsup(a_neg)] = powerint(-asup(a_neg),-ainf(a_neg),abs(b(a_neg)),-1,1);
    end
  end
  b_odd = ( mod(b,2) == 1 );
  if any(b_odd(:))
    a_pos = b_odd & ( ainf >= 0 );
    if any(a_pos)
      [rinf(a_pos) , rsup(a_pos)] = powerint(ainf(a_pos),asup(a_pos),abs(b(a_pos)),-1,1);
    end
    a_neg = b_odd & ( asup <= 0 );
    if any(a_neg)             % careful with rounding
      [rinf(a_neg) , rsup(a_neg)] = powerint(-ainf(a_neg),-asup(a_neg),abs(b(a_neg)),1,-1);
      rinf(a_neg) = -rinf(a_neg);
      rsup(a_neg) = -rsup(a_neg);
    end
    a_0 = b_odd & ~a_pos & ~a_neg;
    if any(a_0)               % careful with rounding
      [rinf(a_0) , rsup(a_0)] = powerint(-ainf(a_0),asup(a_0),abs(b(a_0)),1,1);
      rinf(a_0) = -rinf(a_0);
    end
  end
  
  index0 = ( ainf == 0 & asup == 0);
  if any(index0)
    b = intval(b);
    index = index0 & ( b > 0 );
    rinf(index) = 0;
    rsup(index) = 0;
    index = index0 & ( b == 0 );
    rinf(index) = 1;
    rsup(index) = 1;
  end
  
  index = isnan(ainf) | isnan(ainf) | isnan(b);
  if any(index(:))
    error('NaN occured');
  end
  
end  % function power1

function [r1,r2] = powerint(a1,a2,b,rnd1,rnd2)
% non-negative a1,a2, positive integer b, size(a2)=size(b), result ri=ai.^b with rounding rndi
% Rounding is already switched to upwards by caling functions 'power1' and 'power'. 
 r1 = a1;
 r2 = a2;
 b = b(:)-1;
 while any(b>0)
   index = ( mod(b,2) == 1 );
   if any(index)
     r1(index) = rnd1 * ( (rnd1 * r1(index)) .* a1(index) );  % equivalent to setround(rnd1); r1(index) = r1(index).*a1(index);
     r2(index) = rnd2 * ( (rnd2 *  r2(index)) .* a2(index) ); % equivalent to setround(rnd2); r2(index) = r2(index).*a2(index);   
   end
   b = floor(b/2);
   index = ( b~=0 );
   if any(index)
     a1(index) = rnd1 * ( (rnd1 * a1(index)) .* a1(index) );  % equivalent to setround(rnd1); a1(index) = a1(index).*a1(index);
     a2(index) = rnd2 * ( (rnd2 * a2(index)) .* a2(index) );  % equivalent to setround(rnd2); a2(index) = a2(index).*a2(index);
   end
 end
end  % function powerint