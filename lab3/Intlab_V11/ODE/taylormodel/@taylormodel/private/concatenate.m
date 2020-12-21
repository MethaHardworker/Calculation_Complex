function r = concatenate(a,b)
%CONCATENATE     Taylor model concatenation a o b
%
%   r =  concatenate(a,b)  
% 
% The following is assumed: 
%  
%  1) a.type = 1 = b.type, i.e., a and b are ODE-Taylor models having domain 
%     D := [-1,1]^(d-1) x [t,t+h] and center point (0,...,0,t).
%  2) a is an mxn-Taylor model matrix and b is a Taylor model vector of length d-1.
%
%  The result "r" is, like "a", an mxn-Taylor-model matrix of type 1, such that
%  r+I = (a+J)o(b+K) := a(b+K,t)+J 
%  is the concatenation of (a+J), (b+K) in Taylor model arithmetic, where 
%  I,J,K denote the corresponding error intervals for a,b,r, respectively.
%  Roughly speaking, b_1,....,b_{d-1} are inserted in the d-1 space coordinates
%  of each single Taylor model a(i,j) while the time variable "t" remains unchanged.

% written  04/25/17     F. Buenger

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards in verified mode
end

S_a = size(a);
S_b = size(b);

if min(S_b)~= 1
    error('Second input paraneter must be a vector.')
end
if S_b(1) ~= 1
    b = b'; % transformation to row vector
    S_b = size(b);
end

d = a(1).dim;

if S_b(2) ~= d-1
    error('inconsisitent input dimensions')
end

r = a; % cheap storage preallocation of result

% First generate a matrix M that contains all monomials of all Taylor models 
% components of "a", one beneath the other. 
numel_a = S_a(1)*S_a(2); % numel_a = numel(a)
m = ones(numel_a+1,1);
M = [];
for k = 1:numel_a
  a_ = a(k);  
  M = [M;a_.monomial];     
  m(k+1) = m(k)+length(a_.coefficient); % M(m(k),:) is first monomial of a(k).
                                        % M(m(k+1)-1,:) is last monomial of a(k).
end
Mt = M(:,d);     % Mt contains the exponents of the time variable. 
Mx = M(:,1:d-1); % Mx contains the exponents of the space variables.
bMx = b.^Mx;

% For performance reasons the computation of the rowwise product P := prod(bMx,2) is split into 
% 
%   P = prod(bMx(:,1:2),2) * prod(bMx(:,3:4),2) * ..... 
%
% The benefit is that the rowwise subproducts of two-columns, for example p_12 := prod(bMx(:,1:2),2),  
% contain - with high probability - duplicate lines corresponding to the same pair of exponents. 
% Then, only products of "unique" lines can be computed and the results can be duplicated afterwards.  

for k = 1:2:d-1
    if k < d-1
        [U,iMx,iU] = unique(Mx(:,k:k+1),'rows');    % Find "unique" lines in the two columns k,k+1 of the exponent matrix Mx. U is not needed, only iMx and iU are wanted. 
        p = prod(bMx(iMx,k:k+1),2);                 % Compute the rowwise product of the unique lines ...
        p = p(iU,:);                                % ... and duplicate the results to all lines.
    else
        p = bMx(:,d-1);                             % If the number of columns of bMx is odd, then a single column remains at the end and becomes the final factor.
    end
    if k == 1
        P = p;                                      % Initialize total product P with first factor p for k = 1.
    else
        P = P.*p;                                   % Build up the whole product: prod(bMx,2) = prod(bMx(:,1:2),2) * prod(bMx(:,3:4),2) * .....
    end
end
                                           
for k = 1:numel_a
    a_ = a(k);
    c = a_.coefficient;
    r_ = P(m(k):m(k+1)-1);                          % r_ = prod(b.^M_,2) where M_ = a_.monomial 
    Mt_ = Mt(m(k):m(k+1)-1);                        % Mt_ contains the time exponents.   
    r_ = sdotprod(c,r_,Mt_);                        % Compute r_ := sum c(j) * r_(j) * t^Mt_(j), where j runs from 1 to n := # of monomials of a_ .
    r_.interval = iv_plus(r_.interval,a_.interval); % Sum up error intervals.
    r(k) = r_ ;                                     % Roughly speaking we computed r_+I = (a_+J) o (b+K) = a_(b+K)+J .
end

if rndold ~= 1 
    setround(rndold)
end

end % function concatenate
