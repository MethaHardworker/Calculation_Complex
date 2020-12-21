function r = evaltaylormodel2(a,x)
% EVALTAYLORMODEL2  computes verified enclosures of (image + error)         
%                   of the single Taylor model "a" on all subdomains 
%                   of "a" given in the rows of x. 
%
%   r = evaltaylormodel2(a,x)
%                       
% Speaking roughly, the function computes res(i) := a(x(i,:)-a.center') + a.interval
% for all row indices i of x.  
% The function has the same purpose as the function evaltaylormodel but is implemented
% differently and can be used for comparison.

% written  12/02/15     F. Buenger
% modified 02/10/16     F. Buenger  "intval"-components --> intval-like structures 
 
if numel(a) ~= 1
    error('Only single Taylor models can be evaluated');
end

if isfloat(x)
    len_x = size(x,1);    
else
    len_x = size(x.inf,1);
end

M = a.monomial;                     
M = repmat(M,len_x,1);              % M is duplicated len_x times.

c = a.coefficient;                  
len_c = length(c);
c = repmat(c,len_x,1);              % c is duplicated len_x times.

xs = iv_minus(x,a.center');         % centering of the evaluation points xs := x - a.center'

xs.inf = repelem(xs.inf,len_c,1);   % Duplicate rows of x len_c times.
xs.sup = repelem(xs.sup,len_c,1);

r = iv_intpower(xs,M);              % P := xs.^M . These are all needed monomial powers

r.inf = [r.inf,c];                  
r.sup = [r.sup,c];
r = iv_prod(r,2);                   % Build all products p(i,k) := c(i) * xs(k,:)^M(i,1) * ... * xs(k,n)^M(i,n), i = 1,...,len_c, k = 1,...,len_x, n := size(M,2). 
r.inf = reshape(r.inf,len_c,len_x);
r.sup = reshape(r.sup,len_c,len_x);
r = iv_transpose(iv_sum(r,1));      % Build sums  s(k) := p(1,k) +...+ p(len_c,k), k = 1,...,len_x. 

r = iv_plus(r,a.interval);          % Add error interval.

end  % function evaltaylormodel2