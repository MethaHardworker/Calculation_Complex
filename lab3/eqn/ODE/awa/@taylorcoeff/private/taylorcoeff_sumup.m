function r = taylorcoeff_sumup(a,p)
%TAYLORCOEFF_SUMUP  sums up Taylor coefficients up to order p. 
%
%   r = taylorcoeff_sumup(a,p)

% written  08/07/17     F. Buenger

if nargin < 2 || isempty(p)
    p = length(a(1).inf);
end

n = numel(a);
r.inf = zeros(size(a));
r.sup = r.inf;
for i = 1:n
    x.inf = a(i).inf(p+1:-1:1); % Reverse the order of the first p coefficients so that they are presumably brought in ascending order
    x.sup = a(i).sup(p+1:-1:1); % which is numerically better for summation.
    x = iv_sum(x);
    r.inf(i) = x.inf;
    r.sup(i) = x.sup;
end

end % taylorcoeff_sumup




