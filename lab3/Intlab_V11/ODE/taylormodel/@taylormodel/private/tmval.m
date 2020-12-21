function res = tmval(a,x)
% TMVAL verified evaluation of the image of the polynomial part of 
%       a single Taylor model "a" at x. 
%       The input parameter x can be a point vector for evaluating "a" at 
%       a single point, or an interval (given as intval-like structure)
%       for evaluating "a" over an interval.
%
%   res = tmval(a,X)

% written  5/02/16     F. Buenger

if isfloat(x)   % x is a point vector
    X.inf = x;
    X.sup = x;
else           % x is an interval
    if size(x.inf,2) == 1 % convert column vector to row vector 
      X = iv_transpose(x);
    else
      X = x;        
    end    
end

P = iv_minus(X,a.center.');        % P := X - a.center, centering of X 
P = monomial_power(P,a.monomial);  % P := P.^(a.monomial)
P = iv_prod(P,2);                  % rowwise product of the monomial powers
res = iv_dotprod(a.coefficient,P); % res := (a.coefficient)'*P

end  % function tmval