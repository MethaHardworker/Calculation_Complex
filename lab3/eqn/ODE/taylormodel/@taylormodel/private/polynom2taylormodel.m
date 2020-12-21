function res = polynom2taylormodel(p)
%POLYNOMIAL2TAYLORMODEL   transform p of type 'polynomial' to a taylormodel.
%
%   res = polynom2taylormodel(p)

% written  01/22/16     F. Buenger
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures

if ~isa(p,'polynom') ||  ~( isfloat(p.c) && isreal(p.c) )
    error('p must be a polynomial with real floating-point coefficients.')
end
if ~ ( isfloat(p.c) && isreal(p.c) )
    error('p must have type polynom')    
end


% In the univariate case q.e is simply the polynomial degree.
% In the multivariate case q.e is a nontrivial matrix which contains
% in each row the exponents of a monomial.

if numel(p.e) == 1 % univariate case
    t_0 = 0;                                          % default evaluation point zero
    domain.inf = -1;                                  % default domain [-1,1]
    domain.sup = 1;
    order = p.e;                                      % = degree of p   
    [dummy,monomial,coefficient] = find(fliplr(p.c)); % type 'taylormodel' uses sparse notation while 'polynom' 
                                                      % uses full notation for univariate polynomials 
    monomial = monomial' - 1;                         % polynomial exponents                                             
    coefficient = coefficient';  
else % multivariate case
    order = max(sum(p.e,2));                          % maximal degree
    coefficient = p.c;
    monomial = p.e;
    n = length(monomial(1,:));                        % number of variables, p = p(x1,...,xn)
    t_0 = zeros(n,1);                                 % default evaluation point t_0(i) = 0, i = 1,...,n.
    domain.inf = -ones(n,1);                          % default domain := [-1,1]^n 
    domain.sup = ones(n,1);  
end

res = taylormodelinit(t_0,domain,order,monomial,coefficient);

end % function polynom2taylormodel