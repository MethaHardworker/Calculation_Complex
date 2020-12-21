function p = taylormodel2polynom(a,idx)
%POLYNOMIAL2TAYLORMODEL   transforms a single taylor Model "a" to an interval 
%                         polynomial of INTLAB-type "polynom".
%                          
%   1) p = taylormodel2polynom(a)
%   2) p = taylormodel2polynom(a,idx)  The dependence of all variables with index 
%                                      in "idx" is omitted (default: idx = []).

% written  11/01/16     F. Buenger

e = 1e-30;
if 1+e > 1   % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards in verified mode
end

if ~isa(a,'taylormodel') || numel(a) ~= 1
    error('The input parameter "a" must be a single Taylor model.')
end

if nargin == 1
    idx = [];
else
    if ~isa(idx,'double') || any(round(idx) ~= idx) || size(idx,1) ~= 1 ...
            || any(idx < 0) || any(idx > a.dim) || length(unique(idx)) < length(idx)     
        error('inconsistent parameter "idx"')   
    end
end

n = a.dim;
a_interval = iv2intval(a.interval);
b = a; 

if isfloat(b.coefficient)
    [c0,i0] = get_constant_term(b);
    b.coefficient = intval(b.coefficient);
else 
    i0 = find(all(b.monomial == 0,2),1); % index of the constant term in the monomial and coefficient array
    if isempty(i0)
        i0 = 0;
    end
    b.coefficient = iv2intval(b.coefficient);    
end

if a_interval ~= 0
    if i0 ~= 0
        b.coefficient(i0) = a_interval + b.coefficient(i0); % Shift error interval to constant term.
    else
        b.coefficient = [a_interval;b.coefficient]; % The error interval becomes constant term if a does not have a constant term.
        b.monomial = [zeros(1,n);b.monomial];
    end
end

pv = cell(1,n); % initialization of cell array for variable names

for i = 1:n-1
  pv(i) = {['x',num2str(i)]};    
end

if b.type == 1 % ODE-type, last variable corresponds to time t.
  pv(n) = {'t'}; 
else
  pv(n) = {['x',num2str(n)]};  
end

if ~isempty(idx) % the dependence of all variables with index in "idx" is omitted.
  pv(idx) = [];  % Delete variable names.   
  D = iv2intval(a.domain)-a.center;
  D = D';
  D = D(idx); 
  M = b.monomial(:,idx);
  c = b.coefficient;
  D = repmat(D,size(M,1),1);
  c = c .* prod(D.^M,2);  % Dependency of variables x_i with i in idx is omitted.   
  b.monomial(:,idx) = []; % Delete monomial entries for x_i, i in idx.
  % Coefficients corresponding to the same monomial must be summed up.
  [M,iM] = unique_rows(b.monomial); % Delete double entries and get corresponding indices: U(iU)=M.
  b.monomial = M;
  c_lower = -accumarray(iM,-c.inf); % Recall that rounding is upwards
  c_upper = accumarray(iM,c.sup);
  b.coefficient = infsup(c_lower,c_upper);    
end

pc = b.coefficient;
pe = b.monomial;
p = polynom(pc,pe,pv);

if rndold ~= 1
    setround(rndold)
end

end % function taylormodel2polynom


