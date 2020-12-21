function da = deriv(a,k,idx)
%PDERIV       k-th derivative of polynomial part of Taylor model "a" 
%             with respect to variable with index idx 
%
%   da = deriv(a,k,idx) 
%

% written  11/08/16     F. Buenger
%

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

if nargin < 3
    error('Not enough input parameters')
end
if numel(a) == 0
    error('empty input');
end
if ~( isreal(k) && numel(k) == 1 && k >= 0 && round(k) == k )
    error('Second input parameter k must be a nonnegative integer.');
end
if ~( isreal(idx) && numel(idx) == 1 && idx >= 1 && round(idx) == idx )
    error('Third input parameter idx must be a positive integer.');
end
S_a = size(a);
if length(S_a) ~= 2 
    error('Only two dimensions for Taylor models are allowed.')
end

da = a; % Initialize derivative.
if k == 0 % zero-th derivative => do nothing
    if rndold ~= 1
        setround(rndold)
    end    
    return
end

for i = 1:S_a(1)
    for j = 1:S_a(2)
        da_ = a(i,j);
        n = da_.dim;
        if idx > n
            error('index exceeds Taylor model dimension.')
        end
        E = da_.monomial(:,idx);       % Exponents of variable idx
        I = (E < k);
        
        da_.interval.inf = 0;          % zero error interval
        da_.interval.sup = 0;        
        da_.coefficient(I) = [];  
        da_.monomial(I,:) = [];
        E(I) = [];
        if isempty(da_.coefficient)
            da_.coefficient = 0;       % Return zero Taylor model.
            da_.monomial = zeros(1,n);
            da_.image.inf = 0;
            da_.image.sup = 0;
        else
            da_.monomial(:,idx) = E-k; % Reduce the exponents of variable idx by k.
            
            % Compute coefficients of derivative.
            P = [ repmat(E,1,k) - repmat(0:k-1,length(E),1) , da_.coefficient ];
            setround(-1)
            c.inf = prod(P,2);
            setround(1)
            c.sup = prod(P,2);
            c_mid = (c.inf+c.sup)/2;
            da_.coefficient = c_mid;   % approximate midpoint
            rest = da_; 
            rc.sup = c.sup-c_mid;      % Recall that rounding is upwards.
            rc.inf = -(c_mid-c.inf); 
            rest.coefficient = rc;
            da_.interval = image(rest);
            da_.image = image(da_);
        end        
        da(i,j) = da_;
    end
end

if rndold ~= 1
    setround(rndold)
end

end % function deriv
