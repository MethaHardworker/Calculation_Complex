function r = typeadjust(x,y)
%TYPEADJUST  fast conversion (type adjustment) of floating-point or interval  
%            vector/matrix to Taylor model.  
%
%   r = typeadjust(x,r)
%
%          x: float- or intval- or iv - vector/matrix to be converted 
%          r: dummy Taylor model (type "taylormodel") for fast type cast of x

% written  06/01/18     F. Buenger

dummy = y(1);
dummy.monomial = zeros(1,dummy.dim);
dummy. coefficient = 0;
z.inf = 0; % zero interval
z.sup = 0;
dummy.interval = z; 

if isfloat(x)
    [m,n] = size(x);
    %r(1:m,1:n) = dummy; % just preallocation
    r = repmat(dummy,m,n); % just preallocation
    for i = 1:m*n
        dummy.coefficient = x(i);
        dummy.image.inf = x(i);
        dummy.image.sup = x(i);
        r(i) = dummy;
    end
else
    if isa(x,'intval')
        x = intval2iv(x);
    end
    [m,n] = size(x.inf);
    %r(1:m,1:n) = dummy; % just preallocation
    r = repmat(dummy,m,n); % just preallocation
    [xm,xr] = iv_getmidrad(x);
    for i = 1:m*n
        dummy.coefficient = xm(i);
        dummy.interval.inf = -xr(i);
        dummy.interval.sup = xr(i);        
        dummy.image.inf = x.inf(i);
        dummy.image.sup = x.sup(i);
        r(i) = dummy;
    end
end

end % function typeadjust

    