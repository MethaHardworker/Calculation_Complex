function r = typeadjust(x,y)
%TYPEADJUST  fast conversion (type adjustment) of floating-point or interval  
%            vector/matrix to Taylor coefficient.  
%
%   r = typeadjust(x,r)
%
%          x: float- or intval- or iv - vector/matrix to be converted 
%          r: dummy Taylor coefficient (type "taylorcoeff") for fast type cast of x

% written  02/21/18     F. Buenger

dummy = y(1);
z = zeros(length(dummy.inf)-1,1);

if isfloat(x)
    [m,n] = size(x);
    r(1:m,1:n) = dummy; % fast preallocation
    for i = 1:m*n
        r(i).inf = [x(i);z];
        r(i).sup = [x(i);z];
    end
elseif isa(x,'intval')
    x = struct(x);
    [m,n] = size(x.inf);
    r(1:m,1:n) = dummy; % fast preallocation
    for i = 1:m*n
        r(i).inf = [x.inf(i);z];
        r(i).sup = [x.sup(i);z];
    end
elseif isiv(x)
    [m,n] = size(x.inf);
    r(1:m,1:n) = dummy; % fast preallocation
    for i = 1:m*n
        r(i).inf = [x.inf(i);z];
        r(i).sup = [x.sup(i);z];
    end
end

end % function typeadjust

    