function r = sum(a,dim)
%SUM  Implements  sum(a,dim)  for Taylor coefficients
%
%   r = sum(a,dim)
%
% functionality as Matlab function sum for matrices, parameter dim optional

% written  08/02/17     F. Buenger
% modified 01/15/18     F. Buenger  vectorized version 

global INTLAB_AWA_VARS

e = 1e-30;
if 1+e > 1 % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

status = INTLAB_AWA_VARS.STATUS;
dummy = a(1);
if ~isempty(status)
    treenr = INTLAB_AWA_VARS.TREENR;
    vertexnr = INTLAB_AWA_VARS.VERTEXNR + 1;
    INTLAB_AWA_VARS.VERTEXNR = vertexnr;
    if treenr > 0
        r = INTLAB_AWA_VARS.TREE{treenr}{vertexnr};
    end
end
a = tcoeff2tc(a);

if ~isempty(status) % Compute Taylor coefficients of order = status.  
    k = status +1;
    if nargin < 2 || isempty(dim)
        r_.inf = -sum(-a.inf(:,:,k)); % Recall that rounding is upwards 
        r_.sup = sum(a.sup(:,:,k));
    else
        r_.inf = -sum(-a.inf(:,:,k),dim); % Recall that rounding is upwards 
        r_.sup = sum(a.sup(:,:,k),dim);
    end    
    if status == 0
        r = a; % just storage preallocation for r
        [m,n,p] = size(a.inf);    
        if nargin < 2 || isempty(dim)
            if m == 1 
                r.inf = zeros( [1 1 p]);
            else
                r.inf = zeros( [1 n p]);
            end
        elseif dim == 1
            r.inf = zeros( [1 n p]);
        else
            r.inf = zeros( [m 1 p]);            
        end
        r.sup = r.inf;
    end
    r.inf(:,:,k) = r_.inf;
    r.sup(:,:,k) = r_.sup;
else % Compute Taylor coefficients of all orders.  
    r = a; % just storage preallocation for r   
    if nargin < 2 || isempty(dim)
        if size(a.inf,1) == 1
            dim = 2;  
        else
            dim = 1;
        end
    end
    r.inf = -sum(-a.inf,dim); % Recall that rounding is upwards 
    r.sup = sum(a.sup,dim);
end

if ~isempty(status)
    INTLAB_AWA_VARS.TREE{abs(treenr)}{vertexnr} = r; % Store result.
end
r = tc2tcoeff(r,dummy);

if rndold ~= 1
    setround(rndold)
end

end % function sum