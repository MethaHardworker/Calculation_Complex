function r = iv_hull(varargin)
%IV_HULL  interval hull
%
%   r = hull(a,b,c,...)
%
% It is assumed that all input parameters are either interval-like strucrures or real floats. 

% written  08/10/17     F. Buenger

l = length(varargin);
if l == 2   
    a = varargin{1};
    b = varargin{2};
    if isfloat(a)
        x.inf = a;
        x.sup = a;
        a = x;
    end
    if isfloat(b)
        x.inf = b;
        x.sup = b;
        b = x;
    end  
    if ~isiv(a) || ~isiv(b)
        error('interval hull called with invalid input')    
    end
    if isempty(a) || (isempty(a.inf) && isempty(a.sup)) 
        r = b;
        return
    end
    if isempty(b) || (isempty(a.inf) && isempty(a.sup))
        r = a;
        return
    end 
    if ~ ( isequal(size(a.inf),size(a.sup)) && ...
           isequal(size(b.inf),size(b.sup)) && ...
           (isequal(size(a.inf),size(b.inf)) || numel(a.inf) == 1 || numel(b.inf) == 1 ) )        
       
        error('interval hull called with non-matching dimensions')
    end    
    r.inf = min(a.inf,b.inf,'includenan');
    r.sup = max(a.sup,b.sup,'includenan');    
elseif l > 2
    r = varargin{1};
    for i = 2:l
        r = iv_hull(r,varargin{i});          % calls intval function hull
    end
elseif l == 1
    r = varargin{1};
    if isfloat(r)
        x.inf = r;
        x.sup = r;
        r = x;
    end           
else
    error('hull called without parameter')
end

end % function iv_hull
