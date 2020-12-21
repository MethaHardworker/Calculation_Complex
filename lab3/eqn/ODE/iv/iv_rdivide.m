function c =  iv_rdivide(a,b)
%IV_RDIVIDE  interval structure right division  a ./ b
%
%   c = iv_rdivide(a,b)

% written  02/09/16     F. Buenger
% modified 08/01/17     F. Buenger  parameter 'includenan' in min-, max-calculations added

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

if isfloat(a) 
    if isfloat(b) % case "float / float"
        if any(isinf(b(:)))
            error('division by +/- Inf.')
        elseif any( b(:) == 0 )
            error('division by interval conatining zero')
        end
        c.inf = -((-a)./b); % Recall that rounding is upwards.
        c.sup = a./b;        
    else  % case "float / iv"
        s = [b.inf,b.sup];
        if any(isinf(s(:)))
            error('division by +/- Inf.')
        elseif any(any( and(b.inf <= 0, b.sup >= 0) ))
            error('division by interval conatining zero')
        end
        c.inf = min(-((-a)./b.inf),-((-a)./b.sup),'includenan');
        c.sup = max(a./b.inf,a./b.sup,'includenan');
    end
elseif isfloat(b) % case "iv / float"        
    if any(isinf(b(:)))
        error('division by +/- Inf.')
    elseif any(b(:) == 0)  
        error('division by zero')        
    end
    c.inf = min(-((-a.inf)./b),-((-a.sup)./b),'includenan');
    c.sup = max(a.inf./b,a.sup./b,'includenan');
else % case "iv / iv"
    s = [b.inf,b.sup];
    if any(isinf(s(:))) 
        error('division by +/- Inf.')
    elseif any( and(b.inf(:) <= 0, b.sup(:) >= 0) )
        error('division by interval conatining zero')        
    end
    c.inf = min( min( -((-a.inf)./b.inf) , -((-a.inf)./b.sup) ,'includenan' ) , min( -((-a.sup)./b.sup) , -((-a.sup)./b.inf) , 'includenan' ) , 'includenan' );
    c.sup = max( max( a.inf./b.inf , a.inf./b.sup ,'includenan' ) , max( a.sup./b.sup , a.sup./b.inf , 'includenan' ) ,'includenan' );
end

% s = [c.inf,c.sup];
% if any(isnan(s(:))) || any(isinf(s(:)))
%     error('NaN/Inf occured.')
% end

if rndold ~= 1 
    setround(rndold)
end

end % function iv_rdivide
