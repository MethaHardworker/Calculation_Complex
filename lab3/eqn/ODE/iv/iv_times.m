function c =  iv_times(a,b)
%IV_TIMES  interval structure multiplication  a .* b
%
%   c = iv_times(a,b)

% written   02/09/16     F. Buenger
% modified  08/01/17     F. Buenger  parameter 'includenan' in min-, max-calculations added

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

if isfloat(a) % case "float .* float/iv"
    if isfloat(b) % case "float .* float"
        c.inf = -(-a.*b);
        c.sup = a.*b;
    else % case "float * iv"
        c.inf = min(-(-a.*b.inf),-(-a.*b.sup),'includenan');
        c.sup = max(a.*b.inf,a.*b.sup,'includenan');
%         c.inf = min(-(-a.*b.inf),-(-a.*b.sup));
%         c.sup = max(a.*b.inf,a.*b.sup);
    end
elseif isfloat(b) % case "iv * float"
    c.inf = min(-(a.inf.*(-b)),-(a.sup.*(-b)),'includenan');
    c.sup = max(a.inf.*b,a.sup.*b,'includenan');
%     c.inf = min(-(a.inf.*(-b)),-(a.sup.*(-b)));
%     c.sup = max(a.inf.*b,a.sup.*b);
else % case "iv * iv"
%     c.inf = min( min( -(-a.inf.*b.inf) ,-(-a.inf.*b.sup),'includenan' ) , min( -(-a.sup.*b.sup) , -(-a.sup.*b.inf),'includenan' ) , 'includenan' );
%     c.sup = max( max( a.inf.*b.inf , a.inf.*b.sup ,'includenan') , max( a.sup.*b.sup , a.sup.*b.inf ,'includenan') , 'includenan' );    
    
%     c.inf = min( min( -(-a.inf.*b.inf) ,-(-a.inf.*b.sup)) , min( -(-a.sup.*b.sup) , -(-a.sup.*b.inf)));
%     c.sup = max( max( a.inf.*b.inf , a.inf.*b.sup) , max( a.sup.*b.sup , a.sup.*b.inf));

    p1 = -(-a.inf.*b.inf);
    p2 = -(-a.inf.*b.sup);
    p3 = -(-a.sup.*b.sup);
    p4 = -(-a.sup.*b.inf);

    p5 = a.inf.*b.inf;
    p6 = a.inf.*b.sup;
    p7 = a.sup.*b.sup;
    p8 = a.sup.*b.inf;

    c.inf = min(min(p1,p2),min(p3,p4));
    c.sup = max(max(p5,p6),max(p7,p8));
    if any(any(isnan(p1) | isnan(p2) | isnan(p3) | isnan(p4) | ...
       isnan(p5) | isnan(p6) | isnan(p7) | isnan(p8) ))
        error('NaN occured during multiplication');
    end
end

if rndold ~= 1
    setround(rndold)
end

end % function iv_times
