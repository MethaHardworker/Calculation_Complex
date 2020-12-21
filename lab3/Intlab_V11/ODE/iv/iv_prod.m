function c =  iv_prod(a,dim)
%IV_PROD  interval structure prod
%
%   c = iv_prod(a,dim)

% written  02/10/16     F. Buenger
% modified 05/16/18     F. Buenger   float input
% modified 05/18/18     F. Buenger   correction for multidimensional arrays

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

float_a = isfloat(a);
if float_a
    size_a = size(a);
else
    size_a = size(a.inf);
end
if nargin == 1
    dim = max([find(size_a>1,1),1]);
end

if float_a
    aa = abs(a);
    setround(-1)
    c.inf = prod(aa,dim);
    setround(1)
    c.sup = prod(aa,dim);
    c = iv_times(prod(sign(a),dim),c);
else
    s = 0;
    ainf = a.inf;
    asup = a.sup;
    
    index = ( asup < 0 );
    if any(index(:))                   % make sure asup>=0
        cc = ainf(index);
        ainf(index) = -asup(index);
        asup(index) = -cc;
        s = sum(index,dim);            % sign information
    end
    index = ( abs(ainf) > abs(asup) );
    if any(index(:))                   % make sure asup >= -ainf >= 0
        cc = ainf(index);
        ainf(index) = -asup(index);
        asup(index) = -cc;
        s = s + sum(index,dim);    	   % sign information
    end
    
    setround(-1)
    cinf = prod(abs(ainf),dim);        % correct for ainf>0 subject to sign
    setround(1)
    csup = prod(asup,dim);             % largest absolute value
    index = ( ainf < 0 );
    if any(index(:))                   % treat (proper) zero intervals
        N = NaN(size(ainf));           % ignore ~index
        D = N;
        N(index) = abs(ainf(index));
        D(index) = abs(asup(index));
        q = max(N./D,[],dim);          % NaN is ignored; bound correct because setround(1)
        index0 = isnan(q);             % not affected indices
        if any(index0(:))
            cc = cinf;
            cinf = -( q.*csup );
            cinf(index0) = cc(index0); % recover not affected indices
        else
            cinf = -( q.*csup );
        end
    end
    
    index = ( 2*round(s/2) ~= s );
    if any(index(:))                   % correct sign
        cc = cinf(index);
        cinf(index) = -csup(index);
        csup(index) = -cc;
    end
    
    index = ( isnan(cinf) | isnan(csup) ) & ( ~any(iv_isnan(a),dim) );
    if any(index(:))                   % take care of 0*inf
        cinf(index) = -inf;
        csup(index) = inf;
    end

    c.inf = cinf;
    c.sup = csup;
end

if rndold ~= 1
    setround(rndold)
end

end % function iv_prod
