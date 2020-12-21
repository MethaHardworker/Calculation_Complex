function y = iv_sqr(x)
%SQR  Implements (elementwise)  sqr(x) = x.^2  for interval like structures
%
%   y = sqr(x)
%

% written  08/01/17     F. Buenger

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

y = x;

y.inf = -(-x.inf .* x.inf); % Recall that rounding is upwards.

index0 = ( x.inf <= 0 ) & ( x.sup >= 0 );
if any(index0(:))
    y.inf(index0) = 0;
end

indexneg = ( x.sup < 0 );
if any(indexneg(:))
    y.inf(indexneg) = x.sup(indexneg) .* x.sup(indexneg);
end

y.sup = x.sup .* x.sup; % Recall that rounding is upwards.

if any(index0(:))
    y.sup(index0) = max( y.sup(index0) , x.inf(index0) .* x.inf(index0), 'includenan' );
end

if any(indexneg(:))
    y.sup(indexneg) = x.inf(indexneg) .* x.inf(indexneg);
end

if rndold ~= 1
    setround(rndold)
end

end % function iv_sqr
  