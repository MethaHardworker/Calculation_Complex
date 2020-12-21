function r = rdivide(a,b)
%RDIVIDE  Taylor model right division a ./ b
%
%   r = rdivide(a,b)

% written  09/03/15     F. Buenger  scalar input
% modified 11/24/15     F. Buenger  matrix input (componentwise evaluation)
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures

e = 1e-30;
if 1+e == 1-e                          % fast check for rounding to nearest
    rndold = 0;
else
    rndold = getround;
    setround(0)
end

if (length(size(a)) > 2) || (length(size(b)) > 2)
    error('maximally two dimensions are allowed for type taylormodel');
end

if isa(b,'intval')  
  b = intval2iv(b);
end        
    
if isfloat(b) || isiv(b)  
   c = iv_rdivide(1,b);
   r = a.*c;
else
    r = a.*inv(b);
end

if rndold
    setround(rndold)
end

end % function rdivide