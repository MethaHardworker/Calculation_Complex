function res = float2iv(a)
%FLOAT2IV  convert float to interval structure
%
%   s = float2iv(a)

% written  02/09/16     F. Buenger

res.inf = a;
res.sup = a;

end % function float2iv
