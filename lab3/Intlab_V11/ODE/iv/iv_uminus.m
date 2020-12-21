function a = iv_uminus(a)
%IV_UMINUS  unary minus for interval structure
%
%   s = iv_uminus(a)

% written  08/01/17     F. Buenger
% modified 06/13/18     F. Buenger

if isfloat(a)
    x.inf = -a;
    x.sup = x.inf;
    a = x;
else
    x = a.inf;
    a.inf = -a.sup;
    a.sup = -x;
end

end % function iv_uminus
