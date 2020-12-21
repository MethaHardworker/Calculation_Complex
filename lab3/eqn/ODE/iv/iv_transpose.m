function a =  iv_transpose(a)
%IV_TRANSPOSE  interval structure transpose a'
%
%   c = iv_transpose(a)

% written  02/10/16     F. Buenger

a.inf = a.inf';
a.sup = a.sup';

end % function iv_transpose
