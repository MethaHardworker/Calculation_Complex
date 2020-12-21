function [m,n] =  iv_size(a)
%IV_SIZE  size of interval structure
%
%   s = iv_size(a)

% written  02/09/16     F. Buenger

m = size(a.inf,1);
n = size(a.inf,2);

end % function iv_size
