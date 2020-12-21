function p = invperm(p)
%INVPERM      Inverse permutation,  p(q) = q(p) = 1:n  for length(p)=n
%
%  q = invperm(p);
%

% written  12/01/97     S.M. Rump
% modified 03/10/04     S.M. Rump  improved performance
%

  [dummy p] = sort(p);
