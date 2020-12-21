function res = parity(p)
%PARITY       Parity of a permutation vector or matrix
%
%   res = parity(p)
%
%Given a permutation vector or matrix p, the result is the parity +1 or -1
%

% written  03/09/18     S.M. Rump
%

  [m,n] = size(p);
  if m==n
    res = det(sparse(p));
  elseif ( m~=1 ) && ( n~=1 )
    error('input must be a permutation vector or matrix')
  else
    I = speye(n);
    res = det(I(:,p));
  end
  