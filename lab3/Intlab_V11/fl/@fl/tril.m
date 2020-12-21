function A = tril(A,k)
%TRIL         Implements  tril(a,k)  for fl-type
%
%   C = tril(A,k)
%
% functionality as Matlab function tril for matrices
%

% written  10/21/13     S.M. Rump
% modified 03/10/18     S.M. Rump   global removed
%

  if nargin==1
    k = 0;
  end

  A.value = tril(A.value,k);
