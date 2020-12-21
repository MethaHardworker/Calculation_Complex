function R = pred(A,k)
%PRED         fl-predecessor of fl-type A
%
%   R = pred(A)
%
%On return, R consists of the (elementwise) largest fl-numbers with A > R.
%Similarly,
%
%   R = pred(A,k)
%
%produces the  k-th predecessor. Negative k yields k-th successor.
%

% written  10/31/13     S.M. Rump
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
% modified 03/10/18     S.M. Rump  comment
% modified 03/17/18     S.M. Rump  const eta
%

  if isa(A,'intval')
    error('function predecessor not for fl-type intervals')
  end

  rndold = getround;
  if rndold
    setround(0)
  end

  eta = realmin*eps;

  if nargin==1
    k = 1;
  end

  R = A;
  if k<0
    setround(1)
    for i=1:(-k)
      R = R + eta;
    end
  else
    setround(-1)
    for i=1:k
      R = R - eta;
    end
  end
    
  setround(rndold)
  