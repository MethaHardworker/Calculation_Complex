function res = isfinite(a)
%ISFINITE       Array of 1's for finite components
%
%   res = isfinite(a)
%

% written  08/03/14     S.M. Rump
% modified 04/17/18     S.M. Rump  spell check
%

  res = isfinite(intval(a));
  