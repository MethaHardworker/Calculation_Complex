function k = numvars(p)
%NUMVARS      Number of variables of polynomial
%
%   k = numvars(p)
%
%Careful: number of variables of a constant polynomial (without variables)
%is zero, e.g. removevars(polynom(2)).
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 07/25/18     S.M. Rump  comment
%

  if isempty(p.v)
    k = 0;
  else
    k = size(p.e,2);
  end
