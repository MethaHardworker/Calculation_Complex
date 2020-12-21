function res = affarivars
%AFFARIVARS   Current number of error terms
%
%   res = affarivars
%

% written  03/07/14   S.M. Rump
% modified 04/23/14   S.M. Rump  set/getappdata replaced by global
% modified 12/13/17   S.M. Rump  comment
%

  global INTLAB_CONST
  
  res = INTLAB_CONST.AFFARI;
