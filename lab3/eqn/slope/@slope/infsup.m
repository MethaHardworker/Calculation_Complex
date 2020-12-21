function infsup(u,str)
%INFSUP       Display of slope variables by infimum and supremum
%
%  infsup(u)
%

%Second parameter name for internal purposes: display full structure
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
% modified 08/26/12     S.M. Rump  global variables removed
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
% modified 09/23/15     S.M. Rump  avoid get(0,... in Octave
%

  global INTLAB_CONST

  INTLAB_SLOPE = INTLAB_CONST.SLOPE;

  if INTLAB_CONST.OCTAVE
    loose = false;
  else
    loose = strcmp(get(0,'FormatSpacing'),'loose');
  end

  numvar = size(u.s,2);
  if numvar~=INTLAB_SLOPE.NUMVAR
    warning('**** number of dependent variables in slope expansion inconsitent')
  end

  if nargin==2
    if isequal(str,'internal')
      disp([ 'slope range (internal structure) ' inputname(1) '.r = ' ])
      display(u.r,'',1);
      disp([ 'slope slope ' inputname(1) '.s = ' ])
      display(u.s,'',1);
      return
    else
      error('Invalid call of slope infsup')
    end
  end

  if loose, disp(' '); end
  disp([ 'slope intval center ' inputname(1) '.c = ' ])
  if loose, disp(' '); end
  ur = u.r(:,1);
  infsup( reshape(ur,u.size) , '' , 1 )
  if loose, disp(' '); end

  if loose, disp(' '); end
  disp([ 'slope intval range ' inputname(1) '.r = ' ])
  if loose, disp(' '); end
  ur = u.r(:,INTLAB_SLOPE.NUMVAR+1);
  infsup( reshape(ur,u.size) , '' , 1 )
  if loose, disp(' '); end

  if loose, disp(' '); end
  disp([ 'slope intval slope ' inputname(1) '.s = ' ])
  if loose, disp(' '); end
  infsup( reshape(u.s,[u.size numvar]) , '' , 1 )
  if loose, disp(' '); end
