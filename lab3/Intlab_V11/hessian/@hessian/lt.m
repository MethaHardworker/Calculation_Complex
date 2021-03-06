function res = lt(a,b)
%LT           Implements  a < b  for hessians, compares only a.x and b.x
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/18/14     S.M. Rump  code optimization
%

  if ~isa(a,'hessian')
    if isreal(a) && isreal(b.x)
      res = ( a<b.x );
    else
      res = ( real(a)<real(b.x) ) & ( imag(a)<imag(b.x) );
    end
  elseif ~isa(b,'hessian')
    if isreal(a.x) && isreal(b)
      res = ( a.x<b );
    else
      res = ( real(a.x)<real(b) ) & ( imag(a.x)<imag(b) );
    end
  else
    if isreal(a.x) && isreal(b.x)
      res = ( a.x<b.x );
    else
      res = ( real(a.x)<real(b.x) ) & ( imag(a.x)<imag(b.x) );
    end
  end
