function A = sin(A)
%SIN           Implements  sin(x)  for fl-type (intervals)
%
%  y = sin(x)
%

% written  04/22/14     S.M. Rump
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
%

  if isa(A.value,'intval')
    A = fl(sin(A.value));
  else                      % double input
    
    rndold = getround;
    if rndold
      setround(0)
    end
    
    A = fl(sin(A.value));   % make sure rounding to nearest
    
    if rndold               % restore rounding mode
      setround(rndold)
    end
    
  end
  