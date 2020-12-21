function y = Griewank(x)
%GRIEWANK     Griewank's global optimization test function
%
%   y = Griewank(x)
%
%     n     estimated number of stationary points
%---------------------------------------------------
%     5                    1e13
%    10                    1e25
%    20                    1e51
%    30                    1e77
%    40                    1e103
%    50                    1e129
%  

% written  09/05/15     S.M. Rump
% modified 01/15/16     S.M. Rump  rounding
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
% modified 06/26/18     S.M. Rump  comment
%

  rndold = getround;
  if rndold
    setround(0)
  end
  
  n = size(x,1);
  if isa(x,'intval')
    I = intval((1:n)');
  else
    I = (1:n)';
  end
  y = 1 + sum(sqr(x),1)/4000 - prod(cos(x./sqrt(repmat(I,1,size(x,2)))),1);
  
  if rndold
    setround(rndold)
  end
    