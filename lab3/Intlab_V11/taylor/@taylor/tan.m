function r = tan(a)
%TAN          Taylor tangent  tan(a)
%
%Thanks to George Corliss for providing the Taylor expansion
%

% written  06/03/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
%

  global INTLAB_CONST
  
  rndold = getround;
  if rndold
    setround(0)
  end

  K = INTLAB_CONST.TAYLOR_ORDER;
  
  r = a;
  ct = a.t;
  N = size(a.t,2);
  r.t(1,:) = tan(a.t(1,:));
  ct(1,:) = 1 + r.t(1,:).^2;     % 1+tan(a)^2
  r.t(2,:) = ct(1,:) .* a.t(2,:);
  for j=2:K
    ct(j,:) = sum( r.t(1:j,:).*r.t(j:-1:1,:) , 1 );
    r.t(j+1,:) = sum( ct(1:j,:).*a.t(j+1:-1:2,:).*repmat((j:-1:1)',1,N) , 1 )/j;
  end

  if rndold
    setround(rndold)
  end