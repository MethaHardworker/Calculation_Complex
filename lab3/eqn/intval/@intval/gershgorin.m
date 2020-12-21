function v = gershgorin(A,p)
%GERSHGORIN   Complex interval vector containing eigenvalues of matrix A
%
%   v = gershgorin(A)
%
% mid(v) = diag(mid(A)),  rad(v) computed by Gershgorin circles
%
%The call
%   v = gershgorin(A,true)
%plots the Gershgorin circles.
%

% written  10/16/98     S.M. Rump
% modified 08/07/02     S.M. Rump  abss instead of abs
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 04/04/14     S.M. Rump  function name
% modified 07/30/16     S.M. Rump  rounding upwards 
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
% modified 03/08/17     S.M. Rump  plot added
%

  rndold = getround;

  v.complex = 1;
  v.inf = [];
  v.sup = [];
  v.mid = mid(diag(A));
  setround(1)
  v.rad = sum( mag( A - diag(v.mid) ) , 2 );

  v = class(v,'intval');
  
  if nargin>1
    if p
      K = 100;          % K meshpoints on circles
      f = 0.04;         % factor to increase x/y range
      vinf = inf(v);
      vsup = sup(v);
      ax = [min(real(vinf)) max(real(vsup)) min(imag(vinf)) max(imag(vsup))];
      dx = (ax(2)-ax(1))*f;
      dy = (ax(4)-ax(3))*f;
      ax = [ax(1)-dx ax(2)+dx ax(3)-dy ax(4)+dy];
      axis(ax)
      hold on
      plot(ax(1:2),[0 0])
      plot([0 0],ax(3:4))
      axis equal
      circ = exp(1i*linspace(0,2*pi,K));
      x = real(circ);
      y = imag(circ);
      vmid = v.mid;
      vrad = v.rad;
      for i=1:dim(A)
        plot(real(vmid(i))+vrad(i)*x,imag(vmid(i))+vrad(i)*y);
      end
      hold off
    end
  end
  
  setround(rndold)
