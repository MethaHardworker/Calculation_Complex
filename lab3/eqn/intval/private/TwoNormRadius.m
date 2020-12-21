function N = TwoNormRadius(A,sym)
% upper bound N on the 2-norm of non-negative A
% sym=1  A symmetric, thus norm(A,2)=rho(A)
% rounding mode upwards after execution
  setround(0)           % set rounding to nearest
  % approximation of norm(rad(A)) = rho(rad(A))
  x = ones(size(A,2),1);
  m = 1;
  M = 2;
  iter = 0;
  if sym                % input matrix symmetric
    while ( abs((M-m)/(m+M))>.1 ) && ( iter<10 )
      iter = iter+1;
      y = A*x;
      x = y./x;
      M = max(x);
      m = min(x);
      scale = max(y);
      x = max( y/scale , 1e-12 );
    end
    % upper bound N of norm(rad(A)) = rho(rad(A))
    setround(1)
    y = A*x;
    N = max(y./x);
  else
    At = A';
    while ( abs((M-m)/(m+M))>.1 ) && ( iter<10 )
      iter = iter+1;
      y = At*(A*x);
      x = y./x;
      M = max(x);
      m = min(x);
      scale = max(y);
      x = max( y/scale , 1e-12 );
    end
    % upper bound N of norm(rad(A)) = rho(rad(A))
    setround(1)
    y = At*(A*x);
    N = sqrt(max(y./x));
  end
end  % function TwoNormRadius
