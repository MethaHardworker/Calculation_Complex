function [X,I] = min(X,Y,dim)
%MIN          Smallest element for non-interval gradients
%
%Calling conventions are as in Matlab, i.e. for non-interval data:
%
%   min(X)
%   [Z,I] = min(X)
%   [Z,I] = min(X,[],dim)
%   Z = min(X,Y)
%
%and for up to 2-dimensional interval data:
%
%   min(X)
%   Z = min(X)
%   Z = min(X,[],dim)
%   Z = min(X,Y)
%
%The minimum is taken based on the .x components of gradients. For min(X,Y)
%one parameter may be scalar. The result is a gradient quantity.
%
%For intval (non-gradient) quantities the minimum is defined as
%  min(A,B) := { min(a,b) : a in A, b in B }
%Since the indices of the minimum is not unique, computing I for interval
%data makes no sense.
%For interval data, input must be real.
%

% written  04/04/14     S.M. Rump 
% modified 05/18/14     S.M. Rump  code optimization
% modified 06/29/18     S.M. Rump  interval data
% modified 08/31/18     S.M. Rump  compatible sizes
%

  IV = isintval(X) || isaffari(X);
  
  if (~IV) && ( nargin>=2 )
    IV =  isintval(Y) || isaffari(Y);
  end
  if nargin==3
    if ~isempty(Y)
      error('invalid call')
    end
  end
  
  if IV                       % interval data
    
    if nargout>1
      error('no index computation for interval data')
    end
    
    if nargin~=2              % min(X) or min(X,[],dim)
      [m,n] = size(X.x);
      transp = ( ( nargin==3 ) && ( dim==2 ) ) || ...
        ( ( m==1 ) && ( nargin==1 ) );
      if transp
        X = X';               % now dim=1, i.e. min(X)
        dummy = m;
        m = n;
        n = dummy;
      end
      if m~=1                 % not only one row
        Xx = X.x;
        Xdx = X.dx;
        N = size(Xdx,2);
        if any(isnan(Xx(:))) | any(isnan(Xdx(:)))
          X.x = intval(NaN(1,n));
          X.dx = intval(NaN(1,N));
        else
          X.x = min(Xx,[],1);
%           Xdx(inf(Xx) > sup(X.x),:) = NaN;  % ignore those
          Xdx(bsxfun(@gt,inf(Xx),sup(X.x)),:) = NaN;  % ignore those
          Xdx = Xdx(reshape(permute(reshape(reshape(1:m*n*N,m*n,N)',N,m,n),[2 1 3]),m,n*N));
          X.dx = reshape(intval(min(inf(Xdx),[],1),min(sup(Xdx),[],1),'infsup'),N,n)';
        end
      end
      if transp
        X = X';
      end
    else                      % min(X,Y)
      sX = size(X);
      sY = size(Y);
      if ~isequal(sX,sY)
        if prod(sX)==1
          X = repmat(X,sY);
          sX = sY;
        elseif prod(sY)==1
          Y = repmat(Y,sX);
        else
          error('size not compatible')
        end
      end
      %VVVV  Z = [ X(:) Y(:) ];
      s.type = '()'; s.subs = {':'}; Z = [ subsref(X,s) subsref(Y,s) ];
      %AAAA  Matlab bug fix
      X = reshape(min(Z,[],2),size(Y));
    end
    
  else                        % non-interval data
    
    if nargin==2              % Z = min(X,Y)
      sX = size(X);
      sY = size(Y);
      if ~isequal(sX,sY)
        if prod(sX)==1
          X = repmat(X,sY);
          sX = sY;
        elseif prod(sY)==1
          Y = repmat(Y,sX);
        else
          error('size not compatible')
        end
      end
      %VVVV  Z = [ X(:) Y(:) ];
      s.type = '()'; s.subs = {':'}; Z = [ subsref(X,s) subsref(Y,s) ];
      %AAAA  Matlab bug fix
      [dummy,I] = min(Z.x,[],2);
      K = prod(sX);
      index = (1:K)'+(I-1)*K;
      %VVVV  X = Z(I);
      s.type = '()'; s.subs = {index}; X = subsref(Z,s);
      %AAAA  Matlab bug fix
      X = reshape(X,sX);
    else                      % [Z,I] = min(X)  or  [Z,I] = min(X,[],dim)
      [m,n] = size(X.x);
      if nargin==1
        [dummy,I] = min(X.x);
        if m==1
          %VVVV  X = X(I);
          s.type = '()'; s.subs = {I}; X = subsref(X,s);
          %AAAA  Matlab bug fix
        else
          index = I + (0:n-1)*m;
          %VVVV  X = X(index);
          s.type = '()'; s.subs = {index}; X = subsref(X,s);
          %AAAA  Matlab bug fix
        end
      else                        % nargin==3
        [dummy,I] = min(X.x,[],dim);
        if dim==1
          index = I + (0:n-1)*m;
        else
          index = (I-1)*m + (1:m)';
        end
        %VVVV  X = X(index);
        s.type = '()'; s.subs = {index}; X = subsref(X,s);
        %AAAA  Matlab bug fix
      end
    end
    
  end