function [X,I] = max(X,Y,dim)
%MAX          Largest element for gradients
%
%Calling conventions are as in Matlab, i.e. for non-interval data:
%
%   max(X)
%   [Z,I] = max(X)
%   [Z,I] = max(X,[],dim)
%   Z = max(X,Y)
%
%and for up to 2-dimensional interval data:
%
%   max(X)
%   Z = max(X)
%   Z = max(X,[],dim)
%   Z = max(X,Y)
%
%The maximum is taken based on the .x components of gradients. For max(X,Y)
%one parameter may be scalar. The result is a gradient quantity.
%
%For intval (non-gradient) quantities the maximum is defined as
%  max(A,B) := { max(a,b) : a in A, b in B }
%Since the indices of the maximum is not unique, computing I for interval
%data makes no sense.
%For interval data, input must be real.
%

% written  04/04/14     S.M. Rump 
% modified 05/18/14     S.M. Rump  code optimization
% modified 08/01/16     S.M. Rump  comment
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

    if nargin~=2              % max(X) or max(X,[],dim)
      [m,n] = size(X.x);
      transp = ( ( nargin==3 ) && ( dim==2 ) ) || ...
        ( ( m==1 ) && ( nargin==1 ) );
      if transp
        X = X';               % now dim=1, i.e. max(X)
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
          X.x = max(Xx,[],1);
%           Xdx(sup(Xx) < inf(X.x),:) = NaN;  % ignore those
          Xdx(bsxfun(@lt,sup(Xx),inf(X.x)),:) = NaN;  % ignore those
          Xdx = Xdx(reshape(permute(reshape(reshape(1:m*n*N,m*n,N)',N,m,n),[2 1 3]),m,n*N));
          X.dx = reshape(intval(min(inf(Xdx),[],1),max(sup(Xdx),[],1),'infsup'),N,n)';
        end
      end
      if transp
        X = X';
      end
    else                      % max(X,Y)
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
      X = reshape(max(Z,[],2),size(Y));
    end
            
  else                        % non-interval data
  
    if nargin==2              % Z = max(X,Y)
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
      [dummy,I] = max(Z.x,[],2);
      K = prod(sX);
      index = (1:K)'+(I-1)*K;
      %VVVV  X = Z(I);
      s.type = '()'; s.subs = {index}; X = subsref(Z,s);
      %AAAA  Matlab bug fix
      X = reshape(X,sX);
    else                      % [Z,I] = max(X)  or  [Z,I] = max(X,[],dim)
      [m,n] = size(X.x);
      if nargin==1
        [dummy,I] = max(X.x);
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
        [dummy,I] = max(X.x,[],dim);
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
  