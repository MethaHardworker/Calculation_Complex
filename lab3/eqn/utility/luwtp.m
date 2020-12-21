function [L,U,p,q,growth] = luwtp(A)
%LUWTP        Gaussion elimination with total pivoting
%
%    [L,U,p,q,growth] = luwtp(A);
%
%Input A may be rectangular and double, complex, intval, fl, affari, cpair
%or dd. %In either case the decomposition is performed in the corresponding 
%arithmetic. 
%If no zero pivots occur, L*U = A(p,q) after execution.
%If specified, growth is the growth factor.
%
%To solve a linear system with total pivoting use  solvewtp(A,b).
%

% written  10/17/17     S.M. Rump
% modified 03/18/18     S.M. Rump  cpair, dd
% modified 04/17/18     S.M. Rump  spell check
%

  [m,n] = size(A);
  is_ival = isintval(A);
  p = 1:m;
  q = 1:n;
  growth = 0;

  for i=1:n-1
    if i<m
      [v,indexi] = max(mag(A(i:m,i:n)),[],1);
      [dummy,indexj] = max(v);
      indexi = indexi(indexj);          % pivot element (indexi,indexj)
      indexi = indexi+i-1;
      if indexi~=i                      % row interchange due to pivoting
        A([i indexi],:) = A([indexi i],:);
        p([i indexi]) = p([indexi i]);
      end
      indexj = indexj+i-1;
      if indexj~=i                      % col interchange due to pivoting
        A(:,[i indexj]) = A(:,[indexj i]);
        q([i indexj]) = q([indexj i]);
      end
      if ( is_ival && in(0,A(i,i)) ) || ( A(i,i)==0 )
        L = NaN(m);
        U = NaN(m,n);
        return
      end
      A(i+1:m,i) = A(i+1:m,i)/A(i,i);
      j1 = i+1:m;
      j2 = i+1:n;
      A(j1,j2) = A(j1,j2) - A(j1,i)*A(i,j2);
      growth = max(growth,max(max(mag(A(:)))));
    end
  end
  if m>n
    j = n+1:m;
    A(j,n) = A(j,n)/A(n,n);
    growth = max(growth,max(max(mag(A(:)))));
  end

  mn = min(m,n);
  L = tril(A(1:m,1:mn),-1) + eye(m,mn);
  U = triu(A(1:mn,1:n));
  growth = growth/max(mag(A(:)));
  