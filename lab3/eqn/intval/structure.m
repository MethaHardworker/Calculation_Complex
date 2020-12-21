function [As,flag] = structure(A,Astruct,Afactor)
%STRUCTURE    Compute structure information, used for structured solvers
%
%   [As,flag] = structure(A,Astruct,Afactor)
%
%input   A        real input matrix A
%        Astruct  structure of input matrix A
%        Afactor  (optional) entrywise factors for A, only for user-defined
%                   structures (NaN avoids check for consistency, see below)
%output  As       structure information for structured solvers
%        flag     (optional) set to 1 if A does not comply with Astruct
%                   (without causing error message)
%
%Given a real matrix A (point or interval) with structure such as being
%  symmetric, toeplitz etc., output As characterizes the input structure.
%  Structure is specified numerically in Astruct (and optionally in
%  Afactor), or symbolically (see below).
%If input matrix A does not comply to the specified structure, an error message
%is given unless there are two input parameters. In the latter case, the flag is 
%set to 1 and output As is set to NaN.
%
%The structure is stored in As with field elements .Phi and .p such that
%  (*)     A(:) = As.Phi * As.p  .
%Structure information As is used by algorithms such as structlss, structeig.
%  Any user-defined linear structure satisfying (*) can be stored in As
%  and passed to those algorithms.
%For example, a matrix without structure is be represented by
%          As.Phi = speye(n^2);  As.p = A(:);
%
%Structured elements in the matrix A are numbered columnwise. The matrix
%  As.Phi is sparse, dimension is n^2 rows and k columns, where k is the
%  number of independent parameters.
%
%Careful: needs quite some memory to avoid interpretation overhead;
%  not recommended for dimension >100
%
%The (real) entries of the matrix Astruct depict the anticipated structure,
%  where equal entries are interpreted as dependent, i.e. to be identical
%  entries in an input matrix A. For example, symmetry could be depicted
%  in Astruct by
%
%       (  1   4   7  )
%       (  4   5   8  )
%       (  7   8   9  )   .
%
%This implies A_12=A_21, A_13=A_31, and A_23=A_32.
%
%For example, structure information of a Hilbert matrix afflicted with absolute
%  tolerances 1e-10 can be obtained as follows:
%
%    n = 5;
%    A = hilb(n);
%    Astruct = reshape(1:n^2,n,n);
%    Astruct = triu(Astruct) + triu(Astruct,1)';
%    As = structure(A,Astruct);
%
%Other structures, like toeplitz matrices, can, for example, be generated by
%
%    A = toeplitz(1:n);
%    As = structure(A,A);
%
%If input Afactor is specified, the entries of A are entrywise multiplied by
%  the elements of Astruct. For example,
%
%         (  0   2   2  )
%    A =  (  7  -1   2  )
%         (  7   7   7  )
%
%               (  1   4   4  )                    (  1   1   2  )
%    Astruct =  (  5   2   4  )   and   Afactor =  (  1   1   3  )
%               (  5   5   3  )                    (  2   3   1  )
%
%implies for the structured system matrix AA
%  AA_11=A11     AA_12=A_12      AA_13=2*A_12
%  AA_21=A21     AA_22=A_22      AA_23=3*A_12
%  AA_31=2*A21   AA_32=3*A_21    AA_33=A_33
%This information about AA is stored in As. Note that elements 23 and 33
%  in A are equal, but treated as independent. Furthermore, the matrix
%  A must comply with Astruct before entrywise multiplication with
%  Afactor, otherwise an error message is given.
%
%Input Afactor=NaN avoids check for consistency. This is useful if
%  given entries should be mirrored according to the given structure.
%  For example, let
%     A = midrad(rand(n),1e-5); A = A-A';
%  then A is not skewsymmetric because diagonal elements are nonzero. Now
%     As = structure(A,'skewsymmetric',NaN);
%  mirrors the lower triangle to the upper and sets diagonal elements to zero.
%
%For convenient use the structure can be specified symbolically. For
%  Astruct being one of the strings
%
%     'symmetric'              A = A'
%     'skewsymmetric'          A = -A'
%     'symmetricToeplitz'      A Toeplitz and A = A'
%     'generalToeplitz'        A unsymmetric Toeplitz
%     'persymmetric'           A = E*A'*E for E = flipud(eye(n))
%     'Hankel'                 E*A is Toeplitz for E = flipud(eye(n))
%     'persymmetricHankel'     E*A is symmetric Toeplitz for E = flipud(eye(n))
%     'circulant'              A = P*A for P = circ(0,1,0,...,0)
%
%the corresponding structure is stored in Astruct. In this case
%  Afactor shall not be specified. Note that input matrix A must be real.
%
%For the first example above, the structure could also be computed by
%
%    n = 5;
%    A = midrad( hilb(n) , 1e-8 );
%    As = structure(A,'symmetric')
%
%If only the structure matrix is needed, use
%    As = structure(zeros(n),struct)
%for struct any of the above.
%
%Based on
%   S.M. Rump: Rigorous Sensitivity Analysis for Systems of Linear and
%     Nonlinear Equations. Math. Comput., 54(10):721-736, 1990.
%See also
%  S.M. Rump: Verification methods: Rigorous results using floating-point arithmetic.
%    Acta Numerica, 19:287-449, 2010. 
%

% written  07/26/99     S.M. Rump
% modified 02/20/01     S.M. Rump  new structure
% modified 10/21/01     S.M. Rump  Hankel added
% modified 08/27/12     S.M. Rump  persymmetricHankel added, real input data
% modified 10/04/12     S.M. Rump  checking correct structure
% modified 05/15/14     S.M. Rump  code optimization
% modified 04/17/18     S.M. Rump  spell check
%

  if nargin==2
    Afactor = 1;
  end
  if ( ~isreal(A) ) || ( ~isreal(Afactor) )
    error('Input data must be real')
  end
  checkstruct = ( nargout~=2 );
  n = dim(A);
  flag = 0;

  if ischar(Astruct)         % prespecified structure by string
    if ( ~isequal(Afactor,1) ) & ( ~isnan(Afactor) )
      error('If structure is specified symbolically Afactor shall not be specified.')
    end
    As_factor = 1;
    Astruct = lower(Astruct);
    switch lower(Astruct)
      case 'symmetric'
        pindex = cumsum(diag(cumsum([1 n:-1:2]))+tril(ones(n),-1));
        As_struct = pindex+tril(pindex,-1)';
        As.Phi = sparse( 1:n^2 , As_struct , 1 , n^2 , n*(n+1)/2 );
        As.p = A(pindex~=0);
      case 'skewsymmetric'
        pindex = cumsum(diag(cumsum([0 (n-1):-1:2]),-1)+tril(ones(n),-1));
        As_struct = pindex+tril(pindex,-1)';
        As_factor = tril(ones(n),-1)-triu(ones(n),1);
        v = 1:n^2;
        v(1:n+1:n^2) = [];
        As_struct = As_struct(v);
        As_factor = As_factor(v)';
        As.Phi = tril( sparse( v , As_struct , As_factor , n^2 , n*(n-1)/2 ) , -1 );
        As.p = A(pindex~=0);
      case 'symmetrictoeplitz'
        As_struct = toeplitz(1:n);
        As.Phi = sparse( 1:n^2 , As_struct , 1 , n^2 , n );
        As.p = A(:,1);
      case 'generaltoeplitz'
        As_struct = toeplitz((1:n)',[1 (n+1):(2*n-1)]);
        As.Phi = sparse( 1:n^2 , As_struct , 1 , n^2 , 2*n-1 );
        As.p = [ A(1:n,1) ; A(1,2:n).' ];
      case 'persymmetric'
        I = 1:n;
        Astruct = [ n.*(I-1)-(I-3).*I/2 ; ones(n-1,n)];
        Astruct = fliplr(cumsum(Astruct,1));
        As_struct = fliplr( triu(Astruct) + triu(Astruct,1)' );
        As.Phi = sparse( 1:n^2 , As_struct , 1 , n^2 , n*(n+1)/2 );
        As.p = A(find(fliplr(triu(ones(n)))));
      case 'hankel'
        As_struct = hankel(1:n,n:2*n-1);
        As.Phi = sparse( 1:n^2 , As_struct , 1 , n^2 , 2*n-1 );
        As.p = [ A(1:n,1) ; A(n,2:n).' ];
      case 'persymmetrichankel'
        As_struct = hankel(1:n,n:-1:1);
        As.Phi = sparse( 1:n^2 , As_struct , 1 , n^2 , n );
        As.p = A(1:n,1);
      case 'circulant'
        As_struct = toeplitz(1:n,[1 n:-1:2]);
        As.Phi = sparse( 1:n^2 , As_struct , 1 , n^2 , n );
        As.p = A(1:n,1);
      otherwise
        error('invalid string for structure')
    end

    % Test for correct input
    % take care of skewsymmetric structure (containing zeros)
    %   if checkstruct & ~isequal( reshape(As.Phi*As.p,n,n) , A )
    if ~isnan(Afactor)
      skew = isequal(Astruct,'skewsymmetric');
      if skew
        v = 1:n^2;
        v(1:(n+1):n^2) = [];
        AA = zeros(size(A));
        if isa(A,'intval'), AA=intval(AA); end
        if n==2
          dummy = As.p(As_struct)';     % careful n=2: must be col vector
        else
          dummy = As.p(As_struct);
        end
        AA(v) = dummy.*As_factor;
      else
        AA = As.p(As_struct).*As_factor;
      end
      if ~isequal(AA,A)
        if checkstruct
          error('specified structure and input matrix A do not match')
        else
          As = NaN;
          flag = 1;
        end
      else
        flag = 0;
      end
    end

  else                       % user-defined structure

    % make sure, Astruct consists of integers 1..k
    Astruct = Astruct(:);
    [v,I] = sort(Astruct);
    w = [ 1 ; diff(v)~=0 ];
    k = sum(w);
    u = find(w~=0);
    [z,J] = sort(I(u));
    J(J) = 1:k;
    v = cumsum(w);
    v = J(v);
    Astruct(I) = v;

    % calculate field elements
    if nargin==3
      As_factor = Afactor;
    else
      As_factor = 1;
    end
    [v,I] = sort(Astruct);
    d = [ true ; diff(v)~=0 ];
    As.Phi = sparse( 1:n^2 , Astruct , As_factor , n^2 , k );
    As.p = A(I(d));

    % Test for correct input
    if ~isequal( As.p(Astruct) , A(:) )
      if checkstruct
        error('specified structure and input matrix A do not match')
      else
        flag = 1;
        As = NaN;
      end
    end

  end