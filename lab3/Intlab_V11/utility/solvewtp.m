function x = solvewtp(A,b)
%SOLVEWTP     Solution of square linear system by GE with total pivoting, generic routine
%Solution of linear system using Gaussian elimination with total pivoting
%
%    x = solvewtp(A,b);
%
%A is square, b may be matrix of right hand sides.
%Input A may be double, complex, intval, fl, affari, gradient, hessian,
%cpair or dd.
%In either case the linear system is solved in the corresponding arithmetic. 
%

% written  10/17/17     S.M. Rump
% modified 03/18/18     S.M. Rump  cpair, dd, and correction
% modified 04/17/18     S.M. Rump  spell check
%

  global INTLAB_CONST

  if INTLAB_CONST.OCTAVE
    b = typeadj(b,typeof(A));
  end

  [L,U,p,q] = luwtp(A);

  y = forward(L,b(p,:));
  x = backward(U,y);
  x = x(invperm(q),:);
