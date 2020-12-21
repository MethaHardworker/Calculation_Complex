function B = blunt(A)
%BLUNT   blunting of matrix A
%
%   B = blunt(A)
%
%   input: 
%       A:       nxn - matrix
%
%   output: 
%       B:       nxn - matrix, q-bluntig of A. 
%
%
%   First, the columns of A are sorted according to Euclidean length in deshending order.
%   For this column-reordered matrix As = A*P a QR-factorization is carried out so that As = QR, 
%   where Q is taken such that all diagonal entries of R are nonnegative and P denotes a suitable 
%   permutation matrix.
%
%   Then, B := As + Q*diag(q) = Q*(R+diag(q)) for a vector q > 0 componentwise is the "q-blunting" of A.  
%   Since the diagonal entries of the upper triangular matrix R are nonnegative, B is always regular.
%   The purpose of blunting is that B has a better condition number then As which of course depends on the 
%   blunting factors q_i. 
%
%   These blunting factors q_i are heuristically chosen. We follow [MB_2003] p.33, [MB_2006_a] p.17, [MB_2006_b] p.20:
%   "Chose the blunting factors q_i to be 1E−3 times the length of the longest column vector of the linear matrix." 
%   (Thus all q_i have the same single value which is just called q in the sequel,i.e., the vector q becomes a scalar number.) 
% 
%   The largest Euclidian column norm of A is the length of the first column of As, namely q = norm(As(:,1)). 
%   Thus B = As + q*Q  = Q*(R+q*I) = A*P + q*Q, where I is the identity matrix.
%
%   Finally, the columns of B are rearranged to the original order of A, i.e.: B := B*P^T = A + q*Q*P^T. 
%
%   Note that if A is singular or nearly singular such that q = 0, then B = A, i.e., 
%   A remains unchanged upto column reordering.
%
%   Note also that blunting does not need any verified computations.
%
%
%  [MB_2003]   K. Makino and M. Berz, "Suppression of the wrapping effect by Taylor model - based validated integrators",
%                MSU HEP Report 40910, 2003, mainly p. 20-21 

%  [MB_2006_a] K. Makino and M. Berz, "Suppression of the Wrapping Effect by Taylor Model-based Verified Integrators: Long-term Stabilization by Preconditioning",
%                International Journal of Differential Equations and Applications 10(4) (2005) 353-384. 

%  [MB_2006_b] K. Makino and M. Berz, "Suppression of the Wrapping Effect by Taylor Model-based Verified Integrators: Long-term Stabilization by Shrink Wrapping",
%                International Journal of Differential Equations and Applications 10(4) (2005) 385-403

% written  06/28/17     F. Buenger

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

cond_bound = 1E2; % bound for the nonverified condition number cond(A) of A. 
                  % This bound is chosen quite arbitrarily. Feel free to change that. 
                  % If blunting is switched on and cond(A) > cond_bound, then blunting is executed.

cond_A = cond(A);
if cond_A < cond_bound
    B = A;
    return; 
end

normsA = sqrt(sum(A.*A));  % non-verified computation of the Euclidean norms of the column vectors of A
[normsA,iA] = sort(normsA,'descend');  

As = A(:,iA);   % The columns of A are sorted by Euclidean length in descending order.
[Q,R] = qr(As); % Compute a non-verified QR decomposition As = Q*R.
dummy = sign(diag(R));
dummy(dummy == 0) = 1;
Q = Q.*dummy';  % Q is normalized such that R has nonnegative diagonal. (R is not explicitly needed.)         

% The blunting factors q_i are heuristically chosen. We follow [MB_2003] p.33, [MB_2006_a] p.17, [MB_2006_b] p.20:
%   "Choose the blunting factors q_i to be 1E−3 times the length of the longest column vector of the linear matrix." 
 
q_stretch = 1E-3; % heuristic stretching factor
q_min = 1e-13;    % heuristic lower bound for scaling factor q. Feel free to change that! 
q = max(q_stretch * normsA(1),q_min);
B = As + q*Q;     % B is the q-blunting of A.
B(:,iA) = B;      % B = A + q*Q*P^T final reordering of columns of B. 

% cond(B)/cond(A)  % ratio of condition number of input A and output B. Only for testing!

if rndold ~= 1 
    setround(rndold)
end
end  % function precondition
