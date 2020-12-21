function A = get_linear_terms(a)
%GET_LINEAR_TERMS   returns the matrix A consisting of the linear terms of the Taylor polynomial vector "a".
%
%   A = get_linear_terms(a)
%
% Precisely, A(i,j) equals a(i).coefficient(k) with a(i).monomial(k) = e_j
% where e_j = (0...0 1 0...0) is the j-th standard basis vector of length a(i).dim.
% If e_j is not a monomial of a(i), then A(i,j):=0.
% It is supposed that all components a(i) have the same dimension a(i).dim

% written  11/17/15     F. Buenger

s = size(a);
if s(1) > 1 && s(2) >1
    error('The input parameter must be a Taylor model vector.')
end
m = length(a);
n = a(1).dim;
A = zeros(m,n);                          % Initialize result wit zeros
for i = 1:m
    idx = find(sum(a(i).monomial,2)==1); % Find linear terms of Taylor model a(i).
    if ~isempty(idx)
      M = a(i).monomial(idx,:);          % Extract corresponding monomials.
      c = a(i).coefficient(idx);         % Extract corresponding coefficients.
      A(i,:) = c' * M;                   % Bring coefficients in correct order and insert them in the i-th row of A.
    end
end

end % function get_linear_terms