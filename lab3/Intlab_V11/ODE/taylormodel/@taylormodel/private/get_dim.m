function D = get_dim(A)
%GET_DIM          dimension of Taylor model
%
%    n = dim(A)


% written  11/17/15     F. Buenger
% modified 11/23/15     F. Buenger  matrix input
%

S = size(A);
D = zeros(S);
for i = 1:S(1)
    for j = 1:S(2)
        D(i,j) = A(i,j).dim;
    end
end

end % function get_dim
