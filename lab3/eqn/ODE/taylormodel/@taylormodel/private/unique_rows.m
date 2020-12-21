function [U,iU] =  unique_rows(M)
%UNIQUE_ROWS     direct, faster implementation of [U,~,iU] = unique(M,'rows')
%
%   [U,iU] =  unique_rows(M)

% written  01/12/16     F. Buenger   

% We partially use Kronecker's trick to encode each monomial row of M by one 
% nonnegative integer so that M is encoded in a column vector D of length size(M,1).
% (Note that the Kronecker-trick-encoding requires that M is a nonnegative integer matrix.) 
% Then, for better performance, D instead of M is used to find the unique monomial
% matrix U and the indices iU such that U(iU,:) = M.

b = [1,max(M(:,1:end-1),[],1)+1]; 
b = cumprod(b); 
D = M*b'; % D is the encoding vector of M.

[DU,idx] = sort(D);               % DU = D(idx)
U = M(idx,:);                     % equivalent to U = sortrows(M) 
idx_rev(idx) = 1:length(idx);     % U(idx_rev,:) = M;
I = ( DU(1:end-1) == DU(2:end) ); % Find duplicate rows of DU (these are also the duplicate rows of U).
U(I,:) = [];                      % equivalent to U = unique(M,'rows')
I = [I;0];
iU = flip( size(U,1)+1 -(cumsum(flip(~I))) );
iU = iU(idx_rev);                 % U(iU) = M;  
end % function unique_rows