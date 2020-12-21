function [D,d] =  encode_monomials(M,max_monomial,d)
%ENCODE_MONOMIALS     used for conversion of multivariate polynomials to univariate ones.
%
%   [D,d] =  encode_monomials(M,max_monomial,d)
%
% Kronecker's trick is used to encode each monomial row of M by one nonnegative integer 
% so that M is encoded in a column vector D of length size(M,1).
% Recall that the Kronecker's trick requires that M is a nonnegative integer matrix  
% and that max_monomial > max(M) componentwise.

% written  01/13/16     F. Buenger  
if nargin < 3 || isempty(d)
  d = [1,max_monomial(1:end-1)+1]; 
  d = cumprod(d);
end
D = M*d'+1; % D is the encoding vector of M. The final shift "+1" ensures that all components of D are positive.
            % This is needed because D shall be used as an index set. 
end % function encode_monomials