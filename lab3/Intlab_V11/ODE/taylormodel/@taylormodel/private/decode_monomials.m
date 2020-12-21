function M = decode_monomials(D,d)
%DECODE_MONOMIALS   reverses the encoding of the monomial matrix M by
%  encode_monomials, i.e. M = decode_monomials(encode_nonomials(M)). 
%
%   M = decode_monomials(D,d)

% written  01/13/16     F. Buenger   

n = length(d);          % monomial length
M = zeros(length(D),n); % Initialize result M;
D = D-1;                % Reverse the final shift "+1" done in encode_nonomials

for i = n:-1:2
  r = mod(D,d(i));  
  M(:,i) = (D-r)/d(i);
  D = r;        
end
M(:,1) = D;             % Always d(1) = 1 holds true. Therefore this abbreviates running the previous for-loop downto 1.

end % function decode_monomials