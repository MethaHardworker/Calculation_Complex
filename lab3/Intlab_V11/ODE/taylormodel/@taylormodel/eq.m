function res = eq(a,b)
%EQ  Taylor model logical equality  a == b
%
%  res = eq(a,b)

% written  08/28/15     F. Buenger
% modified 11/23/15     F. Buenger  matrix input (componentwise equality)
% modified 02/10/16     F. Buenger  "intval"-components --> intval-like structures 

if ~isa(a,'taylormodel') || ~isa(b,'taylormodel')
    error('invalid call of "==" for Taylor models')
end

S_a = size(a);
S_b = size(b);

if (length(S_a) > 2) || (length(S_b) > 2)
    error('maximally two dimensions for type taylormodel');
end

if any(S_a ~= S_b)
    error('Matrix dimensions must agree.');
end

res = false(S_a);
for i = 1:S_a(1)
    for j = 1:S_a(2)
        a_ = a(i,j);
        b_ = b(i,j);
        res(i,j) = (    all(a_.dim == b_.dim) ...
                     && all(a_.center == b_.center) ...
                     && all(iv_eq(a_.domain,b_.domain)) ...
                     && a_.order == b_.order ...
                     && all(size(a_.monomial) == size(b_.monomial)) ...
                     && all(min(a_.monomial == b_.monomial)) ...
                     && all(a_.coefficient == b_.coefficient) ...
                     && iv_eq(a_.interval,b_.interval) );
    end 
end 

end % function eq
