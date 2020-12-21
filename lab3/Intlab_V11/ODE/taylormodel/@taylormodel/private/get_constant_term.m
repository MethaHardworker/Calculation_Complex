function [c0,i0] = get_constant_term(a)
%GET_CONSTANT_TERM   returns the constant coefficient c0 of the Taylor polynomial, i.e.,
%                    the coefficient for the zero-monomial, and also its index i0 in the
%                    coefficient array so that c0 = a.coefficient(i0).
%
%   [c0,i0] = get_constant_term(a)

% written  09/03/15     F. Buenger
% modified 11/23/15     F. Buenger  matrix input

S_a = size(a);
if (length(S_a) > 2)
    error('maximally two dimensions for type taylormodel');
end

c0 = zeros(S_a); % Initialize c0.
i0 = zeros(S_a); % Initialize i0.

for i = 1:S_a(1)
    for j = 1:S_a(2)        
        a_ = a(i,j);        
        i0_ = find(all(a_.monomial == 0,2),1); % index of the constant term in the monomial and coefficient array
        if ~isempty(i0_)
            c0_ = a_.coefficient(i0_); 
        else
            c0_ = 0;
            i0_ = 0; 
        end
        c0(i,j) = c0_; 
        i0(i,j) = i0_; 
    end % j    
end % i

end % function get_constant_term

