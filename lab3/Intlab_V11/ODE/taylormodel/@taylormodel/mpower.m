function r = mpower(a,b)
%MPOWER  implements  a^b  for Taylor model matrix "a" and single positive integer b.
%        If "a" is a single Taylor model, then a^b is the same as a.^b.     
%
%   r = mpower(a,b)

% written  11/05/15     F. Buenger  pointwise matrix power
% modified 11/23/15     F. Buenger  matrix power a^b, only for Taylor model matrix 'a' and positive integer 'b', same as a.^b 
% modified 05/24/18     F. Buenger  "real" matrix power a^b

if numel(a) == 1 
    r = a.^b;   % If "a" is a singleton, than a^b is the same as a.^b. 
                % It is supposed that a user wanted to write "a.^b" and forgot the "."
                % which happens very often. With this convention it doesn't matter.
                % Note that in contrast MATLAB causes an error for 2^[3,4].
else
    if ~( isa(a,'taylormodel') && ... % 'a' must be a Taylor model
            isa(b,'double') && isreal(b) && numel(b) == 1 && b == round(b) && b > 0 ...  % 'b' must be a nonnegative integer
            )
        error('incompatible input arguments')
    end
    if b == 2*floor(b/2)
        b = b/2;
        b_is_even = 1;
    else
        b_is_even = 0;
    end
    b = b - 1; % b is at least 1.
    r = a;
    while b > 0
        if mod(b,2) == 1
            r = r * a;
        end
        b = floor(b/2);
        if b ~= 0
            a = a*a;
        end
    end
    if b_is_even
        r = r*r;
    end
end

end % function mpower