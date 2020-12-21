function r = mpower(a,b)
%MPOWER  implements  matrix power a^b only for Taylor coefficient 
%        matrix 'a' and positive integer 'b'.
%
%   r = mpower(a,b)

% written  08/02/17     F. Buenger

if numel(a) == 1
    r = a.^b;
    return;
end

if ~ ( isa(a,'taylorcoeff') && ... % 'a' must be a Taylor coefficient
       isa(b,'double') && isreal(b) && numel(b) == 1 && b == round(b) && b > 0 ...  % 'b' must be a positive integer    
     )
    error('wrong input')
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
        r = r*a;
    end
    b = floor(b/2);
    if b ~= 0
        a = a*a;
    end
end
if b_is_even
    r = r*r;
end

end % function mpower