function a = set_domain(a,idx,iv,x)
%SET_DOMAIN     sets a(i,j).domain(idx) := iv and a(i,j).center(idx) := x 
%               for all i,j and recomputes the image a.image with the new domain.
%               It is assumed that all domains "a(i,j).domain" are equal!
%
%   a = set_domain(a,idx,interval,x)

% written  11/09/15     F. Buenger  
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures

s = size(a);
for i = 1:s(1)
    for j = 1:s(2)
        a(i,j).domain.inf(idx) = iv.inf;
        a(i,j).domain.sup(idx) = iv.sup;
        a(i,j).center(idx) = x;
        if any(a(i,j).monomial(:,idx)) % if a(i,j) depends on the variable with index idx, then the image must be recomputed.
            a(i,j).image = image(a(i,j));
        end
    end
end

end % function set_domain