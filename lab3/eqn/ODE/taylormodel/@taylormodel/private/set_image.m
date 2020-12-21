function a = set_image(a)
%SET_INTERVAL     computes and sets a(i,j).image for all i,j.
%
%   a = set_image(a)

% written  12/15/15     F. Buenger

for i=1:size(a,1)
    for j=1:size(a,2)
        a(i,j).image = image(a(i,j));
    end
end

end % function set_image 