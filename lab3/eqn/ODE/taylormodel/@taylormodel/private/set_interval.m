function a = set_interval(a,interval)
%SET_INTERVAL     sets a(i,j).interval := interval for all i,j.      
%
%   a = set_interval(a,interval)

% written  11/10/15     F. Buenger
% modified 02/12/16     F. Buenger  "intval"-components --> intval-like structures

size_a = size(a);
size_iv = size(interval.inf);
if ~all(size_a == size_iv)
    error('Matrix dimensions must agree.')
end

%iv = struct(interval);

for i=1:size_a(1)
    for j=1:size_a(2)
        a(i,j).interval.inf = interval.inf(i,j);
        a(i,j).interval.sup = interval.sup(i,j);
    end
end
end % function set_interval