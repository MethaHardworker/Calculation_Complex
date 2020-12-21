function res = prod(a,dim)
%PROD  Taylor model product prod(a,dim)
%
%   res = prod(a,dim)
%
% functionality as MATLAB function prod for matrices, parameter dim optional

% written  04/26/17     F. Buenger

S_a = size(a);

no_dim = (nargin == 1 || isempty(dim));

if (S_a(1) == 1 || S_a(2) == 1) && no_dim
    res = a(1);
    for i = 2:length(a)
        res = res .* a(i);
    end
else
    if no_dim
        dim = 1; % default is columnwise multiplication
    end
    switch dim
        case 1 % columnwise product
            res = a(1,:);            % initialize res with first row of a
            for i = 2:S_a(1)
                res = res .* a(i,:); % res := res .* j-th column of a
            end
        case 2 % rowwise product
            res = a(:,1);            % initialize res with first column of a
            for j = 2:S_a(2)
                res = res .* a(:,j); % res := res .* j-th column of a
            end
        otherwise
            error('maximal two dimensions for type taylormodel');
    end
end

end % function prod