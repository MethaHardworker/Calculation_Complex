function p = taylormodel(q,str)
%taylormodel  class constructor for taylormodel
%
%  p = taylormodel(q,str)
%

% written  08/20/15     F. Buenger
% modified 01/11/16     F. Buenger  polynomial -> taylormodel

superiorto('intval');

switch nargin
    case 2
        if isequal(str,'taylormodelinit')         % call by taylormodelinit
            p = class(q,'taylormodel');
            [m,n] = size(q);
            for i = 1:m
                for j = 1:n
                    p(i,j).image = image(p(i,j)); % Compute image of Taylor model p(i,j).
                end
            end
        else
            error('invalid call of constructor taylormodel')
        end
    case 1
        if isa(q,'taylormodel')
            p = q;
        elseif isa(q,'polynom')
            p = polynom2taylormodel(q);
        elseif isfloat(q) && isreal(q)
            p = taylormodel(polynom(q));
        else
            error('invalid call of constructor taylormodel')
        end
    otherwise
        error('invalid call of constructor taylormodel')
end % switch nargin

end % function taylormodel