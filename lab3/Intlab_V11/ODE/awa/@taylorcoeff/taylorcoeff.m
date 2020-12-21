function p = taylorcoeff(q,str)
%TAYLORCOEFF  class constructor for taylorcoeff
%
%  p = taylorcoeff(q,str)

% written  07/30/15     F. Buenger

superiorto('intval');

if nargin == 2 && isequal(str,'taylorcoeffinit')
    p = class(q,'taylorcoeff');
else
    error('invalid call of constructor taylorcoeff')
end

end % function taylorcoeff