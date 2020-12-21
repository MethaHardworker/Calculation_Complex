function disp(x,newline)
%DISP  display function for pop-up windows in debugger
%

% written  08/27/15     F. Buenger

if nargin < 2 || isempty(newline)
    newline = false;
end
display(x,[],newline);

end % function disp