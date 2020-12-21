function display(a,name)
%DISPLAY  command window display of Taylor coefficients
%    
%   display(a,name)

% written  07/31/17     F. Buenger

if nargin < 2 || isempty(name)
    name = inputname(1);
    if isempty(name)
        name = 'ans';
    end
end

disp(['taylorcoeff ' name ' = ' ])

if isempty(a) || isempty(a(1).inf)
    display([]);
    return
end

[m,n] = size(a);

if m == 1 || n == 1
    for i = 1:m*n
        disp(iv2intval(struct(a(i))).');
    end
else
    disp(iv2intval(tcoeff2tc(a)));
end
end  % function display

  
  