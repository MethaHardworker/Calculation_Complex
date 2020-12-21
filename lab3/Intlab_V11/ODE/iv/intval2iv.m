function res =  intval2iv(a)
%INTVAL2IV  convert intval to interval structure
%
%   s = intval2iv(a)

% written  02/09/16     F. Buenger
% modified 02/22/18     F. Buenger  float input is also considered 

if isfloat(a)
    res.inf = a;
    res.sup = a;    
else
    s = struct(a);
    res.inf = s.inf;
    res.sup = s.sup;
end

end % function intval2iv
