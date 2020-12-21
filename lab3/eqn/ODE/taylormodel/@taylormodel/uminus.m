function r = uminus(a)
%UMINUS  Taylor model unary minus -a
%
%   r = uminus(a)

% written  08/24/15     F. Buenger
% modified 11/23/15     F. Buenger  matrix input 
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures

% Remark: Considering the nonverified ODEMODE = 0 and/or record feature write/read mode is not relevant for uminus.

S_a = size(a);

r = a; % Initialize result r with a. 
for i = 1:S_a(1)
    for j = 1:S_a(2)          
        r_ = r(i,j);  
        
        r_.coefficient = -r_.coefficient;
        
        x = r_.interval; 
        r_.interval.inf = -x.sup;
        r_.interval.sup = -x.inf;
        
        x = r_.image;
        r_.image.inf = -x.sup;
        r_.image.sup = -x.inf; 
        
        r(i,j) = r_;
    end
end

end % function uminus