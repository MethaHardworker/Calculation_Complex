function r = get_interval(a)
%GET_INTERVAL   returns interval matrix a(i,j).interval of Taylor model matrix "a".
%
%   r = get_interval(a)

% written  02/04/16     F. Buenger
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures 

S = size(a);
r_inf = zeros(S);
r_sup = r_inf;
for i = 1:S(1)
    for j = 1:S(2)
        x =  a(i,j).interval;
        % [x_inf,x_sup] = getinfsup(x);
        % faster implementation than the previous line in comments        
        if ~(isempty(x.inf) || isempty(x.sup)) 
            r_inf(i,j) = x.inf; 
            r_sup(i,j) = x.sup; 
        else
            r = x; % x is an empty interval
            return;
        end
    end
end
r.inf = r_inf;
r.sup = r_sup;

end % function get_interval
