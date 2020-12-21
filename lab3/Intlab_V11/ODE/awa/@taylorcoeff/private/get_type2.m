function [r,type_a,type_b,x,m_a,n_a,m_b,n_b,m,n,K] = get_type2(a,b,funname)
%GET_TYPE2  determines meta data of input a and b of a binary function "funname".
%
%   [r,type_a,type_b,x,m_a,n_a,m_b,n_b,m,n] = get_type2(a,b)
%
% It is supposed that a,b have type float, iv, intval or tcoeff and that
% a or b has type tcoeff.
%
%  0 : type(b) = tcoeff
%  1 : type(b) = iv
%  2 : type(b) = float

% written  12/11/17     F. Buenger

if isa(a,'taylorcoeff')     
    type_a = 0;
    S_a = size(a);  
    K = length(a(1).inf);
    if isa(b,'taylorcoeff')
        type_b = 0;
        S_b = size(b);
        x = [];
    elseif isfloat(b) && isreal(b)
        type_b = 2;
        S_b = size(b);
        x.inf = b;
        x.sup = b;
    elseif isa(b,'intval')
        type_b = 1;
        x = intval2iv(b); % Convert intval to iv-structure.
        S_b = size(x.inf);
    elseif isiv(b)
        type_b = 1;
        S_b = size(b.inf);
        x = b;
    else
        error('type combination not supported')
    end    
else
    type_b = 0; 
    S_b = size(b);
    K = length(b(1).inf);
    if isa(a,'taylorcoeff')
        type_a = 0;
        S_a = size(a);
        x = [];
    elseif isfloat(a) && isreal(a)
        type_a = 2;
        S_a = size(a);
        x.inf = a;
        x.sup = a;
    elseif isa(a,'intval')
        type_a = 1;
        x = intval2iv(a); % Convert intval to iv-structure.
        S_a = size(a.inf);
    elseif isiv(a)
        type_a = 1;
        S_a = size(a.inf);
        x = a;
    else
        error('type combination not supported')
    end    
end

if (length(S_a) > 2) || (length(S_b) > 2)
    error('Maximally two dimensions are allowed in taylorcoeff arithmetic.');
end
if (max(S_a) == 0 || max(S_b) == 0)
    error('empty input');
end

m_a = S_a(1); 
n_a = S_a(2); 
m_b = S_b(1); 
n_b = S_b(2); 

switch funname
    case 'mtimes'
        if n_a ~= m_b
            error('inconsistent dimensions for matrix multiplication');
        end
        m = m_a;
        n = n_b;
    otherwise
        m = max(m_a,m_b);
        n = max(n_a,n_b);
        m_ = min(m_a,m_b);
        n_ = min(n_a,n_b);
        if (m_a ~= m_b && m_ ~= 1) || (n_a ~= n_b && n_ ~= 1)
            error('inconsistent matrix dimensions');
        end
end
r(1:m,1:n) = a(1); % just storage preallocation 

end % function get_type2








