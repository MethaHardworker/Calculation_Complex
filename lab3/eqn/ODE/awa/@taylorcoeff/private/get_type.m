function [r,type,s,x,m_a,n_a,m_b,n_b,m,n,K] = get_type(a,b)
%GET_TYPE  determines meta data of a and b.
% 
%   [r,type,x,m_a,n_a,m_b,n_b,m,n] = get_type(a,b)
%
% It is supposed that a,b have type float, iv, intval or tcoeff and that 
% a or b has type tcoeff.
%
%  0 : type(b) = tcoeff         
%  1 : type(b) = iv         
%  2 : type(b) = float

% written  12/11/17     F. Buenger

s = ~isa(a,'taylorcoeff'); % switch a and b
if s
    c = a;
    a = b;
    b = c;
end
K = length(a(1).inf);
S_a = size(a);
if isa(b,'taylorcoeff')
    type = 0;
    S_b = size(b);
    x = [];
elseif isfloat(b) && isreal(b)
    type = 2;
    S_b = size(b);
    x.inf = b;
    x.sup = b;
elseif isa(b,'intval')
    type = 1;
    x = intval2iv(b); % Convert intval to iv-structure.
    S_b = size(x.inf);
elseif isiv(b)
    type = 1;
    S_b = size(b.inf);
    x = b;
else
    error('type combination not supported')
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

m = max(m_a,m_b);
n = max(n_a,n_b);

m_ = min(m_a,m_b);
n_ = min(n_a,n_b);

if (m_a ~= m_b && m_ ~= 1) || (n_a ~= n_b && n_ ~= 1)
    error('inconsistent matrix dimensions');
end

r(1:m,1:n) = a(1); % just storage preallocation 

end % function get_type








