function r = mtimes(a,b)
%MTIMES  Taylor model matrix multiplication a * b
%
%   r = mtimes(a,b)

% written  11/23/15     F. Buenger
% modified 05/24/18     F. Buenger,

if isa(a,'taylormodel')
    type_a = 0; % a is Taylor model
    if isa(b,'taylormodel')
        type_b = 0;
    elseif isiv(b)
        type_b = 1;
    elseif isa(b,'intval')
        type_b = 1;
        b = intval2iv(b);
    elseif isfloat(b)
        type_b = 2;
    else
        error('incompatible type')
    end
else
    type_b = 0; % b is Taylor model
    if isiv(a)
        type_a = 1;
    elseif isa(a,'intval')
        type_a = 1;
        a = intval2iv(a);
    elseif isfloat(a)
        type_a = 0;
    else
        error('incompatible type')
    end
end
if type_a == 1
    S_a = size(a.inf);
else
    S_a = size(a);
end
if type_b == 1
    S_b = size(b.inf);
else
    S_b = size(b);
end

if (length(S_a) > 2) || (length(S_b) > 2)
    error('maximally two dimensions for type taylormodel');
end

m = S_a(1);
na = S_a(2);
mb = S_b(1);
n = S_b(2);

if (m == 1 && na == 1) || (mb == 1 && n == 1) % Case "singleton * matrix" or vice versa.
    r = a.*b;                                 % in this case a*b = a.*b
else
    if na ~= mb
        error('matrix sizes not compatible for matrix multiplication');
    end
    k = na;
    switch type_b
        case 0 % Taylor model
            b = repmat(b,m,1);
        case 1 % interval
            b.inf = repmat(b.inf,m,1);
            b.sup = repmat(b.sup,m,1);
        case 2 % float
            b = repmat(b,m,1);
    end    
    switch type_a
        case 0 % Taylor coefficient
            a = a';
            a = a(:);
        case 1 % interval
            a.inf = a.inf';
            a.inf = a.inf(:);
            a.sup = a.sup';
            a.sup = a.sup(:);
            x = a; % switch a and b
            a = b;
            b = x;
        case 2 % float
            a = a';
            a = a(:);
            x = a; % switch a and b
            a = b;
            b = x;
    end
    r = reshape(a.*b,k,m*n);
    r = reshape(sum(r,1),m,n);
    
end
end % function mtimes
  