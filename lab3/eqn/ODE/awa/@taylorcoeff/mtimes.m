function r = mtimes(a,b)
%MTIMES  multiplication  a * b   of Taylor coefficients 
%
%   r = mtimes(a,b)
%
% Multiplication of the (mxl)-matrix a and the (lxn)-matrix b is, for full 
% vectorization without any loops, rewritten by using only pointwise 
% multiplication ".*" and summation "sum" for Taylor coefficients as 
% follows:
%
%   a*b = reshape ( sum( ( [a(1,:),...,a(m,:)]' .* [b;b;...;b] )' ) , m,n) 
% 
% For convenience, the factors 
%
%   a_ := [a(1,:),...,a(m,:)]' and b_ := [b;b;...;b]
%
% in the commutative(!) pointwise product a_ .* b_ are switched, if b is a 
% Taylor coeficient matrix and a is not.

% written  08/02/17     F. Buenger 
% modified 01/17/18     F. Buenger  full vectorization 

global INTLAB_AWA_VARS

e = 1e-30;
if 1+e > 1 % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

status = INTLAB_AWA_VARS.STATUS;

if ~isempty(status)
    k = status + 1;
    treenr = INTLAB_AWA_VARS.TREENR;
    vertexnr = INTLAB_AWA_VARS.VERTEXNR + 1;
    INTLAB_AWA_VARS.VERTEXNR = vertexnr;    
else
    k = 0;
end

if k && treenr > 0
    [type_a,type_b,x_,m,l,n,K] = INTLAB_AWA_VARS.TREE{treenr+2}{vertexnr}{:};
    if l == 1 % => matrix multiplication "*" same as pointwise multiplication ".*"
        r = a.*b;
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = r;
        return;
    end
    r = INTLAB_AWA_VARS.TREE{treenr}{vertexnr};
else  
    [~,type_a,type_b,x_,m,l,~,n,~,~,K] = get_type2(a,b,'mtimes');
    if k && l == 1 % => matrix multiplication "*" same as pointwise multiplication ".*"
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {type_a,type_b,x_,m,l,n,K}; % Store meta data.
        r = a.*b;
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = r; % Store result.
        return        
    end
end

switch type_b
    case 0 % Taylor coefficient
        dummy = b(1);
        b = tcoeff2tc(b);        
        b.inf = repmat(b.inf,m,1);
        b.sup = repmat(b.sup,m,1);
    case 1 % interval
        b = x_;
        b.inf = repmat(b.inf,m,1);    
        b.sup = repmat(b.sup,m,1);    
    case 2 % float
        b = repmat(b,m,1);
end

switch type_a
    case 0 % Taylor coefficient
        dummy = a(1);
        a = tcoeff2tc(a);        
        a.inf = reshape(permute(a.inf,[2 1 3]),m*l,1,[]);
        a.sup = reshape(permute(a.sup,[2 1 3]),m*l,1,[]);
        type = type_b;
    case 1 % interval
        a = x_;
        a.inf = a.inf';
        a.inf = a.inf(:);
        a.sup = a.sup';
        a.sup = a.sup(:);
        type = type_a;
        a_ = a; % switch a and b
        a = b;
        b = a_;
    case 2 % float
        a = a';
        a = a(:);    
        type = type_a;
        a_ = a; % switch a and b
        a = b;
        b = a_;
end

if k % Compute Taylor coefficients of order = status.    
    switch type
        case 0 % right factor "b" is a Taylor coefficient
            a_.inf = a.inf(:,:,1:k);
            a_.sup = a.sup(:,:,1:k);
            b_.inf = b.inf(:,:,k:-1:1);
            b_.sup = b.sup(:,:,k:-1:1);
            r_ = iv_dotprod(a_,b_,3);            
        case {1,2} % right factor "b" is a float or interval
            r_.inf = a.inf(:,:,k);
            r_.sup = a.sup(:,:,k);
            r_ = iv_times(r_,b);
    end   
    if status == 0
        r = a; % just storage preallocation for r
        r.inf = zeros([1,m*n,K]);
        r.sup = r.inf;
    end  
    r_.inf = reshape(r_.inf,l,m*n);
    r_.sup = reshape(r_.sup,l,m*n);
    r.inf(:,:,k) = -sum(-r_.inf,1); % Recall that rounding is upwards.
    r.sup(:,:,k) = sum(r_.sup,1);        
    
else % Compute Taylor coefficients of all orders.    
    
    switch type
        case 0 % right factor "b" is a Taylor coefficient                      
            K2 = K^2;
            
            idx1 = reshape(tril(repmat(1:K,K,1))',K2,1);
            idx1 = idx1(idx1 > 0);  % idx1 = [1, 1 2 , 1 2 3 , 1 2 3 4 , .... , 1 2 ...K]'
            
            idx3 = reshape(tril(repmat((1:K)',1,K))',K2,1);
            idx3 = idx3(idx3 > 0);  % idx3 = [1, 2 2 , 3 3 3 , 4 4 4 4 , .... , K K ...K]'  grouping of Taylor orders  
            
            idx2 = idx3-idx1+1;     % idx2 = [1, 2 1 , 3 2 1 , 4 3 2 1 , .... , K K-1 ... 1]'
            %idx2 = reshape(tril(toeplitz(1:K))',K2,1); 
            %idx2 = idx2(idx2 > 0); % idx2 = [1, 2 1 , 3 2 1 , 4 3 2 1 , .... , K K-1 ... 1]'
            
            a.inf = a.inf(:,:,idx1);
            a.sup = a.sup(:,:,idx1);
            b.inf = b.inf(:,:,idx2);
            b.sup = b.sup(:,:,idx2);
            r = iv_times(a,b);
            
            r_ = r.sup;            
            [x,y,z] = size(r_);                        
            M = full(sparse(1:z,idx3,1)); % Indicator matrix from group indices in idx3.
            %E = eye(n);
            %M = E(idx3,1:K); % Indicator matrix from group indices in idx3.
            r_ = reshape(r_,[],z); % Reshape into a 2d-matrix. "Taylor orders" are grouped along columns 1,...,z. The grouping is given in idx3.
            r.sup = reshape(r_*M,x,y,[]); % Compute result and reshape. Recall that rounding is upwards and that M has 0,1-entries.
            
            r_ = -r.inf; % Take negative to avoid switching the rounding mode to downwards for computing the lower bound.
            r_ = reshape(r_,[],z); % Reshape into a 2d-matrix. "Taylor orders" are grouped along columns 1,...,z. The grouping is given in idx3.
            r.inf = -reshape(r_*M,x,y,[]); % Compute result and reshape. Recall that rounding is upwards and that M has 0,1-entries.                                   
            
        case {1,2} % right factor "b" is a float or interval
            r = iv_times(a,b);
    end    
    
    r.inf = reshape(r.inf,l,m*n,[]);
    r.sup = reshape(r.sup,l,m*n,[]);

    r.inf = -sum(-r.inf,1); % Recall that rounding is upwards 
    r.sup = sum(r.sup,1);
    
end

if k
    if treenr > 0
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = r; % Store result.
    else
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {type_a,type_b,x_,m,l,n,K}; % Store meta data.
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = r; % Store result.
    end
end

r.inf = reshape(r.inf,m,n,[]);
r.sup = reshape(r.sup,m,n,[]);
r = tc2tcoeff(r,dummy);

if rndold ~= 1
    setround(rndold)
end

end % function mtimes
  