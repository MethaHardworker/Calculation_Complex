function r = times(a,b)
%TIMES Taylor coefficient pointwise multiplication  a .* b
%
%   r = times(a,b)  
%
% Corresponds to Lohner's function TKOEFF, case OPER = MAL, file awa.p. 
%
% Taylor coefficient recursion formula, see [L], (1.1.4.3), p. 18:    
%   
%   (a*b)_s = sum_{j=0}^s (a)_j (b)_{s-j}     (convolution)

% written  08/01/17     F. Buenger
% modified 01/17/18     F. Buenger  fully vectorized version added

global INTLAB_AWA_VARS

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
    r = INTLAB_AWA_VARS.TREE{treenr}{vertexnr};
    [type,s,x_,m_a,n_a,m_b,n_b,m,n,K,trigger] = INTLAB_AWA_VARS.TREE{treenr+2}{vertexnr}{:};
else
    [r,type,s,x_,m_a,n_a,m_b,n_b,m,n,K] = get_type(a,b);
    numel_bound = 3; % If status is nonempty and if m*n = numel(a.*b) <= numel_bound, then the non-vectorized multiplication is faster than vectorized multiplication. Feel free to change this heuristic value!
    trigger = k && (m*n <= numel_bound);
    %trigger = true;    % only for testing !!! 
    %trigger = false;   % only for testing !!!
    if ~trigger
        r0.inf = zeros(m,n,K);
        r0.sup = r0.inf;
        r = r0;
    end
end

if s  % switch a and b
    c = a;
    a = b;
    b = c;
end

if trigger
    
    for i = 1:m
        i_a = min(m_a,i);
        i_b = min(m_b,i);
        for j = 1:n
            j_a = min(n_a,j);
            j_b = min(n_b,j);
            a_ = a(i_a,j_a);
            if k < 2 % case status = [] (k = 0) or status = 0 (k = 1)
                r_ = a_;  % just cheap preallocation
            else
                r_ = r(i,j);
            end
            
            if type % tcoeff .* (float or interval)
                if type  == 1 % interval type
                    b_.inf = x_.inf(i_b,j_b);
                    b_.sup = x_.sup(i_b,j_b);
                else % float
                    b_ = b(i_b,j_b);
                end
                if k % Compute Taylor coefficients of order status.
                    x.inf = a_.inf(k);
                    x.sup = a_.sup(k);
                    x = iv_times(x,b_);
                    r_.inf(k) = x.inf;
                    r_.sup(k) = x.sup;
                else % Compute Taylor coefficients of all orders.
                    x.inf = a_.inf;
                    x.sup = a_.sup;
                    x = iv_times(x,b_);
                    r_.inf  = x.inf;
                    r_.sup  = x.sup;
                end
            else % tcoeff .* tcoeff
                b_ = b(i_b,j_b);
                if k % Compute Taylor coefficients of order status.
                    x.inf = a_.inf(1:k);
                    x.sup = a_.sup(1:k);
                    y.inf = b_.inf(k:-1:1);
                    y.sup = b_.sup(k:-1:1);
                    x = iv_dotprod(x,y);
                    r_.inf(k) = x.inf;
                    r_.sup(k) = x.sup;
                else % Compute Taylor coefficients of all orders.
                    for kk = 1:K
                        x.inf = a_.inf(1:kk);
                        x.sup = a_.sup(1:kk);
                        y.inf = b_.inf(kk:-1:1);
                        y.sup = b_.sup(kk:-1:1);
                        x = iv_dotprod(x,y);
                        r_.inf(kk) = x.inf;
                        r_.sup(kk) = x.sup;
                    end
                end
            end
            r(i,j) = r_;
        end
    end
    r_tree = r;
    
else
    
    dummy = a(1);
    a = tcoeff2tc(a);
    if type == 0
        b = tcoeff2tc(b);   
    elseif type == 1
        b = x_;    
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
        r.inf(:,:,k) = r_.inf;
        r.sup(:,:,k) = r_.sup;
        
    else % Compute Taylor coefficients of all orders.
        
        switch type
            case 0 % right factor "b" is a Taylor coefficient
                e = 1e-30;
                if 1+e > 1 % fast check for rounding upwards
                    rndold = 1;
                else
                    rndold = getround;
                    setround(1) % rounding upwards
                end
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
                
                if rndold ~= 1
                    setround(rndold)
                end
                
            case {1,2} % right factor "b" is a float or interval
                r = iv_times(a,b);
        end
        
    end
    
    r_tree = r;  
    r = tc2tcoeff(r_tree,dummy);
    
end

if k
    if treenr > 0
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = r_tree; % Store result.
    else % negative treenr signalizes that the function tree is traversed for the first time and that meta data shall be stored once.
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {type,s,x_,m_a,n_a,m_b,n_b,m,n,K,trigger}; % Store meta data.
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = r_tree; % Store result.
    end
end

end % function times