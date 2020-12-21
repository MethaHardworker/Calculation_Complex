function r = rdivide(a,b)
%RDIVIDE  Taylor coefficient division a ./ b
%
%   r = rdivide(a,b)  
%
% Corresponds to Lohner's function TKOEFF, case OPER = DURCH, file awa.p. 
%
% Taylor coefficient recursion formula, see [L], (1.1.4.3), p. 18:    
%   
%   (a./b)_s = [(a)_s-sum_{j=1}^s (b)_j (a/b)_{s-j}] / b  

% written  08/01/17     F. Buenger
% modified 01/23/18     F. Buenger  fully vectorized version added

global INTLAB_AWA_VARS

status = INTLAB_AWA_VARS.STATUS;

if ~isempty(status)
    k = status + 1;
    L = k;
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
    if ~k
         L = 1:K; % L-1 is the order of the Taylor expansion, i.e., index of highest derivative.
    end
    numel_bound = 2; % If status is nonempty and if m*n = numel(a.*b) < numel_bound, then the non-vectorized division is faster than vectorized multiplication. Feel free to change this heuristic value!
    trigger = k && (m*n <= numel_bound);
    %trigger = true;    % only for testing !!! 
    %trigger = false;   % only for testing !!!
    if ~trigger
        r0.inf = zeros(m,n,K);
        r0.sup = r0.inf;
        r = r0;
    end
end

if trigger
    
    if s  % Switch a and b. (I know that division is not commutative.)
        c = a;
        a = b;
        b = c;
        if type
            zeroiv.inf = 0;
            zeroiv.sup = 0;
        end
    end   
            
    for i = 1:m
        i_a = min(m_a,i);
        i_b = min(m_b,i);
        for j = 1:n
            j_a = min(n_a,j);
            j_b = min(n_b,j);
            a_ = a(i_a,j_a);
            if type == 1 % b is (constant) interval
                b_.inf = x_.inf(i_b,j_b);
                b_.sup = x_.sup(i_b,j_b);
            else
                b_ = b(i_b,j_b);
            end
            if k < 2 % case status = [] (k = 0) or status = 0 (k = 1)
                r_ = a_;  % just cheap preallocation
            else
                r_ = r(i,j);
            end
            
            for l = L
                if ~s && type            % case:  tcoeff ./ (interval or float)
                    u.inf = a_.inf(l);
                    u.sup = a_.sup(l);
                    x = iv_rdivide(u,b_);
                else
                    if type              % case:  (interval or float) ./ tcoeff
                        if l == 1
                            uk = b_;
                        else
                            uk = zeroiv;
                        end
                        v1.inf = a_.inf(1);
                        v1.sup = a_.sup(1);
                        x.inf = a_.inf(2:l);
                        x.sup = a_.sup(2:l);
                    else                 % case:  tcoeff ./ tcoeff
                        uk.inf = a_.inf(l);
                        uk.sup = a_.sup(l);
                        v1.inf = b_.inf(1);
                        v1.sup = b_.sup(1);
                        x.inf = b_.inf(2:l);
                        x.sup = b_.sup(2:l);
                    end
                    if l > 1
                        y.inf = r_.inf(l-1:-1:1);
                        y.sup = r_.sup(l-1:-1:1);                        
                        x = iv_dotprod(x,y);
                        x = iv_minus(uk,x);
                    else
                        x = uk;
                    end
                    x = iv_rdivide(x,v1);
                end
                r_.inf(l) = x.inf;
                r_.sup(l) = x.sup;
            end
            
            r(i,j) = r_;
        end
    end
    r_tree = r;
    
else % vectorized division
    
    if s
        dummy = b(1);
        b = tcoeff2tc(b);
        if type == 1
            a = x_;
        end
        type_a = type;
        type_b = 0;
    else
        dummy = a(1);
        a = tcoeff2tc(a);
        if type == 0
            b = tcoeff2tc(b);
        elseif type == 1
            b = x_;
        end
        type_a = 0;
        type_b = type;
    end        
    
    if k % Compute Taylor coefficients of order = status.        
        
        if type_b                         %  tcoeff / (float or interval)
            a.inf = a.inf(:,:,k);
            a.sup = a.sup(:,:,k);
            w = iv_rdivide(a,b);
        else
            v1.inf = b.inf(:,:,1);
            v1.sup = b.sup(:,:,1);
            if k == 1
                if type_a
                    w = a;
                else
                    w.inf = a.inf(:,:,k);
                    w.sup = a.sup(:,:,k);
                end
            else
                v.inf = b.inf(:,:,2:k);
                v.sup = b.sup(:,:,2:k);
                
                r_.inf = r.inf(:,:,k-1:-1:1);
                r_.sup = r.sup(:,:,k-1:-1:1);
                
                w = iv_dotprod(v,r_,3);
                switch type_a
                    case 0                          %  tcoeff / tcoeff
                        uk.inf = a.inf(:,:,k);
                        uk.sup = a.sup(:,:,k);
                    case 1                          %  interval / tcoeff
                        uk = zeros(m_a,n_a);
                    case 2                          %  float / tcoeff
                        uk = zeros(m_a,n_a);
                end
                w = iv_minus(uk,w);
            end
            w = iv_rdivide(w,v1);
        end
        r.inf(:,:,k) = w.inf;
        r.sup(:,:,k) = w.sup;
        
    else % Compute Taylor coefficients of all orders.
        
        switch type_b            
            case 0
                for k_ = 1:K
                    if type_b                         %  tcoeff / (float or interval)
                        a.inf = a.inf(:,:,k_);
                        a.sup = a.sup(:,:,k_);
                        w = iv_rdivide(a,b);
                    else
                        v1.inf = b.inf(:,:,1);
                        v1.sup = b.sup(:,:,1);
                        if k_ == 1
                            if type_a
                                w = a;
                            else
                                w.inf = a.inf(:,:,k_);
                                w.sup = a.sup(:,:,k_);
                            end
                        else
                            v.inf = b.inf(:,:,2:k_);
                            v.sup = b.sup(:,:,2:k_);
                            
                            r_.inf = r.inf(:,:,k_-1:-1:1);
                            r_.sup = r.sup(:,:,k_-1:-1:1);
                            
                            w = iv_dotprod(v,r_,3);
                            switch type_a
                                case 0                          %  tcoeff / tcoeff
                                    uk.inf = a.inf(:,:,k_);
                                    uk.sup = a.sup(:,:,k_);
                                case 1                          %  interval / tcoeff
                                    uk = zeros(m_a,n_a);
                                case 2                          %  float / tcoeff
                                    uk = zeros(m_a,n_a);
                            end
                            w = iv_minus(uk,w);
                        end
                        w = iv_rdivide(w,v1);
                    end
                    r.inf(:,:,k_) = w.inf;
                    r.sup(:,:,k_) = w.sup;
                end                
            case 1                        %  tcoeff / interval                
                b.inf = repmat (b.inf,[1,1,K]);
                b.sup = repmat (b.sup,[1,1,K]);
                r = iv_rdivide(a,b);                
            case 2                        %  tcoeff / float                
                b = repmat (b,[1,1,K]);
                r = iv_rdivide(a,b);
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

end % function rdivide