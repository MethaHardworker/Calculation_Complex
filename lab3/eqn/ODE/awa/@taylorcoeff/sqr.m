function r = sqr(a)
%SQR  Taylor coefficient square  a.^2
%
%   r = sqr(a)  
%
% Corresponds to Lohner's function TKOEFF, case OPER = QUAD, file awa.p. 

% written  08/01/17     F. Buenger
% modified 01/18/18     F. Buenger  fully vectorized version added
% modified 06/18/18     F. Buenger  correction for order = 1 (K = 2) in case status = []

global INTLAB_AWA_VARS

status = INTLAB_AWA_VARS.STATUS;

if ~isempty(status)
    k = status + 1;
    status_ = status;
    treenr = INTLAB_AWA_VARS.TREENR;
    vertexnr = INTLAB_AWA_VARS.VERTEXNR + 1;
    INTLAB_AWA_VARS.VERTEXNR = vertexnr;    
else
    k = 0;
end

if k && treenr > 0
    r = INTLAB_AWA_VARS.TREE{treenr}{vertexnr};
    [m,n,K,trigger] = INTLAB_AWA_VARS.TREE{treenr+2}{vertexnr}{:};
else
    S_a = size(a);
    if (length(S_a) > 2)
        error('Maximally two dimensions are allowed for type taylorcoeff.');
    end    
    m = S_a(1);
    n = S_a(2);
    K = length(a(1).inf);
    if ~k
        status_ = 0:K-1; 
    end     
    numel_bound = 2; % If status is nonempty and if m*n = numel(a) < numel_bound, then the non-vectorized squaring is faster than vectorized squaring. Feel free to change this heuristic value!
    trigger = k && (m*n < numel_bound);
    %trigger = true;    % only for testing !!!
    %trigger = false;   % only for testing !!!
    if trigger
        r = a;
    else
        r.inf = zeros(m,n,K);
        r.sup = r.inf;
    end
end

if trigger
    
    for i = 1:m
        for j = 1:n
            a_ = a(i,j);
            if k < 2 % case status = [] (k = 0) or status = 0 (k = 1)
                r_ = a_;  % just cheap preallocation
            else
                r_ = r(i,j);
            end
            
            for s = status_
                t = s+1;
                if s == 0
                    x.inf = a_.inf(1);
                    x.sup = a_.sup(1);
                    x = iv_sqr(x);
                else
                    k_ = floor(t/2);
                    idx1 = (0:k_-1).';
                    idx2 = s-idx1;
                    x.inf = a_.inf(idx1+1);
                    x.sup = a_.sup(idx1+1);
                    y.inf = a_.inf(idx2+1);
                    y.sup = a_.sup(idx2+1);
                    x = iv_dotprod(x,y);
                    x.inf = 2*x.inf; % Multiplication by 2 causes no rounding errors.
                    x.sup = 2*x.sup;
                    if s == 2*floor(s/2)  % Check if s is even, i.e., divisible by 2.
                        y.inf = a_.inf(k_+1);
                        y.sup = a_.sup(k_+1);
                        y = iv_sqr(y);
                        x = iv_plus(x,y);
                    end
                end
                r_.inf(t) = x.inf;
                r_.sup(t) = x.sup;
            end
            r(i,j) = r_;
        end
    end
    r_tree = r;
    
else % vectorized version      
    
    dummy = a(1);
    a = tcoeff2tc(a);

    if k % Compute Taylor coefficients of order = status.
        if status == 0
            r = a;
            r_.inf = a.inf(:,:,1);
            r_.sup = a.sup(:,:,1);
            r_ = iv_sqr(r_);
            r.inf(:,:,1) = r_.inf;
            r.sup(:,:,1) = r_.sup;
        else % status >= 1
            k_ = floor(k/2);
            idx1 = (0:k_-1).';
            idx2 = status-idx1;
            r_.inf = a.inf(:,:,idx1+1);
            r_.sup = a.sup(:,:,idx1+1);
            hlp.inf = a.inf(:,:,idx2+1);
            hlp.sup = a.sup(:,:,idx2+1);
            r_ = iv_dotprod(r_,hlp,3);
            r_.inf = 2*r_.inf; % Multiplication by 2 causes no rounding errors.
            r_.sup = 2*r_.sup;
            if status == 2*floor(status/2)  % Check if status is even, i.e., divisible by 2.
                hlp.inf = a.inf(:,:,k_+1);
                hlp.sup = a.sup(:,:,k_+1);
                hlp = iv_sqr(hlp);
                r_ = iv_plus(r_,hlp);
            end
            r.inf(:,:,k) = r_.inf;
            r.sup(:,:,k) = r_.sup;
        end
        
    else % Compute Taylor coefficients of all orders.
        
        if K == 1
            r = iv_sqr(a);
        else % K >= 2
            e = 1e-30;
            if 1+e > 1 % fast check for rounding upwards
                rndold = 1;
            else
                rndold = getround;
                setround(1) % rounding upwards
            end
            k2 = 2:K;
            k_ = floor(k2/2);
            K_ = k_(end);      % = floor(K/2)
            if K > 2
                K__ = k_(end-1);   % = floor((K-1)/2)
            else % K == 2
                K__ = 0;
            end
            N = K_+K__;
            
            % Construct     idx1 = [ [1] [1], [1 2] [1 2], [1 2 3] [1 2 3], ....., [1 2 ... K_] [1 2 ... K_]]                      if K is odd
            %               idx1 = [ [1] [1], [1 2] [1 2], [1 2 3] [1 2 3], ....., [1 2 ... K_-1] [1 2 ... K_-1], [1 2 ... K_] ]   if K is even
            idx = tril(repmat(1:K_,K_,1))';
            idx1 = zeros(K_,2*K_);
            idx1(:,1:2:N) = idx;
            idx1(:,2:2:N) = idx(:,1:K__);
            idx1 = idx1(:);
            idx1 = idx1(idx1>0);
            
            % Construct     idx3 = [ [2] [3], [4 4] [5 5], [6 6 6] [7 7 7], ....., [(K-1)...(K-1)] [K...K]]                      if K is odd
            %               idx3 = [ [2] [3], [4 4] [5 5], [6 6 6] [7 7 7], ....., [(K-2)...(K_-2)] [(K-1)...(K-1)], [K...K] ]   if K is even
            idx = triu(repmat(1:K_,K_,1));
            idx3 = zeros(K_,2*K_);
            idx3(:,1:2:N) = 2*idx;
            idx3(:,2:2:N) = triu(2*idx(:,1:K__)+1);
            idx3 = idx3(:);
            idx3 = idx3(idx3>0);
            
            idx2 = idx3+1-idx1;   % idx2 = [ [2] [3], [4 3] [5 4], [6 5 4] [7 6 5], ....., [(K-1)...(K-K_)] [K...(K+1-K_)] ]                      if K is odd
                                  % idx2 = [ [2] [3], [4 3] [5 4], [6 5 4] [7 6 5], ....., [(K-2)...(K-1-K_)] [(K-1)...(K-K_)], [K...(K+1-K_)]]   if K is even
            
            r0.inf = a.inf(:,:,idx1);
            r0.sup = a.sup(:,:,idx1);
            hlp.inf = a.inf(:,:,idx2);
            hlp.sup = a.sup(:,:,idx2);
            r0 = iv_times(r0,hlp);
            r0.inf = 2*r0.inf; % Multiplication by 2 causes no rounding errors.
            r0.sup = 2*r0.sup;
            
            r_ = r0.sup;
            [x,y,z] = size(r_);
            %        M = full(sparse(1:z,idx3,1)); % Indicator matrix from group indices in idx3.
            M = full(sparse(1:z,idx3-1,1)); % Indicator matrix from group indices in idx3.
            r_ = reshape(r_,[],z); % Reshape into a 2d-matrix. "Taylor orders" are grouped along columns 1,...,z. The grouping is given in idx3.
            r0.sup = reshape(r_*M,x,y,[]); % Compute result and reshape. Recall that rounding is upwards and that M has 0,1-entries.
            
            r_ = -r0.inf; % Take negative to avoid switching the rounding mode to downwards for computing the lower bound.
            r_ = reshape(r_,[],z); % Reshape into a 2d-matrix. "Taylor orders" are grouped along columns 1,...,z. The grouping is given in idx3.
            r0.inf = -reshape(r_*M,x,y,[]); % Compute result and reshape. Recall that rounding is upwards and that M has 0,1-entries.
            
            if K == 2*K_  % Check if K is even, i.e., divisible by 2.
                m = K_;
            else
                m = K_+1;
            end
            hlp.inf = a.inf(:,:,1:m);
            hlp.sup = a.sup(:,:,1:m);
            hlp = iv_sqr(hlp);
            
            r = a; % just storage preallocation for result
            
            r.inf(:,:,1) = hlp.inf(:,:,1); % simple square of function values (= coefficients of Taylor order 0)
            r.sup(:,:,1) = hlp.sup(:,:,1);
            
            r.inf(:,:,2:K) = r0.inf; % Result for Taylor orders > 0 is initialized by r0. For odd Taylor orders > 1 this is already the final result.
            r.sup(:,:,2:K) = r0.sup;
            
            % Even Taylor orders get the additional addend hlp.
            r.inf(:,:,3:2:K) = -(-r.inf(:,:,3:2:K)-hlp.inf(:,:,2:m)); % Recall that rounding is upwards.
            r.sup(:,:,3:2:K) = r.sup(:,:,3:2:K) + hlp.sup(:,:,2:m);
            
            if rndold ~= 1
                setround(rndold)
            end
        end
        
    end
    
    r_tree = r;
    r = tc2tcoeff(r_tree,dummy(1));        
    
end

if k
    if treenr > 0
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = r_tree; % Store result.
    else % negative treenr signalizes that the function tree is traversed for the first time and that meta data shall be stored once.
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {m,n,K,trigger}; % Store meta data.
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = r_tree; % Store result.
    end
end

end % function sqr