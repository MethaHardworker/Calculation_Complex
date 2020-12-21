function r = minus(a,b)
%MINUS  Taylor coefficient subtraction  a - b
%
%   r = minus(a,b)
%
% Taylor coefficient recursion formula, see [L], (1.1.4.3), p. 18:    
%   
%   (a-b)_s = (a)_s - (b)_s   

% written  08/01/17     F. Buenger
% modified 01/05/18     F. Buenger  fully vectorized version added

global INTLAB_AWA_VARS

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
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
    r = INTLAB_AWA_VARS.TREE{treenr}{vertexnr};
    [type,s,x,m_a,n_a,m_b,n_b,m,n,K,trigger] = INTLAB_AWA_VARS.TREE{treenr+2}{vertexnr}{:};
else
    [r,type,s,x,m_a,n_a,m_b,n_b,m,n,K] = get_type(a,b);
    numel_bound = 10; % If m*n = numel(a+b) < numel_bound, the non-vectorized summation is faster than vectorized summation. Feel free to change this heuristic value!
    trigger = (m*n <= numel_bound);
    %trigger = true;    % only for testing !!! 
    % trigger = false;   % only for testing !!!    
    if ~trigger       
        r0.inf = zeros(m,n,K);
        r0.sup = r0.inf;
        r = r0;
    end
end

if s  % Switch a and b. (I know that subtraction is not commutative, but still a-b = -(b-a) holds true.) 
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
            if type % tcoeff - (float or interval)
                b_.inf = x.inf(i_b,j_b);
                b_.sup = x.sup(i_b,j_b);
                if k % Compute Taylor coefficients of order status.
                    if status == 0 
                        if s
                            r_.sup(1) = b_.sup - a_.inf(1); % Recall that rounding is upwards.
                            r_.inf(1) = -(a_.sup(1) - b_.inf); 
                        else
                            r_.sup(1) = a_.sup(1) - b_.inf; % Recall that rounding is upwards.
                            r_.inf(1) = -(b_.sup - a_.inf(1));
                        end
                    else 
                        if s
                            r_.sup(k) = -a_.inf(k);
                            r_.inf(k) = -a_.sup(k);
                        else
                            r_.sup(k) = a_.sup(k);
                            r_.inf(k) = a_.inf(k);
                        end
                    end
                else % Compute Taylor coefficients of all orders.
                    if s
                        r_.sup = -a_.inf;
                        r_.inf = -a_.sup;
                        r_.sup(1) = b_.sup - a_.inf(1);    % Recall that rounding is upwards. 
                        r_.inf(1) = -(a_.sup(1) - b_.inf);
                    else
                        r_.sup(1) = a_.sup(1) - b_.inf; % Recall that rounding is upwards.
                        r_.inf(1) = -(b_.sup - a_.inf(1));
                    end
                end
            else % tcoeff - tcoeff
                b_ = b(i_b,j_b);
                if k % Compute Taylor coefficients of order status.
                    r_.sup(k) = a_.sup(k) - b_.inf(k); % Recall that rounding is upwards.
                    r_.inf(k) = -(b_.sup(k) - a_.inf(k));
                else % Compute Taylor coefficients of all orders.
                    r_.sup = a_.sup - b_.inf; % Recall that rounding is upwards.
                    r_.inf = -(b_.sup - a_.inf);
                end
            end
            r(i,j) = r_;
        end
    end   
    r_tree = r;
    
else % fully vectorized subtraction
    
    dummy = a(1);
    a = tcoeff2tc(a);
    if type == 0
        b = tcoeff2tc(b);  
    end
    
    if k % Compute Taylor coefficients of order = status
        switch type
            case 0 % right addend "b" is a Taylor coefficient
                r_inf = -( b.sup(:,:,k) - a.inf(:,:,k) );  % Recall that rounding is upwards.
                r_sup = a.sup(:,:,k) - b.inf(:,:,k);
                r.inf(:,:,k) = r_inf;
                r.sup(:,:,k) = r_sup;
            case 1 % right addend "b" is an interval
                if status == 0
                    r_inf = -( b.sup - a.inf(:,:,k) ); % Recall that rounding is upwards.
                    r_sup = a.sup(:,:,k) - b.inf;
                    r.inf(:,:,1) = r_inf;
                    r.sup(:,:,1) = r_sup;
                else
                    r.inf(:,:,k) = repmat(a.inf(:,:,k),m/m_a,n/n_a);
                    r.sup(:,:,k) = repmat(a.sup(:,:,k),m/m_a,n/n_a);
                end
            case 2 % right addend "b" is a float
                if status == 0
                    r_inf = -( b - a.inf(:,:,k) ); % Recall that rounding is upwards.
                    r_sup = a.sup(:,:,k) - b;
                    r.inf(:,:,1) = r_inf;
                    r.sup(:,:,1) = r_sup;
                else
                    r.inf(:,:,k) = repmat(a.inf(:,:,k),m/m_a,n/n_a);
                    r.sup(:,:,k) = repmat(a.sup(:,:,k),m/m_a,n/n_a);
                end
        end
    else % Compute Taylor coefficients of all orders
        switch type
            case 0 % right addend "b" is a Taylor coeffficient
                r.inf = -( b.sup - a.inf );  % Recall that rounding is upwards.
                r.sup = a.sup - b.inf;
            case 1 % right addend "b" is an interval
                r_inf = -( b.sup - a.inf(:,:,1) ); % Recall that rounding is upwards.
                r_sup = a.sup(:,:,1) - b.inf;
                r.inf(:,:,1) = r_inf;
                r.sup(:,:,1) = r_sup;
                r.inf(:,:,2:K) = repmat(a.inf(:,:,2:K),m/m_a,n/n_a,1);
                r.sup(:,:,2:K) = repmat(a.sup(:,:,2:K),m/m_a,n/n_a,1);
            case 2 % right addend "b" is a float
                r_inf = -( b - a.inf(:,:,1) ); % Recall that rounding is upwards.
                r_sup = a.sup(:,:,1) - b;
                r.inf(:,:,1) = r_inf;
                r.sup(:,:,1) = r_sup;
                r.inf(:,:,2:K) = repmat(a.inf(:,:,2:K),m/m_a,n/n_a,1);
                r.sup(:,:,2:K) = repmat(a.sup(:,:,2:K),m/m_a,n/n_a,1);
        end
    end
    
    r_tree = r;
    if s 
        if k
            r_tree_inf = r_tree.inf(:,:,k);
            r_tree.inf(:,:,k) = -r_tree.sup(:,:,k);
            r_tree.sup(:,:,k) = -r_tree_inf;
        else
            r_tree_inf = r_tree.inf;
            r_tree.inf = -r_tree.sup;
            r_tree.sup = -r_tree_inf;
        end
    end
    r = tc2tcoeff(r_tree,dummy);  
    
end

if k
    if treenr > 0        
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = r_tree; % Store result.
    else % negative treenr signalizes that the function tree is traversed for the first time and that meta data shall be stored once.
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {type,s,x,m_a,n_a,m_b,n_b,m,n,K,trigger}; % Store meta data.
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = r_tree; % Store result.
    end
end

if rndold ~= 1
    setround(rndold)
end

end % function minus
