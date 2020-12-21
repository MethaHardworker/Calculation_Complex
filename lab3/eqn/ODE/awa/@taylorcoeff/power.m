function r = power(a,b)
%POWER  Taylor coefficient power a .^ b
%
%   r = power(a,b)  
%
% For real floating-point or interval b the implementation corresponds to 
% [AWA] function TKOEFF, case OPER = HOCH, file awa.p. 
% For Taylor coefficient b the result is t = exp(b.*log(a)).
% 
% Different to [AWA], for non-negative integer b, we implemented a.^b by  
% squaring and multiplication since the general recursion formula used by 
% [AWA] becomes undefined or ill-conditioned if some function value interval 
% [a(i,j).inf(1),a(i,j).sup(1)] includes zero or is "close" to zero.
%
% Taylor coefficient recursion formula, see [L], (1.1.4.10), p. 20:    
%   
%   (a^b)_s = 1/(s*a) * sum_{j=0}^{s-1} (b*(s-j)-j) * (a)_{s-j} * (a^b)_j

% written  08/02/17     F. Buenger
% modified 01/26/18     F. Buenger  vectorized version added
% modified 02/01/18     F. Buenger  nonnegative integer power implemented by squaring and multiplication 

global INTLAB_AWA_VARS

status = INTLAB_AWA_VARS.STATUS;

if ~isempty(status)
    status_ = status;
    k = status + 1;
    treenr = INTLAB_AWA_VARS.TREENR;
    vertexnr = INTLAB_AWA_VARS.VERTEXNR + 1;
    INTLAB_AWA_VARS.VERTEXNR = vertexnr;    
else
    k = 0;
end

if k && treenr > 0 
    
    [type_a,type_b,x,m_a,n_a,m_b,n_b,m,n,K,trigger,flag_int_power] = INTLAB_AWA_VARS.TREE{treenr+2}{vertexnr}{:};
    if ~type_b %  b is tcoeff => a.^b := exp(b.*log(a))
        if type_a == 1 
            r = exp(b.*log(x));
        else
            r = exp(b.*log(a));
        end
        return;
    elseif type_b == 2 && flag_int_power
        r = int_power(a,b);
        return;
    end
    r = INTLAB_AWA_VARS.TREE{treenr}{vertexnr};
    
else  % The "else-part" is only executed once for a function tree. Thus, its performance is unimportant.
    
    [r,type_a,type_b,x,m_a,n_a,m_b,n_b,m,n,K] = get_type2(a,b,'power');
    
    if ~type_b          %  b is tcoeff => a.^b := exp(b.*log(a))        
        if type_a == 1
            x = iv2intval(x);
            a = x;
        end
        if k 
            INTLAB_AWA_VARS.TREE{abs(treenr)}{vertexnr} = {[]};
            INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {type_a,type_b,x,m_a,n_a,m_b,n_b,m,n,K,false,false}; % Store meta data.        
        end
        r = exp(b.*log(a));                     
        return;  
    elseif type_b == 2 % b is float.
        flag_int_power = ( m_b*n_b == 1 && b >= 0 && b == round(b) );   % b is non-negative integer
        %flag_int_power = false; % only for testing !!!
        if flag_int_power
            if k % only k = 0 and k = 1 are possible since treenr < 0 implies k = 1. 
                INTLAB_AWA_VARS.TREE{abs(treenr)}{vertexnr} = {[]};
                INTLAB_AWA_VARS.TREE{abs(treenr)+2}{vertexnr} = {type_a,type_b,x,m_a,n_a,m_b,n_b,m,n,K,false,flag_int_power}; % Store meta data.           
            end
            r = int_power(a,b);                     
            return; 
        end        
    end
    % b is interval or not a non-negative integer
    flag_int_power = false;
    numel_bound = 1;
    trigger = (m*n <= numel_bound);
    %trigger = true;    % only for testing !!!
    %trigger = false;   % only for testing !!!
    
    if ~trigger
        r0.inf = zeros(m,n,K);
        r0.sup = r0.inf;
        r = r0;
    end
    if ~k
        status_ = 0:K-1; 
    end
    
end

if trigger
    
    for i = 1:m
        i_a = min(m_a,i);
        i_b = min(m_b,i);
        
        for j = 1:n
            r_ = r(i,j);
            j_a = min(n_a,j);
            j_b = min(n_b,j);
            a_ = a(i_a,j_a);
            if type_b == 1                     % b is intval
                b_.inf = x.inf(i_b,j_b);
                b_.sup = x.sup(i_b,j_b);                
            else                               % b is float
                b_ = b(i_b,j_b);
            end
            a_0.inf = a_.inf(1);
            a_0.sup = a_.sup(1);
            
            for s = status_
                t = s+1;
                if s == 0
                    if type_b == 1
                        w = iv2intval(a_0).^iv2intval(b_) ; 
                    else
                        w = iv2intval(a_0).^b_ ;                         
                    end
                    if isnan(w) || isinf(w)
                        error('Inf/NaN occured in exponentiation')
                    end
                    w = intval2iv(w);
                else
                    idx1 = (0:s-1).';
                    idx2 = s-idx1;
                    u.inf = r_.inf(idx1+1);
                    u.sup = r_.sup(idx1+1);
                    v.inf = a_.inf(idx2+1);
                    v.sup = a_.sup(idx2+1);
                    c = iv_minus(iv_times(idx2,b_),idx1);        % c = idx2.*b_-idx1
                    u = iv_times(c,u);                           % u = c.*u
                    w = iv_dotprod(u,v);                         % w = sum(u.*v)
                    w = iv_rdivide(w,iv_times(s,a_0));           % w = w/(s.*a_0)
                end
                r_.inf(t) = w.inf;
                r_.sup(t) = w.sup;
            end
            
            r(i,j) = r_;
        end
    end
    r_tree = r;
    
else
    
    a_ = tcoeff2tc(a);    
    if type_b == 1                     % b is intval                        
        b_ = x; 
    else                               % b is float
        b_ = b;
    end
    a_0.inf = a_.inf(:,:,1);
    a_0.sup = a_.sup(:,:,1);
    for s = status_
        t = s+1;
        if s == 0
            if type_b == 1                
                a0.inf = repmat(a_0.inf,m/m_a,n/n_a);
                a0.sup = repmat(a_0.sup,m/m_a,n/n_a);
                b0.inf = repmat(b_.inf,m/m_b,n/n_b);
                b0.sup = repmat(b_.sup,m/m_b,n/n_b);                
                w = iv2intval(a0).^iv2intval(b0);                  
            else
                a0.inf = repmat(a_0.inf,m/m_a,n/n_a);
                a0.sup = repmat(a_0.sup,m/m_a,n/n_a);
                b0 = repmat(b_,m/m_b,n/n_b);
                w = iv2intval(a0).^b0;  
            end
            if any(isnan(w(:))) || any(isinf(w(:)))
                error('Inf/NaN occured in exponentiation')
            end
            w = intval2iv(w);
        else
            idx1 = (0:s-1).';
            idx2 = s-idx1;
            u.inf = r.inf(:,:,idx1+1);
            u.sup = r.sup(:,:,idx1+1);
            v.inf = a_.inf(:,:,idx2+1);
            v.sup = a_.sup(:,:,idx2+1);
            idx1_ = reshape(idx1,1,1,[]);
            idx2_ = reshape(idx2,1,1,[]);
            c = iv_minus(iv_times(idx2_,b_),idx1_);  % c = idx2.*b_-idx1
            u = iv_times(c,u);                       % u = c.*u
            w = iv_dotprod(u,v,3);                   % w = sum(u.*v)
            w = iv_rdivide(w,iv_times(s,a_0));       % w = w/(s.*a_0)
        end
        r.inf(:,:,t) = w.inf;
        r.sup(:,:,t) = w.sup;
    end    
    r_tree = r;
    r = tc2tcoeff(r,a(1));
    
end

if k
    if treenr > 0
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = r_tree; % Store result
    else
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {type_a,type_b,x,m_a,n_a,m_b,n_b,m,n,K,trigger,flag_int_power}; % Store meta data.
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = r_tree; % Store result
    end
end

end % function power

%----------------------------------r = set_exception(a,b,type_b,m_a,n_a,m_b,n_b,K)---------------------------------------%

function r = int_power(a,b)
if b == 0
    r = a.*0 + 1;
else
    start = true;
    while b
        if mod(b,2) == 1
            if start
                r = a;
                start = false;
            else
                r = r.*a;
            end
        end
        b = floor(b/2);
        if b ~= 0
            a = sqr(a);
        end
    end
end
end % function int_power
