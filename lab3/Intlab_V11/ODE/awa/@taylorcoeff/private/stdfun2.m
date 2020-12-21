function r = stdfun2(a,b,fun,funname)
%STDFUN2  collection of standard functions "fun" in
%         {tan,cot,tanh,coth,asin,acos,atan,acoth,asinh,acosh,atanh,acoth}.
%
%   r = stdfun2(a,b,fun,funname)  
%
% Computes fun(a) in Taylor coefficient arithmetic.
% Corresponds to Lohner's function TKOEFF,
% case OPER = {TANG,COTANG,TANGH,COTANGH,ARCSI,ARCCO,
%              ARCTG,ARCCTG,ARSIH,ARCOH,ARTGH,ARCTGH},
% file awa.p. 
%
% Example for Taylor coefficient recursion formula, 
% see [L], (1.1.4.10) last equation, p. 20:    
% 
%   fun = atan 
%   b = 1+a^2   
%   (atan(a))_s = [(a)_s - 1/s * sum_{j=1}^{s-1} j*(atan(a))_j * (b)_{s-j}] / b


% written  08/01/17     F. Buenger
% modified 01/26/18     F. Buenger  vectorized version added

global INTLAB_AWA_VARS

status = INTLAB_AWA_VARS.STATUS;

if ~isempty(status)
    status_ = status;
    k = status + 1;
    treenr = INTLAB_AWA_VARS.TREENR;
    vertexnr = INTLAB_AWA_VARS.VERTEXNR + 1;
    INTLAB_AWA_VARS.VERTEXNR = vertexnr;   % Each standard function in this collection is stored as a consecutive pair in the tree.
                                           % Thus INTLAB_AWA_VARS.VERTEXNR is increased by two. 
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
    numel_bound = 1; 
    trigger = (m*n <= numel_bound);
    %trigger = true;    % only for testing !!! 
    %trigger = false;   % only for testing !!!    
    if trigger
        r = a; % just preallocation
    else
        r.inf = zeros(m,n,K); % just preallocation
        r.sup = r.inf;
    end
end

% Example for Taylor coefficient recursion formula, see [L], (1.1.4.10) last equation, p. 20:    
% 
%             fun = atan 
%             b = 1+a^2   
%             (atan(a))_s = [ (a)_s - 1/s * sum_{j=1}^{s-1} j*(atan(a))_j * b_{s-j} ] / b

if trigger
    
    for i = 1:m
        for j = 1:n
            r_ = r(i,j);
            a_ = a(i,j);
            b_ = b(i,j);
            b_1.inf = b_.inf(1);
            b_1.sup = b_.sup(1);            
            for s = status_
                t = s+1;
                if s == 0
                    u.inf = a_.inf(1);
                    u.sup = a_.sup(1);
                    x = intval2iv(fun(iv2intval(u)));
                else
                    if s == 1
                        x.inf = 0;
                        x.sup = 0;
                    else
                        idx1 = (1:s-1).';
                        idx2 = s-idx1;
                        u.inf = r_.inf(idx1+1);
                        u.sup = r_.sup(idx1+1);
                        u = iv_times(idx1,u);
                        v.inf = b_.inf(idx2+1);
                        v.sup = b_.sup(idx2+1);
                        x = iv_dotprod(u,v);
                    end
                    x = iv_rdivide(x,s); % x = x/status
                    a_k.inf = a_.inf(t);
                    a_k.sup = a_.sup(t);
                    switch funname
                        case {'cot','coth','acos','acot'}
                            x =  iv_minus(iv_uminus(a_k),x);  % x = -a_k - x
                        otherwise
                            x = iv_minus(a_k,x);              % x =  a_k - x
                    end
                    x =  iv_rdivide(x,b_1);                   % x = x / b_1
                end
                r_.inf(t) = x.inf;
                r_.sup(t) = x.sup;
            end
            r(i,j) = r_;
        end
    end
    r_tree = r;
    
else
    
    a_ = tcoeff2tc(a);
    b_ = tcoeff2tc(b);
    b_1.inf = b_.inf(:,:,1);
    b_1.sup = b_.sup(:,:,1);
    for s = status_
        t = s+1;
        if s == 0
            u.inf = a_.inf(:,:,1);
            u.sup = a_.sup(:,:,1);
            x = intval2iv(fun(iv2intval(u)));
        else
            if s == 1
                x.inf = 0;
                x.sup = 0;
            else
                idx1 = (1:s-1).';
                idx2 = s-idx1;
                u.inf = r.inf(:,:,idx1+1);
                u.sup = r.sup(:,:,idx1+1);
                u = iv_times(reshape(idx1,1,1,[]),u); % This simply means u = idx1 .* u. 
                v.inf = b_.inf(:,:,idx2+1);
                v.sup = b_.sup(:,:,idx2+1);
                x = iv_dotprod(u,v,3);
            end
            x = iv_rdivide(x,s); % x = x/status
            a_k.inf = a_.inf(:,:,t);
            a_k.sup = a_.sup(:,:,t);
            switch funname
                case {'cot','coth','acos','acot'}
                    x =  iv_minus(iv_uminus(a_k),x);  % x = -a_k - x
                otherwise
                    x = iv_minus(a_k,x);              % x =  a_k - x
            end
            x =  iv_rdivide(x,b_1);                   % x = x / b_1
        end
        r.inf(:,:,t) = x.inf;
        r.sup(:,:,t) = x.sup;
    end
    r_tree = r;
    r = tc2tcoeff(r_tree,a(1));
        
end

if k
    if treenr > 0
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = r_tree; % Store result.
    else % negative treenr signalizes that the function tree is traversed for the first time and that meta data shall be stored once.
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {m,n,K,trigger}; % Store meta data.
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = r_tree; % Store result.
    end
end

end % function stdfun2