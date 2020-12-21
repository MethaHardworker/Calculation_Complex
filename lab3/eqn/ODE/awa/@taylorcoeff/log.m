function r = log(a)
%LOG  natural logarithm for Taylor coefficient
%
%   r = log(a)  
%
% Corresponds to Lohner's function TKOEFF, case OPER = LOGA, file awa.p. 
%
% Taylor coefficient recursion formula, see [L], (1.1.4.10), p. 20:    
%   
%   (a)_s = [ (a)_s - 1/s * sum_{j=1}^{s-1} (s-j) * (a)_s * (log(a))_{s-j} ] / a

% written  08/01/17     F. Buenger
% modified 08/11/17     F. Buenger  new treemode INTLAB_AWA_VARS.TREE_OFF


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
        r = a;
    else
        r.inf = zeros(m,n,K);
        r.sup = r.inf;
    end
end

if trigger
    
    for i = 1:m
        for j = 1:n
            r_ = r(i,j);
            a_ = a(i,j);
            a_1.inf = a_.inf(1);
            a_1.sup = a_.sup(1);
            for s = status_
                t = s+1;
                if s == 0
                    u.inf = a_.inf(1);
                    u.sup = a_.sup(1);
                    w = intval2iv(log(iv2intval(u)));
                else
                    if s == 1
                        w = 0;
                    else
                        idx1 = (1:s-1).';
                        idx2 = s-idx1;
                        u.inf = a_.inf(idx1+1);
                        u.sup = a_.sup(idx1+1);
                        v.inf = r_.inf(idx2+1);
                        v.sup = r_.sup(idx2+1);
                        v = iv_times(idx2,v);
                        w = iv_dotprod(u,v);
                    end
                    a_k.inf = a_.inf(t);
                    a_k.sup = a_.sup(t);                    
                    w = iv_rdivide(w,-s);
                    w = iv_plus(a_k,w);
                    w = iv_rdivide(w,a_1);
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
    a_1.inf = a_.inf(:,:,1);
    a_1.sup = a_.sup(:,:,1);
    for s = status_
        t = s+1;
        if s == 0
            u.inf = a_.inf(:,:,1);
            u.sup = a_.sup(:,:,1);
            w = intval2iv(log(iv2intval(u)));
        else
            if s == 1
                w = 0;
            else
                idx1 = (1:s-1).';
                idx2 = s-idx1;
                u.inf = a_.inf(:,:,idx1+1);
                u.sup = a_.sup(:,:,idx1+1);
                v.inf = r.inf(:,:,idx2+1);
                v.sup = r.sup(:,:,idx2+1);
                v = iv_times(reshape(idx2,1,1,[]),v); % This simply means v = idx2 .* v .
                w = iv_dotprod(u,v,3);
            end
            a_k.inf = a_.inf(:,:,t);
            a_k.sup = a_.sup(:,:,t);
            w = iv_rdivide(w,-s);
            w = iv_plus(a_k,w);
            w = iv_rdivide(w,a_1);
        end
        r.inf(:,:,t) = w.inf;
        r.sup(:,:,t) = w.sup;
    end        
    r_tree = r;
    r = tc2tcoeff(r,a(1));
    
end

if k
    if treenr > 0
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = r_tree; % Store result.
    else % negative treenr signalizes that the function tree is traversed for the first time and that meta data shall be stored once.
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {m,n,K,trigger}; % Store meta data.
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = r_tree; % Store result.
    end
end

end % function log