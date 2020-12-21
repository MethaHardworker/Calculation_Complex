function r = stdfun1(a,fun1,fun2,funname1,funname2)
%STDFUN1  collection of standard functions fun1 in {sin,cos,sinh,cosh}
%                    
%   r = stdfun1(a,fun1,fun2,funname1,funname2)  
%
% Computes fun1(a) in Taylor coefficient arithmetic.  
% Corresponds to Lohner's function TKOEFF,
% case OPER = {SINUS,COSIN,SINUSH,COSINH}, file awa.p. 
%
% Example for Taylor coefficient recursion formula, 
% see [L], (1.1.4.10), p. 20:    
% 
%   fun1 = sin
%   fun2 = cos
%    
%   (sin(a))_s =  1/s * sum_{j=0}^{s-1} (j+1)*(cos(a))_{s-j-1} * (a)_{j+1}
%   (cos(a))_s = -1/s * sum_{j=0}^{s-1} (j+1)*(sin(a))_{s-j-1} * (a)_{j+1}

% written  08/01/17     F. Buenger
% modified 01/26/18     F. Buenger  vectorized version added

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
    [r,r2] = INTLAB_AWA_VARS.TREE{treenr}{vertexnr}{:};
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
        r = a;  % just preallocation 
    else
        r.inf = zeros(m,n,K); % just preallocation
        r.sup = r.inf;
    end
    r2 = r; % just preallocation   
end

if trigger
    
    for i = 1:m
        for j = 1:n
            r_ = r(i,j);
            r2_ = r2(i,j);
            a_ = a(i,j);            
            for s = status_
                t = s+1;
                if s == 0
                    u.inf = a_.inf(1);
                    u.sup = a_.sup(1);
                    x = intval2iv(fun1(iv2intval(u)));
                    y = intval2iv(fun2(iv2intval(u)));
                else
                    idx1 = (1:s).';
                    idx2 = s-idx1;
                    u.inf = a_.inf(idx1+1);
                    u.sup = a_.sup(idx1+1);
                    u = iv_times(idx1,u);
                    v.inf = r2_.inf(idx2+1);
                    v.sup = r2_.sup(idx2+1);
                    w.inf = r_.inf(idx2+1);
                    w.sup = r_.sup(idx2+1);
                    x = iv_dotprod(u,v);
                    x = iv_rdivide(x,s);
                    y = iv_dotprod(u,w);
                    y = iv_rdivide(y,s);
                    switch funname1
                        case 'sin'
                            % y = iv_uminus(y);
                            y_inf = y.inf;
                            y.inf = -y.sup;
                            y.sup = -y_inf;
                        case 'cos'
                            % x = iv_uminus(x);
                            x_inf = x.inf;
                            x.inf = -x.sup;
                            x.sup = -x_inf;
                    end
                end
                r_.inf(t) = x.inf;
                r_.sup(t) = x.sup;
                r2_.inf(t) = y.inf;
                r2_.sup(t) = y.sup;
            end
            r(i,j) = r_;
            r2(i,j) = r2_;
        end
    end
    r_tree = r;
    
else
    
    a_ = tcoeff2tc(a);
    for s = status_
        t = s+1;
        if s == 0
            u.inf = a_.inf(:,:,1);
            u.sup = a_.sup(:,:,1);
            x = intval2iv(fun1(iv2intval(u)));
            y = intval2iv(fun2(iv2intval(u)));
        else
            idx1 = (1:s).';
            idx2 = s-idx1;
            u.inf = a_.inf(:,:,idx1+1);
            u.sup = a_.sup(:,:,idx1+1);
            u = iv_times(reshape(idx1,1,1,[]),u); % This simply means u = idx1 .* u 
            v.inf = r2.inf(:,:,idx2+1);
            v.sup = r2.sup(:,:,idx2+1);
            w.inf = r.inf(:,:,idx2+1);
            w.sup = r.sup(:,:,idx2+1);
            x = iv_dotprod(u,v,3);
            x = iv_rdivide(x,s);
            y = iv_dotprod(u,w,3);
            y = iv_rdivide(y,s);
            switch funname1
                case 'sin'
                    % y = iv_uminus(y);
                    y_inf = y.inf;
                    y.inf = -y.sup;
                    y.sup = -y_inf;
                case 'cos'
                    %x = iv_uminus(x);
                    x_inf = x.inf;
                    x.inf = -x.sup;
                    x.sup = -x_inf;
            end
        end
        r.inf(:,:,t) = x.inf;
        r.sup(:,:,t) = x.sup;
        r2.inf(:,:,t) = y.inf;
        r2.sup(:,:,t) = y.sup;
    end
    r_tree = r;
    r = tc2tcoeff(r_tree,a(1));
    
end

if k
    if treenr > 0
        INTLAB_AWA_VARS.TREE{treenr}{vertexnr} = {r_tree,r2}; % Store result.
    else % negative treenr signalizes that the function tree is traversed for the first time and that meta data shall be stored once.
        INTLAB_AWA_VARS.TREE{-treenr+2}{vertexnr} = {m,n,K,trigger}; % Store meta data.
        INTLAB_AWA_VARS.TREE{-treenr}{vertexnr} = {r_tree,r2}; % Store result.
    end
end

end % function stdfun1