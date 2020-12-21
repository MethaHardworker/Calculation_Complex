function r = LDB(a)
%LDB         compute the image of taylormodel-like structures
%            using a L(inear)D(ominated)B(ounder) method.
%
%            See [M] K. Makino, Rigorous Analysis of Nonlinear Motion in Particle Accelerators,
%                     Dissertation, Michigan State University, 1998, p. 128-130
%                [E] I. Eble, Über Taylormodelle, Dissertation, Universität Karlsruhe, 2007, p. 5,6
%                [N] A. Neumaier, Taylor Forms -- Use and Limits, Reliable Computing 9: 43-79, 2003
%                    Theorem 9.1, p. 58
%
%
%   r = LDB(a) 
%
% Let p(t) be the multivariate polynomial given by a.monomial and a.coefficient.
% Set I := a.domain and z := a.center.
% For determining the image of the Taylor model a, the LDB method tries to enclose
% values t_min in I_min and t_max in I_max such that f(t):=p(t-z) attains its global minimum
% and maximum at t_min and t_max, respectively. Here I_min and I_max are sought subdomains of I.
% Clearly, the image of "a" is contained in the interval [inf(f(I_min)),sup(f(I_max))].
%
% Our implementation follows [N], Theorem 9.1, p. 58. We have the following minor comment to this theorem:
%
% (1) On the right hand-side of formula (9.3) a supremum seems to be missing, i.e., it must read
%
%     (9.3)    |mid(B_i)| diam(X_i) > delta := diam(A) + sup_{x in X} diam(B)^T |x-z| .
%
%     (Here capital A, B, and X stand for Neumaier's bold marked interval \bf{a} and interval vectors \bf{b} and \bf{x}.)
%
% Theorem 9.1 assumes a given interval A and a given column interval vector B such that condition
%
%      (9.1)    f(X) is a subset of A + B^T(X-z)
%
% holds true.
%
% Our implementation simply takes B as the linear part of p(t).
% Accordingly, we take A as the naive interval image of dummy(t):=p(t-t_0)-B^T(t-t_0).
% Note that if B is a point-vector (diam(B) = 0), then (9.3) becomes
%
%     (9.3')    B_i * diam(X_i) > delta := diam(A) .
%
% Note also that our LDB method may improve naive interval evaluation only if the Taylor model a
% has a nonzero linear part B.

% written  10/26/15     F. Buenger

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

S_a = size(a);
if (length(S_a) > 2)
    error('maximally two dimensions for type taylormodel');
end

c_improve = 0.99; % LDB improvement steps are repeated only if the diameter of at least one component of the interval vectors 
                  % I_min and I_max, which contain a point at which the minimum and the maximum is attained, respectively, 
                  % decreased by 1 percent, i.e.: diam_new < c_improve * diam_old. 
                  % The constant is chosen quite arbitrarily. Feel free to change that!
                  
max_iter = 3; % maximum number of LDB iterations. The number is chosen quite arbitrarily. Feel free to change that!

% Initialize result interval (matrix).
r.inf = zeros(S_a) ;
r.sup = r.inf; 

for i = 1:S_a(1)
    for j = 1:S_a(2)
        a_ = a(i,j); 
        coef = a_.coefficient;
        M = a_.monomial;
        L = (sum(M,2) == 1); % indicator for monomials corresponding to linear terms of a_
        [c_mid,c_rad] = iv_getmidrad(coef);
        idx_lin_pos = find(L .* (c_mid > 0)); % Determine indices of all positive midpoints of linear terms of a_ .
        idx_lin_neg = find(L .* (c_mid < 0)); % Determine indices of all negative midpoints of linear terms of a_ .
        M_pos = a_.monomial((idx_lin_pos),:);
        M_neg = a_.monomial((idx_lin_neg),:);
        is_float =  isfloat(coef);    
                
        if isempty(idx_lin_pos) &&  isempty(idx_lin_neg) % LDB not applicable, use naive method
            hlp = image(a_,0);
            r.inf(i,j) = hlp.inf;
            r.sup(i,j) = hlp.sup;
            continue;
        end
            
        if ~isempty(idx_lin_pos)
            B_pos = c_mid(idx_lin_pos)' * M_pos; % linear terms with positive midpoints. 
                                                 % The ordering corresponds to that of a_.domain.
            if is_float
                B = B_pos';
            else
                B.inf = (coef.inf(idx_lin_pos)' * M_pos).';
                B.sup = (coef.sup(idx_lin_pos)' * M_pos).';
            end                                                                                                  
            diam_B_pos = 2*c_rad(idx_lin_pos)' * M_pos; % corresponding diameters
            diam_B = diam_B_pos;
            B_pos_idx = find(B_pos > 0);                % indices of non-zero entries in the order of a_.domain.
            B_pos = B_pos(B_pos_idx);
        else
            diam_B = zeros(1,a_.dim);
            if is_float
                B = diam_B';
            else
                B.inf = diam_B';
                B.sup = B.inf;                
            end
        end
        if ~isempty(idx_lin_neg)
            B_neg = c_mid(idx_lin_neg)' * M_neg; % linear terms with negative midpoints.
                                                 % The ordering corresponds to that of a.domain.
            if is_float
                B = B + B_neg'; % Note that no rounding errors occur since idx_lin_pos and idx_lin_neg are disjoint.
            else
                B.inf = B.inf + (coef.inf(idx_lin_neg)' * M_neg).'; % Note that no rounding errors occur since idx_lin_pos and idx_lin_neg are disjoint.
                B.sup = B.sup + (coef.sup(idx_lin_neg)' * M_neg).'; % Note that no rounding errors occur since idx_lin_pos and idx_lin_neg are disjoint.
            end                                                                                                                                                   
            diam_B_neg = 2*c_rad(idx_lin_neg)' * M_neg; % corresponding diameters
            diam_B = diam_B + diam_B_neg;               % Note that no rounding errors occur since idx_lin_pos and idx_lin_neg are disjoint.
            B_neg_idx = find(B_neg < 0);                % indices of non-zero entries in the order of a_.domain.
            B_neg = B_neg(B_neg_idx);          
        end
        
        % The interval A in formula (9.1) is computed by building the naive interval image of
        % dummy(t):=a_(t-t_0)-B^T(t-t_0) over the whole domain a_.domain. 
        % Here B^T(t-t_0) denotes the linear part of a_ having nonzero coefficient midpoints, which is eliminated by subtraction.
        dummy = a_;
        if is_float % case a_.coefficient is a point vector
            % Eliminate linear part with nonzero coefficient midpoints.
            dummy.coefficient(idx_lin_pos) = 0;  %  
            dummy.coefficient(idx_lin_neg) = 0;  %  
        else % case a_.coefficient is an interval vector
            % Eliminate linear part with nonzero coefficient midpoints.
            dummy.coefficient.inf(idx_lin_pos) = 0;  
            dummy.coefficient.sup(idx_lin_pos) = 0;  
            dummy.coefficient.inf(idx_lin_neg) = 0;   
            dummy.coefficient.sup(idx_lin_neg) = 0;   
        end        
        
        % Initialize variables A_min and A_max which correspond to the interval A in formula (9.1) 
        % depending on whether the minimum or maximum is sought.
        % In the beginning A_min and A_max are both taken as the naive interval image of "dummy",
        % which differs from that of a_ by the contribution of the linear terms having nonzero coefficient midpoints.

        A_min = image(dummy,0);
        A_max = A_min; 

        I_min = a_.domain; % Initialize the domain that contains t_min with the whole domain a_.domain.
        I_max = I_min;     % Initialize the domain that contains t_max also with the whole domain a_.domain.
        c = a_.center;
        z_min = c;         % Initialize non-verified center of I_min
        z_max = c;         % Initialize  non-verified center of I_max
        improve = true;
        
        diam_I_min_old = iv_diam(I_min);         
        diam_I_max_old = diam_I_min_old; % = iv_diam(I_max) at the beginning         
        
        iter = 1;
%       while iter <= max_iter && improve            
        while iter <= max_iter && improve            
            if iter > 1
                A_min = tmval(dummy,I_min);      % Update A_min by considering only the reduced domain I_min.                
                z_min = (I_min.sup-I_min.inf)/2; % non-verified center of I_min (Any point in I_min would be sufficient.)
                
                % f(x) in A_min + B'*(x-c) = (A_min + B'(z_min-c)) + B'*(x-z_min) 
                % => Update A_min := A_min + B'(z_min-c))
                A_min = iv_plus(A_min,iv_dotprod(B,iv_minus(z_min,c))); 
                
                A_max = tmval(dummy,I_max);      % Update A_max by considering only the reduced domain I_max.
                z_max = (I_max.sup-I_max.inf)/2; % non-verified center of I_min (Any point in I_min would be sufficient.)
                
                % f(x) in A_max + B'*(x-c) = (A_max + B'(z_max-c)) + B'*(x-z_max) 
                % => Update A_max := A_max + B'(z_max-c))
                A_max = iv_plus(A_max,iv_dotprod(B,iv_minus(z_max,c))); 
            end
            % (I) update I_min
            diam_A_min = iv_diam(A_min);
            % Define delta according to formula (9.3) of Neumaier's article [N].
            X = iv_abs(iv_minus(I_min,z_min)); % interval computation of |I_min-z|
            delta = diam_A_min+diam_B*X.sup; % diam_B,X >= 0 and recall that rounding is upwards. 
            
            % (I.1) The following corresponds to the first line of (9.4), p. 58, of Neumaier's article [N].
            if ~isempty(idx_lin_pos)
                D = delta ./ B_pos'; % Recall that rounding is upwards and that delta,B_pos >= 0.
                J = I_min.inf(B_pos_idx);
                %I_min(B_pos_idx) = infsup( J , min(sup(J+D),sup(I_min(B_pos_idx))) );
                I_min.inf(B_pos_idx) = J;
                I_min.sup(B_pos_idx) = min(J+D,I_min.sup(B_pos_idx)); % Recall that rounding is upwards.
            end
            
            % (I.2) The following corresponds to the second line of (9.4), p.58, of Neumaier's article [N].
            if  ~isempty(idx_lin_neg)
                D = delta ./ B_neg'; % Recall that rounding is upwards and that delta, B_neg >= 0.
                J = I_min.sup(B_neg_idx);
                % I_min(B_neg_idx) = infsup( max(inf(J+D),inf(I_min(B_neg_idx))) , J);
                I_min.inf(B_neg_idx) = max(-(-J-D),I_min.inf(B_neg_idx)); % Recall that rounding is upwards.
                I_min.sup(B_neg_idx) = J;
            end

            % (II) update I_max
            diam_A_max = iv_diam(A_max);
            % Define delta according to formula (9.3) of Neumaier's article [N].
            X = iv_abs(iv_minus(I_max,z_max)); % interval computation of |I_max-z|
            delta = diam_A_max+diam_B*X.sup; % diam_B,X >= 0 and recall that rounding is upwards. 
            
            % (II.1) The following corresponds to the first line of (9.5), p.58, of Neumaier's article [N].
            if  ~isempty(idx_lin_neg)
                D = delta ./ B_neg';
                J = I_max.inf(B_neg_idx);
                %I_max(B_neg_idx) = infsup( J , min(sup(J-D),sup(I_max(B_neg_idx))) );
                I_max.inf(B_neg_idx) = J;
                I_max.sup(B_neg_idx) = min(J-D,I_max.sup(B_neg_idx));% Recall that rounding is upwards.               
            end
            
            % (II.2) The following corresponds to the second line of (9.5), p.58, of Neumaier's article [N].
            if ~isempty(idx_lin_pos)
                D = delta ./ B_pos';
                J = I_max.sup(B_pos_idx);
                %I_max(B_pos_idx) = infsup( max(inf(J-D),inf(I_max(B_pos_idx))) , J);
                I_max.inf(B_pos_idx) = max(-(-J+D),I_max.inf(B_pos_idx)) ; % Recall that rounding is upwards. 
                I_max.sup(B_pos_idx) = J;
            end
            
            diam_I_min = iv_diam(I_min);
            diam_I_max = iv_diam(I_max);
            
            improve = any( [diam_I_min;diam_I_max] < c_improve*[diam_I_min_old;diam_I_max_old] );
            
            diam_I_min_old = diam_I_min;
            diam_I_max_old = diam_I_max; 
            iter = iter + 1;
        end % while
        
        % Compute enclosure of the minimum of a_ which is assumed on the domain I_min according to the LDB-iteration.
        image_min = tmval(a_,I_min); % enclosure of min(p(t-t_0))
        
        % Compute enclosure of the maximum of a_ which is assumed on the domain I_max according to the LDB-iteration.        
        image_max = tmval(a_,I_max); % enclosure of max(p(t-t_0))
        
        % The image of a_ is enclosed in the interval [image_min.inf,image_max.sup].        
        r.inf(i,j) = image_min.inf;
        r.sup(i,j) = image_max.sup;
    end % j
end %i

if rndold ~= 1 
    setround(rndold)
end

end % function LDB