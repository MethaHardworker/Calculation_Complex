function r = integral(a,k)
%INTEGRAL        Taylor model integration with respect to one Variable t_k
%
%   r = integral(a,k)
%
% a : Taylor model. Maximally two dimensions are allowed, i.e., "a" can be
%     a matrix but not a higher dimensional object.
%
% k : Index of the variable with respect to which shall be integrated.
%     Thus k must be a positive integer  with 1 <= k <= a(i,j).dim
%     for all i<= size(a,1), j <= size(a,2).
%

% written  09/10/15     F. Buenger
% modified 11/23/15     F. Buenger  matrix input (pointwise integration)
% modified 12/09/15     F. Buenger  new parameter 'mode'
% modified 12/14/15     F. Buenger  parameter 'mode' moved to INTLAB_ODE_VARS.ODEMODE
% modified 01/20/16     F. Buenger  record feature
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures 

global INTLAB_ODE_VARS
ODEMODE = INTLAB_ODE_VARS.ODEMODE;
RECMODE = INTLAB_ODE_VARS.RECMODE;

if RECMODE ~= 2
    sparsity_tol = get_sparsity_tol; % interval coefficients c with |mid(c)| < sparsity_tol are set to zero to keep the Taylor polynomial sparse.
                                     % The error caused by that will be transferred to the error interval r.interval
end

r = a; % initialize the result r with a
S_a = size(a);

for i = 1:S_a(1)
    for j = 1:S_a(2)
        a_ = a(i,j);
        if RECMODE ~= 2
            r_ = a_;
            % integration of Taylor polynomial with respect to the single Variable t_k
            M = a_.monomial;
            M(:,k) = M(:,k)+1; 
            K = M(:,k); % In particular, we have K > 0 componentwise , since M >= 0. Thus division by K is safe.
            c = r_.coefficient;
            c_mid = c ./ K; % non-verified computation of coefficients of integrated Taylor model
            if ODEMODE == 1 % verified computation of coefficients of integrated Taylor model
                % Recall that rounding is upwards in verified mode and that K > 0 componentwise.
                c_iv.inf = -((-c)./K); 
                c_iv.sup = c./K; 
            end
            
            row = and( abs(c_mid) >= sparsity_tol , sum(M,2) <= r_.order );
            % Find row indices i for which |c_mid(i)| >= sparsity_tol and degree(r.monomial(i)) <= r.order
            if ~any(row)
                r_.monomial = zeros(1,r_.dim); % Empty polynomials are avoided!
                r_.coefficient = 0;            
            else
                r_.coefficient = c_mid(row); % Store the corresponding coefficients as those of the result r_ .
                r_.monomial = M(row,:);      % Store also the corresponding monomials.
            end
            
            if ODEMODE == 1 % verified computation ; computation of error interval and image of the integrated Taylor model.
                % Determine error interval r.interval according to [E], p.110,111.
                rest = r_;
                rest.order = r_.order+1; % Due to integration the order of rest increases by 1.
                rest.monomial = M;
               
                % Center "relevant" coefficients intervals around zero. Recall that rounding is upwards in verified mode.
                c_iv.inf(row) = -(c_mid(row) - c_iv.inf(row)); 
                c_iv.sup(row) = c_iv.sup(row) - c_mid(row); 
                
                rest.coefficient = c_iv;
                rest.image = image(rest);
                                
                % Compute error interval according to [E], p.111:
                %    r_.interval = rest.image + a_.interval * (a_.domain(k)-a_.center(k))                  
                
                % First compute interval bounds of a_.domain(k)-a_.center(k). Recall that rounding is upwards. 
                diff_dom_k.inf = -(a_.center(k) - a_.domain.inf(k)); 
                diff_dom_k.sup = a_.domain.sup(k) - a_.center(k);
                
                r_.interval = iv_plus( rest.image , iv_times(a_.interval,diff_dom_k) );
                
                if RECMODE == 1 % record write mode
                   % Store rest.image,  a_dom_k
                   push_reclist('rest.image', rest.image);
                   push_reclist('diff_dom_k', diff_dom_k); 
                end
            else % In the non-verified mode r_.interval and r_.image are set to the empty interval.
                r_.interval = INTLAB_ODE_VARS.EMPTYIV;
                r_.image = INTLAB_ODE_VARS.EMPTYIV;
            end
        else % RECMODE == 2, record read mode, only the error interval component is computed !!!
            r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;    % Initialize the component result r_ .
            rest_image = pull_reclist('rest.image'); % Read rest_image (former rest.image) from record list.
            diff_dom_k = pull_reclist('diff_dom_k'); % Read diff_dom_k = a_.domain(k)-a_.center(k) from record list.
            r_.interval = iv_plus( rest_image , iv_times(a_.interval,diff_dom_k) );
        end % RECMODE ~=2
        r(i,j) = r_;
    end % j
end % i

end % function integral

