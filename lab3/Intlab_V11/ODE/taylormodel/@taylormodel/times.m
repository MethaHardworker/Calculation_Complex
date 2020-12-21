function r = times(a,b)
%TIMES  Taylor model multiplication a .* b
%
%   r =  times(a,b)  

% written  08/31/15  F. Buenger
% modified 11/23/15  F. Buenger  componentwise multiplication of matrices
% modified 12/10/15  F. Buenger  check rounding 'upwards' instead of 'to nearest', code optimization
% modified 12/14/15  F. Buenger  switch between verified/non-verified Taylor model arithmetic
% modified 01/19/16  F. Buenger  record feature
% modified 02/11/16  F. Buenger  "intval"-components --> intval-like structures 
% modified 05/09/17  F. Buenger  switch between sparse and full convolution   

global INTLAB_ODE_VARS

ODEMODE = INTLAB_ODE_VARS.ODEMODE;
RECMODE = INTLAB_ODE_VARS.RECMODE;
SP = INTLAB_ODE_VARS.SPARSE;

e = 1e-30;
if ODEMODE == 1     % verified mode for Taylor model arithmetic
    if 1+e > 1      % fast check for rounding upwards
        rndold = 1;
    else
        rndold = getround;
        setround(1) % rounding upwards in verified mode
    end
else                % non-verified mode for Taylor model arithmetic
    if 1+e == 1-e   % fast check for rounding to nearest
        rndold = 0;
    else
        rndold = getround;
        setround(0) % rounding to nearest in non-verified mode
    end
end

% Case "taylormodel x float (intval)"

if ~isa(b,'taylormodel')
    % interchange a and b for uniform treatment of the case "... + taylormodel"
    c = b;
    b = a;
    a = c;
end

if isfloat(a) && isreal(a)
    a = float2iv(a);   % Convert float to iv in verified mode.
elseif isa(a,'intval')
    a = intval2iv(a);  % Convert intval to iv for common treatment.
end
isiv_a = isiv(a);

if ~(isiv_a || isa(a,'taylormodel'))
   error('incompatible factors')  % Send error massage for "x * taylormodel" if x is not float/intval/iv/taylormodel.
end

if isiv_a
    S_a = size(a.inf);
else
    S_a = size(a);
end
S_b = size(b);
if (length(S_a) > 2) || (length(S_b) > 2)
    error('Maximally two dimensions are allowed for type taylormodel.');
end
if (max(S_a) == 0 || max(S_b) == 0)
    error('empty input');
end

m = max(S_a(1),S_b(1));
n = max(S_a(2),S_b(2));
ismat_a = (max(S_a) > 1); % Recognize if a is a nontrivial matrix.
ismat_b = (max(S_b) > 1); % Recognize if b is a nontrivial matrix.
if ismat_a && ismat_b && any(S_a ~= S_b)
    if ~all( S_a == 1 | S_b == 1 | S_a == S_b)
        error('incompatible dimension sizes')        
    end
    if S_a(1) < m
      if isiv_a
        a = iv_repmat(a,m,1);
      else
        a = repmat(a,m,1);
      end
    end
    if S_a(2) < n
      if isiv_a
        a = iv_repmat(a,1,n);            
      else
        a = repmat(a,1,n);  
      end
    end
    if S_b(1) < m
      b = repmat(b,m,1);
    end
    if S_b(2) < n
      b = repmat(b,1,n);  
    end
end

if RECMODE ~= 2
    sparsity_tol = get_sparsity_tol; % interval coefficients c with |mid(c)| < sparsity_tol are set to zero to keep the Taylor polynomial sparse.
                                     % The error caused by that will be transferred to the error interval r.interval.
end

r(1:m,1:n) = b(1); % Initialize result r as an (mxn)-matrix of Taylor models.
for i = 1:m
    for j = 1:n
        if ismat_a % a is a nontrivial matrix.
            if isiv_a
                a_.inf = a.inf(i,j);
                a_.sup = a.sup(i,j);
            else
                a_ = a(i,j);
            end
        else
            a_ = a;
        end
        if ismat_b % b is a nontrivial matrix.
            b_ = b(i,j);
        else
            b_ = b; 
        end
        if RECMODE ~= 2
            r_ = b_;  % Initialize the result r with b_.        
            if isiv_a % Case "iv * taylormodel" 
                if ODEMODE == 1 % verified mode
                    a_image = a_;
                    a_interval = INTLAB_ODE_VARS.ZEROIV; % zero error interval [0,0]
                    c_ = iv_times(a_,b_.coefficient);    % Multiplication of float/intval a_ and taylormodel b_ means multiplying all coefficients of b_ with a_.
                    c_lower = c_.inf;
                    c_upper = c_.sup;                   
                end
                c_mid = (0.5.*(a_.inf+a_.sup)) .*  b_.coefficient; % An exact midpoint is not needed here.
                U = b_.monomial;   
            else % Case "taylormodel .* taylormodel"
                if a_.dim ~= b_.dim
                    error('Dimensions must agree.')
                end
                if a_.order ~= b_.order
                    error('Orders of Taylor models must agree.')
                end
                if any(a_.center ~= b_.center)
                    error('Center points must be equal.')
                end
                if ~iv_eq(a_.domain,b_.domain)
                    error('Domains of factors must be equal.')
                end                
                max_monomial = max(a_.monomial,[],1)+max(b_.monomial,[],1); % Exponents for encoding multivariate polynomials according to Kronecker's trick.
                [D_a,d_a] = encode_monomials(a_.monomial,max_monomial);
                D_b = encode_monomials(b_.monomial,[],d_a);                  
                if SP
                    a_poly = sparse(D_a,1,a_.coefficient);
                    b_poly = sparse(D_b,1,b_.coefficient);
                else 
                    a_poly = full(sparse(D_a,1,a_.coefficient)); % Weird but fast creation of a full univariate polynomial representation of a_
                                                                 % with coefficients in correct order according to their degrees given in D_a.
                    b_poly = full(sparse(D_b,1,b_.coefficient)); % same for b_
                end                                
                if ODEMODE == 1 % verified mode
                    a_image = a_.image;
                    a_interval = a_.interval;                    
                    if SP 
                        c_upper = sconv(a_poly,b_poly);
                        c_lower = -sconv(-a_poly,b_poly);
                        c_upper = full(c_upper);
                        c_lower = full(c_lower);                                       
                    else 
                        c_upper = conv(a_poly,b_poly);
                        c_lower = -conv(-a_poly,b_poly);
                    end                    
                    % Suppress common zero entries.
                    D_unique = find(abs(c_lower)+abs(c_upper)); % Encoded monomials of the nonzero interval coefficients of the product polynomial
                    c_lower = c_lower(D_unique);                % according compressed lower coefficient bounds, see below command line: "U = decode_monomials(D_unique,d_a);" 
                    c_upper = c_upper(D_unique);                % according compressed upper coefficient bounds
                    c_mid = 0.5 * (c_lower+c_upper);            % Approximate coefficient midpoints of the product polynomial
                else % non-verified mode                                        
                    if SP 
                        c_mid = sconv(a_poly,b_poly);           % some approximate polynomial product, not necessarily a "midpoint".                    
                    else    
                        c_mid = conv(a_poly,b_poly);            % Compute some approximate polynomial product, not necessarily a "midpoint".
                    end                                        
                    [D_unique,dummy,c_mid] = find(c_mid);
                end
                U = decode_monomials(D_unique,d_a);             % U is the decoded matrix of unique monomials in correct order corresponding to c_mid, c_lower, c_upper.
            end     
            
            % Find row indices i for which |c_mid(i)| >= sparsity_tol and degree(r.monomial(i)) <= r_.order. 
            % The contribution of the terms of all other row indices will be moved to the error interval. See "rest" below.
            row = and( (abs(c_mid) >= sparsity_tol) , (sum(U,2) <= r_.order) );
                                
            if ~any(row)
                r_.monomial = zeros(1,r_.dim); % Empty Taylor models shall be avoided.
                r_.coefficient = 0;
            else
                r_.coefficient = c_mid(row); % Store the corresponding coefficients as those of the result r_ = a_ * b_ .
                r_.monomial = U(row,:);      % Store also the corresponding monomials.
                if ODEMODE == 1
                    % Recall that rounding is upwards.
                    c_lower(row) = -(r_.coefficient-c_lower(row)); % Compute centered lower bound for non-sparse coefficients.
                    c_upper(row) = c_upper(row)-r_.coefficient;    % Compute centered upper bound for non-sparse coefficients.
                end
            end
            
            % Compute image and error interval.
            if ODEMODE == 1 % verified mode.
                % Determine error interval r_.interval of a*b according to p.107,108 of [E].
                if ~isempty(U)
                    rest = r_;
                    if isa(a,'taylormodel')
                        rest.order = 2*r_.order; % The product of two polynomials of max degree n is a polynomial of max degree 2n.
                    end
                    rest.monomial = U;               
                    coeff.inf = c_lower;
                    coeff.sup = c_upper;
                    rest.coefficient = coeff; 
                    rest.image = image(rest);
                else
                    rest.image = INTLAB_ODE_VARS.ZEROIV; % rest.image := zero interval [0,0];
                end
                
                % Computation of 
                %
                %   I_ab := rest.image + a_interval .* (b_.image + b_.interval) + a_image .* b_.interval;
                %   I_ba := rest.image + b_.interval .* (a_image + a_interval) + b_.image .* a_interval;
                %
                %   r_.interval := intersect(I_ab,I_ba);  see [E], p.107,108 of [E]
                                
                r_.interval = times_interval(a_interval,a_image,b_.interval,b_.image,rest.image);
                r_.image = image(r_);
                if RECMODE == 1  % record write mode
                    push_reclist('rest_image',rest.image); % Write rest.image to record list.
                    push_reclist('a_image',a_image);       % Write a_image to record list.
                    push_reclist('b_image',b_.image);      % Write b_.image to record list.
                end
            else  % In non-verified mode r_.interval and r_.image are set to the empty interval for clarity.
                r_.interval = INTLAB_ODE_VARS.EMPTYIV;
                r_.image = INTLAB_ODE_VARS.EMPTYIV;
            end            
        else % RECMODE == 2, record read mode, only the error interval component is computed !!!
            r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;    % Initialize the component result r_ .        
            rest_image = pull_reclist('rest_image'); % Read rest_image (former rest.image) from record list.
            a_image = pull_reclist('a_image');       % Read a_image from record list.
            b_image = pull_reclist('b_image');       % Read b_image (former b_.image) from record list.
            if isiv_a        
                a_interval = INTLAB_ODE_VARS.ZEROIV; % trivial error interval [0,0] 
            else
                a_interval = a_.interval;
            end
            r_.interval = times_interval(a_interval,a_image,b_.interval,b_image,rest_image);
        end % RECMODE ~= 2
        r(i,j) = r_;
    end % j
end % i

if (ODEMODE == 1 && rndold ~= 1) || (ODEMODE == 0 && rndold ~= 0)
    setround(rndold)
end

end % function times
