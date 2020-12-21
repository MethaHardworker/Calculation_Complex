function r = sqr(a)
%SQR Taylor model (elementwise) square sqr(a)
%
%   r = sqr(a)
%
%   The optional parameter "sp" (sparse) switches between using the MATLAB function
%   "conv", which only accepts non-sparse input (sp = false), and the private  
%   function "sconv", which can have sparse or non-sparse input (sp = true).
%   The default is sp = true. (It's all about performance.)

% written  11/05/15     F. Buenger
% modified 12/10/15     F. Buenger  check rounding 'upwards' instead of 'to nearest', code optimization
% modified 12/18/15     F. Buenger  switch between verified/non-verified Taylor model arithmetic
% modified 01/20/16     F. Buenger  record feature
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures 
% modified 05/09/17     F. Buenger  switch between sparse and full convolution   


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
else % non-verified mode for Taylor model arithmetic
    if 1+e == 1-e   % fast check for rounding upwards
        rndold = 0;
    else
        rndold = getround;
        setround(0) % rounding to nearest in non-verified mode
    end
end

if RECMODE ~= 2
    sparsity_tol = get_sparsity_tol; % interval coefficients c with |mid(c)| < sparsity_tol are set to zero to keep the Taylor polynomial sparse.
                                     % The error caused by that will be transferred to the error interval r.interval
end

% The computation of the polynomial part Q of the Taylor model [Q,J]:=a^2=[P,I]^2 is (more or less) the same
% as that of a*a = [P,I]*[P,I] computed by the function "times.m".
% The main difference is a slightly tighter inclusion of the error interval J.

S_a = size(a);
r = a; % Initialize result r with a. (just storage preallocation)

for i = 1:S_a(1)
    for j = 1:S_a(2)
        a_ = a(i,j);
        if RECMODE ~=2
            r_ = r(i,j); % = a(i,j)
            max_monomial = 2*max(a_.monomial,[],1); % Exponents for encoding multivariate polynomials according to Kronecker's trick.
            [D_a,d_a] = encode_monomials(a_.monomial,max_monomial);
                        
            if SP             
                a_poly = sparse(D_a,1,a_.coefficient); 
            else
                a_poly = full(sparse(D_a,1,a_.coefficient)); % Weird but fast creation of a full univariate polynomial representation of a_
            end
            
            if ODEMODE == 1 % verified mode
                                
                if SP 
                    c_upper = sconv(a_poly,a_poly);
                    c_lower = -sconv(-a_poly,a_poly);
                    c_upper = full(c_upper);
                    c_lower = full(c_lower);
                else
                    c_upper = conv(a_poly,a_poly);
                    c_lower = -conv(-a_poly,a_poly);
                end                
                
                % Suppress common zero entries.
                D_unique = find(abs(c_lower)+abs(c_upper)); % Encoded monomials of the "nonzero" coefficients of the product polynomial
                c_lower = c_lower(D_unique);                % according compressed lower coefficient bounds
                c_upper = c_upper(D_unique);                % according compressed upper coefficient bounds
                c_mid = 0.5 * (c_lower+c_upper);            % approximate coefficient midpoints of the product polynomial
                U = decode_monomials(D_unique,d_a);         % U is the decoded matrix of unique monomials in correct order corresponding to c_mid, c_lower, c_upper.
                row = and((abs(c_mid) >= sparsity_tol),(sum(U,2) <= r_.order));
                % Find row indices i for which |c_mid(i)| >= sparsity_tol and degree(r.monomial(i)) <= r_.order
                if ~any(row) % avoid empty Taylor models
                    r_.monomial = zeros(1,r_.dim);
                    r_.coefficient = 0;
                else
                    r_.coefficient = c_mid(row); % Store the corresponding coefficients as those of the result r_ = a_ * b_ .
                    r_.monomial = U(row,:);      % Store also the corresponding monomials.
                    % Recall that rounding is upwards.
                    c_lower(row) = -(r_.coefficient-c_lower(row)); % Compute centered lower bound for non-sparse coefficients.
                    c_upper(row) = c_upper(row)-r_.coefficient;    % Compute centered upper bound for non-sparse coefficients.
                end
                
                % Determine error interval r_.interval of a(i,j)^2 according to p.107,108 of [E].
                if ~isempty(U)
                    rest = r_;
                    rest.order = 2*r_.order; % The product of two polynomials of max degree n is a polynomial of max degree 2n.
                    rest.monomial = U;
                    coeff.inf = c_lower;
                    coeff.sup = c_upper;
                    rest.coefficient = coeff;
                    rest.image = image(rest);
                else
                    rest.image = INTLAB_ODE_VARS.ZEROIV; % rest.image := zero interval [0,0];
                end
                
                % Formula (2.20) in [E], p.17, for y_1=y_2:=y simplifies to
                %
                % (2.20')   I_n,y^2 = W(P_{>n}) + (I_n,y)^2 + 2*W(P_n,y)*I_n,y
                %
                % which is easily deduced by [E],(2.19) with \tilde{y}_1=\tilde{y}_2.
                % Therefore, the first to formulas in [E], p.108, reduce to one formula:
                %
                % \tilde{I}_n,y^2 = I + (I_n,y)^2 + 2*I_{P_n,y}*I_n,y    where I := rest.image.
                
                r_.interval = sqr_interval(rest.image,a_.interval,a_.image); % Compute r_.interval := rest.image + (a_.interval).^2 + 2*a_.image.*a_.interval;                
                r_.image = image(r_);                                        % Finally compute the image of the result r_.                   
                if RECMODE == 1  % record write mode
                    push_reclist('rest.image',rest.image); % Write rest_image to record list.
                    push_reclist('a_.image',a_.image);     % Write a_.image to record list.
                end
            else  % non-verified mode
                
                if SP 
                    c_mid = sconv(a_poly,a_poly); % Compute some approximate polynomial product, not necessarily a "midpoint".
                else
                    c_mid = conv(a_poly,a_poly);  % Compute some approximate polynomial product, not necessarily a "midpoint".
                end
                
                [D_unique,dummy,c_mid] = find(c_mid);           
                U = decode_monomials(D_unique,d_a); % U is the decoded matrix of unique monomials in correct order corresponding to c_mid.                
                row = and( abs(c_mid) >= sparsity_tol , sum(U,2) <= r_.order );
                % Find row indices i for which |c_mid(i)| >= sparsity_tol and degree(r_.monomial(i)) <= r_.order.
                if ~any(row)
                    r_.monomial = zeros(1,r_.dim); % Avoid empty result by setting the constant term to zero.
                    r_.coefficient = 0;
                else
                    r_.coefficient = c_mid(row) ;  % Store the corresponding coefficients as those of the result r_ = a(i,j)^2.
                    r_.monomial = U(row,:);        % Store also the corresponding monomials.
                end
                r_.interval = INTLAB_ODE_VARS.EMPTYIV;
                r_.image = INTLAB_ODE_VARS.EMPTYIV;
            end            
        else % RECMODE == 2, record read mode, only the error interval component is computed!
            r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;                       % Initialize the component result r_.
            rest_image = pull_reclist('rest.image');                    % Read rest_image (former rest.image) from record list.
            a_image = pull_reclist('a_.image');                         % Read a_image (former a_.image) from record list.
            r_.interval = sqr_interval(rest_image,a_.interval,a_image); % Compute r_.interval := rest_image + (a_.interval).^2 + 2*a_image.*a_.interval;
        end 
        r(i,j) = r_;
    end 
end 

if (ODEMODE == 1 && rndold ~= 1) || (ODEMODE == 0 && rndold ~= 0)
    setround(rndold)
end
end % function sqr