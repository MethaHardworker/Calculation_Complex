function r = plus(a,b)
%PLUS  Taylor model plus a + b
%
%   r = plus(a,b)  

% written  08/24/15     F. Buenger
% modified 11/23/15     F. Buenger  summation of matrices
% modified 12/10/15     F. Buenger  check rounding 'upwards' instead of 'to nearest', code optimization
% modified 12/14/15     F. Buenger  switch between verified/non-verified Taylor model arithmetic
% modified 01/19/16     F. Buenger  record feature
% modified 02/10/16     F. Buenger  "intval"-components --> intval-like structures 

global INTLAB_ODE_VARS

ODEMODE = INTLAB_ODE_VARS.ODEMODE;
RECMODE = INTLAB_ODE_VARS.RECMODE;

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

% Case "taylormodel + float/intval/iv"
if ~isa(b,'taylormodel')
    % Interchange a and b for uniform treatment of the case "float/intval/iv + taylormodel".
    c = b;
    b = a;
    a = c;    
end

isfloat_a = isfloat(a) && isreal(a);
if isfloat_a
  if ODEMODE == 1  
    a = float2iv(a); % Convert float to iv in verified mode.
  end
elseif isa(a,'intval')
  a = intval2iv(a);  % Convert intval to iv for common treatment. 
end
isiv_a = isiv(a);

if ~(isfloat_a || isiv_a || isa(a,'taylormodel'))
   error('incompatible addends')  % Send error massage for "x + taylormodel" if x is not float/intval/iv/taylormodel.
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
                                     % The error caused by that will be transferred to the error interval r.interval
end

r(1:m,1:n) = b(1); % Initialize result r as an (mxn)-matrix of taylormodels.
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
            if isfloat_a || isiv_a                                         % case "float/iv + taylormodel"               
                r_ = b_;                                                   % Initialize the component result r_ with b_ .
                [c0,i0] = get_constant_term(b_);                           % Find index of zero-monomial of  taylor model b_ .
                
                if ODEMODE == 1 % verified mode. Recall that rounding is upwards.
                    c_ = iv_plus(a_,c0);                                   % interval evaluation of c_ = a_ + b_;
                                                                           % Recall that float was transformed to iv in verified mode.
                    c_mid = 0.5*(c_.inf+c_.sup);                           % approximate midpoint of c
                    if abs(c_mid) >= sparsity_tol
                        if i0 == 0                                         % b_ has no constant term 0 .
                            r_.monomial = [zeros(1,r_.dim);r_.monomial];   % Insert zero monomial.
                            r_.coefficient = [c_mid;r_.coefficient];       % Set coefficient for constant term.
                        else
                            r_.coefficient(i0) = c_mid;                    % Set constant term.
                        end
                        c_diff = iv_minus(c_,c_mid);                       % interval evaluation of c_ - c_mid
                                                                           % This is the error update, since c_mid >= sparsity_tol.
                        r_.image = image(r_);                              % Adapt image.
                    else
                        if i0 ~= 0
                            r_ = subtract_constant_term(r_,i0);            % Delete constant term (and compute new range).
                        else                                               % The constant term was already zero. Thus, only the new range is computed.
                            r_.image = image(r_);
                        end
                        c_diff = c_;                                       % Since c_mid < sparsity_tol, the whole interval coefficient c_ is moved to the error interval.                        
                    end                    
                    r_.interval = iv_plus(b_.interval,c_diff);             % interval evaluation of r_interval = b_.interval + (c_-c_mid)
                    
                    if RECMODE == 1                                        % record write mode
                        push_reclist('c_diff', c_diff);                    % store c_diff
                    end
                else % non-verified mode
                    % Determine an approximate midpoint of new constant term (c0 = 0 if i0 is empty)
                    if isfloat_a
                        c_mid = a_ + c0;                        
                    else % isiv_a
                        c_mid = 0.5*(a_.inf+a_.sup) + c0;
                    end
                    if(abs(c_mid) >= sparsity_tol)
                        if i0 == 0                                         % b_ has no constant term 0 .
                            r_.monomial = [zeros(1,r_.dim);r_.monomial];   % Insert zero monomial.
                            r_.coefficient = [c_mid;r_.coefficient];       % Set coefficient for constant term.
                        else
                            r_.coefficient(i0) = c_mid;                    % Set constant term.
                        end
                    else
                        if i0 ~= 0
                            if length(r_.coefficient) > 1                  % Avoid empty Taylor models.
                                r_.coefficient(i0) = [];
                                r_.monomial(i0,:) = [];
                            else
                                r_.coefficient(i0) = 0;
                            end
                        end
                    end
                    r_.interval = INTLAB_ODE_VARS.EMPTYIV;
                    r_.image = INTLAB_ODE_VARS.EMPTYIV;
                end
            else % case "taylormodel + taylormodel"
                if a_.dim ~= b_.dim
                    error('Dimensions of addends must be equal.')
                end
                if a_.order ~= b_.order
                    error('Orders of Taylor expansions must be equal.')
                end
                if any(a_.center ~= b_.center)
                    error('Center points must be equal.')
                end                                
                if ~iv_eq(a_.domain,b_.domain)
                    error('Domains of addends must be equal.')
                end

                r_ = a_;                                                   % Initialize result r_ with a_ .
                M = [a_.monomial;b_.monomial];                             % Concatenate monomial matrices of a_ and b_ .
                %[U,iM,iU] = unique(M,'rows');                             % Delete double entries and get corresponding indices: U(iU)=M, M(iM)=U.
                [U,iU] = unique_rows(M);                                   % Delete double entries and get corresponding indices: U(iU)=M.
                c = [a_.coefficient;b_.coefficient];
                if ODEMODE == 1 % verified mode. Recall that rounding is switched to upwards.
                    c_lower = -accumarray(iU,-c);                          % Accumulate coefficient vector c, by summing up entries according to index iU,
                                                                           % that is summing up all coefficients having the same monomial, i.e.,
                                                                           % the same index in iU, and store that sum at this index in c_lower,
                                                                           % see MATLAB-documentation of the function accumarray.
                    c_upper = accumarray(iU,c);                            % Same accumulation as before for upper interval bound.
                    c_mid = 0.5 * (c_lower + c_upper);                     % Estimate for interval midpoint. 
                    row = (abs(c_mid) >= sparsity_tol);                    % Find row indices i for which |c_mid(i)| >= sparsity_tol.
                    if ~any(row)
                        r_.monomial = zeros(1,r_.dim);
                        r_.coefficient = 0;
                    else
                        r_.coefficient = c_mid(row) ;                      % Store the corresponding coefficients as those of the result r=a*b.
                        r_.monomial = U(row,:);                            % Store also the corresponding monomials.
                        c_lower(row) = -(r_.coefficient-c_lower(row));     % Compute centered lower bound for non-sparse coefficients.
                                                                           % Recall that rounding is upwards.
                        c_upper(row) = c_upper(row)-r_.coefficient;        % Compute centered upper bound for non-sparse coefficients.
                    end
                    
                    % Determine error interval r.interval := I + a.interval + b.interval
                    % where the additional interval I := rest.image is computed according to (4.6), p.105, of [E].
                    if ~isempty(U)
                        rest = r_;
                        rest.monomial = U;
                        coeff.inf = c_lower;
                        coeff.sup = c_upper;
                        rest.coefficient = coeff;
                        rest.image = image(rest);
                        r_.interval = iv_plus( rest.image, iv_plus(a_.interval,b_.interval) ); % r_.interval = rest.image + a_.interval + b_.interval
                        r_.image = image(r_);
                    else
                        rest.image = INTLAB_ODE_VARS.ZEROIV;               % rest.image := zero interval [0,0];
                    end
                    
                    if RECMODE == 1 % record write mode
                        push_reclist('rest.image', rest.image);            % store rest.image
                    end
                else % non-verified mode
                    c_mid = accumarray(iU,c);                              % c_mid is some point in "intval(a.coefficient) + intval(b.coefficient)".
                                                                           % Recall that rounding is to nearest in non-verified mode.
                    row = (abs(c_mid) >= sparsity_tol);                    % Find row indices i for which |c_mid(i)| >= sparsity_tol.
                    if ~any(row)
                        r_.monomial = zeros(1,r_.dim);
                        r_.coefficient = 0;
                    else
                        r_.coefficient = c_mid(row) ;                      % Store the corresponding coefficients as those of the result r = a*b.
                        r_.monomial = U(row,:);                            % Store also the corresponding monomials.
                    end
                    % Just for clarity set .interval- and .image-component of r_ to empty intervals in non-verified mode.
                    r_.interval = INTLAB_ODE_VARS.EMPTYIV; 
                    r_.image = INTLAB_ODE_VARS.EMPTYIV;
                end
            end 
        else % RECMODE == 2, record read mode, only the error interval component is computed !
            r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;                          % Initialize the component result r_ .
            if isiv_a
                c_diff = pull_reclist('c_diff');                           % Read c_lower_diff from record list.
                r_.interval = iv_plus(b_.interval,c_diff);                 % interval evaluation of r_interval = b_.interval + (c_-c_mid)
            else % case "taylor model + taylor model"
                rest_image = pull_reclist('rest.image');                   % Read former rest.image from record.                
                r_.interval = iv_plus( rest_image, iv_plus(a_.interval,b_.interval) );                
            end 
        end 
        r(i,j) = r_;
    end 
end

if (ODEMODE == 1 && rndold ~= 1) || (ODEMODE == 0 && rndold ~= 0)
    setround(rndold)
end

end % function plus

