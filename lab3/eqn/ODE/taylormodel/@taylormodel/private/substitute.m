function r = substitute(a,t_curr,t_next)
%SUBSTITUTE  substitution of the time component in Taylor models of type 1.
%            Substitute the time variable [(n+1)-th variable] in all a(i) by t_next.
%
%   r = substitute(a,t_curr,t_next)

% written  11/13/15     F. Buenger
% modified 12/10/15     F. Buenger  check rounding 'upwards' instead of 'to nearest', code optimization
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures

% Remark: This function is never called during the execution of the integral operator function integral_operator.m.
%         Therefore it is not included in the record feature.

e = 1e-30;
if 1+e > 1   % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

S_a = size(a);
r = a; % Initialize result r with a. (Just preallocation of storage)

n = a(1).dim; 

% Determine interval enclosure of h := t_next - t_curr. Recall that rounding is upwards.
% Note that h >= 0, since  t_next >= t_curr.
h.inf = -(t_curr - t_next);
h.sup = t_next - t_curr;

e = (0:a(1).order)'; % Column vector of all potential exponents that may occur in monomials of a(i,j).

% Determine interval H := h.^e and note carefully that e>=0, and h >= 0 as t_next >= t_curr.
setround(-1); % Switch to rounding downwards.
H.inf = (h.inf).^e;  
setround(1);  % Switch back to rounding upwards.
H.sup = (h.sup).^e;

sparsity_tol = get_sparsity_tol; % Interval coefficients c with |mid(c)| < sparsity_tol are set to zero to keep the Taylor polynomial sparse
                                 % the error caused by that will be transferred to the error interval r.interval.

for i = 1:S_a(1)
    for j = 1:S_a(2)
        a_ = a(i,j); %                 
        
        % Compute c = a_.coefficient .* h.^(a_.monomial(:,n)). This means: h is substituted in the polynomial part of a_.
        N = a_.monomial(:,n); 
        H_.inf = H.inf(N+1);
        H_.sup = H.sup(N+1);
        c = iv_times(a_.coefficient, H_); % c = a_.coefficient .* h.^(a_.monomial(:,n))   
        
        % Clear Exponents of time variable
        M = a_.monomial;
        M(:,n) = 0 ; % Clear dependence on time variable which has index n.
        
        %[U,iM,iU] = unique(M,'rows'); % Delete possible duplicate entries and get corresponding indices: U(iU) = M, M(iM) = U.
        [U,iU] = unique_rows(M); % Delete duplicate entries and get corresponding indices: U(iU) = M.
        
        % Accumulate coefficients of duplicate monomials.
        % Recall that rounding is upwards. 
        c_lower = -accumarray(iU,-c.inf); % Accumulate coefficient vector c.inf by summing up entries according to index iU.
        c_upper = accumarray(iU,c.sup);   % Same accumulation as before for upper interval bound.
        
        c_mid = 0.5 * (c_lower + c_upper);  % Calculate approximate midpoint of c. It is not needed that c_mid is the exact midpoint. 
        row = (abs(c_mid) >= sparsity_tol); % Find row indices i for which |c_mid(i)| >= sparsity_tol.
        
        r_ = a_; % Initialize component result with a.
        if ~any(row) % All coefficients are below the sparcity tolerance. Therefore, the polynomial part of the result r_ is the zero polynomial.
            r_.monomial = zeros(1,r_.dim);
            r_.coefficient = 0;  % Prevent result from being an "empty" Taylor model.
        else         % Not all coefficients are below the sparsity tolerance.
            r_.coefficient = c_mid(row) ;                  % Store the 'large' coefficients as those of the result r_.
            r_.monomial = U(row,:);                        % Store also the corresponding monomials.
            c_lower(row) = -(r_.coefficient-c_lower(row)); % Compute centered lower bound for non-sparse coefficients. 
                                                           % Recall that rounding is upwards.
            c_upper(row) = c_upper(row)-r_.coefficient;    % Compute centered upper bound for non-sparse coefficients.                                                               
                                                           % Recall that rounding is upwards.                                                                      
        end
        
        % Accumulate rounding errors and sparsity terms and add them to the error interval r_.interval.
        rest = r_;
        rest.monomial = U;
        coeff.inf = c_lower;
        coeff.sup = c_upper;
        rest.coefficient = coeff;
        r_.interval = iv_plus(a_.interval,image(rest)); % r_.interval :=  a_.interval + image(rest)
        
        % Finally compute the image of r_.
        r_.image = image(r_); 
        
        r(i,j) = r_;
    end % j
end % i
    

if rndold ~= 1
    setround(rndold)
end
end % function substitute