function r = rescale(a,s)
%RESCALE  Taylor model rescaling r(x,t) = a(x./s,t), where the column vector s 
%         contains the positive scaling factors for the space variables x1,...,x_n. 
%
%   r = rescale(a,s)

% written  06/14/17     F. Buenger

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end
sparsity_tol = get_sparsity_tol;
S_a = size(a);
n = a(1).dim-1; 
q = iv_rdivide(1,s'); % q = 1./s > 0 since s > 0 
r = a;                % Initialize result r with a. Just cheap storage preallocation.
for i = 1:S_a(1)
    for j = 1:S_a(2)
        a_ = a(i,j);
        r_ = r(i,j);                     % = a(i,j)
        M = a_.monomial;
        N = M(:,1:n);                    % Cut off powers for time variable so that N only contains the powers for the space variables.  
        c = a_.coefficient; 
        u.sup = prod(q.sup.^N,2);        % Recall that q.sup > 0, N >= 0 and that rounding is upwards,
        setround(-1)
        u.inf = prod(q.inf.^N,2);        % Recall that q.inf  > 0 and M >= 0.       
        setround(1)
        v = iv_times(u,c);               % This is the rescaling of the Taylor model coefficients: v := ((1./s').^M).*c.
        c_lower = v.inf;                 % according compressed lower coefficient bounds
        c_upper = v.sup;                 % according compressed upper coefficient bounds
        c_mid = 0.5 * (c_lower+c_upper); % approximate coefficient midpoints of the product polynomial
        
        row = (abs(c_mid) >= sparsity_tol);
        % Find row indices i for which |c_mid(i)| >= sparsity_tol and degree(r.monomial(i)) <= r_.order
        if ~any(row)
            r_.monomial = zeros(1,r_.dim);
            r_.coefficient = 0;
        else
            r_.coefficient = c_mid(row);                   % Store the corresponding coefficients as those of the result r_ = a_ * b_
            r_.monomial = M(row,:);                        % Store also the corresponding monomials.
                                                           % Recall that rounding is upwards.
            c_lower(row) = -(r_.coefficient-c_lower(row)); % Compute centered lower bound for non-sparse coefficients.
            c_upper(row) = c_upper(row)-r_.coefficient;    % Compute centered upper bound for non-sparse coefficients.
        end
        % Comute error interval
        rest = r_;
        rest.monomial = M;
        coeff.inf = c_lower;
        coeff.sup = c_upper;
        rest.coefficient = coeff;
        rest.image = image(rest);        
        r_.interval = iv_plus(rest.image,a_.interval);        
        
        r_.image = image(r_); % Finally compute the image of the result r_.
        r(i,j) = r_;        
    end % j
end %i

if rndold ~= 1 
    setround(rndold)
end
end % function rescale