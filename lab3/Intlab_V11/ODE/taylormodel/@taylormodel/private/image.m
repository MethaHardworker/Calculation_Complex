function res = image(a,dummy)
%IMAGE         image of taylormodel-like structures
%
%   res = image(a)
%
% The input parameter a must contain a taylormodel-like structure.
% It must not necessarly be a Taylor model.
% For example, the component array a.coefficient may consist of intervals.

% written   08/31/15     F. Buenger
% modified  10/28/15     F. Buenger  LDB method implemented
% modified  11/23/15     F. Buenger  matrix input (componentwise evaluation)
% modified  12/10/15     F. Buenger  check rounding 'upwards' instead of 'to nearest', code optimization
% modified  02/10/16     F. Buenger  "intval"-components --> intval-like structures  

global INTLAB_ODE_OPTIONS

% Computation of the image by using LDB (linear dominated bounder)
if nargin == 1 && ~isempty(INTLAB_ODE_OPTIONS) &&  strcmp(INTLAB_ODE_OPTIONS.bounder,'LDB')
  res = LDB(a);
  return;  
end

% naive (but fast) interval computation of the image
e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

S_a = size(a);
res_inf = zeros(S_a);
res_sup = res_inf;

for i = 1:S_a(1)
    for j = 1:S_a(2)
        a_ = a(i,j);                
        if a_.type == 1 % ODE-Taylor models have a domain D := [-1,1]^(n-1) x [t,t+h] and center point (0,...,0,t).
                        % This is used to compute the image more efficiently.
            n = a_.dim;
            M = a_.monomial(:,1:n-1); % exponents for [-1,1]x...x[-1,1], the (n-1)-times Cartesian product of [-1,1]
            N = a_.monomial(:,n);     % exponents for time domain [0,h]
            c = a_.coefficient;
            if isfloat(c)
                c_inf = c;
                c_sup = c;
            else
                c_inf = c.inf;
                c_sup = c.sup;
            end
            i0 = all(a_.monomial == 0 , 2); % Get index of constant term (monomial = 0...0).
            if any(i0)                      % a_ has non-zero constant term.
                c0_inf = c_inf(i0);         % Store upper bound of c0.
                c0_sup = c_sup(i0);         % Store lower bound of c0.
                c_inf(i0) = 0;              % Set coefficient of constant term in c_inf to zero. The information is stored in c0_inf and will be added separately.
                c_sup(i0) = 0;              % Set coefficient of constant term in c_sup to zero. The information is stored in c0_sup and will be added separately.
            else                            % a_ has zero constant term.
                c0_inf = 0;
                c0_sup = 0;
            end

            idx_odd = any( 2*round(M/2)-M , 2); % Indicator matrix for monomials M(i,:) containing odd exponents.
            % In this case [-1,1].^M(i,:) = [-1,1]. Otherwise [-1,1].^M(i,:) = [0,1] if M(:,i) is not the zero vector.
            % If M(i,:) is the zero vector, then [-1,1].^M(i,:) ) = [1,1].
            
            % In the following we compute a lower bound u of [-1,1]^M.*c with u(i0):=0 if the constant term index i0 is not empty.
            % Note carefully that u <= 0 componentwise.
            % Note also, if M(i,:) is the zero vector for some i<>i0 and c_inf(i)>0, then u(i) = 0 < c_inf(i) = inf(prod([-1,1].^M(i,:))*c(i)).
            % But since N(i) <> 0 in this case we have inf( prod([-1,1].^M(i,:))*c(i)*[0,h]^N(i) ) = inf(c_inf(i)*[0,h]) = 0 = u(i)*h.
            mu = max(abs(c_inf(idx_odd)),abs(c_sup(idx_odd)));
            u = min(c_inf,0);
            u(idx_odd) = -mu;
            
            % In the following we compute an upper bound v of [-1,1]^M.*c with v(i0):=0 if the constant term index i0 is not empty.
            % Note carefully that v >= 0 componentwise.
            % Like before, if M(i,:) is the zero vector for some i<>i0 and c_sup(i)<0, then v(i) = 0 > c_sup(i) = sup(prod([-1,1].^M(i,:))*c(i)).
            % Again N(i) <> 0 in this case causes sup( prod([-1,1].^M(i,:))*c(i)*[0,h]^N(i) ) = sup(c_sup(i)*[0,h]) = 0 = v(i)*h.
            v = max(c_sup,0);
            v(idx_odd) = mu;
            
            % Recall that the rounding mode is switched to upwards!
            h = a_.domain.sup(n) - a_.center(n); % Compute time step size h rounded upwards.
            
            H = [1;cumprod(h*ones(a_.order,1))]; 
            %H = h.^((0:a_.order).'); % slower than previous cumprod, but of higher accuracy in rounding upwards mode
            
            H = H(N+1);              % four times faster implementation of H = h.^N;
            y = v'*H + c0_sup;       % y is an upper bound for the range of a_. Here we use that v>=0 and H>=0 componentwise.
            x = -((-u')*H - c0_inf); % x is a lower bound for the range of a_. Here we use that u<=0 and H>=0 componentwise.
            res_inf(i,j) = x;
            res_sup(i,j) = y;
            
% Subdivision of each domain interval [-1,1] into [-1,0] and [0,1]. This gives a slightly better inclusion of the image which results 
% in slightly tighter inclusion of the ode solution. On the other hand the computation is more expensive. 
            
%             M_odd = 2*round(M/2)-M; % indicator matrix for odd entries of M
%             p = 2^(n-1); 
%             D = double(dec2bin(0:p-1)-'0').'; % D is an (n-1)x2^(n-1)-matrix containing in its columns all vectors d of length n-1
%               % with components 0 or 1. Each column vector d stands for a specific subinterval (row) vector z=z_d of [-1,1]^(n-1):=[-1,1]x...x[-1,1],
%               % namely z(i) = [-1,0] if d(i) = 1, and z(i) = [0,1] if d(i) = 0. Note that all z_d build a partition of [-1,1]^(n-1).
%             I = M_odd*D; 
%             I_odd = logical(2*round(I/2)-I); % Let d be a column of D , J be the corresponding column of I_odd, and z:=z_d as above.
%               % Moreover let w:=M(k,1:n-1). First suppose that w is not the zero vector.
%               % Then prod(z.^w) = [-1,0] if and only if J(k) = 1, and prod(z.^w) = [0,1] if and only if J(k) = 0.  
%               % The special case where w is the zero vector, i.e. prod(z.^w) = [1,1], leads to J(k) = 0.          
%             U = repmat(min(c_inf,0),1,p);
%             W = repmat(min(-c_sup,0),1,p); 
%             U(I_odd) = W(I_odd);
%             
%             V = repmat(max(c_sup,0),1,p);
%             W = repmat(max(-c_inf,0),1,p); 
%             V(I_odd) = W(I_odd);
% 
%             h = a_.domain.sup(n) - a_.center(n); % Compute time step size h rounded upwards.
%             
%             H = [1;cumprod(h*ones(a_.order,1))];
%             H = H(N+1); % four times faster implementation of H = h.^N;
%             
%             Y = V'*H + c0_sup; 
%             X = -((-U')*H - c0_inf); 
%             
%             res_inf(i,j) = min(X);
%             res_sup(i,j) = max(Y);
                            
        elseif a_.type == 0 % default type               
            n = a_.dim;
            M = a_.monomial;
            coeff = a_.coefficient;
            order = a_.order;
            if isfloat(coeff)
                c.inf = coeff;
                c.sup = coeff;
            else
                c = coeff;
            end
            i0 = all(a_.monomial == 0 , 2); % Get index of constant term (monomial = 0...0).
            if any(i0)                      % a_ has non-zero constant term.
                c0_inf = c.inf(i0);         % Store upper bound of c0.
                c0_sup = c.sup(i0);         % Store lower bound of c0.
                c.inf(i0) = 0;              % Set coefficient of constant term in c_inf to zero. The information is stored in c0_inf and will be added separately.
                c.sup(i0) = 0;              % Set coefficient of constant term in c_sup to zero. The information is stored in c0_sup and will be added separately.
            else                            % a_ has zero constant term.
                c0_inf = 0;
                c0_sup = 0;
            end
            idx_even = ( 2*round(M/2) - M == 0 ); % Indicator matrix for even exponents M(i,j)
            
            % Determine nonnegative float vectors h_inf >= 0, h_sup >= 0 such that the centered interval vector a_.domain-a_.center 
            % is contained in [-h_inf,h_sup]. 
            % Recall that rounding is switched to upwards.
            h_inf = a_.center - a_.domain.inf; 
            h_sup = a_.domain.sup - a_.center;  
            
            H_inf = cumprod([ones(1,n); repmat(h_inf',order,1)]); % [ h_inf.^0 ; h_inf.^1 ; h_inf.^2 ;...; h_inf.^order]
            H_sup = cumprod([ones(1,n); repmat(h_sup',order,1)]); % [ h_sup.^0 ; h_sup.^1 ; h_sup.^2 ;...; h_sup.^order]            
            H_max = max(H_inf,H_sup); 

            idx_shift = (0:n-1)*(order+1);       % index shift for linear indexing of H_sup, H_inf.  
            M_shift = M+1+idx_shift;
            HM_max = H_max(M_shift); 

            % a) [-h_inf,h_sup].^k = [-h_inf.^k,h_sup.^k] for odd k > 0
            %      Note that here [-h_inf,h_sup].^k := intervalhull( {x.^k ; x in [-h_inf,h_sup]} ).          
            %
            % b) [-h_inf,h_sup].^k = [0,max(h_inf,h_sup).^k] for even k > 0

            
            HM_sup = H_sup(M_shift);             % Compute upper bound in a) which is h_sup.^k.
            HM_sup(idx_even) = HM_max(idx_even); % Compute upper bound in b) which is max(h_inf,h_sup).^k
            
            HM_inf = H_inf(M_shift);             % Compute absolute value of lower bound in a) which is h_inf.^k.
            HM_inf(idx_even) = 0;                % Set lower bound in b) which is zero.
               
            % Define P as an enclosure of "(a.domain-a.center)'.^a.monomial" 
            % (accept constant term monomial, if existent).
            P.inf = -HM_inf;
            P.sup = HM_sup;
            
            P = iv_prod(P,2);                    % interval rowwise prod()
            r = iv_dotprod(c,P);                 % Multiply the coefficients with the h-powers and sum everything up to obtain the image (shifted minus the constant term).
            
            % Recall that [c0_inf,c0_sup] encloses the constant term of the polynomial part, which was set to zero 
            % in the coefficient vector c in the beginning. This has to be adapted in the final result as follows:
            res_inf(i,j) = -((-r.inf)-c0_inf);   % Recall that rounding is upwards.  
            res_sup(i,j) = r.sup + c0_sup;
            
        end  % a_.type == 1
    end % j
end % i

res.inf = res_inf;  
res.sup = res_sup;

if rndold ~= 1
    setround(rndold)
end
end % function image