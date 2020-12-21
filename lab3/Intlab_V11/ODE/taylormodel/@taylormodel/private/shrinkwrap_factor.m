function [q,err] = shrinkwrap_factor(f,A)
%SHRINKWRAP_FACTOR   computes the shrink wrap factor for the Taylor model vector f.
%
%   res = shrinkwrap_factor(a,A)
%
% The implementation is based on 
% [Bue] F. Buenger, "Shrink wrapping for Taylor models revisited", Numerical Algorithms 78(4), pp. 1001-1017, 2018

% written  07/04/16     F. Buenger

global INTLAB_ODE_VARS
global INTLAB_ODE_OPTIONS

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

q_max = 1.01;      % bound for the shrink wrap factors. Feel free to change this constant!
q_tol = 1E-12;     % relative "error tolerance" for shrink wrap factor iteration. Feel free to change this constant! 
s_buffer = 1+2E-1; % safety margin for estimated shrink wrap factors. Feel free to change this constant! 
iter_max = 3;      % maximum number of iterations for estimating the shrink wrap factors. Feel free to change that value!

err = 0; % no error 
n = f(1).dim-1;
                  
if INTLAB_ODE_OPTIONS.blunting
     A = blunt(A);
end

R = inv(A); % Compute an approximate inverse R of A.

I = get_interval(f); % error interval of f

r = max(abs(I.inf),abs(I.sup)); % I is contained in the interval vector [-r,r] which symmetric to zero.
rs = abs(R)*r;                  % Since rounding is upwards we have rs=fl(abs(R)*r)>=abs(R)*r.
                                % rs.*Q contains the transformed interval R*(r.*Q) and therefore also R*I

% Compute a verified inclusion df of the Jacobian of f.   

df(1:n,1:n) = f(1); % storage preallocation for df.
for i = 1:n
    f_ = f(i);
    for j = 1:n
        M_j = f_.monomial(:,j); % exponents of variable x_j
        idx = (M_j > 0);
        dfi_dxj = f_; % Initialize partial derivative dfi/dxj.
        if any(idx)
            M_j = M_j(idx);
            % Differentiate polynomial part of f_ with respect to x_j:  df_/dxj = (c * x_j^k )' = (k*c) * x_j^(k-1) .
            dfi_dxj.monomial = f_.monomial(idx,:);
            dfi_dxj.monomial(:,j) = M_j-1; % Differentiation with respect to x_j decreases the nonzero exponents of x_j by one.
            dfi_dxj.coefficient = iv_times(f_.coefficient(idx),M_j);
        else % f(i) does not depend on x_j
            dfi_dxj.monomial = zeros(1,n+1);
            dfi_dxj.coefficient = 0;            
        end
        % Error intervals and images are not relevant in this context and are set to empty intervals.
        dfi_dxj.interval = INTLAB_ODE_VARS.EMPTYIV; 
        dfi_dxj.image = INTLAB_ODE_VARS.EMPTYIV;
        dfi_dxj.type = 0;        % Switch to general, non-ode type. This is necessary since later on the domain [-1,1]^n will be enlarged. 
        dfi_dxj.center(n+1) = 0; % Center of time coordinate, which is superfluous in this context, is set to zero.
        df(i,j) = dfi_dxj;        
    end    
end

% Compute a verified inclusion df of the Jacobian of g = R*f-id. 
% This is done by computing a verified inclusion of R*df-eye(n).                
dg = df; % storage preallocation for dg.
for i = 1:n
  % Compute R*df.  
  for j = 1:n      
     dgi_dxj = df(1,j); 
     dgi_dxj.coefficient = iv_times(R(i,1),dgi_dxj.coefficient); % dgi_dxj :=  R(i,1) * df(1,j)
     for k = 2:n
         dgi_dxj = axpy(R(i,k),df(k,j),dgi_dxj);  % dgi_dxj := dgi_dxj + R_ik * dfk_dxj
     end     
     % Compute R*df-I.
     if i == j % 
       idx = find(all(dgi_dxj.monomial == 0,2),1);                         % index of constant term 
       dgi_dxj.coefficient.sup(idx) = dgi_dxj.coefficient.sup(idx) - 1;    % Recall that rounding is upwards
       dgi_dxj.coefficient.inf(idx) = -(1 - dgi_dxj.coefficient.inf(idx)); % verified lower bound for dgi_dxj.coefficient.inf(idx) - 1
     end
     dg(i,j) = dgi_dxj;
  end  
end
                         
% Calculate an approximate (non-verified) shrink wrap factor
q = 1+rs; % starting value for non-verified, approximate shrink wrap factor    
improve = true;
iter = 0; 

% estimation of shrink wrap factors q
while improve && iter < iter_max
    s = zeros(size(rs));
    q_old = q; 
    for i = 1:n
        for j = 1:n
          dgi_dxj = dg(i,j); 
          % Set domain of dgi_dxj to [-q,q];
          dgi_dxj.domain.inf = [-q;0]; % The trailing zero is the lower bound for the time coordinate which is superfluous in this context.
          dgi_dxj.domain.sup = [q;0];  % The trailing zero is the upper bound for the time coordinate which is superfluous in this context.
          dgi_dxj.image = image(dgi_dxj);        
          im_maxabs = max(abs(dgi_dxj.image.inf),abs(dgi_dxj.image.sup)); % upper bound for partial derivative dgi/dx_j in absolute value on [-q,q].
          s(i) = s(i) + im_maxabs * (q(j)-1);
        end
        q(i) = 1 + rs(i) + s(i); % The newly computed shrink wrap factor q(i) is already used for computing q(i+1),...,q(n).
        if q(i) > q_max % too large shrink wrap factor 
            q = 0; 
            err = 1;
            return;
        end
    end
    
%     q = 1 + rs + s; % The newly computed shrink wrap factor q(i) is already used for computing q(i+1),...,q(n).
%     if any(q > q_max) % too large shrink wrap factor
%         q = 0;
%         err = 1;
%         return;
%     end

    improve = any((q-q_old)./q > q_tol);
    iter = iter + 1;
end

s = s_buffer * s;    % The estimated vector s is slightly increased (safety margin).
q = 1 + rs + s;      % Recall that rounding is upwards so that q >= 1+rs+s in exact arithmetic. 
dg_bound = zeros(n); % Initialization of the array which will contain the verified upper bounds for |dgi/dxj| on [-q,q].

% verification of the shrink wrap factors q
for i = 1:n
    for j = 1:n
        dgi_dxj = dg(i,j);
        % Set domain of dgi_dxj to [-q,q];
        dgi_dxj.domain.inf = [-q;0];
        dgi_dxj.domain.sup = [q;0];
        dgi_dxj.image = image(dgi_dxj);
        dg_bound(i,j) = max(abs(dgi_dxj.image.inf),abs(dgi_dxj.image.sup)); % upper bound for partial derivative dgi/dx_j in absolute value on [-q,q].
    end    
    z = dg_bound(i,:)*(q-1);
    if z > s(i)  % Recall that rounding is upwards so that dg_bound(i,:)*(q-1) <= s(i) in floating-point arithmetic implies
                 % dg_bound(i,:)*(q-1) <= s(i) in exact arithmetic. Note also that fl(q-1) >= rs+s.
        q = 0;
        err = 2; % shrink wrap factors could not be verified.
        return;
    end
end

if rndold ~= 1 
    setround(rndold)
end

end % function shrinkwrap_factor

