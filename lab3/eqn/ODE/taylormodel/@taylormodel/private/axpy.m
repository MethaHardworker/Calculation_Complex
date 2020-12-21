function r = axpy(alpha,a,b)
%axpy  computes alpha*a+b, where alpha is a real number or an interval
%      and a and b are Taylor models or pseudo Taylor models 
%      with interval coefficients.
%
%   r = fma(alpha,a,b)

% written  07/12/16     F. Buenger

global INTLAB_ODE_VARS

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end

r = a; % Initialize result r with a.

M = [a.monomial;b.monomial]; % Concatenate monomial matrices of a and b.
[U,iU] = unique_rows(M); % Delete double entries and get corresponding indices: U(iU)=M.

% Set interval coefficients for alpha * a + b, recall that rounding is upwards 
coef = iv_times(alpha,a.coefficient); 
coef.inf = [coef.inf;b.coefficient.inf];
coef.sup = [coef.sup;b.coefficient.sup];
coef.inf = -accumarray(iU,-coef.inf); % Sum up all coefficients having the same monomial. Recall that rounding is upwards
coef.sup = accumarray(iU,coef.sup); 

r.coefficient = coef; 
r.monomial = U; 

% Error intervals and images are not relevant in this context and are set to empty intervals.
r.interval = INTLAB_ODE_VARS.EMPTYIV; 
r.image = INTLAB_ODE_VARS.EMPTYIV;

if rndold ~= 1
    setround(rndold)
end

end % function axpy

