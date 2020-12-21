function set_sqrt_factor_(n)
%SET_SQRT_FACTOR_ computes f(i):=(-1)^(i-1) (2i-3)!!/(i!*2^i), i=0,...,n+1,
%                and stores the result globally in INTLAB_ODE_SQRT_FACTOR.
%
%   set_sqrt_factor_(n)

% written  11/24/15     F. Buenger
% modified 08/08/16     F. Buenger  "intval"-components --> intval-like structures 

global  INTLAB_ODE_SQRT_FACTOR

iv_3.inf = 3; 
iv_3.sup = 3;

v = [1,1/2,zeros(1,n)];
F.inf = v;
F.sup = v; 

for j = 3:n+2    
    i = j-1;
    hlp.inf = F.inf(i);
    hlp.sup = F.sup(i);
    
    hlp = iv_times( hlp, iv_minus( iv_rdivide( iv_3 , 2*i ) , 1) ) ; % f(i) := f(i-1)*(iv_3/(2*i)-1)  =>  f(i) = (-1)^(i-1) (2i-3)!!/(i!*2^i), i = 0,...,n+2
                                                                     % F(j) := f(i) since Matlab indices start with 1 and not with zero
    F.inf(j) = hlp.inf;
    F.sup(j) = hlp.sup;    
end
INTLAB_ODE_SQRT_FACTOR.IV = F;   % for verified computations
INTLAB_ODE_SQRT_FACTOR.FLP = iv_getmidrad(F); % for fast non-verified computations

end % function set_sqrt_factor_