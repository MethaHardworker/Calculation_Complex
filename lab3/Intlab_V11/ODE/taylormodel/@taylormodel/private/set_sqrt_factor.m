function set_sqrt_factor(n)
%SET_SQRT_FACTOR computes f(i):=(-1)^(i-1) (2i-3)!!/(i!*2^i), i=0,...,n+1,
%                and stores the result in the global constant INTLAB_ODE_TCOEFF
%
%   set_sqrt_factor(n)

% written  11/24/15     F. Buenger
% modified 08/08/16     F. Buenger  "intval"-components --> intval-like structures 
% written  06/08/18     F. Buenger  new global constant INTLAB_ODE_TCOEFF

global INTLAB_ODE_TCOEFF

iv_3.inf = 3; 
iv_3.sup = 3;

x.inf = [1;1/2;zeros(n,1)];
x.sup = x.inf; 

dummy.inf = 0;
dummy.sup = 0;

for j = 3:n+2    
    i = j-1;
    dummy.inf = x.inf(i);
    dummy.sup = x.sup(i);
    
    dummy = iv_times( dummy, iv_minus( iv_rdivide( iv_3 , 2*i ) , 1) ) ; % f(i) := f(i-1)*(iv_3/(2*i)-1)  =>  f(i) = (-1)^(i-1) (2i-3)!!/(i!*2^i), i = 0,...,n+2
                                                                         % x(j) := f(i) since MATLAB indices start with 1 and not with zero
    x.inf(j) = dummy.inf;
    x.sup(j) = dummy.sup;    
end
INTLAB_ODE_TCOEFF.SQRT_IV = x;                % for verified computations
INTLAB_ODE_TCOEFF.SQRT_FLP = iv_getmidrad(x); % for fast non-verified computations

end % function set_sqrt_factor