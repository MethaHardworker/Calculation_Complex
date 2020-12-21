function set_inverse_factorial(n)
%SET_INVERSE_FACTORIAL  computes 1./factorial(0:n+1) and stores the 
%                       result in the global constant INTLAB_ODE_TCOEFF.                         
%
%   set_inverse_factorial(n)

% written  09/09/15     F. Buenger
% written  06/08/18     F. Buenger  intval --> interval-like structures

global INTLAB_ODE_TCOEFF

x.inf = [1;1;zeros(n,1)]; % 0! = 1, 1! = 1, the rest is initialized with zeros  
x.sup = x.inf;
dummy.inf = 0;
dummy.sup = 0;

for i = 3:n+2
   dummy.inf = x.inf(i-1);  
   dummy.sup = x.sup(i-1);  
   dummy = iv_rdivide(dummy,(i-1)); % x(i) = 1/i!
   x.inf(i) = dummy.inf;
   x.sup(i) = dummy.sup;
end

INTLAB_ODE_TCOEFF.INVFACT_IV = x;                % for verified computations  
INTLAB_ODE_TCOEFF.INVFACT_FLP = iv_getmidrad(x); % for non-verfied computations 

end % function set_inverse_factorial