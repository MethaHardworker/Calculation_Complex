function set_inverse_factorial_(n)
%SET_INVERSE_FACTORIAL_  computes 1./factorial(0:n+1) with INTLAB and stores the 
%                        result globally in INTLAB_ODE_INVERSE_FACTORIAL.                         
%
%   set_inverse_factoria_l(n)


% written  09/09/15     F. Buenger

global  INTLAB_ODE_INVERSE_FACTORIAL

v = [1,1,zeros(1,n)];
F = intval(v,v,'infsup');
for i = 3:n+2
   F(i) = F(i-1)/(i-1); % F(i) = 1/i!
end
INTLAB_ODE_INVERSE_FACTORIAL.INTVAL = F;   % for verified computations  
INTLAB_ODE_INVERSE_FACTORIAL.FLP = mid(F); % For fast non-verified computations 

end % function set_inverse_factorial_