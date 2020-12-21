function J = lorenz_jac(t,y)
% LORENZ_JAC   Jacobi matrix function of ode-function for Lorenz attractor:
%
%     dy(1) = sigma.*(y(2)-y(1));
%     dy(2) = (rho-y(3)).*y(1)-y(2);
%     dy(3) = y(1).*y(2)-beta.*y(3);

% written  07/24/18   F. Buenger

sigma = 10;
beta_ = iv2intval(iv_rdivide(-8,3)); % beta_ := -beta := -8/3 in verified arithmetic
rho = 28;

J(1:3,1:3) = y(1); % cheap preallocation

J(1,1) = typeadjust(-sigma,y);
J(1,2) = typeadjust(sigma,y);
J(1,3) = typeadjust(0,y);

J(2,1) = rho-y(3);
J(2,2) = typeadjust(-1,y);
J(2,3) = -y(1);

J(3,1) = y(2);
J(3,2) = y(1);
J(3,3) = typeadjust(beta_,y);

end % function lorenz_jac

