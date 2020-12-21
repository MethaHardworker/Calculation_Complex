function J = vdp_jac(t,y)
%VDP_JAC     Jacobi matrix for Van der Pol's - differential equation

% written  04/19/18     F. Buenger

J = typeadjust([0,1;0,0],y);
J(2,1) = -2.*y(1).*y(2)-1;
J(2,2) = 1-sqr(y(1));

end % function vdp_jac

