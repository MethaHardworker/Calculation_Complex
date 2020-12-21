function J = r3body_jac(t,y)
% R3BODY_JAC  Jacobi matrix of the following ODE function. 
%
% An example of a nonstiff system is the system of equations describing
% the motion of a rigid body without external forces
%
% y_1' = y_2*y_3           y_1(0) = 0
% y_2' = -y_1*y_3          y_2(0) = 1
% y_3' = -0.51*y_1*y_2     y_3(0) = 1

% written  08/29/17   F. Buenger

persistent c;

if isempty(c) 
    c = -intval('0.51');
end

J = typeadjust(zeros(3),y);

J(1,2) = y(3);
J(1,3) = y(2);

J(2,1) = -y(3);
J(2,3) = -y(1);

J(3,1) = c.*y(2);
J(3,2) = c.*y(1);

end % r3body_jac

