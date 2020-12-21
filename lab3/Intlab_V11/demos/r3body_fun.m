function dydt = r3body_fun(t,y)
%R3BODY_FUN   
%
% An example of a nonstiff system is the system of equations describing
% the motion of a rigid body without external forces:
%
% y_1' = y_2*y_3           y_1(0) = 0
% y_2' = -y_1*y_3          y_2(0) = 1
% y_3' = -0.51*y_1*y_2     y_3(0) = 1

% written  04/20/18   F. Buenger

persistent c;

if isempty(c) 
    c = -intval('0.51');
end
  
dydt = y;
dydt(1) = y(2) .* y(3);
dydt(2) = -y(1) .* y(3);
dydt(3) = c .* y(1) .* y(2);

end % function r3body_fun

