function dydt = vdp_fun(t,y)
%VDP_FUN     ODE-function for Van der Pol's - differential equation

% written  04/19/18     F. Buenger

dydt = [y(2);(1-sqr(y(1))).*y(2)-y(1)];

end % function vdp_fun

