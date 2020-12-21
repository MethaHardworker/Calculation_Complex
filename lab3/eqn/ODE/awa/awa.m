function [T,Y,Z] = awa(odefun,jacobimat,tspan,y0,options)
%AWA  implements AWA in MATLAB/INTLAB.  
% 
%   r = awa(odefun,jacobimat,tspan,y0,options)
%
% The function implements the well-known original Pascal program AWA for 
% enclosing the solutions of ordinary initial value problems written by 
% Rudolf Lohner, Institute for Applied Mathematics, Univ. of Karlsruhe.
%
% This MATLAB/INTLAB implementation was written by Florian Buenger,
% Institute for Reliable Computing, Hamburg University of Technology.
%
% input parameters: 
%              
% 1) odefun: This is, as for all MATLAB ODE solvers like ode45, a function
%           handle for the right-hand side f of the ODE y' = f(t,y).
%           The implementation of odefun is due to the user and
%           the output must be a column vector! The implementation
%           can contain INTLAB interval arithmetic. Thus, if, for example
%           f(t,y) = cos(pi*y)+t, then odefun must be implemented 
%           as odefun = @(t,y) cos(intval('pi')*y) + t, since pi is not
%           a floating-point number.
% 
% 2) jacobimat: This is a function handle to the Jacobi matrix J of f
%           with respect to y, i.e., J_{i,j}(t,y) = df_i/dy_j(t,y), 
%           i,j = 1,...,n, where n is the dimension of the ODE.
%           The implementation of the function jacobimat is - like the
%           implementation of odefun - due to the user.
%
% 3) tspan = [t_start,t_end] specifies the integration interval like for 
%           MATLAB ODE solvers. It is allowed that t_end < t_start.
%           Note that t_start and t_end must be floating-point numbers! 
%           In general this is no severe restriction and is either
%           fulfilled anyway or can easily be achieved by a shift and/or
%           scaling of the time variable in the ODE function. 
%           (This may lead to interval parameters in the ODE function.)
%
% 4) y0:    This is the initial value vector, i.e. y(t_start) = y0, which 
%           can also be an interval vector, i.e., y0 = [a,b],
%           a_i <= y_i(t_start) <= b_i, i = 1,...,n.  
%           
% 5) options: This structure has components that are parameters of awa. 
%             The function awaset can be used to create "options" similar 
%             to the MATLAB function "odeset". 
%             The following options can be set:
% 
%    - 'order'    order of Taylor expansion up to which the solution is 
%                 computed in each time step
%    - 'h0'       initial step size (h0 := 0 means that AWA chooses the 
%                 initial step size automatically.) 
%    - 'h_min'    minimum step size (only used for estimating the total 
%                 number of integration steps)
%    - 'EvalMeth' evaluation method (corresponds to AWA variable E_ART)
%                 The following methods can be chosen:
%                 0: (interval vector). After each integration step the 
%                    solution is enclosed by an interval vector (axis 
%                    parallel box) which is used as initial value set for 
%                    the next integration step.
%                 1: (parallelepiped). In each integration step, the so-
%                    lution is enclosed into a parallelepiped which is used 
%                    as initial value set for the next integration step. 
%                 2: (QR decomposition). In each integration step, the 
%                    solution is enclosed into a not necessarily axis 
%                    parallel box which is used as initial value set for 
%                    the next integration step. (In general this works 
%                    quite well for badly conditioned fundamental systems. 
%                    The box is implicitly given by the orthogonal matrix 
%                    Q of a QR-decomposition.)
%                 3: intersection of the inclusion methods 0 and 1.
%                 4: intersection of the inclusion methods 0 and 2.
%    - 'AbsTol'   absolute tolerance for local error 
%                 (corresponds to AWA variable E_A)
%    - 'RelTol'   relative tolerance for local error 
%                 (corresponds to AWA variable E_R) 
%  
% output parameters:
%
% 1) T: Array of time grid points with T(1) = t0 and T(m) = tf where 
%       m := length(T)
% 2) Y: Enclosure of solution at time grid points, y_i(T(j)) is contained  
%       in Y(j,i) = [Y(j,i).inf,Y(j,i).sup], j = 1,...,m, i = 1,...,n.
% 3) Z: Enclosusure of Taylor coefficients. If t is an intermediate time
%       point in (t_j, t_{j+1}), then y_i(t) is contained in the interval
%       evaluation of
%       
%       sum_{k=0}^p [Z.inf(j,i,k),Z.sup(j,i,k)]*([t-t_j]/[t_{j+1}-t_j])^k, 
%      
%       where p:= order.
%
% The implementation is based on:
%    
% [L]   R. Lohner, Einschliessung der Loesung gewoehnlicher Anfangs- und
%       Randwertaufgaben und Anwendungen, Diss. Univ. Karlsruhe, 1988.
%
% [AWA] R. Lohner, Pascal implementation according to [L], available from 
%       http://www2.math.uni-wuppertal.de/~xsc/xsc/pxsc_software.html#awa  

% written  07/18/17     F. Buenger
% modified 04/09/18     F. Buenger, additional output parameter Z for continuous solution inclusion

% taylorcoeffinit(); % <--- This line moved to "startintlab.m". 
dummy.inf = 0;
dummy.sup = 0;
dummy = taylorcoeffinit(dummy); 
y0 = intval2iv(intval(y0)); % Convert y0 to interval structure with components y0.inf and y0.sup.

if nargin < 5
    options = [];
end                                         
[T,Y,Z] = awaTC(dummy,odefun,jacobimat,tspan,y0,options);
Y = iv2intval(Y);
Z = iv2intval(Z);

end % function awa