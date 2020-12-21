function [T,Y,Yr] = verifyode(odefun,tspan,y0,options)
% VERIFYODE   verified ODE-solver for initial vaue problems  
%               
%   y'(t) = odefun(t,y(t))  with t in (tspan(1),tspan(end)) 
%   y(tspan(1)) = y0.
%
%   call:  [T,Y,Yr] = verifyode(odefun,tspan,y0,options)
%
% This verified ODE-solver is based on Taylor models.   
% The four input parameters "odefun, tspan, y0, options" have the same meaning as for other (non-verified) 
% MATLAB ODE-solvers, e.g. ode45. 
% 
% Input parameters:  
%
% 1) "odefun" 
%     For a given ODE y′ = f(t,y), "odefun" is a function handle to f: R^n x [t0,tf] -> R^n. 
%     Here t is the time variable ranging between t0 := tspan(1) and tf := tspan(end), and n is the dimension of the ODE. 
%     This means f(t,y) returns the column vector of length n containing the components f(1),...,f(n). 
%     In addition, the call f(t,y,i) returns the single components f(i).
%
%     Important note: It is assumed that each component function f_i fulfills the following:
%        a) f_i consists only of arithmetic operations +,-,*,/ and of standard functions
%           like exp,log, sin, cos, tan,... ,
%           see folder @taylormodel for the complete list of implemented standard functions.
%        b) f_i is (m+1)-times continuously differentiable, where m is the user specified order of the 
%           used Taylor models (= maximum degree of the multivariate polynomial part).
%        c) f_i is implemented rigorously. For example, non-floating-point constants like pi
%           must be implemented by intval('pi'), see @intval/intval. 
%  
% 2) "tspan" contains the starting point and the endpoint of integration, i.e. tspan = [t0,tf].
%     (Contrary to MATLAB ode-solvers, tspan is not allowed to contain additional intermediate points.)  
%     Both t0 and tf must be floating point numbers. In particular this means that if, for example, 
%     a non-floating-point starting time t0, say t0 = 0.1, is necessary, then the ode-function odefun   
%     must be changed/rescaled by the user so that the starting time t0 is moved to a floating-point number. 
%     In general, such a rescaling is easily done. In many/most practical cases we simply have t0 = 0. 
%
% 3) y0 = (y0(1);...;y0(n)) contains the starting values. This can be a float or an interval vector.  
%    (A row vector will automatically be transformed into a column vector.)
%
% 4) "options" contains additional parameters. Analogously to other MATLAB ode-solvers, "options"  
%     can be set by using the function verifyodeset: options = verifyodeset('name1',value1,'name2',value2,...).
%     The following parameters can be used (see also documentation of verifyodeset.m) : 
%
%     - 'order'         degree of Taylor polynomials 
%     - 'sparsity_tol'  threshold amount for coefficients c of a multivariate polynomial. 
%                       If |c| < sparsity_tol, then the coefficient is removed from the polynomial and the (small) error 
%                       caused by this is transferred to the error interval of the Taylor model.                        
%     - 'loc_err_tol'   tolerance for local error for step size control
%     - 'h_min'         minimum time step size
%     - 'h0'            initial step size. If h0 = 0, then the initial step size is computed automatically.
%     - 'shrinkwrap'    flag for "shrink wrapping", 0 = off, 1 = on 
%     - 'precondition'  Preconditioning of Taylor models, 
%                         0: off
%                         1: QR preconditioning
%                         2: parallelepiped preconditioning   
%     - 'blunting'      Blunting of ill-conditioned matrices during shrink wrapping or parallelepiped preconditioning, 0 = off, 1 = on 
%     - 'bounder'       method for bounding the image of the polynomial part of a Taylor model.
%                         'NAIVE': use interval arithmetic directly (default)
%                         'LDB'  : Linear Dominated Bounder, see the documentation of private function "LDB" 
%
% Output parameters : 
%
% 1) T  : Column vector of time points T = [T(1);...;T(k+1)], where [T(i),T(i+1)], i = 1,...,k  
%         are the integration intervals with T(1) = t0 and T(k+1) = tf.  
% 2) Y  : Solution array. Each row Y(i,:) consists of the n Taylor models Y(i,1),....,Y(i,n) that 
%         enclose the solutions y_1,...,y_n of the ODE on the time interval [T(i),T(i+1)], i = 1,...,k.
%         If preconditioning is chosen (see 4)), then the yl_j := Y(i,j), j=1,...,n, 
%         are the left (time-dependent) Taylor models.   
% 3) Yr : The corresponding (time-independent) right Taylor models are yr_j := Yr(i,j), j=1,...,n. 
%         The ODE-flow on the time interval [T(i),T(i+1)] is given by the concatenation of left and right Taylor models: 
%         yl_j(yr_1(x),...,yr_n(x),t), j=1,...,n, where x = (x_1,...,x_n) \in [-1,1]^n describes the space variables 
%         and t in [T(i),T(i+1)] is the time variable.          
%        
%
% The implementation is based on: 
%  
%  [E]  I. Eble, "Über Taylor-Modelle", Dissertation at Karlsruhe Institute of Technology, 2007 (written in German),
%         Riot, C++-implementation, http://www.math.kit.edu/ianm1/~ingo.eble/de
%  [M]  K. Makino, "Rigorous analysis of nonlinear motion in particle accelerators", 
%         Dissertation at Michigan State University, 1998
%  [MB] K. Makino and M. Berz, "Suppression of the wrapping effect by Taylor model - based validated integrators",
%         MSU HEP Report 40910, 2003 
%  [NJN] M. Neher, K.R. Jackson, N.S. Nedialkov, "On Taylor model based integration of ODEs", 
%          SIAM J. Numer. Anal. 45(1), pp. 236-262, 2007
%  [Bue] F. Buenger, "Shrink wrapping for Taylor models revisited", Numerical Algorithms 78(4), pp. 1001-1017, 2018

% written  11/04/15     F. Buenger
% modified 01/18/16     F. Buenger  record feature                                     
% modified 02/04/16     F. Buenger  code-move to "verifyode_tm" in @taylormodel  
% modified 05/15/17     F. Buenger  preconditioning of Taylor models   
% modified 06/29/17     F. Buenger  blunting   

global INTLAB_ODE_VARS

% taylormodelinit; % <-- This line moved to "startintlab.m". 
if nargin < 4
    options  = [];
end

[T,Y,Yr] = verifyodeTM(INTLAB_ODE_VARS.ZEROTAYLORMODEL,odefun,tspan,y0,options);

end % function verifyode
