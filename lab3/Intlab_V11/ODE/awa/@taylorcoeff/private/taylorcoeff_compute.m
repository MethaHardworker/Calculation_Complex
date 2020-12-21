function r = taylorcoeff_compute(odefun,y0,hdk,order,max_status,treenr,y,t)

%TAYLORCOEFF_COMPUTE    computes generalized Taylor coefficients. 
%
%   r = taylorcoeff_compute(fun,u0,hdk,order,max_status,treenr,treemode)
%
% Generalized Taylor coefficients c of the solution y(t) of the ODE
%                                                      
%   y' = odefun(y), y(t0) = y0
% 
% at (not explicitly given) start time t0 are computed up to order "max_status". 
% Precisely, the generalized Taylor coefficients c_i(k), k = 1,...,max_status+1 
% of the i-th solution component y_i(t), i= 1,...,n, are defined as follows: 
%
%   a) c_i(1) = y0_i,  
%   b) c_i(k) = (hdk(1)*...* hdk(k-1)) * y_i^(k-1)(t0), k = 2,..., max_status+1. 
%
% Here y_i^(j)(t0) denotes the j-th derivative of y_i at t0.  
% The typical setting is hdk(k) = h/k ("hdk" shortens "h divided by k"),
% so that (hdk(1)*...* hdk(k-1)) = h^(k-1)/(k-1)! and 
%
%       c_i(k) = h^(k-1)/(k-1)! y_i^(k-1)(t0)
%
% becomes the k-th addend of the Taylor expansion 
%
%       y_i(t0+h) = sum_{k=0}^infty h^k/k! y_i^(k)(t0).  
% 
% The generalized, scaled Taylor coefficients c_i are stored 
% in the component r(i).coeff of the i-th array entry r(i) wich is 
% of type "taylorcoeff". 
%
% The function mainly corresponds to what is done by the [AWA] functions 
% "F" and "TKOEFF", see file awa.p.

% written  08/07/17     F. Buenger

global INTLAB_AWA_VARS

if isfloat(y0)
    % Convert y0 to point interval (vector or matrix).
    x.inf = y0;
    x.sup = y0;
    y0 = x;
end

% h :=  [hdk.inf(1),hdk.sup(1)] is an inclusion for the actual 
% time step size. The function id(t) = t has the time derivatives
%
%       t0, 1, 0, 0, ....     at t = t0. 
%
% Thus its scaled generalized Taylor coefficients are  
%
%       t0, h, 0, 0, .... .
%
% The input parameter t, which is of type taylorcoeff, is expected to contain
% t0 already as its first entry in the component t.coeff. 
% Moreover it is expected that t(k) = 0 for k >= 3. Thus, only the 
% second generalized Taylor coefficient, which is h, has to be adapted.

t.inf(2) = hdk.inf(1); % t(2) := h
t.sup(2) = hdk.sup(1);

% Create Taylor coefficients y with initial interval vector y0 as coefficients of order zero.
n = numel(y0.inf);  % dimension of the ode
if isempty(y)  
    y.inf = zeros(order+1,1);
    y.sup = y.inf;
    %y = repmat(y,n,1);
    y = repmat(y,size(y0.inf));
    for i = 1:n
        y(i).inf(1) = y0.inf(i);  
        y(i).sup(1) = y0.sup(i);
    end  
    y = taylorcoeffinit(y); % taylorcoeff constructor
else % Reuse given Taylor coefficients y (just for performance reasons).  
    z = zeros(order,1);
    for i = 1:n
        y(i).inf = [y0.inf(i);z];  
        y(i).sup = [y0.sup(i);z];
    end    
end

% The ode function odefun is internally evaluated by MATLAB as a "tree" 
% consisting of basic arithmetic functions like +,-,*,/,^ and standard
% functions like sin, cos, exp, log, atan, ... .
% The vertices of that evaluation tree correspond to intermediate results. 
% The Taylor coefficients of them are globally stored, so that they can be 
% reused in the iterative computatation of generalized Taylor coefficients 
% of orders 1,2,..., max_status.

INTLAB_AWA_VARS.TREENR = treenr; % The main ode function y'= f(t,y) has tree number = 1.
                                 % The internally created ode function for the matrix variation equation
                                 % Y' = (df/du)*Y has the tree number = 2.

% Iterative computation of the generalized Taylor coefficients up to order max_status
for status = 0:max_status-1 
    INTLAB_AWA_VARS.VERTEXNR = 0;
    INTLAB_AWA_VARS.STATUS = status;        
    dy = odefun(t,y); % This is the main point where the generalized Taylor coefficients 
                      % are computed by simply evaluating the ode function odefun 
                      % with input t,y which have type taylorcoeff. 
                      % The operator concept causes that all computations are done
                      % in taylorcoeff arithmetic which is implemented in the folder 
                      % @taylorcoeff. This arithmetic implements the recursion equations
                      % for generalized Taylor coefficients like in [AWA] function TKOEFF, file awa.p.
                      % (See also [AWA] file fknoten.p, where "accompanying" functions for certain 
                      % standard functions like tan, cot,... are stated.)
    hdk_.inf = hdk.inf(status+1);                 
    hdk_.sup = hdk.sup(status+1);                 
    y = shift(y,dy,status,hdk_); % This simply "shifts" the generalized Taylor coefficient of order status of dy 
                                 % to that of order status+1 of y and multiplies it by hdk(status+1). 
                                 % This simply resembles the relation y' = odefun(t,y).
    if treenr < 0
        treenr = -treenr;
        INTLAB_AWA_VARS.TREENR = treenr;
    end                            
end

r = y;

end % taylorcoeff_compute


