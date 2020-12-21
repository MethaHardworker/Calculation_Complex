function [T,Y,Z] = awaTC(dummy,odefun,jacobimat,tspan,y0,options)
%AWATC  implements AWA in MATLAB/INTLAB.  
% 
%   r = awa(dummy,odefun,jacobimat,tspan,y0,options)
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
% 1) dummy: This is a variable of type taylorcoeff which is only used for  
%           calling awaTC in the folder @taylorcoeff. The function 
%           awaTC (shortening of "awa Taylor coefficients") is placed in  
%           that folder for performance reasons, since direct, fast access 
%           to components of variables of type taylorcoeff is only possible 
%           for functions inside that folder. Apart from that "dummy" has 
%           no further meaning and is not used in the function at all.
%              
% 2) odefun: This is, as for all MATLAB ODE solvers like ode45, a function
%           handle for the right-hand side f of the ODE y' = f(t,y).
%           The implementation of odefun is due to the user and
%           the output must be a column vector! The implementation
%           can contain INTLAB interval arithmetic. Thus, if, for example
%           f(t,y) = cos(pi*y)+t, then odefun must be implemented 
%           as odefun = @(t,y) cos(intval('pi')*y) + t, since pi is not
%           a floating-point number.
% 
% 3) jacobimat: This is a function handle to the Jacobi matrix J of f
%           with respect to y, i.e., J_{i,j}(t,y) = df_i/dy_j(t,y), 
%           i,j = 1,...,n, where n is the dimension of the ODE.
%           The implementation of the function jacobimat is - like the
%           implementation of odefun - due to the user.
%
% 4) tspan = [t_start,t_end] specifies the integration interval like for 
%           MATLAB ODE solvers. It is allowed that t_end < t_start.
%           Note that t_start and t_end must be floating-point numbers! 
%           In general this is no severe restriction and is either
%           fulfilled anyway or can easily be achieved by a shift and/or
%           scaling of the time variable in the ODE function. 
%           (This may lead to interval parameters in the ODE function.)
%
% 5) y0:    This is the initial value vector, i.e. y(t_start) = y0, which 
%           can also be an interval vector, i.e., y0 = [a,b],
%           a_i <= y_i(t_start) <= b_i, i = 1,...,n.  
%           
% 6) options: This structure has components that are parameters of awa. 
%           The function awaset can be used to create "options" similar 
%           to the MATLAB function "odeset". 
%           The following options can be set:
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
% 2) Y: Enclosure of solution at time grid points, y_i(T(j)) is contained in 
%       Y(j,i) = [Y(j,i).inf,Y(j,i).sup], j = 1,...,m, i = 1,...,n.
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
%       Randwertaufgaben und Anwendungen, Diss. Univ. Karlsruhe, 1988
%
% [AWA] R. Lohner, Pascal implementation according to [L], available from 
%       http://www2.math.uni-wuppertal.de/~xsc/xsc/pxsc_software.html#awa  

% written  07/18/17     F. Buenger
% modified 04/09/18     F. Buenger, additional output parameter Z for continuous solution inclusion

global INTLAB_AWA_OPTIONS; % global structure that contains parameters specified in "options" and default values otherwise
global INTLAB_AWA_VARS;    % structure for internal global variables that cannot be modified by the user

e = 1e-30;
if 1+e > 1                 % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1)            % switch to rounding upwards
end

INTLAB_AWA_VARS.eps_abort = 1e20;     % Abort integration if low accuracy.    
INTLAB_AWA_VARS.mach_eps  = 1e-16;    % approx. machine accuracy (IEEE)      
INTLAB_AWA_VARS.TREE = {[];[];[];[]}; % Initialize all function trees in global memory.
treenr_1 = 1;                         % ODE function gets tree number 1.
treenr_2 = 2;                         % Jacobi function gets tree number 2.

awaopt2glob(options); % Write all options specified by the user to global structure INTLAB_AWA_OPTIONS.
                      % Parameters that were not specified by the user get default values.   

%--------------------------------------------------------------------------                      
% Invent shortenings for the components of INTLAB_AWA_OPTIONS.
%--------------------------------------------------------------------------                      
h0 = INTLAB_AWA_OPTIONS.h0;           % initial step size
h_min = INTLAB_AWA_OPTIONS.h_min;     % minimum step size  
p = INTLAB_AWA_OPTIONS.order;         % order of Taylor expansion
e_meth = INTLAB_AWA_OPTIONS.EvalMeth; % Corresponds to [AWA] variable E_ART. (The German word "Art" means "type" or "method".)
e_a = INTLAB_AWA_OPTIONS.AbsTol;      % absolute error tolerance
e_r = INTLAB_AWA_OPTIONS.RelTol;      % relative error tolerance 
                                                                                       
%--------------------------------------------------------------------------                      
% Check y0. 
%
%   It is assumed that y0 is an interval-like structure with components 
%   y0.inf and y0.sup. The initial value interval y0 corresponds to the 
%   [AWA] variable "ANFANG" (German for "BEGINNING"). 
%--------------------------------------------------------------------------                      
s = size(y0.inf); 
if (length(s) ~= 2) || (s(1) ~= 1 && s(2) ~= 1) 
    error('invalid call of verifyode: y0 must be a vector.')
end
if s(2) > s(1)              % The row interval vector y0 is transformed into a column vector.
  y0 = iv_transpose(y0);
end

yapp0 = iv_getmidrad(y0);   % midpoint of y0, corresponds to [AWA] variable "UAPP0"
ydev0 = iv_minus(y0,yapp0); % deviation from midpoint corresponds to [AWA] variable "ANF" (shortening of "ANFANG", German for "BEGINNING")
n = length(y0.inf);         % dimension of the ODE-System y' = f(t,y), f: [t0,tf] x R^n -> R^n

%--------------------------------------------------------------------------                      
% Check tspan.
%--------------------------------------------------------------------------                      
s = size(tspan); 
if (length(s) ~= 2) || (s(1) ~= 1 && s(2) ~= 1) || ~isa(tspan,'double') || ~isreal(tspan) ... 
   || any(isinf(tspan)) || any(isnan(tspan))
    error('invalid call of verifyode: improper parameter tspan.')
else
    tstart = tspan(1); % starting time for integration    
    tend = tspan(2);   % end time for integration
    if tstart == tend
        error('Start- and endpoint of integration should be different. Otherwise nothing has to be done!')
    end    
end

t = tstart; % t corresponds to [AWA] variable tt, see file awa.p, procedure MAIN, and file awaraout.p, procedure EINLESEN1.
% Convert t to type "Taylor coefficient"; variable t_tcoeff. 
t_.inf = [t;1;zeros(p-1,1)];
t_.sup = t_.inf;
t_tcoeff = taylorcoeffinit(t_); % Taylor coefficient version of time t.

y0_abs = iv_abs(y0);
y0_diam = iv_diam(y0);
x = y0_diam./y0_abs.sup; % Note that y0_abs.sup(i) = 0 implies y0_diam(i) = 0 so that x = NaN in this case. 
heps = max(max(x,INTLAB_AWA_VARS.mach_eps,'omitnan'));
eps_ = heps*INTLAB_AWA_VARS.eps_abort; % corresponds to [AWA] variable "eps". The notation is slightly changed because eps = 2^-52 is a MATLAB constant. 
                                       % This heuristic constant is only used as an abort criterion if the accuracy becomes too bad during integration.                                

l1 = intval2iv( odefun( intval(t) , iv2intval(y0) ) ); % This is the initial slope interval containing all possible derivatives y'(t) where t = tstart.
                                                       % Corresponds to [AWA] variable L1, see procedure MAIN, file awa.p. 
if size(l1.inf,2) > 1
    error('odefun must return a column vector.');
end

len = zeros(1,n); % Storage preallocation of a vector "len" that will contain the weighted Euclidean column lengths of a matrix M for which a  
                  % non-verified QR decomposition shall be computed. Prior to the QR decomposition the columns of M are sorted in descending order
                  % with respect to these lengths.

%--------------------------------------------------------------------------                      
% Initialize some further values, see [AWA], procedure MAIN, file awa.p, 
% call of procedure "EINLESEN2" (German for "READ2"), file awarout.p 
%--------------------------------------------------------------------------                      
rest0.inf = zeros(n,1);  % rest0 is initialized as the point interval zero. 
rest0.sup = rest0.inf;   % column vector of length n
rest1 = rest0;    
rest = rest0;
abbmat = eye(n);         % nxn identity matrix 
restmat = abbmat;

%--------------------------------------------------------------------------                      
% All preparations up to this point mainly correspond to what is done in 
% [AWA] before entering the REPEAT-loop in procedure MAIN, file awa.p,
% where the integration starts.
%--------------------------------------------------------------------------                      

%--------------------------------------------------------------------------                      
% Prepare output parameters 
%--------------------------------------------------------------------------                      
len_Y = ceil(abs(tend-tstart)/h_min); % len_Y is a pessimistic estimate for the number of time steps. 
                          % It is only used for memory preallocation of the output variables Y and T.
                          % This is the only point where h_min is used. If the step size becomes smaller than h_min,
                          % then the only consequence is that the preallocation was not optimal so that the subsequent 
                          % result arrays T and Y may grow during integration which may cause bad performance.  
T = zeros(len_Y+1,1);     % preallocation of time step result array T
Y.inf = zeros(len_Y,n);   % preallocation of interval array. The interval [Y.inf(i),Y.sup(i)] will 
Y.sup = Y.inf;            % include the sought solution y of the given ODE at time grid point T(i),
                          % i = 1,..., n_t := number of of time grid points.
                          % In most cases n_t << len_Y will finally hold true. Then, at the end of the function,
                          % when n_t is known, the lengths of T and Y will be reduced to n_t. 
T(1) = tstart;            % Write start time to result array T.    
Y.inf(1,1:n) = y0.inf';   % Write initial values to result interval array Y.
Y.sup(1,1:n) = y0.sup';

% Define Y0 := Y(t=0) = I (with I := (nxn) identity matrix) to be the
% initial value matrix of the matrix variation equation Y' = df/du * Y.
Y0.inf = eye(n); 
Y0.sup = Y0.inf;
Y_tcoeff = []; % initialize the Taylor coefficient version of the sought solution Y as empty.

C.inf = zeros(n,p+1); % The interval matrix C is just a help-variable that temporarily stores generalized 
C.sup = C.inf;        % Taylor coefficients of the solution y = (y_1,...,y_n) up to order p.

%--------------------------------------------------------------------------                      
% The following code mainly corresponds to what is done in [AWA], 
% procedure "ZEITSCHRITT" (German for "TIME STEP"), file awa.p, 
% where the stepwise integration is done.  
%--------------------------------------------------------------------------                      

eps_rel = 0.1;       % relative epsilon inflation, corresponds to [AWA] factor VFAKT = 0.1, see file awa.p, procedure ZEITSCHRITT. 
eps_abs = e_a + e_r; % absolute epsilon inflation, compare with [AWA] epsilon inflation in lines 501, 502 of file awa.p. 
time_steps = 0;      % counter for executed time steps
direction = sign(tend-tstart); % direction of integration. If tend > tstart, then direction = 1, and if tend < tstart, then direction = -1. 

P = (1:p+1).';
h = 1;

hdk.inf = h./P;  % Initialize factors for generalized Taylor coefficients in non-verified arithmetic.
                 % Remark: Here [AWA] performs verified computation: HDK[I]:= INTVAL(1)/(I+1), line 468, file awa.p 
hdk.sup = hdk.inf;

% Compute approximate, non-verified generalized Taylor coefficients of y up to order p.
% The expansion point is yapp0 = mid(y0) and the scaling factors for the taylor coefficients are given in the array hdk. 
% The result is only used for automatic, heuristic step size control so that non-verified computation is sufficient.
% The function "taylorcoeff_compute" corresponds to the [AWA] function "TKOEFF".
% The negative input parameter "-treenr_1" signalizes "first call" during which some meta data are stored in order to increase performance for subsequent calls.
y_tcoeff = taylorcoeff_compute(odefun,yapp0,hdk,p,p,-treenr_1,[],t_tcoeff); 

z_local.inf = zeros(n,1);
z_local.sup = z_local.inf;
for i = 1:n
    z_local.inf(i) = y_tcoeff(i).inf(p+1); % Only the last (p-th), non-verified generalized Taylor coefficient is used for somehow estimating 
    z_local.sup(i) = y_tcoeff(i).sup(p+1); % the local error for the first integration step which will be used for a heuristic step size control.
end
Z.inf = zeros(len_Y,n,p+1); % just preallocation of output array Z
Z.sup = Z.inf;
abort = false;              % The abort flag becomes true if the accuracy becomes too bad during integration.
h_in = direction*abs(h0);   % initial step size

%--------------------------------------------------------------------------
% Define some special constant intervals 
%--------------------------------------------------------------------------
iv_01.inf = 0;   % iv_01 = [0,1].                                          
iv_01.sup = 1;                                                  

iv_m11.inf = -1; % iv_m11 = [-1,1]                                          
iv_m11.sup = 1;                                                  

iv_0q.inf = 0;   % iv_0q = [0,0.25]
iv_0q.sup = 1/4; % floating-point number, no rounding error

%**************************************************************************
%**************************************************************************
%
% main time integration loop
%
% This corresponds to the REPEAT-loop in [AWA] procedure "ZEITSCHRITT" 
% (German word for "TIME STEP"), file awa.p.
%
%**************************************************************************
%**************************************************************************
while ~( t >= tend && direction > 0 || t <= tend && direction < 0 || abort )

%--------------------------------------------------------------------------
% automatic step size control    
%     
%   The following heuristic choice of the (approximate) time step size h in 
%   relation to an estimate of the local error is somehow standard. 
%   Choosing the step size is always a bit like reading tea leaves. 
%   Feel free to implement a different (better) choice of the step size.
%--------------------------------------------------------------------------

    norm_z_local = iv_vecnorm(z_local,inf); % infinity norm of the local error of previous time step
    
    if h_in == 0 && norm_z_local == 0
        h_in = tend - tstart;
    end    
    if time_steps > 0        
        err = e_a + e_r * iv_vecnorm(ytau_local,inf); % y_tau_local is an inclusion for y on the whole integration interval of the previous time step.
        norm_diam_z_local =  norm(iv_diam(z_local),inf);         
        if norm_diam_z_local > 0 
            % In [AWA], "norm_diam_z_local > 0" is not checked! 
            % This causes that, for example, the one-dimensional ODE
            %      y'(t) = sqr(t), y0 = 0, t in [0,1],
            % solved with order 12, crashes since norm_diam_z_local = 0 
            % holds true after the first time step so that division by zero 
            % occurs in the next statement below this comment.
            % (This corresponds to line 482 in 'awa.p' at which [AWA] crashes.)
            %
            % In Matlab, division by zero would result in h_in = Inf 
            % in this case, i.e., infinite large step size.
            % The if-statement causes that the step size stays the same as 
            % in the previous time step if norm_diam_z_local = 0.
            
            h_in = ( abs(h)^(p+1) * err / norm_diam_z_local )^(1/p); 
        end 
    else
        err = e_a + e_r * iv_vecnorm(y0,inf);
        h = norm_z_local;
        if h_in == 0 && h > 0
            h_in = (err/(p*h))^(1/(p-1));
        end
    end;

    h = direction*abs(h_in); % new step size             
    t_next = t+h;            % next time grid point, corresponds to [AWA] variable tt_local, see file awa.p, procedure ZEITSCHRITT
                             
    % h is changed such that |h| <= |t_next-t| = |h_exact|
    if direction > 0
        h = -(t-t_next);   % Recall that rounding is upwards.
        tau.inf = t;       % integration interval: tau = [t,t_next],
        tau.sup = t_next;  % corresponds to [AWA] variable tau = t + [0,1]*h, file awa.p, line 495, 
                           % but it is not exactly the same, see the REMARK below.
    else % direction < 0
        h = t_next-t;      % Recall that rounding is upwards.        
        tau.inf = t_next;  % integration interval: tau = [t_next,t],
        tau.sup = t;       % corresponds to [AWA] variable tau = tt + [0,1]*h, file awa.p, line 495, 
                           % but it is not exactly the same, see the REMARK below.
    end
    diff_tau = iv_minus(tau,t); % diff_tau := tau-t is used in the Picard iteration.
                                % It corresponds to (but is not equal to) 
                                % [AWA] [0,1]*h,see the REMARK below.
                                  
% REMARK: Note again carefully that, possibly due to rounding errors
%         in the summation t_next := t+h, h is only an approximation for 
%         the implicitly given exact step size h_exact = t_next-t which is  
%         the signed diameter of the actual integration interval tau.
%         Therefore, the approximate step size h will only be used for  
%         non-verified computations and never(!) for verified computations.
%
%         We prefer to work with explicitly given floating-point time grid
%         points t, t_next,... at which the function value inclusions are
%         computed and which are collected and returned in the output
%         parameter T and Y, respectively. Step sizes are not returned,
%         so they can remain implicit.
%
%         If in AWA, due to rounding errors, tt_local := tt + h is not equal 
%         to tt + h, then a warning is displayed, see file awa.p, line 535:
%         'step size h and grid point t contaminated by roundoff errors !'
%
%         We avoid that by our slightly different setting.
                                                
% Compute a non-verified approximate solution yapp for y(t_next).    
    
    hdk.inf = h./P; % Set new, non-verified factors for generalized Taylor coefficients. 
                    % Remark: Here [AWA] performs verified computation: HDK[I]:= INTVAL(1)/(I+1), line 493, file awa.p
    hdk.sup = hdk.inf;
    
    y_tcoeff = taylorcoeff_compute(odefun,yapp0,hdk,p,p-1,treenr_1,y_tcoeff,t_tcoeff); % Create generalized Taylor coefficients up to order p-1 with new factors hdk.
    yapp = taylorcoeff_sumup(y_tcoeff,p-1); % This corresponds to [AWA] function call "uapp = naeherung(uapp0,p-1)" in procedure "ZEITSCHRITT", file awa.p, line 496.
                                            % yapp is a non-verified approximation for y(t+h).                                      
                                            
%--------------------------------------------------------------------------
% Compute a rough, verified inclusion [y] of y(t) on the whole interval 
% tau = [t,t_next] (tau = [t_next,t] if direction = -1) as follows:
%  
% Choose some (appropriate) interval [y^0] that includes y0 and iterate
%
%   a) epsilon inflation: 
%
%         [y^0] := [1-eps_rel,1+eps_rel]*[y^0] + [-eps_abs,eps_abs]   
%
%   b) simple form of Picard iteration:
%        
%         [y^1] := y0 + diff_tau * odefun(tau,[y^0])                 
%        
%      (Recall that diff_tau includes [0,t_next-t] respectively 
%       [t_next-t,0] if direction = -1.) 
%
%  until 
%
%    c)  [y^1] is contained in [y^0].     
%
% Step c) is an inclusion property for applying Banach's fixed point  
% theorem to the integral operator
%
%       K(y)(s) := y0 + \int_t^{s} odefun(x,y(x)) dx, s in [t,t_next], 
%
% see [L], Satz (1.2.2.8), p. 29.
%--------------------------------------------------------------------------
                                      
    ytau = iv_hull(y0,yapp); % ytau is the interval hull of y0 and yapp. Clearly, [y^0] := ytau includes y0. 
                        
    %--- first iteration without epsilon inflation
    ftau = intval2iv(odefun(iv2intval(tau),iv2intval(ytau))); % ftau := odefun(tau,ytau)       evaluated in INTLAB interval arithmetic. INTLAB's object orientation is needed here for evaluating the ODE function.
    ytau = iv_plus( y0 , iv_times(diff_tau,ftau) );           % ytau := y0 + diff_tau * ftau   evaluated in fast simplified, non-object-oriented interval arithmetic. 
    
    %--- second iteration with epsilon inflation
    ytau = blow(ytau,eps_rel,eps_abs);                        % Compare with [AWA] epsilon inflation in lines 501, 502 of file awa.p.  
    ftau = intval2iv(odefun(iv2intval(tau),iv2intval(ytau))); % ftau := odefun(tau,ytau)       evaluated in INTLAB interval arithmetic.
                                                              
    % Using ytau and ftau, a reduced step size h_ will now be determined 
    % such that the inclusion property c) stated above still holds true 
    % for h_*[0,1] instead of diff_tau. This is always possible since y0 
    % is contained in the interior of [y^0],  see the explanation in [L] 
    % at the bottom of page 41, and [L], Satz (1.2.3.8), p. 29.
    % 
    % For "direction" > 0, the reduced step size h_ is aimed at fulfilling
    % the following inclusion property componentwise:
    %
    %   d) [y^1] := y0 + [0,h_] * ftau  \subseteq ytau =: [y^0].
    % 
    % Roughly speaking, d) is now "resolved" componentwise for h_.
    %  
    % The following comments refer to "direction" > 0; 
    % the case "direction" < 0 is analogous.     
        
    if direction > 0
        idx1 = (ftau.inf < 0);
        idx2 = (ftau.sup > 0);
        
        dummy = ytau.inf(idx1) - y0.inf(idx1);      % Recall that rounding is upwards and that ytau.inf <= y0.inf componentwise since ytau contains y0 by construction.
                                                    % Thus dummy <= 0 and |dummy| <= y0.inf(idx1)-ytau.inf(idx1).
        dummy = -( -dummy ./ ftau.inf(idx1) );      % By choice of idx1 we have ftau.inf(idx1) < 0 so that dummy >= 0. Rounding upwards yields:
                                                    % 
                                                    %      dummy                                  <= ( ytau.inf(idx1) - y0.inf(idx1) ) ./ ftau.inf(idx1) 
                                                    % <=>  dummy.*ftau.inf(idx1)                  >= ytau.inf(idx1)-y0.inf(idx1)   [Recall that ftau.inf(idx1) < 0]
                                                    % <=>  y0.inf(idx1) + dummy.*ftau.inf(idx1)   >= ytau.inf(idx1)                                                 
        h1 = min(dummy,[],'includenan');            %  =>  y0.inf(idx1) + inf(h1 * ftau(idx1))    >= ytau.inf(idx1)     (*)
        
        dummy = -( y0.sup(idx2) - ytau.sup(idx2) ); % Rounding upwards and ytau.sup >= y0.sup since (as y0 is contained in ytau by construction) yield:
                                                    % 0 <= dummy <= ytau.sup(idx2) - y0.sup(idx2).
        dummy = -( -dummy ./ ftau.sup(idx2) );      % Rounding upwards and ftau.sup(idx2) > 0 by choice of idx2 yield:
                                                    % 0 <= dummy <= (ytau.sup(idx2) - y0.sup(idx2)) ./ ftau.sup(idx2).
                                                    % Thus y0.sup(idx2) + dummy.* ftau.sup(idx2) <= ytau.sup(idx2) which implies 
        h2 = min(dummy,[],'includenan');            %      y0.sup(idx2) + sup(h2*ftau(idx2))     <= ytau.sup(idx2)      (**)  
                       
        h_ = min([h h1 h2],[],'includenan');        % Recall that 0 <= h <= h_exact by construction.
                                                    % From (*) and (**) it follows that the desired inclusion property  
                                                    %   d) y0 + [0,h_] * ftau  \subseteq ytau    
                                                    % holds true. 
                               
       t_next = -(-t-h_);                           % Rounding upwards and h_ > 0 imply h_exact := t_next-t <= t+h_ in exact arithmetic.                                                    
       t_next = min(t_next,tend);                   % t_next <= t_end must always hold true. 
       h = -(t-t_next);                             % Since rounding is upwards the new approximate step size h fulfills again |h| <= |h_exact| <= |h_|.
       tau.sup = t_next;                            % Update tau.       
       
    else % The case "direction" < 0 proceeds analogously.
        
        idx1 = (ftau.inf < 0);
        idx2 = (ftau.sup > 0);
        
        dummy = -(y0.sup(idx1) - ytau.sup(idx1));   % Recall that rounding is upwards.
        dummy = dummy ./ ftau.inf(idx1);
        h1 = max(dummy,[],'includenan');
        
        dummy = ytau.inf(idx2) - y0.inf(idx2);      % Recall that rounding is upwards.
        dummy = dummy ./ ftau.sup(idx2);
        h2 = max(dummy,[],'includenan');
        
        h_ = max([h h1 h2],[],'includenan');        % Recall that |h| <= |h_exact|.  
                                                    % Inclusion property  d) y0 + [0,h_] * ftau \subseteq ytau   holds true.
                                             
        t_next = t+h_;                              % Rounding upwards and h_ < 0 imply 0 >= h_exact := t_next-t >= h_  in exact arithmetic.                                                    
        t_next = max(t_next,tend);                  % t_next >= t_end must always hold true. 
        h = t_next-t;                               % Since rounding is upwards, the new approximate step size h fulfills again |h| <= |h_exact| <= |h_|.
        tau.inf = t_next;                           % Update tau.
    end  
        
    diff_tau = iv_minus(tau,t);                     % Update diff_tau.                 
        
    % --- final, third Picard iteration
    % Since the inclusion property d) holds true, further iterations improve
    % the inclusion interval. Only one such improvement step is carried out.
                     
    ftau = odefun(iv2intval(tau),iv2intval(ytau));
    ftau = intval2iv(ftau);
    ytau = iv_plus( y0 , iv_times(diff_tau,ftau) ); % ytau := y0 + diff_tau * ftau  is the sought rough, verified inclusion of the solution y(t) on the whole interval [t,t_next]. 
    
%--------------------------------------------------------------------------
% Now the Frechet derivative A_j := (I+h*Phi'(y_j)) of the ODE method 
%
%   y_{j+1} := y_j + h * Phi(y_j)   with    y_k := y(t_k) 
%
% is computed where I is the nxn-identity matrix. This corresponds to the 
% result of [AWA] function TFAKT, file awa.p. The computation is done by 
% solving a second ODE which Lohner calls "Variationsgleichung", 
% see [L], (2.2.10), p. 42: 
%
% (ODE_2)   Y' = (df/dy) * Y,   Y0 := Y(t_j) := I.    
%
% Here df/dy := (df_i/dy_j)_{1 <= i,j <= n} denotes the Jacobi matrix 
% w.r.t. y_1,...,y_n of the ODE function f of the original ODE  
%
% (ODE_1)   y' = f(t,y). 
%
% Note carefully that df/dy =: g(t,y(t)) is a matrix function in t 
% given by the input function handle jacobimat.
% Thus, if the Taylor coefficients of y up to order p-1 are known at 
% t = t_j, then, by using the Taylor coefficient arithmetic implemented for 
% the data type "taylorcoeff", Taylor coefficients of df/dy at t = t_j are 
% derived up to the same order by a single function call of jacobimat.
%
% Afterwards, Taylor coefficient arithmetic can be applied to (ODE_2) in 
% order to compute Taylor coefficients of the solution Y at t_j. 
% Their sum up to order p-1 is the sought Frechet derivative A_j. 
% See [L], p. 42f for detailed explanation and proofs.
%--------------------------------------------------------------------------
    
    h_ = iv_minus(t_next,t); % interval inclusion of the the exact step size h_exact = t_next - t
    hdk = iv_rdivide(h_,P); % verified computation of generalized Taylor coefficient factors hdk := h_./P 

    % Compute a verified inclusion of generalized Taylor coefficients of the solution y(t) with start interval vector y0 
    % and new Taylor coefficient factors hdk up to order p-1.    
    y_tcoeff = taylorcoeff_compute(odefun,y0,hdk,p,p-1,treenr_1,y_tcoeff,t_tcoeff);
    Z_ = y_tcoeff; % Z_ contains the first p-1 Taylor coefficients of the output array entry in Z for the current time step.
                   % The final p-th Taylor coefficient, which encloses the local error on the whole time interval, is appended later. 
    
    INTLAB_AWA_VARS.STATUS = []; % Compute Taylor coefficients of all orders.
    % Define Taylor coefficient version of time t for actual step size 
    t_ = t_tcoeff; 
    t_.inf(2) = hdk.inf(1); % generalized Taylor coefficient of order 1:  t_(2) := h = 1 * h,  since d/dt t = 1 
    t_.sup(2) = hdk.sup(1);
    
    df_tcoeff = jacobimat(t_,y_tcoeff); % Compute verified generalized Taylor coefficients df_tcoeff of the Jacobi matrix df/dy 
                                        % up to order p-1 by a single function call of jacobimat.    
    fun = @(t,Y)(df_tcoeff * Y);        % This is the ODE function of (ODE_2) defined as a function handle.    
    
    Y_tcoeff = taylorcoeff_compute(fun,Y0,hdk,p,p-1,-treenr_2,Y_tcoeff,t_tcoeff); % Verified computation of generalized Taylor coefficients of the solution Y of (ODE2) at t = t_j up to order p-1.
                                                                                  % The negative input parameter "-treenr_2" signalizes "first call" during which some meta data are stored in order 
                                                                                  % to increase performance for subsequent calls.
    iphi = taylorcoeff_sumup(Y_tcoeff,p-1); % Sum up the generalized Taylor coefficients Y_tcoeff, i.e., build the Taylor series, up to order p-1.  
                                            % As described above this sum is the sought Frechet-derivative A_j = (I+h*Phi'(y_j)) of the ODE method.
                                            % This corresponds to the [AWA] statement IPHI := TFAKT; in line 549 of procedure "ZEITSCHRITT", file awa.p. 
                                           
%--------------------------------------------------------------------------                                            
%  Estimation of global error, part I, see [L], section 2.3, p. 46 et seq.                                           
%--------------------------------------------------------------------------                                            
                                            
    ima1 = iv_mtimes(iphi,abbmat); % ima1 = iphi*abbmat. The interval matrix ima1 corresponds to [A_j]*C_j in [L], (2.3.10), p. 52. 
    abbmat = iv_getmidrad(ima1);   % matrix that approximately describes the mapping of the initial value vector.    
                                   % The point-matrix abbmat corresponds to C_{j+1} in [L], (2.3.10), p. 52.
    % Compute new approximation yapp of y(t_next) for new t_next. 
    % This corresponds to [AWA] uapp:= naeherung(uapp0,p-1), file awa.p, line 545.
    y_tcoeff = taylorcoeff_compute(odefun,yapp0,hdk,p,p-1,treenr_1,y_tcoeff,t_tcoeff); 
    yapp = taylorcoeff_sumup(y_tcoeff,p-1); 
        
    % Now 5 evaluation methods for the global error are distinguished.
    switch e_meth 
        case 0 % interval vector enclosure
            rest0 = iv_mtimes(iphi,rest0); % rest0 = iphi * rest0;
            yhat = iv_plus( iv_plus( iv_mtimes(ima1,ydev0) , rest0 ) , yapp ); % yhat = ima1*ydev0 + rest0 + yapp, 
                                                                               % yhat corresponds to [AWA] variable "UDACH". (Read "U"-"DACH" which is German for "U"-"HAT")
        case {1,2} % parallelepiped or QR decomposition
            ima0  = iv_mtimes(iphi,restmat); % ima0  = iphi*restmat;
            yhat = iv_plus( iv_plus( iv_mtimes(ima1,ydev0) , iv_mtimes(ima0,rest1) ) , yapp ); % yhat = ima1*ydev0 + ima0*rest1 + yapp            
        case {3,4} % intersection of 0 and 1, or 0 and 2, respectively.
            rest0 = iv_mtimes(iphi,rest0); % rest0 = iphi*rest0;
            ima0  = iv_mtimes(iphi,restmat); % ima0  = iphi*restmat;
            yhat = iv_plus( iv_plus( iv_mtimes(ima1,ydev0) , iv_intersect(rest0,iv_mtimes(ima0,rest1)) ) , yapp ); % yhat = ima1*ydev0 + iv_intersect(rest0,ima0*rest1) + yapp;
    end;
    
%-------------------------------------------------------------------------
% Calculate enclosures z for the local error, and y1 for the solution
% at t = t_next. 
%
% In [AWA] this is done in the procedure CALC_Z_U1, file awa.p .
%
% Using inclusions y0 and y1( = yhat + z) at t_j and t_{j+1}, respectively, 
% and the slope interval l1, which contains all derivatives between t_j and 
% t_{j+1}, a better inclusion ytau of the solution on the whole time 
% interval tau = [tj,t_{j+1}] is computed by considering the intersection 
% of the two cones C1 := y0+(t-tj)*L1 and C2 := y1+(t-t(j+1))*L1, t in tau.
%
% The better inclusion of ytau is derived by using cubically convergent 
% remainder forms stated in section 4 of
%
% [LC] Cornelius, H., Lohner, R.: Computing the range of values of real 
%        functions with accuracy higher than second order, 
%        Computing 33, 331-347, 1984.
%
% Using the improved ytau, a new slope interval l1 as well as a new 
% interval enclosure z for the local error are computed, which is used for 
% a new enclosure of y1 := yhat + z;   
%--------------------------------------------------------------------------
    tau_tcoeff = t_tcoeff;
    tau_tcoeff.inf(1) = tau.inf;
    tau_tcoeff.sup(1) = tau.sup; % taylorcoeff-version of integration time interval tau = [t,t_next].
    tc = taylorcoeff_compute(odefun,ytau,hdk,p,p,treenr_1,y_tcoeff,tau_tcoeff); % Compute verified Taylor coefficient enclosure on the whole interval tau.
    for i = 1:n % Write generalized Taylor coefficient (GTC) into a matrix C with n rows and p+1 columns.
                % The entries are C_{i,j} := "(j-1)-th GTC of i-th component function", where 1 <= i <= n, 1 <= j <= p+1.  
        C.inf(i,:) = tc(i).inf.'; 
        C.sup(i,:) = tc(i).sup.'; 
    end
    l1.inf = C.inf(:,2);    % The interval vector l1 contains all scaled first derivatives y'(t)*h_ = f(y(t))*h_, t in tau = [t,t_next] (if h > 0), tau = [t_next,t] (if h < 0).
    l1.sup = C.sup(:,2);    %  
    l1 = iv_rdivide(l1,h_); % Division by h_ = intval(t_next)-t implies that l1 contains all derivatives  y'(t) = f(y(t)).
                            % Since h_ is a thin interval, this seems acceptable. In [AWA] this corresponds to "L1K := T[K,1]/H", procedure "CALC_Z_U1", file awa.p.
                            % This is a minor drawback of the concept of generalized Taylor coefficients where the true Taylor coefficients are already scaled with powers of h.
                            % Here this scaling must be reversed to get an inclusion of the true first derivative.
    
    z.inf = C.inf(:,p+1);   % An enclosure z of the local error is given as generalized Taylor coefficient of order p in the p+1-th column of C,
    z.sup = C.sup(:,p+1);   % see [L], Equation (2.1.2), p. 37.
    
    y1 = iv_plus(yhat,z);   % y1 = yhat + z is an enclosure of the solution at t_next.
    y2tau.inf = C.inf(:,3); % Scaled second derivatives are contained in third column of C.
    y2tau.sup = C.sup(:,3); 
        
    idx1 = or(l1.sup <= 0,l1.inf >= 0); % Indices of components y(i) where the derivative does not change its sign. 
                                        % => y(i) is monotonous on tau. => The image of y(i) on tau is contained in the convex hull of y0(i) and y1(1).  
    in0_y2tau = iv_in(0,y2tau);
    idx2 = and(~idx1,in0_y2tau);        % Indices of non-monotonous components where curvature becomes zero.
    idx3 = and(~idx1,~in0_y2tau);       % Indices of non-monotonous components where curvature has constant sign.
    
    % Case 1: monotonous components i   
    %         =>   ytau(i) := hull(y0(i),y1(i))     
    hlp1 = iv_hull(y0,y1);    
    ytau.inf(idx1) = hlp1.inf(idx1);
    ytau.sup(idx1) = hlp1.sup(idx1);
    
    % Case 2: non-monotonous components i where the curvature becomes zero 
    %         => ytau(i) := hull(y0,y1) - y2tau*[0,0.25]  
    %         This is the degenerate form of the cubically convergent remainder formula (44), p. 342, of [LC], for m := 0.    
    hlp2 = iv_minus(hlp1,iv_times(y2tau,iv_0q)); % hlp2 = hull(y0,y1) - y2tau*[0,0.25]
    ytau.inf(idx2) = hlp2.inf(idx2);
    ytau.sup(idx2) = hlp2.sup(idx2);
    
    % Case 3: non-monotonous components i where the curvature has constant sign
    %         => ytau(i) = 0.5*(y0(i)+y1(i)) + 0.25*( y2tau_mid(i) * (([-1,1]+y0y1m)^2-y0y1m^2-1) - (y2tau(i)-y2tau_mid(i))*[0,1] ), see below for the definitions of y0y1m.
    %         This is the cubically convergent remainder formula (44), p. 342, of [LC].    
    %
    % preparation
    y2tau_mid = direction*iv_getmidrad(y2tau);    % y2tau_mid = direction*mid(y2tau);    
    y0_idx3.inf = y0.inf(idx3);
    y0_idx3.sup = y0.sup(idx3);
    y1_idx3.inf = y1.inf(idx3);
    y1_idx3.sup = y1.sup(idx3);
    y2tau_idx3.inf = y2tau.inf(idx3);
    y2tau_idx3.sup = y2tau.sup(idx3);    
    y0y1m = iv_rdivide(iv_minus(y1_idx3,y0_idx3),y2tau_mid(idx3)); % y0y1m = direction*(y1(k)-y0(k))./y2tau_mid;    
    %
    % Compute  ytau(idx3) = 0.5*(y0(idx3)+y1(idx3)) + 0.25*( y2tau_mid(idx3) * (([-1,1]+y0y1m)^2-y0y1m^2-1) - (y2tau(idx3)-y2tau_mid(idx3))*[0,1] );    
    x3 = iv_times( 0.5 , iv_plus( y0_idx3 , y1_idx3 ) );           % x3 = 0.5 *( y0(idx3)+y1(idx3) );
    x4 = iv_sqr( iv_plus(iv_m11,y0y1m) );                          % x4 = ([-1,1]+y0y1m)^2
    x5 = iv_minus( iv_minus( x4 , iv_sqr(y0y1m) ) , 1 );           % x5 = ([-1,1]+y0y1m)^2 - y0y1M^2 - 1 
    x6 = iv_times(y2tau_mid(idx3),x5);                             % x6 = y2tau_mid(idx3)*( ([-1,1]+y0y1m)^2 - y0y1M^2 - 1 )
    x7 = iv_times( iv_minus(y2tau_idx3,y2tau_mid(idx3)) , iv_01 ); % x7 = (y2tau(idx3)-y2tau_mid(idx3))*[0,1]
    hlp3 = iv_plus( x3 , iv_times( 0.25 , iv_minus(x6,x7) ) );        
    ytau.inf(idx3) = hlp3.inf;
    ytau.sup(idx3) = hlp3.sup;
        
    % The new, better enclosure ytau is now used to compute better enclosures l1, z, y1.
    tc = taylorcoeff_compute(odefun,ytau,hdk,p,p,treenr_1,y_tcoeff,tau_tcoeff); % Compute verified Taylor coefficient enclosure on the whole time interval tau = [t,t_next].    
    for i = 1:n % Write generalized Taylor coefficient into matrix C again.
        C.inf(i,:) = tc(i).inf.'; 
        C.sup(i,:) = tc(i).sup.'; 
    end
    l1.inf = C.inf(:,2);    % scaled first derivatives y'(t)*h_ = f(y(t))*h_.
    l1.sup = C.sup(:,2);    
    l1 = iv_rdivide(l1,h_); % better enclosure of the interval slope (on the whole time interval tau) 
    
    z.inf = C.inf(:,p+1);   % better enclosure of the local error (on the whole time interval tau) 
    z.sup = C.sup(:,p+1);  
    
    y1 = iv_plus(yhat,z);   % better enclosure of the solution at the time point t_next
        
%--- end of calculation of z and y1    
      
%--------------------------------------------------------------------------                                            
%  Estimation of global error, part II, see [L], sections 2.3 and 2.4.,
%  p. 46 et seq.                                           
%--------------------------------------------------------------------------   

    yapp0 = iv_getmidrad(iv_plus(yapp,z)); % yapp0 = mid(yapp + z), new, improved approximation at t_next 
    z_local = z;                           % Update the local error z_local. This is stored as final coefficient of the Taylor expansion in the output array Z. 
                                           % Furthermore, z_local is used for heuristic step size control in the beginning of the time integration loop.     
    ytau_local = ytau;                     % Update ytau_local which is only used for heuristic step size control in the beginning of the time integration loop, see variable "err". 
    
    z = iv_plus( iv_plus( iv_minus(yapp,yapp0) , z ) , iv_mtimes( iv_minus(ima1,abbmat) , ydev0 ) );  % z = yapp - yapp0 + z + (ima1-abbmat)*ydev0;
    z = iv_hull(z,0);
    
    if time_steps == 0 % first time step 
        rest0 = z;
        rest1 = z;
    else
        if e_meth == 0
            rest0 = iv_plus(rest0,z); %  rest0 = rest0 + z;  
        else
            phi = iv_getmidrad(ima0);
            if e_meth == 2 || e_meth == 4
                % The following code corresponds to [AWA] function call "QR_ZERL(PHI,REST1);", see files awa.p and awarout.p.
                rest1_diam = iv_diam(rest1);
                for i = 1:n
                    len(i) = rest1_diam(i) * norm(phi(:,i)); % (non-verified) Euclidean length of i-th column of phi weighted with rest1_diam(i).                                                  
                end
                [~,idx] = sort(len,'descend');
                phi = phi(:,idx);  % The matrix phi is sorted in descending order with respect to weighted column lengths.    
                [phi,~] = qr(phi); % Perform a (non-verified) QR decomposition of phi and store the orthogonal part Q in phi. 
            end
            restmat = phi;
            iphi = intval2iv(inv(intval(phi))); % enclosure of phi^{-1}
            rest1 = iv_plus( iv_mtimes( iv_mtimes(iphi,ima0) , rest1 ) , iv_mtimes(iphi,z) ); %  rest1 = ( iphi * ima0 ) * rest1 + iphi * z;
            if e_meth > 2
                rest0 = iv_intersect( iv_plus(rest0,z) , iv_mtimes(restmat,rest1) );          %  rest0 = intersect( rest0 + z , restmat * rest1 );
                rest1 = iv_intersect( rest1 , iv_mtimes(iphi,rest0) );                        %  rest1 = intersect( rest1 , iphi * rest0 );
            end;
        end;
    end;
    
    switch e_meth
        case 0
            rest = rest0;
        case {1,2}
            rest = iv_mtimes(restmat,rest1); %  rest = restmat * rest1;
        case {3,4}
            rest = iv_intersect( iv_mtimes(restmat,rest1) , rest0 ); % rest = intersect( (restmat * rest1) , rest0 );
    end;
    y1 = iv_intersect(y1 ,iv_plus( iv_plus( rest , iv_mtimes(abbmat,ydev0) ) , yapp0 ) ); % y1 = intersect( y1 , ( rest + abbmat * ydev0 + yapp0 ) );

    % Initializations for next time step
    t = t_next;     
    t_tcoeff.inf(1) = t; % Adapt type taylorcoeff version of t. 
    t_tcoeff.sup(1) = t;    

    y0 = y1;
    
    time_steps = time_steps + 1;         % Increase time step counter.
    T(time_steps+1) = t;                 % Write t to output array T.
    Y.inf(time_steps+1,:) = y0.inf';     % Write y0 to result array Y.
    Y.sup(time_steps+1,:) = y0.sup';             
    for i = 1:n
        for p_ = 1:p
            Z.inf(time_steps,i,p_) = Z_(i).inf(p_);  % Write Taylor polynomial Z_ to result array Z.
            Z.sup(time_steps,i,p_) = Z_(i).sup(p_);  
        end
        Z.inf(time_steps,i,p+1) = z_local.inf(i);    % Append enclosure of local error as p-th coefficient to Taylor polynomial Z_.
        Z.sup(time_steps,i,p+1) = z_local.sup(i);
    end
    
    abort = false;
    heps = 1e308;
    y0_diam = iv_diam(y0);
    y0_abs = iv_abs(y0); 
    y0_in0 = iv_in(0,y0);
    for i = 1:n
        act_eps = y0_diam(i);
        if ~y0_in0(i) 
            heps = act_eps / y0_abs.inf(i);
        end
        if heps < act_eps
            act_eps = heps;
        end
        abort = abort || ( act_eps > eps_ );
    end;
    
    if abort
        disp(' ')
        disp('Accuracy too bad: integration aborted !');
        if y0.inf ~= y0.sup
            disp('                  Try again with different error tolerances and/or');
            disp('                  larger order p and/or reduce size of initial value set.');
        else
            disp('                  Try again with different error tolerances.');
        end
        disp(' ');
    end;
end
%**************************************************************************
%**************************************************************************
%
%   End of main integration loop
%
%**************************************************************************
%**************************************************************************

% Shrink output arrays T, Y to the true number of time steps. 
idx = 1:time_steps+1;
T = T(idx);                 
Y.inf = Y.inf(idx,:);  
Y.sup = Y.sup(idx,:);  

if rndold ~= 1
    setround(rndold)
end

end % function awa
