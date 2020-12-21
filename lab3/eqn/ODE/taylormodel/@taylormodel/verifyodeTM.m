function [T,Y,Yr] = verifyodeTM(dummy_tm,odefun,tspan,y0,options)
% VERIFYODETM  verified ODE-solver for the initial vaue problem  
%               
%   y'(t) = odefun(t,y(t))  with t in (tspan(1),tspan(end)) 
%   y(tspan(1)) = y0.
%
%   call:  [T,Y,Yr] = verifyodeTM(dummy_tm,odefun,tspan,y0,options)
%
% This verified ODE-solver is based on Taylor models.   
% The five input parameters "odefun, tspan, y0, options" have the same meaning as for other (non-verfied) 
% MATLAB ODE-solvers, for example, ode45. 
% 
% Input parameters:  
%
% 1) "dummy" is a dummy - Taylor model to reach the folder @taylormodel. It has no other purpose.  
% 2) "odefun" 
%     For a given ODE y′ = f(t,y), "odefun" is a function handle to f: R^n x [t0,tf] -> R^n. 
%     Here t is the time variable ranging between t0 := tspan(1) and tf := tspan(end), and n is the dimension of the ODE. 
%     This means f(t,y) returns the column vector of length n containing the components f(1),...,f(n). 
%     In addition, the call f(t,y,i) returns the single components f(i).
%
%     Important note: It is assumed that each component function f_i fulfills the following:
%        a) f_i consists only of arithmetic operations +,-,*,/ and of standard functions
%           like exp,log, sin, cos, tan,... ,
%           see @taylormodel for the complete list of implemented standard functions.
%        b) f_i is (m+1)-times continuously differentiable, where m is the user specified order of the 
%           used Taylor models (= maximum degree of the multivariate polynomial part).
%        c) f_i is implemented rigorously. For example, non-floating-point constants like pi
%           must be implemented by intval('pi'), see @intval/intval. 
%  
% 3) "tspan" contains the starting point and the endpoint of integration, i.e. tspan = [t0,tf].
%     (Contrary to MATLAB ode-solvers, tspan is not allowed to contain additional intermediate points.)  
%     Both t0 and tf must be floating point numbers. In particular this means that if, for example, 
%     a non-floating-point starting time t0, say t0 = 0.1, is necessary, then the ode-function odefun   
%     must be changed/rescaled by the user so that the starting time t0 is moved to a floating-point number. 
%     In general, such a rescaling is easily done. In many/most practical cases we simply have t0 = 0. 
%     
%
% 4) y0 = (y0(1);...;y0(n)) contains the starting values. This can be a float or an interval vector.  
%    (A row vector will automatically be transformed into a column vector.)
%
% 5) "options" contains additional parameters. Analogously to other MATLAB ode-solvers, "options"  
%     can be set by using the function verifyodeset: options = verifyodeset('name1',value1,'name2',value2,...).
%     The following parameters can be transferred (see also documentation of odeset.m) : 
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
%                         0 = off
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
%          MSU HEP Report 40910, 2003 
%  [NJN] M. Neher, K.R. Jackson, N.S. Nedialkov, "On Taylor model based integration of ODEs", 
%          SIAM J. Numer. Anal. 45(1), pp. 236-262, 2007
%  [Bue] F. Buenger, "Shrink wrapping for Taylor models revisited", Numerical Algorithms 78(4), pp. 1001-1017, 2018

% written  11/04/15     F. Buenger
% modified 01/18/16     F. Buenger  record feature, see [E], Section 4.4.3.1, p.112.   
% modified 05/15/17     F. Buenger  preconditioning of Taylor models   
% modified 06/29/17     F. Buenger  blunting   
% modified 07/24/17     F. Buenger  stepwise order increase in Picard iteration    
% modified 07/07/18     F. Buenger  initial stepwise h0 can be chosen automatically    

e = 1e-30;
if 1+e > 1   % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

global INTLAB_ODE_OPTIONS % global copy of ODE-Options that can be set by the user in the 'options' parameter
global INTLAB_ODE_VARS    % additional, program internal, global variables that cannot be set by the user  
global INTLAB_CONST       % global INTLAB constants, only used for for prgress window 

finishup = onCleanup(@() cleanup); % Call cleanup function if programm is terminated.

%**************************************************************************
%**************************************************************************
%
% preparations
%
%**************************************************************************
%**************************************************************************

%-------------------------------------------------------------------------- 
% Initialize global parameter INTLAB_ODE_VARS.ODEMODE .
%-------------------------------------------------------------------------- 

% 1: verified computation in Taylor model arithmetic
% 0: non-verified computation in Taylor model arithmetic  
%    In particular, error intervals and images are not computed and set to INTLAB_ODE_VARS.EMPTYIV.
%    Compares to NO_REMAINDER in Riot, see [E], Section 4.4.3.2, p.112-113, and code on p.121 for description.
VERIFIED = true;      % constant for verified ODEMODE
NON_VERIFIED = false; % constant for non-verified ODEMODE                            
INTLAB_ODE_VARS.ODEMODE = VERIFIED; 

%-------------------------------------------------------------------------- 
% Initialize global parameter INTLAB_ODE_VARS.RECMODE .
%--------------------------------------------------------------------------   

% The record mode controls the so called "record feature", 
% see [E], section 4.4.3.1, p.112, and code on p.120-127 for description.
% It enhances the performance drastically!
%
% Three values are distinguished: 
%   0 : Record feature off (default). This compares to NORMAL_MODE in Riot, see [E], p.121,127. 
%   1 : Record feature on, write mode. This compares to REMAINDER_REC in Riot, see [E], p.122. 
%   2 : Record feature on, read mode. This compares to REMAINDER_PLAY in Riot, see [E], p.123,126.
REC_OFF = 0;   % constant for record feature off
REC_WRITE = 1; % constant for record feature on, write mode 
REC_READ = 2;  % constant for record feature on, read mode                              
INTLAB_ODE_VARS.RECMODE = REC_OFF;                                
  
%-------------------------------------------------------------------------- 
% Check y0.
%-------------------------------------------------------------------------- 

y0 = intval(y0);
s = size(y0); 
if (length(s) ~= 2) || (s(1) ~= 1 && s(2) ~= 1) 
    error('invalid call of verifyode: y0 must be a vector.')
end
if s(2) > s(1)  % The row vector y0 is transformed to a column vector.
  y0 = y0';    
end
n = length(y0); % dimension of the ODE-System y' = f(y,t), f: R^n x [t0,tf] -> R^n

%-------------------------------------------------------------------------- 
% Check tspan.
%-------------------------------------------------------------------------- 

s = size(tspan); 
if (length(s) ~= 2) || (s(1) ~= 1 && s(2) ~= 1) || ~isa(tspan,'double') || ~isreal(tspan) || tspan(1) >= tspan(2) 
    error('invalid call of verifyode: improper parameter tspan.')
else
    t0 = tspan(1);    % starting time for integration    
    t_end = tspan(2); % t_end = tf , end time for integration
end

%-------------------------------------------------------------------------- 
% Set global parameter INTLAB_ODE_VARS.SPARSE.
%-------------------------------------------------------------------------- 

if n >= 4 % switch on sparse convolution for higher dimensional ODEs. 
          % The dimension bound 4 is chosen quite arbitrarily. Feel free to change that!
    INTLAB_ODE_VARS.SPARSE = true;  
else
    INTLAB_ODE_VARS.SPARSE = false;
end

%-------------------------------------------------------------------------- 
% Transfer "options" to global structure INTLAB_ODE_OPTIONS
% All parameters not specified in "options" get default values.
%-------------------------------------------------------------------------- 

if nargin < 5
    options = [];
end
opt2glob(options); 

order = INTLAB_ODE_OPTIONS.order;  
loc_err_tol = INTLAB_ODE_OPTIONS.loc_err_tol;  
h_min = INTLAB_ODE_OPTIONS.h_min;  

%-------------------------------------------------------------------------- 
% Create an "identity" Taylor model base_tm with polynomial part 
%
%   base_tm(x_1,...,x_n,t) = (x_1,...,x_n), x_i in [-1,1], i = 1,...,n.
%
% Note that base_tm does not depend on the time variable "t"  
% which formally is the (n+1)-th input parameter. 
%-------------------------------------------------------------------------- 
                                              
% All domains of the space-like variables x_1,...,x_n are [-1,1]. 
% The time variable t gets the dummy domain [0,0], which, as for all other
% Taylor models, will be adapted in each integration step by [t_i,t_{i+1}].
% Thus, the whole domain for (x_1,...,x_n,t) is [-1,1]^n x [0,0].
v.sup = [ones(1,n),0]; 
v.inf = -v.sup;
domain = mat2cell(repmat(v,n,1),ones(1,n),1); 
center = mat2cell( zeros(n,n+1) , ones(1,n) , n+1 ); % The center points of all domains for x_1,...,x_n, and t are set as zero.
type = ones(n,1);                                    % Taylor models with type = 1 must have a standard domain [-1,1]x...x[-1,1]x[t,t+h], t, h >= 0, 
                                                     % and center point (0,...,0,t). 
monomial = num2cell([eye(n),zeros(n,1)],2);          % cell array: monomial{i} = (0,...,0,1,0,...,0) where the 1 is at position i = 1,...,n. 
                                                     % monomial{i} has length n+1 so that the last entry is always zero, i.e., no time dependence.   
coefficient = num2cell(ones(n,1));                   % cell column vector with 1 in each component, i.e., coefficient{i} = 1 for i = 1,...,n   

% Create n elementary Taylor models base_tm(i) with polynomial part p_i(x_1,...,x_n,t) = x_i , i = 1,...,n. 
% All have the same order, domain = [-1,1]^n x[0,0], center point center = (0,....,0) and trivial error interval [0,0].
% These elementary Taylor models serve in some sense as "basis" polynomials to build up others easily. 
base_tm = taylormodelinit(center,domain,order*ones(n,1),monomial,coefficient,[],type); 

%-------------------------------------------------------------------------- 
% Transform initial value set/interval y0 to a Taylor model y0_tm. 
%-------------------------------------------------------------------------- 

[mid0,rad0] = iv_getmidrad(y0); % midpoints and radii of initial value intervals
y0_tm = base_tm .* rad0 + mid0; % Transform initial values y0, which are of type float or intval, to type Taylor model

%-------------------------------------------------------------------------- 
% Initialize time integration interval [t_curr,t_curr+h_curr] 
% and corresponding time Taylor model t.
%-------------------------------------------------------------------------- 

t_curr = t0;                                % t_curr is the current time grid point for integration. It is initialized with the starting time t0. 
h_curr = INTLAB_ODE_OPTIONS.h0;             % h_curr is the current step size of the actual time integration step. 
                                            % It is initalized with the initial step size h0.
if h_curr == 0                              % If h_curr = h0 = 0, then choose the initial step size automatically.
    h_curr = choose_h0(odefun,t0,y0,order,loc_err_tol);
    if h_min > h_curr
        h_min = h_curr/10;                  % h_min is automatically reduce to 0.1*h0. The factor 0.1 is quite arbitrary. Feel free to change that! 
        disp(['Minimum step size is reduced to h_min = ',num2str(h_min)]) % Inform the user that h_min was automatically reduced.
    end
end

center = [zeros(n,1);t_curr];               % Time evaluation starts from the lower bound of the time interval [t_curr,t_next] .
clear domain;
domain.inf = [-ones(n,1);t_curr];           % [-1,1]^n x [t_curr , t_curr] is the initial time domain.
domain.sup = [ones(n,1);t_curr];
t_monomial = [zeros(1,n+1);[zeros(1,n),1]]; % Set exponents for time Taylor model p(x_1,...,x_n,t) = t_curr + (t-t_curr), t in [t_curr, t_next]. 
                                            % The time Taylor model is independent of all space-like variables x_1,...,x_n, 
                                            % which formally occur with zero exponents in t_monomial.  
coefficient = [t_curr;1];                   % coefficients of the time Taylor polynomial p(x_1,...,x_n,t) = t_curr + 1 * (t-t_curr)
t = taylormodelinit(center,domain,order,t_monomial,coefficient,[],1); % Create time Taylor polynomial.

repeat_step = false;                        % repeat_step == true means: Repeat integration step because of accuracy loss. 

%-------------------------------------------------------------------------- 
% Prepare output parameters. 
%-------------------------------------------------------------------------- 
N_block = 1e6;                                      % maximum block size for result array preallocation. If the number of time steps becomes larger, then result arrays are resized.
N_max = ceil((t_end-t0)/h_min);                     % maximum number of time steps according to minimum step size h_min and integration period [t0,t_end]  
N = min(N_block,N_max);                             % first block size for result preallocation 
T_block = zeros(N+1,1);
Y_block(1:N,1:n) = INTLAB_ODE_VARS.ZEROTAYLORMODEL; % memory preallocation for Taylor model result array Y_
T = T_block;                                        % preallocation of time step result array T;
Y(1:N,1:n) = Y_block;                              % memory preallocation for Taylor model result array Y_ 

%-------------------------------------------------------------------------- 
% Do some initializations if preconditioning is switched on.
%-------------------------------------------------------------------------- 
Yr = [];
if INTLAB_ODE_OPTIONS.precondition  
    yr = base_tm;       % The right Taylor model is initialized as the "identity" Taylor model,i.e.,
                        % the polynomail part is yr(x_1,...,x_n,t) = (x_1,...,x_n) and the error interval part is zero for each of the n components.
    lr_return = true;   % If lr_return = true, then left and right Taylor models are returned as Y.left = Y_ and Y.rigth = Yr.
    %lr_return = false; % Only for Testing! If lr_return = false, then the concatenated Taylor model Y = Y_ o Yr is returned if preconditioning is switched on.
    if lr_return
        Yr = Y;        % storage preallocation for right Taylor models if preconditioning is switched on
    end
end

%--------------------------------------------------------------------------
% Set some constants. Their values are heuristic. Feel free to change them!
%-------------------------------------------------------------------------- 

iter_1_max = 3;   % maximum number of iterations for trying to reach error interval inclusion 
                  % (= number of applications of the integral operator induced by the odefun until its 
                  % contraction property is reached to apply Banach's fixed-point theorem).
iter_2_max = 3;   % maximum number of improvement steps.               
c_improve = 0.99; % factor to measure improvement, chosen as 99 percent, see [E], p.126, 
                  % and Riot, file dglalg.cpp, function "terminate". 

% Constants for heuristic, automatic step size control:                 
                  
c_dec = 0.7;      % If error inclusion after Picard iteration was not successful with current step_size h_curr,
                  % then it is decreased to h_curr := c_dec * h_curr.
                  % c_dec = 0.7 is used in [E], Riot.         
                  
c_max = 1.1;      % maximal factor for step size increase: h_next <= c_max * h_curr
                  % The value c_max = 1.1 is taken in [E], p.130, code line 366.  
c_min = 0.7;      % minimum factor for step size increase: h_next => c_min * h_curr 
                  % Taking a minimum value c_min is not considered in [E].
                  
c_fac = 0.1^(1/(order+1));  % scaling factor, see [E], p.129, (4.9) where c_fac = 0.1^(1/(order+1)) is chosen
%c_fac = 0.8;               % The following values for c_fac are stated in 
                            % [HNW] Hairer, Norsett, Wanner, "Solving Ordinary Differential Equations I", p. 168: 
                            % c_fac = 0.8, 0.9, 0.25^(1/(order+1)), 0.38^(1/(order+1))

%-------------------------------------------------------------------------- 
% Initialize some counters which serve for program analysis.
%-------------------------------------------------------------------------- 

% print_report = true; % If this flag is "true", then the subsequent counters are displayed at the end of integration. 
print_report = false;

step_cntr = 0;     % counts steps of integration
incl_cntr = 0;     % total number of error inclusion steps after the Picard iterartion 
impr_cntr = 0;     % total number of improvement steps after error inclusion  
repeat_cntr = 0;   % total number of time step repetitions
dec_cntr_incl = 0; % counts how often the step size was decreased by the factor c_dec during error inclusion.
inc_cntr = 0;      % counts how often the step size increased 
dec_cntr = 0;      % counts how often the step size decreased by local error control
h_min_cntr = 0;    % counts how often the minimal step size was reached

%--------------------------------------------------------------------------
% miscellaneous settings
%-------------------------------------------------------------------------- 

% Define zero interval row vector, used for initializing the error interval 
% after the Picard iteration. 
start_iv.inf = zeros(n,1);
start_iv.sup = start_iv.inf;                                             

disp_progress = false;           % flag: "display progress window", initialized with false
INTLAB_CONST.POPUPWINDOW = [];
TIMER = clock;                   % After a default waiting time of t_wait = 5 seconds a progress window will be
                                 % displayed showing the actual integration time t and stepsize h.

%**************************************************************************
%**************************************************************************
%
%  end of preparations 
%
%**************************************************************************
%**************************************************************************


%**************************************************************************
%**************************************************************************
%
% main time integration loop
%
%**************************************************************************
%**************************************************************************

while t_curr < t_end
    step_cntr = step_cntr+1; % counts steps of integration
    
    % Set "t_next" and check if the step was too large.
    t_next = t_curr + h_curr;     % [t_curr , t_next] is the integration interval of the actual integration step.
                                  % Note that its diameter may differ from h_curr since the sum t_curr + h_curr might be rounded.
                                  % Thus, h_curr is just a non-verified approximate of the true step size t_next-t_curr which
                                  % is only implicitly given.
    if t_next >= t_end            % If t_next > t_end, then t_next is set to t_end.
        t_next = t_end;
        h_curr = t_next - t_curr; % Calculate new non-verified, approximate step size. 
    end
    % Adjust the current time interval and the time Taylor model "t". 
    time_interval.inf = t_curr;
    time_interval.sup = t_next;   % time_interval := [t_curr , t_next] 
    t.domain.inf(n+1) = t_curr;
    t.domain.sup(n+1) = t_next;   % [-1,1]^n x [t_curr , t_next] is the domain for the current integration step.   
    t.center(n+1) = t_curr;       % Time evaluation starts from the lower bound of [t_curr , t_next].
    t.coefficient(1) = t_curr;    % Adjust constant term.
    t.image = time_interval;
        
    % Recall that y0_tm is a Taylor model - indicated by the suffix "_tm" - 
    % the range of which represents the initial value set for the actual 
    % integration step. Its time domain is now updated to [t_curr,t_next] 
    % and t_curr is taken as "center" point/evaluation point.
    % Note that all polynomials of y0_tm (like all polynomials of base_tm) 
    % do not depend on the time variable, which stands at position n+1. 
    % Thus all monomials/exponents have the form (m_1,m_2,...,m_n,0).
        
    y0_tm = set_domain(y0_tm,n+1,time_interval,t_curr); % Replace y0_tm.domain(n+1) by "time_interval" and the  corresponding center point
                                                        % y0_tm.center(n+1) by "t_curr" which is the left bound of "time_interval". 

%-------------------------------------------------------------------------- 
% Picard iteration,      see  [E], Section 3.4.1, p.51-56, Algorithm 1,
%                                  Section 4.5.1, p.120,121.
%-------------------------------------------------------------------------- 

    INTLAB_ODE_VARS.ODEMODE = NON_VERIFIED; % non-verified computation of approximate Taylor polynomial via Picard iteration
    setround(0);
    
    if ~repeat_step
        t_ = t;
        y = y0_tm;  % Initialize the non-verified, approximate Taylor model solution y.
        % Non-verified calculation of the Taylor polynomial of the approximate solution y with Picard iteration.
        for k = 1:order
            % In the k-th iteration step, the polynomial part of y coincides with the Taylor expansion of the exact 
            % solution up to order k, see [E], Satz 11, p.52, p.121, code line 75, and explanation on p.120.
            
            % Evaluate the Picard operator:  
            %
            % (*)  y(eta,t) = eta + \int_{t_curr}^t f(s,y(s)) ds,     
            %
            % see [E], formula (3.10), p.45.
                                   
            % Stepwise increasing the orders of the Taylor models t,y,y0_tm is only done for performance reasons.
            t_ = set_order(t_,k);      
            y = set_order(y,k);
            y0_tm_k = set_order(y0_tm,k);
            
            z = odefun(t_,y);              % z := f(y,t) = odefun(t,y). For MATLAB-ode-solvers like ode45, odefun has the time variable
                                           % as its first parameter. We keep this convention even though internally the time
                                           % variable of a Taylor model stands at the last position n+1 and not at the first position.
            y = y0_tm_k + integral(z,n+1); % [E], formula (3.10), p.45; non-verified integration of all components z(i) 
                                           % with respect to the time variable which is the (n+1)-th variable.
                                           % Recall that the Taylor model y0_tm represents the initial value set  
                                           % and therefore replaces "eta" in formula (*) above.                                                      
        end        
    else
      % Reuse old solution with correct polynomial part but adapt the time domain to new time_interval = [t_curr,t_next] 
      % and the time center point to new t_curr.  
      y = set_domain(y,n+1,time_interval,t_curr);           
    end
    
%-------------------------------------------------------------------------- 
% error inclusion,       see [E] Section 3.4.2, p.56-61, Algorithm 2, 
%                                improved Algorithm 3, p.66,                
%                                Section 4.5.1, p.122-126.
%-------------------------------------------------------------------------- 

    INTLAB_ODE_VARS.ODEMODE = VERIFIED; % Switch to verified computation in Taylor model arithmetic.
    setround(1);                        % Switch to rounding upwards.
    
    inclusion = false;                  % Flag: inclusion ( for Banach fixed-point theorem ) obtained ?
    
    y_ = set_interval(y,start_iv);      % Set all error intervals of y to zero interval (no error),
                                   
    y_ = set_image(y_);                 % Compute image for subsequent verified Taylor model arithmetic       
    
    while ~inclusion        
        % Calculate the inner approximation.
        % Evaluate the integral operator K( P + [0,0] ) - P.
        % The calculated inner approximation will be stored in the interval vector "curr_iv".        
                
        hlp_ = integral_operator(odefun,y0_tm,t,y_,y_,[],REC_WRITE);       % Evaluate ode-integral operator in record write mode.       
        %full_hlp_ = integral_operator(odefun,y0_tm,t,y_,y_,[],REC_OFF);   % for comparison, record off mode      
        
        curr_iv = iv_plus(get_image(hlp_),get_interval(hlp_));             % curr_iv = get_image(hlp_) + get_interval(hlp_)
                                                                           % This implements the left side of [E], p.122, formula in line 3.
        % Calculate enclosures of the solution.
        iter_1 = 0;        
        DIAM = max(iv_diam(curr_iv));
        inclusion_ = false(n,1); % Initialize logical array which indicates componentwise inclusion .
        while ~inclusion && (iter_1 < iter_1_max) && DIAM < 0.1
            incl_cntr = incl_cntr + 1; % Count total number of inclusion steps after Picard iteration. (Just for reporting) 
            iter_1 = iter_1 + 1;
            % Replace the remainder part with the inflated inner approximation I_k.
            curr_iv = inflate(curr_iv,inclusion_); % Epsilon-Inflation.
            y = set_interval(y_,curr_iv); % Define current inclusion candidate curr_iv as error interval for y.
            %  Evaluate the integral operator: K( P + I_k ) - P =: J_k 
            inclusion = true;
           
            for i = 1:n               
              hlp_i = integral_operator(odefun,y0_tm,t,y,y_,i,REC_READ);   % Evaluate ode-integral operator K( P + I_k ) - P in record read mode.               
              % full_hlp = integral_operator(odefun,y0_tm,t,y,y_,i,REC_OFF); % for comparison, record off mode      
              
              curr_iv_i = iv_plus(hlp_(i).image,hlp_i.interval);           % curr_iv_i = hlp_(i).image + hlp_i.interval 
                                                                           % This implements the left side of [E], p.122, formula in line 3.                                                                 
              curr_iv.inf(i) = curr_iv_i.inf;
              curr_iv.sup(i) = curr_iv_i.sup;
              inclusion_(i) = iv_in(curr_iv_i,y(i).interval);
              inclusion = inclusion && inclusion_(i);                      % Check if there is inclusion.  
              
              % If inclusion holds in the current component use the calculated (better) enclosure in the following calculations.                
              % if inclusion 
              if inclusion_(i)
                  % Use the new, better enclosure for the i-th component to obtain faster convergence for the remaining components i+1,...,n.
                  y(i).interval = curr_iv_i;                                                                          
              end
            end
        end 
        
        % If failure, then decrease the step size.
        if ~inclusion
            % No inclusion has been reached, because the time interval of the integration step was too big. 
            % => Decrease the step size by factor c_dec. 
            h_curr = c_dec * h_curr;   
            dec_cntr_incl = dec_cntr_incl + 1;                             % Count how often the step size was decreased during error inclusion. (Just for reporting)
            if h_curr < h_min 
                error('Step size undergoes minimal step size. Inclusion is not possible.'); % See [E], Section 3.5.1, p.65-67.
            end
            t_next = t_curr + h_curr;                                      % Adapt t_next according to new step size. As before, h_curr is just a non-verified
                                                                           % approximate of the true step size  t_next - t_curr, which is only implicitly given, 
                                                                           % since the sum t_curr + h_curr might be rounded.
            time_interval.inf = t_curr;
            time_interval.sup = t_next;                                    % Change the data of the current time interval.
            
            % Adapt time domain of active Taylor models t, y0_tm, y, and y_:     
            % Set the time domain to time_interval = [t_curr,t_next] wich is at the last position n+1 in the domain interval vector.
            % The time center point remains t_curr.            
            t = set_domain(t,n+1,time_interval,t_curr); 
            y0_tm = set_domain( y0_tm,n+1,time_interval,t_curr);
            y = set_domain(y,n+1,time_interval,t_curr);  
            y_ = set_interval(y,start_iv);                                 % Solution with all error intervals set to zero.               
        end        
    end % while ~inclusion
    
%-------------------------------------------------------------------------- 
% Improve the error inclusion by repeated application of the 
% integral operator to the current solution,  
% 
% see  [E]  Sec.3.5.1, p.65-68; Algorithm 4, 
%           Section 4.5.1, p.126-127.
%-------------------------------------------------------------------------- 

    improve = true;   
    % improve = false; % Only for testing !
    iter_2 = 0;

    while improve && iter_2 < iter_2_max
    % while improve % Improve as long as significant without iteration bound.
        
        impr_cntr = impr_cntr + 1; % Count total number of improvement steps after error inclusion. (Just for reporting)
        iter_2 = iter_2 + 1;
        improve = false;
        
        for i = 1:n 
            curr_iv = y(i).interval;                                       % Save old remainder interval for a later comparison.   
            
            hlp = integral_operator(odefun,y0_tm,t,y,y_,i,REC_READ);       % Evaluate ode-integral operator in record read mode. 
            % full_hlp = integral_operator(odefun,y0_tm,t,y,y_,i,REC_OFF); % for comparison, record off mode    
            
            better_enc = iv_plus(hlp_(i).image,hlp.interval);              % better_enc = hlp_(i).image + hlp.interval                       
            y(i).interval = better_enc;                                    %  Set new error interval of y(i).       
            
            improve_i = (iv_diam(better_enc)/iv_diam(curr_iv)<c_improve ); % Check if the improvement was significant enough to continue the while loop.
            improve = ( improve || improve_i );              
        end
    end 
           
%-------------------------------------------------------------------------- 
% Calculate the solution set at t_next and estimate the local error growth,
% see [E] Section 4.5.1, p.128.
%-------------------------------------------------------------------------- 

    % The following variables I_curr, I_next, eps_, loc_err are only(!) used for the heuristic automatic step size control.

    hlp = substitute(y,t_curr,t_next); % Substitute the time variable [(n+1)-th variable] in all y(i) by t_next.       
    I_next = get_interval(hlp);        % error interval of solution at t_next
    I_curr = get_interval(y0_tm);      % error interval of solution at t_curr
    
    eps_ = max(I_next.sup-I_curr.sup,I_curr.inf- I_next.inf); % non-verified computation of error interval growth from t_curr to t_next, see [E], p.128, line 2                                                                          
    loc_err = max(eps_);               % local error growth estimate, see [E], p.128, (4.8)

%--------------------------------------------------------------------------
% automatic step size control, see 
% [E]   Section 4.5.1, p.129-131,                                  
% [HNW] Hairer, Norsett, Wanner, "Solving Ordinary Differential Equations I". 
%--------------------------------------------------------------------------

    % Choose a heuristic measure of the local error.
    % Feel free to make a different/better choice! (Its a bit reading tea leaves.)

    err = loc_err/loc_err_tol;
    c = (0.1*1/err)^(1/(order+1));          % Determine heuristic factor for new step size. This choice comes from in [E], p.129, (4.9).                                    

    % Alternative choices for err and c
    
%     err = norm(eps_/loc_err_tol)/sqrt(n); % Compare with [HNW], p.168, (4.11).    
%     c = (1/err)^(1/(order+1));
%     c = min(c_max,max(c_min,c_fac*c));    % The factor is bounded by c_min and c_max from below and above, respectively.
%                                           % Compare with [HNW], p.168, (4.11)-(4.13).                                                          
                                       
    % Correct the step size and make sure the step size is not too small.
    h_curr = max( c * h_curr , h_min );        
            
    if  err > 1 && h_curr > h_min      % Repeat integration step.  
        repeat_step = true;            % Set repeat flag.    
        repeat_cntr = repeat_cntr + 1; % Increase step repetition counter. (Just for reporting)     
        step_cntr = step_cntr - 1;     % Decrease integration step counter.
    else                               
        repeat_step = false;           
        if c < 1
            dec_cntr = dec_cntr + 1;   % Count how often the step size was decreased. (Just for reporting) 
        elseif c > 1
            inc_cntr = inc_cntr + 1;   % Count how often the step size was decreased. (Just for reporting)
        end
    end    
    
    if h_curr == h_min
        h_min_cntr = h_min_cntr + 1 ;  % Count how often the minimum step size was reached. (Just for reporting)
    end
    
%--------------------------------------------------------------------------
% - shrink wrapping:  The implementation is according to [Bue].
%                     
%               different to : [E]  Section 3.5.2, p.67-78; Algorithms 5,6, 
%                              [E]  Section 4.5.1, p.131-132
%                              [MB] 
%
% - preconditioning: The implementation is based on [NJN] and [MB].                        
%--------------------------------------------------------------------------

    % If the time step needs not to be repeated then calculate the 
    % solution set at 't_next' and at the intermediate grid points 
    % between 't_curr' and 't_next', do shrink wrapping or preconditioning  
    % and determine the initial value set for the next integration step.
    
    if ~repeat_step        
        if step_cntr > N
           N_ = min(N_max-N,N_block);   % new block size for preallocation 
           N = N + N_;                  % new, extended array size 
           T = [T;T_block(1:N_)];       % Resize T to length N+1. 
           Y = [Y;Y_block(1:N_)];      % Resize Y_ to length N. 
           if INTLAB_ODE_OPTIONS.precondition && lr_return
               Yr = [Yr;Y_block(1:N_)]; % Resize Yr to length N.
           end
        end  
        
        T(step_cntr) = t_curr;          % Save current time step in output array T.
        
%         % If preconditioning is switched on and if the scaling vector S is non-empty, then y is the scaled integrated left Taylor model,
%         % which after integration must be rescaled to y(x,t) := y(x./S,t).
%         if  INTLAB_ODE_OPTIONS.precondition && ~isempty(S) % Only in the first integration step S is empty.
%             y = rescale(y,S);          % y := y(x./S,t)          
%             %hlp = rescale(hlp,S);     % hlp := hlp(x./S,t)   
%             hlp = substitute(y,t_curr,t_next);
%         end
        
        Y(step_cntr,:) = y'; % Save the Taylor model inclusion for the current time interval as a row vector in return array Y. 
                              % If preconditioning is switched on, then save only the left Taylor model for the current time interval as a row vector in the return array Y.
                              % See the description of the output parameters [T,Y] at the beginning of this function.  
                                                                    
        if INTLAB_ODE_OPTIONS.precondition             
            % First rescale the integrated scaled left Taylor model y(x,t) by substituting each space-like variable x_i by x_i/S_i,
            % where S = (S_1,...,S_n)^T are the scaling constants for scaling each component of the image of the right Taylor model yr
            % to the standard domain [-1,1], i.e. the concatenation S o yr(x,t) has its image in [-1,1] componentwise.
            
            yr = set_domain(yr,n+1,time_interval,t_curr); % Sadly the time domain must always be adapted even though yr is time independent.              

            if lr_return
                Yr(step_cntr,:) = yr';  % Save the right Taylor models yr as one row vector in the return array Y.right = Yr. 
            else
                ys = concatenate(y,yr); % y is the integrated left Taylor model. The Taylor model inclusion ys of the ODE-flow for the current time interval
                                        % is the concatenation of y and the actual right Taylor model yr, i.e., ys = y o yr.                                                                                                                  
                Y(step_cntr,:) = ys';  % Save the Taylor model inclusion for the current time interval as a row vector in the return array Y.  
            end
        end
                                
        % Calculate the solution set at t_next which serves as initial value set for the next integration step.        
        y = hlp; % Recall that hlp = substitute(y,n+1,t_next).         
        
        % Apply shrink wrapping to Taylor model y.
        %
        % Remark: If preconditioning is switched on, then y is the integrated left Taylor model.
        %         Since shrink wrapping changes the dependencies on the space-like variables, it cannot be applied to a left part of a Taylor model factorization
        %         but only to the right part which is done inside the function precondition if both, preconditioning and shrink wrapping, are switched on.
        
        if INTLAB_ODE_OPTIONS.shrinkwrap && ~INTLAB_ODE_OPTIONS.precondition && t_next < t_end 
            base_tm = set_domain(base_tm,n+1,time_interval,t_curr); % Sadly the time domain must always be adapted even though base_tm is time independent. 
            y = shrinkwrap(y,base_tm);                              % Try shrink wrapping. If this fails, then other, more simple strategies are performed.
        end         
        
        % Apply preconditioning 
        if INTLAB_ODE_OPTIONS.precondition && t_next < t_end 
            base_tm = set_domain(base_tm,n+1,time_interval,t_curr); % Sadly the time domain must always be adapted.
            [y,yr] = precondition(y,yr,base_tm);                    % Compute new left and right Taylor models for the next integration step.
        end 
        
        y0_tm = y;       % y0_tm represents an enclosure of the initial values for the next integration step.

        if disp_progress 
            str = ['t = ',num2str(t_curr,'%.4e'),'   h = ',num2str(h_curr,'%.2e')];
            progress(str,314);            
        else
            disp_progress = ~(INTLAB_CONST.MINTIME == inf || isempty(INTLAB_CONST.MINTIME)) && ...
                             (etime(clock,TIMER) > INTLAB_CONST.MINTIME);            
        end 
        
        t_curr = t_next; % Set t_curr for next time step.      
        
    end 
end 

%**************************************************************************
%**************************************************************************
%
%  end of main time integration loop 
%
%**************************************************************************
%**************************************************************************

T(step_cntr+1) = t_end;     % Store also endpoint of integration in T. 
T((step_cntr+2):N+1) = [];  % Delete trailing unused preallocated entries T.
Y((step_cntr+1):N,:) = []; % Delete trailing unused preallocated rows in Y_.
if INTLAB_ODE_OPTIONS.precondition
    if lr_return % In this case left and right Taylor models are returned.                  
        Yr((step_cntr+1):N,:) = []; % Delete trailing unused preallocated rows in Yr.
    else
        Yr = [];
    end
end

if print_report
    % output of internal program counters which only serve for program analysis
    disp(' ')
    disp(['Number of integration steps : ' num2str(step_cntr)])
    disp(['          inclusion steps   : ' num2str(incl_cntr)])
    disp(['          improvement steps : ' num2str(impr_cntr)])
    disp(['          step repetitions  : ' num2str(repeat_cntr)])
    disp(['          h-decreases incl. : ' num2str(dec_cntr_incl)])
    disp(['          h-decreases       : ' num2str(dec_cntr)])
    disp(['          h-increases       : ' num2str(inc_cntr)])
    disp(['          h_min hits        : ' num2str(h_min_cntr)])
    disp(' ')
end

if rndold ~= 1 % Set old rounding mode with which verifyode was called. 
               % If this was already "rounding upwards", nothing has to be done.
    setround(rndold);
end

end % function verifyode

%--------------------------------------------------------------------------

function cleanup()

global INTLAB_ODE_VARS     

% Clean up global setting in INTLAB_ODE_VARS:
INTLAB_ODE_VARS.RECMODE = 0;    % Switch recording feature off.
INTLAB_ODE_VARS.RECLIST = [];   % Clear record list.
INTLAB_ODE_VARS.RECIDX = [];    % Clear record index.
INTLAB_ODE_VARS.ODEMODE = 1;    % Switch to verified mode (default).
INTLAB_ODE_VARS.SPARSE = false; % Switch to full convolution (default).          

progress(-1); % Close progress window.

end % function cleanup
