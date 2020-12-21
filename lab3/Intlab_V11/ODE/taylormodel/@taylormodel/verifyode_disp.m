function [ts,y] = verifyode_disp(T,Y,Yr,y0,t,print,fmt)
% VERIFYODE_DISP  computes and displays the inclusions of an ODE at specified time grid points.             
%
%   res = verifyode_disp(T,Y,Yr,y0,t,print)
%                       
%   It is assumed that [T,Y,Yr] is the output of a call  
%  
%       [T,Y,Yr] = verifyode(odefun,tspan,y0,options)
%
%   for solving an ODE given by the function handle odefun
%   on a time interval tspan with interval initial value y0 
%   (type intval). 
% 
%  t: Optional parameter. If t is empty (default), then inclusions at all   
%     automatically determined time grid points T are displayed, i.e., t := T is taken.
%     If t is not empty, then inclusions and diameters at all time points t(k) are displayed.
%     It is assumed that t is sorted in ascending order. Only entries of t that are 
%     contained in the integration interval [T(1),T(end)] are considered.
%
%   print: Optional parameter. Inclusion intervals and their diameters are 
%          displayed. Default = true;
%
%   fmt: number format, e.g. 'longe'
%   y: output argument of type intval. y(k,i) is an enclosure of y_i(ts(k)) 
%   ts: time grid points corresponding to y sorted in ascending order.

% written  05/24/17     F. Buenger
% modified 05/28/18     F. Buenger, redesign, number format according to MATLAB command "format" 
% modified 08/29/18     F. Buenger, INTLAB display formats

global INTLAB_CONST

fmt_old = getformat;
if nargin >= 7 && ~isempty(fmt)
    eval(['format ' fmt]);
end

if nargin < 6 || isempty(print) 
    print = true;
end

n = length(y0); % length of initial conditions = dimension of ODE system.
y0 = intval(y0);
K_ = length(T); % number of time grid points (for initial condition)   

if nargin < 5 || isempty(t)
    t = T;
else
    t = t( and(t >= T(1),t <= T(K_)) ) ;
    t = sort(t);
end

K = length(t); % number of time grid points   

% Initialize result y with zeros.
y = intval(zeros(K,n));

% Save initial conditions y0 to first row of y.
if size(y0,1) > 1
    y(1,:) = y0';
else
    y(1,:) = y0;
end

k_ = 2; % [T(k_-1),T(k_)] is the "current" time interval. 
set_taylor_model = true;

% Initialize evaluation domain.
D.inf = [-ones(1,n),0];  
D.sup = [ones(1,n),0];
i_ = []; 

for k = 1:K
    
    % Adapt domain D for evaluating Taylor model Y(k_-1,:) at time t = t(k).
    % The time domain is set to the  point interval [D.inf(n+1),D.sup(n+1)] = [t(k),t(k)]. 
    D.inf(n+1) = t(k);
    D.sup(n+1) = t(k);  
    
    % Find time interval [T(k_-1),T(k_)] that contains t(k).
    % The corresponding Taylor model vector is Y(k_-1,:) (and also Yr(k_-1,:) if Yr is specified).     
    while t(k) > T(k_)
        k_ = k_ + 1;
        set_taylor_model = true;
    end
    
    if ~isempty(Yr) % right taylor models from preconditioning are provided
        Y_ = substitute(Y(k_-1,:),T(k_-1),t(k)); % First substitute the time variable in the left Taylor model Y(k_-1,:) by t(k)
        Y_ = concatenate(Y_,Yr(k_-1,:));         % concatenation Y_ o Yr(k_-1,:)   
    elseif set_taylor_model
    % Set Taylor model vector Y_ that is responsible for the inclusion on the time interval [T(k_-1),T(k_)].
        Y_ = Y(k_-1,:);
        set_taylor_model = false;
    end
       
    if print
        t_str = strtrim(evalc('disp(t(k))'));
        % Note that t_str might be rounded!
        % Thus, the "display version" t_str of t(k) may not
        % exactly represent the floating-point number t(k).
        % We explicitly did not(!) want to display t(k)
        % artificially as an interval in such a case.
        
        disp(['t = ',t_str]); % Display time point at which an inclusion shall be printed.
    end    
    
    % Display inclusions and their diameters.
    for i = 1:n
        if t(k) == T(1)
            y_ = y0(i);
        else
            y_ = iv2intval(evaltaylormodel(Y_(i),D)); % Evaluate Taylor model Y(k_-1,i) at time t = t(k).
        end
        y(k,i) = y_; % Save inclusion interval to result array.
        if print
            if n > 1 
                i_ = ['_',num2str(i)];
            end            
            switch INTLAB_CONST.INTVAL_DISPLAY
                case 'DisplayMidrad'
                    y_str = midrad(y_(:),[],[]);
                    y_str_diam = []; % In midrad notation the diameter must not bedisplayed separately. 
                otherwise
                    y_str = infsup(y_(:),[],[]);
                    y_str_diam = ['    d([y',i_,']) = ',num2str(diam(y_),'%.2e')];
            end
            y_str_exp = strtrim(y_str.exp);
            y_string = strtrim(y_str.str);
            if ~isempty(y_str_exp)
                y_string = [y_str_exp,' ', y_string];
            end               
            % Display inclusion of component y(i) at time t(k), and also its diameter.
            disp(['    [y',i_,'] = ',y_string,y_str_diam]); 
        end
    end 
    if print
        disp(' ')
    end
end
ts = t(:);
eval(['format ' fmt_old]);

end  % function verifyode_disp