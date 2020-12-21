function res = integral_operator(odefun,y0,t,y,ys,i,recmode )
% INTEGRAL_OPERATOR  Evaluation of the integral operator  
%       
%    res = integral_operator(odefun,y0,t,y,ys,i,recmode )
%
% K(y0,t,y,ys) := ( y0 + quad(odefun(t,y),n+1) ) - ys;  n:=y(1).dim
%
%  odefun : function handle to the ode-function
%      y0 : Taylor model for initial values
%      t  : time  Taylor model
%      y  : y solution Taylor model for the ode y' = odefun(t,y),y(t=t0) = y0
%      ys : Taylor model shift
%      i  : optional parameter. If i is specified, then only the i-th component of K is evaluated
% recmode : record mode, controls the so called "record feature". See [E], section 4.4.3.1, p.112, and code on p.120-127 for description.
%           Three values are distinguished: 
%             0 :  Record feature off (default). This compares to NORMAL_MODE in Riot, see [E], p.121,127. 
%             1 :  Record feature on, write mode. This compares to REMAINDER_REC in Riot, see [E], p.122. 
%             2 :  Record feature on, read mode. This compares to REMAINDER_PLAY in Riot, see [E], p.123,126. 


% written  01/19/16     F. Buenger

global INTLAB_ODE_VARS             % Additional, program internal, global variables that cannot be set by the user. 

n = y(1).dim-1;                    % order of the ode system
if nargin < 7 || isempty(recmode)
    recmode = 0; 
end
if nargin < 6 || isempty(i)
    i = 0;                         % Compute all components of the integral operator.
end
INTLAB_ODE_VARS.RECMODE = recmode; % publish recmode globally.  
if recmode == 1                    % record write mode
  INTLAB_ODE_VARS.RECLIST = [];    % Initialize record list as empty.
elseif recmode == 2 && i == 1      % Record read mode always starts with i == 1, i.e. first component of integral operator.
  INTLAB_ODE_VARS.RECIDX = 1;      % Start reading from the first entry of INTLAB_ODE_VARS.RECLIST.  
end

if i == 0                          % Evaluate all components of the integral operator.
                                   % i == 0 should never happen for recmode == 2 !!!   
    if recmode == 1                % record-write-mode
        res = y;                   % Initialize res.
        for k = 1:n
            z = odefun(t,y,k);     % z := odefun(t,y,k), only the k-th component of odefun is evaluated.
            res(k) = ( y0(k) + integral(z,n+1) ) - ys(k);
        end
        % In record-write-mode the integral operator is evaluated in a loop k = 1,...,n
        % because the subsequent record-read-mode must(!!!) evaluate it that way. 
        % Therefore, the recording is done here in the same fashion so that user "errors"
        % in the implementation of the ode function odefun, namely distinct function evaluation
        % in the vector mode "odefun(t,y)" and in the component mode "odefun(t,y,k)",
        % are excluded.
    else
        z = odefun(t,y); 
        res = ( y0 + integral(z,n+1) ) - ys;
    end
else
    z = odefun(t,y,i);             % Only the i-th component of odefun is evaluated.
    res = ( y0(i) + integral(z,n+1) ) - ys(i);
end

INTLAB_ODE_VARS.RECMODE = 0;       % Switch record feature off

end % function integral_operator