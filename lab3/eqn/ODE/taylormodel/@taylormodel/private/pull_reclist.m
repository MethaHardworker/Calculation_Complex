function res =  pull_reclist(name)
% PULL_RECLIST   returns the value component of the record 
%               INTLAB_ODE_VARS.RECLIST(INTLAB_ODE_VARS.RECIDX)       
%               and increases INTLAB_ODE_VARS.RECIDX by one.
%
%   res =  pull_reclist(name)

% written  01/19/16     F. Buenger   
global INTLAB_ODE_VARS

record =  INTLAB_ODE_VARS.RECLIST(INTLAB_ODE_VARS.RECIDX);
if ~strcmp(name,record.name)
    error('record error'); % soft check
end
res = record.value;
INTLAB_ODE_VARS.RECIDX = INTLAB_ODE_VARS.RECIDX + 1;

end % function pull_reclist