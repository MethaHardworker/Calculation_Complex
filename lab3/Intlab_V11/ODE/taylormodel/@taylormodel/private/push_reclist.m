function push_reclist(name,value)
% PUSH_RECLIST   appends a record to INTLAB_ODE_VARS.RECLIST.      
%
%   push_reclist(name,value)    


% written  01/19/16     F. Buenger   
global INTLAB_ODE_VARS
  
record.name = name; % Storing a record name only serves for better readability and understanding.
record.value = value; 
INTLAB_ODE_VARS.RECLIST =  [INTLAB_ODE_VARS.RECLIST; record];

end % function push_reclist