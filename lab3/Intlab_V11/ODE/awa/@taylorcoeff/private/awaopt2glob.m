function awaopt2glob(options)
%AWAOPT2GLOB  transfer parameters from options to the global structure AWA_OPTIONS
%             and set all parameters which are not specified in options by default values. 
%  
%    awaopt2glob(options)

% written  07/22/17     F. Buenger

global INTLAB_AWA_OPTIONS

order_default = 10;        % default order of Taylor expansion
h0_default = 0;            % AWA determines initial step size automatically
h_min_default = 1E-5;
EvalMeth_default = 4;
AbsTol_default = 1E-16;   
RelTol_default = 1E-16;
 
INTLAB_AWA_OPTIONS.h0 = awaget(options,'h0',h0_default);  
INTLAB_AWA_OPTIONS.h_min = awaget(options,'h_min',h_min_default);  
INTLAB_AWA_OPTIONS.order = awaget(options,'order',order_default); 
INTLAB_AWA_OPTIONS.EvalMeth = awaget(options,'EvalMeth',EvalMeth_default);
INTLAB_AWA_OPTIONS.AbsTol = awaget(options,'AbsTol',AbsTol_default); 
INTLAB_AWA_OPTIONS.RelTol = awaget(options,'RelTol',RelTol_default); 

end % function awaopt2glob

