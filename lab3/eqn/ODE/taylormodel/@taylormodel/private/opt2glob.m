function opt2glob( options )
%OPT2GLOB    transfer parameters from options to the global structure INTLAB_ODE_OPTIONS
%            and set all parameters which are not specified in options by default values. 
%  
%    opt2glob( options )

% written  11/06/15     F. Buenger

global INTLAB_ODE_OPTIONS  

order_default = 12; % multivariate polynomials of degree at most 12
sparsity_tol_default = 1E-20;
loc_err_tol_default = 1E-10;
h_min_default = 1E-4;
h0_default = 0;             % h0 = 0 means that the initial stepsize will be chosen automatically.
bounder_default = 'NAIVE' ; % default := 'NAIVE' : naive interval evaluation
                            % At the moment the only alternative is 'LDB' : linear dominated bounder.
shrinkwrap_default = true ; % default := true, this means shrink wrapping on.    
precondition_default = 0 ;  % default := 0, this means that preconditioning is switched off.                             
blunting_default = false;   % default := false, this means that blunting is switched off.

INTLAB_ODE_OPTIONS.order = verifyodeget(options,'order',order_default);   
INTLAB_ODE_OPTIONS.sparsity_tol = verifyodeget(options,'sparsity_tol',sparsity_tol_default); 
INTLAB_ODE_OPTIONS.loc_err_tol = verifyodeget(options,'loc_err_tol',loc_err_tol_default); 
INTLAB_ODE_OPTIONS.h_min = verifyodeget(options,'h_min',h_min_default);  
INTLAB_ODE_OPTIONS.h0 = verifyodeget(options,'h0',h0_default);  
INTLAB_ODE_OPTIONS.bounder = verifyodeget(options,'bounder',bounder_default); 
INTLAB_ODE_OPTIONS.shrinkwrap = verifyodeget(options,'shrinkwrap',shrinkwrap_default); 
INTLAB_ODE_OPTIONS.precondition = verifyodeget(options,'precondition',precondition_default); 
INTLAB_ODE_OPTIONS.blunting = verifyodeget(options,'blunting',blunting_default); 

if INTLAB_ODE_OPTIONS.h_min > INTLAB_ODE_OPTIONS.h0 && INTLAB_ODE_OPTIONS.h0 > 0
    error('h0 must be greater or equal h_min.')
end

end % function opt2glob

