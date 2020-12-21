function r = get_sparsity_tol()
%GET_SPARSITY_TOL    sparcity tolerance for coefficients of Taylor models.  
%       
%    r = get_sparsity_tol()
%
%                    Taylor coefficients c with |c| < sparsity_tol are set
%                    to zero and the corresponding error is transferred to
%                    the error interval of the Taylor model.

% written  26/08/15     F. Buenger

global INTLAB_ODE_OPTIONS

if ~isfield(INTLAB_ODE_OPTIONS,'sparsity_tol') 
  opt2glob([]); % define INTLAB_ODE_OPTIONS with default values  
end
r = INTLAB_ODE_OPTIONS.sparsity_tol;     

end % function get_sparsity_tol

