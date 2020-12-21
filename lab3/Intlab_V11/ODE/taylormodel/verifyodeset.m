function options = verifyodeset(varargin)
%verifyodeset      sets option structure for verified ode-solver verifyode,
%                  analogous to MATLAB function odeset.
%                  options = verifyodeset('name1',value1,'name2',value2,...)
% 
%   options = verifyodeset(varargin)
%
% The following options can be set:
% 
%     - 'order'         degree of Taylor polynomials 
%     - 'sparsity_tol'  threshold amount for coefficients c of a multivariate polynomial. 
%                       If |c|<sparsity_tol then the coefficient is removed from the polynomial and the (small) error 
%                       caused by this is transferred to the error interval of the Taylor model.                        
%     - 'loc_err_tol'   tolerance bound for local error
%     - 'h_min'         minimum time step size 
%     - 'h0'            initial step size
%     - 'bounder'       LDB, NAIVE
%     - 'shrinkwrap'    flag for 'shrink wrapping', 1 = on, 0 = off
%     - 'precondition'  flag for preconditioning, 0 = off, 1: QR method, 2: parallelepiped method

% written  11/05/15     F. Buenger

if odd(nargin) || nargin == 0 || ~ischar(cell2mat(varargin(1:2:nargin)))
    error('invalid call of verifyodeset.')
end

for k=1:2:nargin
    name  = varargin{k};   % option name
    value = varargin{k+1}; % corresponding option value
    
    if strcmp(name,'order')
        if isa(value,'double') && isreal(value) && numel(value) == 1 && value > 0 && value == round(value) % check that value is a positive integer
            options.order = value;
        else
            error(['invalid call of verifyodeset: The option "order" must be a positive integer.'])
        end
    elseif strcmp(name,'sparsity_tol')
        max_sparsity_tol =  1E-5;  % maximum sparsity tolerance 1E-5. Feel free to change that!
        if  isa(value,'double') && isreal(value) && numel(value) == 1 && value > 0 && value <= max_sparsity_tol % check that value is a positive real number <1
            options.sparsity_tol = value;
        else
            error(['invalid call of verifyodeset: The option "sparsity_tol" must be a positive real number <= ',num2str(max_sparsity_tol),'.'])
        end
    elseif strcmp(name,'loc_err_tol')
        max_loc_err_tol =  0.1; % maximum tolerance = 10% . Feel free to change that!
        if  isa(value,'double') && isreal(value) && numel(value) == 1 && value > 0 && value < max_loc_err_tol % check that value is a positive real number < 1
            options.loc_err_tol = value;
        else
            error(['invalid call of verifyodeset: The option "loc_err_tol" must be a positive real number less than ',num2str(max_loc_err_tol),'.'])
        end
    elseif strcmp(name,'h_min')
        if isa(value,'double') && isreal(value) && numel(value) == 1 && value > 0 % check that value is a positive real number
            options.h_min = value;
        else
            error('invalid call of verifyodeset: The option "h_min" must be a positive real number.')
        end
    elseif strcmp(name,'h0')
        if isa(value,'double') && isreal(value) && numel(value) == 1 && value >= 0 % check that value is a positive real number >= 0.
            options.h0 = value;                                                      % h0 = 0 means that the initial stepsize will be chosen automatically.
        else
            error('invalid call of verifyodeset: The option "h0" must be a nonnegative real number >= h_min .')
        end
    elseif strcmp(name,'bounder')
        if strcmp(value,'LDB') || strcmp(value,'NAIVE')
            options.bounder = value;
        else
            error('invalid call of verifyodeset: invalid bounder');
        end
    elseif strcmp(name,'shrinkwrap')
        options.shrinkwrap = (value == true);
    elseif strcmp(name,'precondition')
        PREC_OFF = 0; % no preconditioning
        PREC_QR  = 1; % QR preconditioning
        PREC_PE  = 2; % parallelepiped preconditioning
        
        if value == PREC_OFF || value == PREC_QR || value == PREC_PE
            options.precondition = value;
        else
            error(['invalid call of verifyodeset.' newline ...
                '  The option "precondition" must be one of the following numbers:' newline ...
                '    0: off' newline ...
                '    1: QR preconditioning.' newline ...
                '    2: parallelepiped preconditioning']);
        end
    elseif strcmp(name,'blunting')
        options.blunting = (value == true);
    else
        error(['invalid call of verifyodeset: unknown option ',name])
    end
end

if isfield(options,'h0') && isfield(options,'h_min') && options.h_min > options.h0 && options.h0 > 0
    error('h0 must be greater or equal h_min.')
end

end % function verifyodeset