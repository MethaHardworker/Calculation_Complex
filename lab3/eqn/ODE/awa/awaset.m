function options = awaset(varargin)
%AWASET  sets option structure for verified ode-solver awa,
%        analogous to MATLAB function odeset.
%        options = awaset('name1',value1,'name2',value2,...)
% 
%   options = awaset(varargin)
%
% The following options can be set:
% 
%    - 'order'    order of Taylor expansion up to which the solution is 
%                 computed in each time step
%    - 'h0'       initial step size
%    - 'h_min'    minimum step size
%    - 'EvalMeth' evaluation method (corresponds to AWA variable E_ART)
%                 The following methods can be chosen:
%                 0: (interval vector). After each integration step the 
%                    solution is enclosed by an interval vector (axis 
%                    parallel box) which is used as initial value set for 
%                    the next integration step.
%                 1: (parallelepiped). In each integration step, the 
%                    solution is enclosed into a parallelepiped which is used 
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

% written  07/22/17     F. Buenger
% modified 04/26/17     F. Buenger   lower case field names 'evalmeth','abstol','reltol' are also accepted.

if odd(nargin) || nargin == 0 || ~ischar(cell2mat(varargin(1:2:nargin))) 
    error('invalid call of awaset.')
end

for k=1:2:nargin
    name  = varargin{k}; % option name
    value = varargin{k+1}; % corresponding option value
    
    if strcmp(name,'order') %
        if  isa(value,'double') && isreal(value) && numel(value) == 1 && value > 0 && value == round(value)  % check that value is a positive integer
            if value < 2
                error('The Taylor expension order must be at least 2.')
            end
            options.order = value;
        else
            error(['invalid call of awaset: The option "order" must be a positive integer.'])
        end                        
    elseif strcmp(name,'h0')
        if isa(value,'double') && isreal(value) && numel(value) == 1 && value >= 0  % check that value is a positive real number > 0
            options.h0 = value;
        else
            error('invalid call of awaset: The option "h0" must be a positive real number.')
        end
    elseif strcmp(name,'h_min')
        if isa(value,'double') && isreal(value) && numel(value) == 1 && value > 0  % check that value is a positive real number > 0
            options.h_min = value;
        else
            error('invalid call of awaset: The option "h_min" must be a positive real number.')
        end
    elseif strcmp(name,'EvalMeth') || strcmp(name,'evalmeth')
        if (value == 0 || value == 1 || value == 2 || value == 3 || value == 4)
            options.EvalMeth = value;
        else
            error('invalid call of awaset: The option "EvalMeth" must be 0, 1, 2, 3, or 4.')
        end
    elseif strcmp(name,'AbsTol') || strcmp(name,'abstol')
        if isa(value,'double') && isreal(value) && numel(value) == 1 && value > 0
            options.AbsTol = value;
        else
            error('invalid call of awaset: The option "AbsTol" must be a small positive number.')
        end
    elseif strcmp(name,'RelTol') || strcmp(name,'reltol')
        if isa(value,'double') && isreal(value) && numel(value) == 1 && value > 0
            options.RelTol = value;
        else
            error('invalid call of awaset: The option "RelTol" must be a small positive number.')
        end
    else
        error(['invalid call of awaset: unknown option ',name])
    end
end

% if isfield(options,'h0') && isfield(options,'h_min') && options.h_min > options.h0
%     error('h0 must be greater or equal h_min.')
% end

end % function awaset