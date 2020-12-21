function r = subtract_constant_term(a,i0)
%SUBTRACT_CONSTANT_TERM   deletes (sets to zero) the constant term c0 := a.coefficient(i0) of a Taylor model.
%                         This is a short/fast alternative for computing a-c0.
%
%   r = subtract_constant_term(a,i0)


% written  09/09/15     F. Buenger
% modified 11/24/15     F. Buenger  matrix input
% modified 12/18/15     F. Buenger  switch between verified/non-verified Taylor model arithmetic 
% modified 01/20/16     F. Buenger  record feature
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures

global INTLAB_ODE_VARS

ODEMODE = INTLAB_ODE_VARS.ODEMODE;
RECMODE = INTLAB_ODE_VARS.RECMODE;

S_a = size(a);

if RECMODE ~= 2
    if nargin == 1
        [c0,i0] = get_constant_term(a);
    else
        S_i0 = size(i0);
        if (length(S_a) > 2) || any(S_i0 ~= S_a)
            error('Input parameters are not consistent.')
        end
    end
end

r = a; % Initialize result r with a. (Just preallocation of memory)
% Remark: In RECMODE == 2, the function could immediately return
%         after this line but for clearity all components r(i,j).cmp
%         with '.cmp' distinct from '.interval' are cleared.

for i = 1:S_a(1)
    for j = 1:S_a(2)
        r_ = r(i,j); 
        if RECMODE ~= 2
            i0_ = i0(i,j);
            if i0_ ~= 0                        % i0 == 0 means that the coefficient of the constant term of r_ = a(i,j) is already zero so that nothing has to be done.
                if length(r_.coefficient) == 1 % If 'r_' only contains the constant monomial,then this is not deleted to avoid empty Taylor models.
                    r_.coefficient(1) = 0;     % The coefficient of the constant term is simply set to zero.
                else
                    idx = [1:i0-1,i0+1:length(r_.coefficient)];
                    r_.monomial = r_.monomial(idx,:);     % Delete zero- monomial of constant term.
                    r_.coefficient = r_.coefficient(idx); % Delete constant term.
                end
                
                if ODEMODE == 1 % verified mode.
                    r_.image = image(r_);
                    % r.image = r.image - intval(c0); % faster but possibly broader image interval
                    % Note that the subtraction a-c0 does not change the error interval a.interval, see [E], p.105, (4.6).
                else % non-verified mode
                    r_.interval = INTLAB_ODE_VARS.EMPTYIV;
                    r_.image = INTLAB_ODE_VARS.EMPTYIV;
                end
            end
        else % RECMODE == 2, record read mode, only the error interval component is computed / of interest !!!
            r_interval = r_.interval;             % Save error interval before clearing all other components
            r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL; % Unnecessarily, just for clearity all components are cleared in RECMODE == 2.
            r_.interval = r_interval;             % Only the interval component is of interest which remains unchanged.
        end
        r(i,j) = r_;
    end
end
       
end % function subtract_constant_term

