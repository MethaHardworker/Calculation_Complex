function r = get_bounder()
%GET_BOUNDER  returns the specified bounder method for including the image
%             of a Taylor model.
%
%  r = get_bounder()
%
%             0 : 'NAIVE' interval evaluation method
%             1 : 'LDB'   linar dominated bounder

% written  11/25/15     F. Buenger

global INTLAB_ODE_OPTIONS
if ~isfield(INTLAB_ODE_OPTIONS,'bounder')
    opt2glob([]); % define INTLAB_ODE_OPTIONS with default values
end
if strcmp(INTLAB_ODE_OPTIONS.bounder,'NAIVE')
    r = 0;
elseif strcmp(INTLAB_ODE_OPTIONS.bounder,'LDB')
    r = 1;
else
    error('wrong input');
end

end % function get_bounder

