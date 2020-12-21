function a = set_order(a,order)
%SET_ORDER     sets a.order to "order" and deletes all monomials and 
%              corresponding coefficients of higher order.
%              The components a.interval and a.image are set to the
%              empty intervals INTLAB_ODE_VARS.EMPTYIV. 
%
%   a = set_order(a,order)

% written  02/05/16     F. Buenger

global INTLAB_ODE_VARS

for i = 1:numel(a)
    a_ = a(i);
    a_.order = order;
    idx = (sum(a_.monomial,2) > order); % Determine indices of monomials of higher order.  
    a_.monomial(idx,:) = [];            % Delete monomials of higher order.
    a_.coefficient(idx) = [];           % Delete corresponding coefficients of higher order.
    a_.interval = INTLAB_ODE_VARS.EMPTYIV;
    a_.image = INTLAB_ODE_VARS.EMPTYIV;    
    a(i) = a_;
end
end % function set_order