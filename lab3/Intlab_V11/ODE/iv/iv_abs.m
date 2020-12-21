function r =  iv_abs(a)
%IV_ABS  absoulute value for interval like structures
%
%   r = iv_abs(a)

% written  04/28/16     F. Buenger
% modified 08/01/17     F. Buenger  parameter 'includenan' in min-, max-calculations added
% modified 04/13/18     F. Buenger  idx0 for intervals containing 0

r.inf = min(abs(a.inf),abs(a.sup),'includenan');
r.sup = max(abs(a.inf),abs(a.sup),'includenan');

idx0 = and(a.inf < 0, a.sup > 0);
r.inf(idx0) = 0;

end % function iv_abs
