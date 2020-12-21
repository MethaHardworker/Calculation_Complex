function a = log10(a)
%LOG10  Taylor model decimal logarithm log10(a)
%
%   a = log10(a)

% written  12/15/15     F. Buenger

global INTLAB_CONST

INTLAB_STDFCTS_LOG10_ = INTLAB_CONST.STDFCTS_LOG10_;

c.inf = INTLAB_STDFCTS_LOG10_.INF;
c.sup = INTLAB_STDFCTS_LOG10_.SUP;

a = log(a) .* c;

end % function log10
