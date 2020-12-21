function a = log10(a)
%LOG10  decimal logarithm for Taylor coefficients
%
%   r = log10(a)

% written  08/02/17     F. Buenger

global INTLAB_CONST

INTLAB_STDFCTS_LOG10_ = INTLAB_CONST.STDFCTS_LOG10_;

c.inf = INTLAB_STDFCTS_LOG10_.INF;
c.sup = INTLAB_STDFCTS_LOG10_.SUP;

a = log(a) .* c;

end % function log10
