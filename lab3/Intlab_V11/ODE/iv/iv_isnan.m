function r = iv_isnan(a)
%IV_ISNAN  Array of 1's for NaN components
%
%   r = iv_isnan(a)

% written  05/02/16     F. Buenger

r = isnan(a.inf) | isnan(a.sup) ;

end % function iv_isnan
