function r =  iv_repmat(a,m,n)
%IV_REPMAT  repmat for interval structures
%
%   c = iv_repmat(a,m,n)

% written  02/10/16     F. Buenger

r = struct('inf',repmat(a.inf,m,n),'sup',repmat(a.sup,m,n));

end % function iv_repmat
