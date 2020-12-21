function r =  iv_intersect(a,b)
%IV_INTERSECT  intersection of interval structures with 
%              maximally two dimensions
%
%   s = iv_intersect(a,b)

% written  05/04/16     F. Buenger

u = max(a.inf,b.inf,'includenan');
v = min(a.sup,b.sup,'includenan');

if any(any(u > v))
  r = []; % empty intersection
else
  %r = struct('inf',u,'sup',v); 
  r.inf = u;
  r.sup = v;
end  

end % function iv_intersect
