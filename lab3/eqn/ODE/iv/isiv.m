function res =  isiv(a)
%ISIV  interval structure check
%
%   c = isiv(a)

% written  02/09/16     F. Buenger

res = false;

%if isstruct(a) && size(fieldnames(a),1) == 2 && isfield(a,'inf') && isfield(a,'sup')    
%if all(isfield(a,{'inf','sup'}))      
if isstruct(a) && isfield(a,'inf') && isfield(a,'sup')      
   s_inf = size(a.inf);
   s_sup = size(a.sup);
   if length(s_inf) > 2 || length(s_sup) > 2
       error('maximally two dimensions for interval structures')
   end
   if s_inf(1) == s_sup(1) && s_inf(2) == s_sup(2) 
      b = [a.inf, a.sup];
      if isfloat(b) && isreal(b)
          res = true;
      end
   end
end

end % function isiv
