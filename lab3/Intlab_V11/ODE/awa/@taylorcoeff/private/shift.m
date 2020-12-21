function y = shift(y,dy,status,hdk)
%SHIFT  shift of coefficients: 
%
%    y = shift(y,dy,status,hdk)
%
%          k = status + 1
%          y(k+1) = dy(k) * hdk  

% written  08/02/17     F. Buenger
    
k = status + 1;

if nargin < 4 || isempty(hdk)
    h = 1;
    hdk = iv_rdivide(h,k); % "hdk" shortens "h divided by k".
end

for i = 1:numel(y)
    y_ = y(i);
    dy_ = dy(i);
    x.inf = dy_.inf(k);
    x.sup = dy_.sup(k);
    x = iv_times(hdk,x);         
    y_.inf(k+1) = x.inf;  % u.coeff(k+1) = du.coeff(k) * (h/k) 
    y_.sup(k+1) = x.sup;      
    y(i) = y_;
end

end % function shift
