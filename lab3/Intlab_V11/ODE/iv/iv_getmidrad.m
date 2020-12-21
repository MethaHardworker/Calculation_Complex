function [c_mid,c_rad] = iv_getmidrad(a)
%IV_GETMIDRAD  return aproximate midpoint "c_mid" and radius "c_rad" of interval structure "a"
%              such that mathematically c_mid-c_rad <= a.inf and a.sup <= c_mid+c_rad
%              componentwise
%
%   [c_mid,c_rad] = iv_getmidrad(a,b)

% written   02/10/16     F. Buenger
% modified  08/01/17     F. Buenger  parameter 'includenan' in min-, max-calculations added

if isfloat(a)
    c_mid = a;
    c_rad = 0*a;
else    
    e = 1e-30;    
    if 1+e > 1      % fast check for rounding upwards
        rndold = 1;
    else
        rndold = getround;
        setround(1) % rounding upwards
    end
    
    c_mid = (a.inf+a.sup)/2;  % approximate midpoint of a.                               
    
    % if any(isnan(c_mid(:))) || any(isinf(c_mid(:)))
    %     error('NAN/INF occured.')
    % end
     
    c_rad = c_mid - a.inf; % Since rounding is upwards, c_mid > c_mid_exact := (a.inf+a.sup)/2. 
                           % Thus, c_mid - a.inf >= a.sup - c.mid
    
    if rndold ~= 1
        setround(rndold)
    end
    
end

end % function iv_getmidrad
