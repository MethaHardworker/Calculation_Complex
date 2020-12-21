function r = tc2tcoeff(a,dummy)
%TC2TCOEFF  converts tc-structure to generalized Taylor coefficient
%
%    r = tc2tcoeff(a)

% written  12/22/17     F. Buenger

if isempty(a)
    r = [];
else
    [m,n,k] = size(a.inf);
    if nargin < 2
        dummy.inf = zeros(k,1);
        dummy.sup = dummy.inf;
        dummy = taylorcoeffinit(dummy);
    end
    r(1:m,1:n) = dummy;  % just storage preallocation
        
    %trigger = true;
    trigger = false;
    
    if trigger

%         x = reshape(permute(a.inf,[3,1,2]),k,m*n);
%         x = num2cell(x,1);
%         [r(:).inf] = x{:};
%         x = reshape(permute(a.sup,[3,1,2]),k,m*n);
%         x = num2cell(x,1);
%         [r(:).sup] = x{:};      
        
        % condensed version of the previous lines in comments 
        x = num2cell(reshape(permute([a.inf,a.sup],[3,1,2]),k,2*m*n),1);  % num2cell is - at the moment - not a Matlab Built-in function and therefore quiet slow. 
                                                                          % Maybe Matlab will provide something faster in a future release. 
        [r(:).inf, r(:).sup] = x{:};  

    else
        
        a_inf = reshape(permute(a.inf,[3,1,2]),k,m*n);
        a_sup = reshape(permute(a.sup,[3,1,2]),k,m*n);
        for i = 1:m*n
            r(i).inf = a_inf(:,i);
            r(i).sup = a_sup(:,i);
        end
        
    end
end

end % function tc2tcoeff
