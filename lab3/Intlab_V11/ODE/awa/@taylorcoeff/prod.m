function r = prod(a,dim)
%PROD  Implements  prod(a,dim)  for Taylor coefficients
%
%   r = prod(a,dim)
%
% functionality as Matlab function prod for matrices, parameter dim optional

% written  08/02/17     F. Buenger


[m,n] = size(a);
if nargin < 2
    if m == 1
        dim = 2;
    else
        dim = 1;
    end
end

% trigger = true;    % only for testing !!!
% trigger = false;   % only for testing !!!

% if trigger  

    switch dim
        case 1 % columnwise product
            r = a(1,:); % initialize r with first row of a
            for i = 2:m
                r = r .* a(i,:); % r := r .* i-th row of a
            end
        case 2 % rowwise product
            r = a(:,1); % initialize r with first column of a
            for j = 2:n
                r = r .* a(:,j);  % r := r .* j-th column of a
            end
        otherwise
            error('maximal two dimension for type taylorcoeff');
    end  
    
% else 
%     r = a;
%     switch dim
%         case 1 % columnwise product
%             while m > 1
%                 p = floor(m/2);
%                 rest = (m - 2*p);
%                 idx1 = 1:p;
%                 idx2 = (p+1):(m-rest);
%                 r = a(idx1,:).*a(idx2,:);
%                 if rest
%                     a = [r;a(m,:)];
%                     m = p+1;
%                 else
%                     a = r;
%                     m = p;
%                 end
%             end
%         case 2 % rowwise product
%             while n > 1
%                 p = floor(n/2);
%                 rest = (n - 2*p);
%                 idx1 = 1:p;
%                 idx2 = (p+1):(n-rest);
%                 r = a(:,idx1).*a(:,idx2);
%                 if rest
%                     a = [r,a(:,n)];
%                     n = p+1;
%                 else
%                     a = r;
%                     n = p;
%                 end
%             end
%         otherwise
%             error('maximal two dimension for type taylorcoeff');
%     end  
% end

end % function prod