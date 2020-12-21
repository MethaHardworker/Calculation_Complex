function C = sconv(A,B)
% SCONV  convolution of two column vectors A and B. 
%        The result is always sparse while A and B can be sparse or full.
%        This is in contrast to the MATLAB functions conv and conv2  which 
%        do not accept sparse input.
%
% C = sconv(A,B)

% written  05/09/17     F. Buenger

m = length(A);
n = length(B);

%[i,~,a] = find(A);      % The "~"-notation is not downward-compatible for older MATLAB versions.
[i,dummy,a] = find(A);   % Therefore a "dummy"-variable is stated which is not used in the sequel.  
if isempty(i) % A is sparse zero
    C = A;    % => return sparse zero
    return;
end

%[j,~,b] = find(B);      % See above.
[j,dummy,b] = find(B);
if isempty(j) % B is sparse zero
    C = B;    % => return sparse zero
    return;
end

K = i+j'-1;
C = a*b.'; % If rounding is upwards (downwards), then C is a componentwise upper (lower) bound for a*b'.

if numel(K) < 1E6
    C = sparse(K(:),1,C(:),m+n-1,1); % Entries of C(:) with same entries in K(:) are automatically summed up with respect to the actual rounding mode.
                                     % Thus, if rounding is upwards (downwards), then C is a verified componentwise upper (lower) bound for the convolution of A and B.
else                                 
    [K,I] = sort(K(:)); % The performance of "sparse()" is increased by previous sorting of the indices in K(:).   
    C = sparse(K,1,C(I),m+n-1,1); 
    %C = accumarray(K,C(I),[m+n-1 1],[],[],true);    
end

end % function sconv

