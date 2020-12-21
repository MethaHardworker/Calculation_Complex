function c =  iv_mtimes(a,b)
%IV_MTIMES  interval structure matrix multiplication  a * b
%
%   c = iv_mtimes(a,b)

% written  02/09/16     F. Buenger

e = 1e-30;

if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

%trigger = false; % only for testing !!!
%trigger = true;  % only for testing !!!

% if trigger
%     
%     if isfloat(a)
%         m = size(a,1);
%         if isfloat(b)  % case "float * float"
%             n = size(b,2);
%             c.inf = zeros(m,n);
%             c.sup = c.inf;
%             for  i = 1:m
%                 for j=1:n
%                     d = iv_dotprod(a(i,:)',b(:,j));
%                     c.inf(i,j) = d.inf;
%                     c.sup(i,j) = d.sup;
%                 end
%             end
%         else          % case "float * iv"
%             n = size(b.inf,2);
%             c.inf = zeros(m,n);
%             c.sup = c.inf;
%             for  i = 1:m
%                 for j=1:n
%                     b_j.inf = b.inf(:,j);
%                     b_j.sup = b.sup(:,j);
%                     d = iv_dotprod(a(i,:)',b_j);
%                     c.inf(i,j) = d.inf;
%                     c.sup(i,j) = d.sup;
%                 end
%             end
%         end
%     elseif isfloat(b)  % case "iv * float"
%         m = size(a.inf,1);
%         n = size(b,2);
%         c.inf = zeros(m,n);
%         c.sup = c.inf;
%         for  i = 1:m
%             a_i.inf = a.inf(i,:)';
%             a_i.sup = a.sup(i,:)';
%             for j=1:n
%                 d = iv_dotprod(a_i,b(:,j));
%                 c.inf(i,j) = d.inf;
%                 c.sup(i,j) = d.sup;
%             end
%         end
%     else               % case "iv * iv"
%         m = size(a.inf,1);
%         n = size(b.inf,2);
%         c.inf = zeros(m,n);
%         c.sup = c.inf;
%         for  i = 1:m
%             a_i.inf = a.inf(i,:)';
%             a_i.sup = a.sup(i,:)';
%             for j=1:n
%                 b_j.inf = b.inf(:,j);
%                 b_j.sup = b.sup(:,j);
%                 d = iv_dotprod(a_i,b_j);
%                 c.inf(i,j) = d.inf;
%                 c.sup(i,j) = d.sup;
%             end
%         end
%     end
% 
%     
% else
    
    % a*b = reshape ( sum( ( [a(1,:),...,a(m,:)]' .* [b;b;...;b] )' ) , m,n)
    
    if isfloat(a)
        m = size(a,1);
        l = size(a,2);
        a_ = a';
        a_ = a_(:);
        if isfloat(b)  % case "float * float"
            b_ = repmat(b,m,1);
            n = size(b,2);
        else           % case "float * iv"
            b_.inf = repmat(b.inf,m,1);
            b_.sup = repmat(b.sup,m,1);
            n = size(b.inf,2);
        end
    else
        m = size(a.inf,1);
        l = size(a.inf,2);
        a_.inf = a.inf';
        a_.inf = a_.inf(:);
        a_.sup = a.sup';
        a_.sup = a_.sup(:);
        if isfloat(b)  % case "iv * float"
            b_ = repmat(b,m,1);
            n = size(b,2);
        else               % case "iv * iv"
            b_.inf = repmat(b.inf,m,1);
            b_.sup = repmat(b.sup,m,1);
            n = size(b.inf,2);
        end
    end
    c = iv_times(a_,b_);
    c.inf = reshape(c.inf,l,m*n);
    c.sup = reshape(c.sup,l,m*n);
    c = iv_sum(c);
    c.inf = reshape(c.inf,m,n);
    c.sup = reshape(c.sup,m,n);

% end

if rndold ~= 1 
    setround(rndold)
end

end % function iv_mtimes
