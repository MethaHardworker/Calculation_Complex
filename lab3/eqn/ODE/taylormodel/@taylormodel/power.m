function r = power(a,b)
%POWER  Taylor model (elementwise) power a .^ b
%
%   r = power(a,b)
%
% The main purpose of this function is to compute a.^b for Taylor model matrix a and integer matrix b.
% The following situations can be distinguished:  
%
%  1) (mxn)-matrix .^ scalar          ->  r(i,j) = a(i,j)^b                       
%  2) scalar .^ (mxn)-matrix          ->  r(i,j) = a^b(i,j)                               
%  3) (mxn)-matrix .^ (mxn)-matrix    ->  r(i,j) = a(i,j)^b(i,j)                                 
%  4) (1xn)-vector .^ (mxn)-matrix    ->  r(i,j) = a(j)^b(i,j)       
%  5) (mx1)-vector .^ (mxn)-matrix    ->  r(i,j) = a(i)^b(i,j)      
%  6) (mxn)-matrix .^ (1xn)-vector    ->  r(i,j) = a(i,j)^b(j)         
%  7) (mxn)-matrix .^ (mx1)-vector    ->  r(i,j) = a(i,j)^b(i)            
%  8) (1xn)-vector .^ (mx1)-vector    ->  r(i,j) = a(j)^b(i)                      
%  9) (mx1)-vector .^ (1xn)-vector    ->  r(i,j) = a(i)^b(j)            
%
% i=1,...,m, j=1,...,n. 
% In particular, the output r is an (mxn)-Taylor model matrix in all cases. 
%
% If b is not an integer matrix, then r = exp( b .* log(a) ) is returned.

% written  11/05/15     F. Buenger  scalar input
% modified 11/23/15     F. Buenger  matrix input (componentwise exponentiation)
% modified 12/18/15     F. Buenger  switch between verified/non-verified Taylor model arithmetic 
% modified 01/20/16     F. Buenger  record feature
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures 
% modified 04/25/17     F. Buenger  a.^b for integer matrix b 
% modified 06/01/18     F. Buenger  special implementation for single interval or non-integer float b 

global INTLAB_ODE_VARS
ODEMODE = INTLAB_ODE_VARS.ODEMODE;

TRIGGER = true;
%TRIGGER = false;

isiv_a = isiv(a);
isiv_b = isiv(b);

if isiv_a
    S_a = size(a.inf);
else
    S_a = size(a);    
end
if isiv_b
    S_b = size(b.inf);
else
    S_b = size(b);    
end


if max(S_a) == 0 || max(S_b) == 0 || length(S_a) > 2 || length(S_b) > 2
    error('wrong input')    
end
% a and b have compatible dimension sizes, if for each dimension d=1,2
% either S_a(d) == 1 or S_b(d) == 1 or S_a(d) == S_b(d).
if ~all( S_a == 1 | S_b == 1 | S_a == S_b)
    error('incompatible dimension sizes')        
end

intpower = isfloat(b) && all(b(:) == round(b(:))) ; % indicates that all entries of b are integer powers

if intpower
    m = max(S_a(1),S_b(1));
    n = max(S_a(2),S_b(2));
    c = pcase(S_a,S_b);    
    r(1:m,1:n) = a(1); % just cheap storage preallocation for result r    
    switch c
        case 1 %  1) (mxn)-matrix .^ scalar          ->  r(i,j) = a(i,j)^b     
            for i = 1:m
                for j = 1:n
                   r(i,j) = vec_power(a(i,j),b);  
                end
            end
        case 2 %  2) scalar .^ (mxn)-matrix          ->  r(i,j) = a^b(i,j)  
            r = reshape(vec_power(a,b(:)),m,n);
        case 3 %  3) (mxn)-matrix .^ (mxn)-matrix    ->  r(i,j) = a(i,j)^b(i,j)  
            for i = 1:m
                for j = 1:n
                   r(i,j) = vec_power(a(i,j),b(i,j));  
                end
            end        
        case 4 %  4) (1xn)-vector .^ (mxn)-matrix    ->  r(i,j) = a(j)^b(i,j)  
           for j = 1:n
             r(:,j) =  vec_power(a(j),b(:,j)); 
           end
        case 5 %  5) (mx1)-vector .^ (mxn)-matrix    ->  r(i,j) = a(i)^b(i,j)   
           for i = 1:m               
             r(i,:) =  vec_power(a(i),b(i,:)); 
           end        
        case 6 %  6) (mxn)-matrix .^ (1xn)-vector    ->  r(i,j) = a(i,j)^b(j)  
           for i = 1:m     
             for j = 1:n
                r(i,j) =  vec_power(a(i,j),b(j)); 
             end
           end                    
        case 7 %  7) (mxn)-matrix .^ (mx1)-vector    ->  r(i,j) = a(i,j)^b(i)   
           for i = 1:m     
             for j = 1:n
                r(i,j) =  vec_power(a(i,j),b(i)); 
             end
           end                    
        case 8 %  8) (1xn)-vector .^ (mx1)-vector    ->  r(i,j) = a(j)^b(i)   
             for j = 1:n
                r(:,j) =  vec_power(a(j),b); 
             end
        case 9 %  9) (mx1)-vector .^ (1xn)-vector    ->  r(i,j) = a(i)^b(j) 
             for i = 1:m
                r(i,:) =  vec_power(a(i),b); 
             end
    end        
else 
    if isfloat(b) || isa(b,'intval') 
        if TRIGGER
            r = pow(a,b,S_a,S_b);        
        else
            r = pow_(a,b,S_a,S_b);        
        end
    elseif isiv_b
        if TRIGGER
            r = pow(a,iv2intval(b),S_a,S_b);
        else
            r = pow_(a,iv2intval(b),S_a,S_b);
        end
    elseif isiv_a
        r = exp( b .* log(iv2intval(a)) );
    elseif isfloat(a) && ODEMODE == 1
        r = exp( b .* log(intval(a)) );
    else
        r = exp( b .* log(a) );
    end
end 

end % function power

%--------------------------------------------------------------------------
% function r = vec_power(a,b)
%--------------------------------------------------------------------------

function r = vec_power(a,b)

global INTLAB_ODE_VARS
RECMODE = INTLAB_ODE_VARS.RECMODE;
% computes a^b for a single Taylor model a and integer vector b.
if isempty(b)
    r = [];
    return;
end
S_b = size(b);
r(1:S_b(1),1:S_b(2)) = a; % storage preallocation for result array

idx_0 = (b == 0);
if any(idx_0)
    r_ = a;
    if RECMODE ~= 2
        r_.monomial = zeros(1,r_.dim);        % monomial for constant term
        r_.coefficient = 1;                   % a_^0 = 1
        r_.interval = INTLAB_ODE_VARS.ZEROIV; % no error
        r_.image.inf = 1;                     % The image is the point interval [1,1].
        r_.image.sup = 1;
    else % RECMODE == 2
        r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL; % In particular r_.interval = intval(0)
    end
    r(idx_0) = r_;
end

b_max = max(b); % storage preallocation for array of intermediate values
p(1:b_max) = a; % Until now a = a^1 is the only known power.
e = 1;
bs = unique(sort(b(b>1))); % Only powers > 1 must still be computed.
while ~isempty(bs)
    bs_ = bs(1);
    c = floor(bs_/2);
    if bs_ == 2*c           % If bs_ is even, then a^bs_ = (a^c)^2 is a square and has nonnegative image.
        i = find(e == c,1); % Other possible factorization a^bs_ = a^x * a^y, with x ~= y are discarded 
                            % since such a product does not ensure the positivity of the resulting image.
                            % Also the error interval is calculated less sharp.
        if ~isempty(i)
           r_ =  sqr(p(c)); 
           p(bs_) = r_;
        end
    %end                                % Uncomment these two lines and comment the subsequent "else"-line if 
    %if  bs_ ~= 2*c || isempty(i)       % a^bs_ = a^x * a^y, x ~= y, shall also be tried for even bs_ .   
    else
        %[X,Y] = meshgrid(e);
        %E = X+Y;
        E = e'+e;
        [i,j] = find(E == bs_,1); % The odd power a^bs_ can be computed as a^p(e(i))*a^p(e(j)) by only one multiplication
                                  % since both factors a^p(i)*a^p(j) are already known.
        if ~isempty(i)
            r_ =  p(e(i)).*p(e(j));
            p(bs_) = r_;
        end
    end
    
    if isempty(i)
        r_ = a;
        a_ = a;
        b_ = bs_;
        % Check if b_ is even to ensure that the result is nonnegative.
        if b_ == 2*floor(b_/2)
            b_ = b_/2;
            b_is_even = 1;
        else
            b_is_even = 0;
        end
        b_ = b_ - 1;
        u = 1; % current exponent for which r_ = a^u holds true
        v = 1; % current exponent for which a_ = a^v holds true
        while b_ > 0
            if mod(b_,2) == 1
                u = u + v;
                idx = (e == u); % Check if a^u was already computed
                if any(idx)
                    r_ = p(u);
                else
                    r_ = r_ .* a_;
                    p(u) = r_;
                    e = [e u];
                    r(b == u) = r_;
                    bs(bs == u) = [];
                end            
            end
            b_ = floor(b_/2);
            if b_ ~= 0
                v = 2*v;
                idx = (e == v); % Check if a^u was already computed
                if any(idx)
                    a_ = p(v);
                else
                    a_ = sqr(a_);
                    p(v) = a_;
                    e = [e v];
                    r(b == v) = a_;
                    bs(bs == v) = [];
                end
            end
        end
        if b_is_even
            r_ = sqr(r_);
            u = 2*u;
            p(u) = r_;
            e = [e u];
        end
    end % isempty(i)
    
    e = [e bs_];
    r(b == bs_) = r_;
    bs(bs == bs_) = [];
end

idx_neg = (b < 0);
if any(idx_neg)
    a_ = 1./a;
    b_ = -b(idx_neg);
    r_ = vec_power(a_,b_);
    r(idx_neg) = r_;
end
end % function vec_power

%--------------------------------------------------------------------------
% function pcase
%--------------------------------------------------------------------------

function c = pcase(S_a,S_b)
% Distinction between the following cases for the dimension sizes S_a and S_b 
% of the input parameters a and b:

%  1) (mxn)-matrix .^ scalar          ->  r(i,j) = a(i,j)^b          
    if max(S_b) == 1 
        c = 1;
        return
    end
%  2) scalar .^ (mxn)-matrix          ->  r(i,j) = a^b(i,j)       
    if max(S_a) == 1 
        c = 2;
        return
    end
%  3) (mxn)-matrix .^ (mxn)-matrix    ->  r(i,j) = a(i,j)^b(i,j)     
    if min(S_a) > 1 && min(S_b) > 1 
        c = 3;
        return        
    end
%  4) (1xn)-vector .^ (mxn)-matrix    ->  r(i,j) = a(j)^b(i,j)   
    if S_a(1) == 1 && S_a(2) == S_b(2)
        c = 4;
        return
    end
%  5) (mx1)-vector .^ (mxn)-matrix    ->  r(i,j) = a(i)^b(i,j)    
    if S_a(2) == 1 && S_a(1) == S_b(1)
        c = 5;
        return
    end
%  6) (mxn)-matrix .^ (1xn)-vector    ->  r(i,j) = a(i,j)^b(j)     
    if S_b(1) == 1 && S_a(2) == S_b(2)
        c = 6;
        return
    end
%  7) (mxn)-matrix .^ (mx1)-vector    ->  r(i,j) = a(i,j)^b(i)       
    if  S_b(2) == 1 && S_a(1) == S_b(1)
        c = 7;
        return
    end
%  8) (1xn)-vector .^ (mx1)-vector    ->  r(i,j) = a(j)^b(i)                      
    if  S_a(1) == 1 && S_b(2) == 1 
        c = 8;
        return
    end
%  9) (mx1)-vector .^ (1xn)-vector    ->  r(i,j) = a(i)^b(j)       
    if  S_a(2) == 1 && S_b(1) == 1 
        c = 9;
        return
    end
end % function pcase

%--------------------------------------------------------------------------
% function pow 
%--------------------------------------------------------------------------
function r = pow(a,b,S_a,S_b)
%--------------------------------------------------------------------------
% Implementation of a.^b for Taylor model "a" and real or interval 
% singleton b.
%
% The code is very similar to that of the well-commented function stdfun 
% in the @taylormodel/private folder which computes r = fun(a) for a unary 
% standard function fun and a Taylor model "a". 
%
% The difference is that a.^b is a binary and not unary function
% although b is a fixed singleton. The unary mimicry fun(x) := x.^b,
% b fixed, has Taylor coefficients depending on b, see the factors B in the 
% following code, wherefore the computation is slightly different 
% to that of stdfun. Since I did not want to parameterize stdfun, 
% a separate implementation is given here.  
%--------------------------------------------------------------------------

global INTLAB_ODE_VARS

ODEMODE = INTLAB_ODE_VARS.ODEMODE;
RECMODE = INTLAB_ODE_VARS.RECMODE;
HORNER_SCHEME = 1; 
USE_H_SQUARES = 2; 
POLY_EVAL_METHOD = HORNER_SCHEME;    % The Taylor polynomial is evaluated by Horner's scheme  
%POLY_EVAL_METHOD = USE_H_SQUARES;   % The Taylor polynomial p(x) = c0 + c1 * h + c2*h^2 + ... + cn * h^n
                                     % is evaluated by computing the powers h,h^2,...,h^n first where
                                     % even powers h^(2k) are obtained by squaring h^k. This diminishes 
                                     % the error interval a little bit compared to Horner's scheme,
                                     % where dependencies of h powers are not taken into account.
                                     % The evaluation in mode USE_H_SQUARES is a bit more expensive than Horner's scheme, 
                                     % since extra multiplications of the constants ci with the h powers must be performed.
N = [a.order];
n_max = max(N(:));
m_a = S_a(1); 
n_a = S_a(2);
m_b = S_b(1); 
n_b = S_b(2);
m_r = max(m_a,m_b);
n_r = max(n_a,n_b);
r(1:m_r,1:n_r) = a(1); % just cheap storage preallocation for result r

if ODEMODE == 1 && isfloat(b)
    b = intval(b);    
end

b_is_intval = isa(b,'intval');

if RECMODE ~= 2
    if b_is_intval
        B.inf = ones(n_max+2,1);
        B.sup = B.inf;
    else
        B = ones(n_max+2,1);
    end
end

i_a_old = 0;
j_a_old = 0;
i_b_old = 0;
j_b_old = 0;

for i = 1:m_r
    for j = 1:n_r
        i_a = min(i,m_a);
        j_a = min(j,n_a);
        i_b = min(i,m_b);
        j_b = min(j,n_b);
        
        new_a = (i_a ~= i_a_old || j_a ~= j_a_old); % new exponent, i.e., new component b(i_b,j_b)
        new_b = (i_b ~= i_b_old || j_b ~= j_b_old); % new base, i.e., new component a(i_a,j_a)
        
        if new_a  
            a_ = a(i_a,j_a);
        end  
        if new_b
            b_ = b(i_b,j_b);
            if RECMODE ~= 2
                if b_is_intval
                    b_iv = intval2iv(b_);
                    for k = 1:n_max+1
                        B_.inf = B.inf(k);
                        B_.sup = B.sup(k);
                        B_ = iv_times(B_,iv_rdivide(iv_plus(b_iv,1-k),k)); % B_ := B_*((b-k+1)/k);
                        B.inf(k+1) = B_.inf; % B(1) = 1, B(k+1) = prod_{j=1}^{k} (b-j+1)/j), k = 1,...,n
                        B.sup(k+1) = B_.sup;
                    end
                else
                    for k = 1:n_max+1
                        B(k+1) = B(k)*((b_+(1-k))/k); % B(1) = 1, B(k) = prod_{j=1}^{k} (b-j+1)/j), k = 1,...,n+1
                    end
                end
            end
        end
                
        if RECMODE ~= 2
            if new_a  
                n = a_.order;
                [a0,i0] = get_constant_term(a_);
                h = subtract_constant_term(a_,i0); % h := a_ - a0
                a_image = a_.image;            
                h_image = h.image;            
            end 
            % The Taylor polynomial iv2intval(p(x) of degree n of the function f(x) := x^b at x_0 := c
            % evaluated at x := a_ reads:
            %
            % p(x) := f(a0) + f'(a0) * h + f''(a0)/2 * h^2 + ... + f^(n)(a0)/n! h^n 
            %       = sum_{k=0}^n c_k * h^k 
            %
            %  where  h := x-a0 = a-a0
            %       c_k := f^(k)(a0)/k! = (prod_{j=1}^{k} (b-j+1)/j)*a0^(b-k) , k = 1,...,n 
            %       c_0 := f(a0)
            
            % computation of the Taylor coefficients C(k) := c_k, k = 0,1,...,n 
            K = (0:n)';
            if b_is_intval 
                C = intval2iv(a0.^(b_-K)); % Since b is an interval in this case, the whole computation is done in interval arithmetic. 
                B_.inf = B.inf(1:n+1);
                B_.sup = B.sup(1:n+1);
                C = iv_times(B_,C);                 
            else % b is a float and ODEMODE == 0        
                C = B(1:n+1) .* a0.^(b_-K);
            end
            if RECMODE == 1 % record write mode
                % Note that RECMODE == 1 implies ODEMODE == 1 so that this must not be checked.
                % Store n,a0,i0,a_.image in record list.
                push_reclist('n', n);
                push_reclist('a0', a0);
                push_reclist('i0', i0);
                push_reclist('C', C);
                push_reclist('B', B);
                push_reclist('a_.image', a_.image);
                push_reclist('h.image', h.image);
            end
        else % RECMODE == 2, record read mode, only the error interval component is computed !!!
            % Read n,c,i0,a_.image from record list                
            n = pull_reclist('n');
            a0 = pull_reclist('a0');
            i0 = pull_reclist('i0');
            C = pull_reclist('C');
            B = pull_reclist('B');
            a_image = pull_reclist('a_.image');
            h_image = pull_reclist('h.image');
            if new_a  
                h = subtract_constant_term(a_,i0); % h := a_ - a0                                        
            end
        end 
        
        if new_a
            if POLY_EVAL_METHOD == USE_H_SQUARES
                h_(1:n) = h;
                k_is_even = false;
                for k = 2:n
                    k_is_even = ~k_is_even;
                    if k_is_even
                        h_(k) = sqr(h_(k/2));
                    else
                        h_(k) = h_(k-1).*h;
                    end
                end
            end
        end
         
        if ODEMODE == 1                                    % verified mode
            C_.inf = C.inf(n+1);
            C_.sup = C.sup(n+1);
            if POLY_EVAL_METHOD == HORNER_SCHEME
                rs = C_;                                   % The Horner scheme starts with the highest Taylor coefficient C(n+1).
                for k = n:-1:1
                    C_.inf = C.inf(k);
                    C_.sup = C.sup(k);
                    rs = C_ + h .* rs;                     % Since h_ is a Taylor model, summation and multiplication are done in Taylor model arithmetic.
                end
            else                                           % POLY_EVAL_METHOD == USE_H_POWERS
                rs = C_.*h_(n);                            % Start with rs := C(n+1)*h^n.
                for k = n-1:-1:1
                    C_.inf = C.inf(k+1);
                    C_.sup = C.sup(k+1);
                    rs = rs + C_.*h_(k);                   % rs := rs + C(k+1)*h^k , k = n-1,n-2,...,1. Numerically, backward summation might be a little bit more stable than forward summation. 
                end  
                C_.inf = C.inf(1);
                C_.sup = C.sup(1);                
                rs = rs + C_;                              % Finally, add constant Taylor coefficient (h power is h^0 = 1). 
            end
            R = iv_plus(a_image,a_.interval);              % R := a_.image + a_.interval    , range of a_
            %H = iv_minus(R,a0);                           % H := R - a0                    , enclosure of s-a0 for all s in R
            H = iv_plus(h_image,h.interval);               % H := R - a0 = Range(h) = h.image + h.interval. This is numerically better than computing R-a0.                        
            H_pow = iv_intpower(H,n+1);                    % H_pow := H.^(n+1)
            xi = R;
            xi.inf = min(R.inf,a0); 
            xi.sup = max(R.sup,a0);                        % xi := convex hull of R and a0  , enclosure of the "intermediate" vector for estimating the Taylor series remainder
            xi_pow = intval2iv(iv2intval(xi).^(b_-(n+1))); % xi_pow = xi^(b-(n+1))
            B_.inf = B.inf(n+2);
            B_.sup = B.sup(n+2);
            I =  iv_times(iv_times(H_pow,B_),xi_pow);      % I = H_pow*B(n+2)*xi^(b-(n+1))  , enclosure of the remainder term of the Taylor expansion
            if RECMODE ~= 2
                r_ = rs;
                r_.interval = iv_plus(r_.interval,I);      % r_.interval := r_.interval + I , update error interval of result r_
                r_.image = image(r_);
            else
                r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;      % Unnecessarily, just for clearity, the component result r_ is cleared in RECMODE == 2.
                r_.interval = iv_plus(rs.interval,I);      % Only error intervals are of interest in RECMODE == 2.
            end
        else                                               % non-verified mode
            if POLY_EVAL_METHOD == HORNER_SCHEME
                r_ = C(n+1);                               % The Horner scheme starts with the highest Taylor coefficient C(n+1).
                for k = n:-1:1
                    r_ = C(k) + h .* r_;                   % Since h_ is a Taylor model, summation and multiplication are done in Taylor model arithmetic.
                end
            else                                           % POLY_EVAL_METHOD == USE_H_POWERS
                r_ = C(n+1).*h_(n);                        % Start with r_ := C(n+1)*h^n.
                for k = n-1:-1:1
                    r_ = r_ + C(k+1).*h_(k);               % Since h_ is a Taylor model, summation and multiplication are done in Taylor model arithmetic.
                end  
                r_ = r_ + C(1);                            % Finally, add constant Taylor coefficient (h power is h^0 = 1).               
            end          
            r_.interval = INTLAB_ODE_VARS.EMPTYIV;
            r_.image = INTLAB_ODE_VARS.EMPTYIV;
        end        
        r(i,j) = r_;
        j_a_old = j_a;
        j_b_old = j_b;
    end 
    i_a_old = i_a;
    i_b_old = i_b;
end 

end % function pow


%--------------------------------------------------------------------------
% function pow_ implements a^b for Taylor model "a" and 
% real or interval singleton b. 
%
% Alternative implementation to subfunction pow
%--------------------------------------------------------------------------

function r = pow_(a,b,S_a,S_b)

global INTLAB_ODE_VARS

ODEMODE = INTLAB_ODE_VARS.ODEMODE;
RECMODE = INTLAB_ODE_VARS.RECMODE;

m_a = S_a(1); 
n_a = S_a(2);
m_b = S_b(1); 
n_b = S_b(2);
m_r = max(m_a,m_b);
n_r = max(n_a,n_b);
r(1:m_r,1:n_r) = a(1); % just cheap storage preallocation for result r

for i = 1:m_r
    for j = 1:n_r
        i_a = min(i,m_a);
        j_a = min(j,n_a);
        i_b = min(i,m_b);
        j_b = min(j,n_b);
                
        a_ = a(i_a,j_a);
        b_ = b(i_b,j_b);
        
        if RECMODE ~= 2
            n = a_.order;
            [c,i0] = get_constant_term(a_);
            if ODEMODE == 1  % verified mode
                c = intval(c);
                D = intval(ones(n+2,1));
            else
                D = ones(n+2,1);
            end
            
            % Taylor polynomial p(x) of degree n of the function f(x) := x^b at x_0 := c
            % evaluated for x := a_. 
            %
            % p(x) := f(c) + f'(c)(x-c) + f''(c)(x-c)^2/2 + ... + f^(n)(c)(x-c)^n / n! 
            %       = sum_{k=0}^n c_k(x-c)^k 
            %
            %  with c_k := f^(k)(c)/k! = (prod_{j=1}^{k} (b-j+1)/j)*c^(b-k) , k = 0,1,...,n 
            %
            
            h = subtract_constant_term(a_,i0); % fast computation of a_-c
            
            % verified computation of the Taylor coefficients C(k) := c_k, k = 0,1,...,n 
            K = (0:n)';
            C = c.^(b_-K);
            %D = intval(ones(n+2,1));
            for k = 1:n
                D(k+1) = D(k)*((b_-k+1)/k);
                C(k+1) = D(k+1)*C(k+1);  
            end
            D(n+2) = D(n+1)*((b_-n)/(n+1));
            
            if RECMODE == 1 % record write mode
                % Note that RECMODE == 1 implies ODEMODE == 1 so that this has not to be checked.
                % Store n,c,i0,a_.image in record list
                push_reclist('n', n);
                push_reclist('c', c);
                push_reclist('i0', i0);
                push_reclist('C', C);
                push_reclist('D', D);
                push_reclist('a_.image', a_.image);
            end
            % Remark: The previous push_reclist statements in RECMODE == 1 are located in front of the subsequent for loop
            %         because in it Taylor model arithmetic is done which also writes to the record list
            %         in a later read mode that linear order of recording must be kept.
            p = C(n+1);
            for k = n:-1:1 % Horner scheme evaluation
                p = C(k)+h.*p;
            end
            r_ = p; %Initialize the result component with p.            
            if ODEMODE == 1 % verified mode.
                I_y = iv2intval(iv_plus(a_.image,a_.interval)); % see [E], p.30
                h_ = I_y-c;
                I_ry =  h_^(n+1) * D(n+2) * I_y^(b_-(n+1));                                
                r_.interval = iv_plus(p.interval,intval2iv(I_ry));  % p.interval + I_ry, see [E], p.31
                r_.image = image(r_);
            else % non-verified mode
                r_.interval = INTLAB_ODE_VARS.EMPTYIV;
                r_.image = INTLAB_ODE_VARS.EMPTYIV;
            end
        else % RECMODE == 2, record read mode, only the error interval component is computed !!!
            % Read n,c,i0,a_.image from record list
            n = pull_reclist('n');
            c = pull_reclist('c');
            i0 = pull_reclist('i0');
            C = pull_reclist('C');
            D = pull_reclist('D');
            a_image = pull_reclist('a_.image');
            
            h = subtract_constant_term(a_,i0); % fast computation of a-c
                                               % Remark: In RECMODE == 2 actually h = a_ would also be o.k.since only interval components 
                                               % are of interest/considered and for h := subtract_constant_term(a_,i0) we have h.interval = a_.interval.             
            p = C(n+1);
            for k = n:-1:1 % Horner scheme evaluation
                p = C(k)+h.*p;
            end
            
            I_y = iv2intval(iv_plus(a_image,a_.interval)); % See [E], p.30.
            h_ = I_y-c;
            I_ry =  h_^(n+1) * D(n+2) * I_y^(b_-(n+1));                                
            
            r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;              % Unnecessarily, just for clarity the component result r_ is cleared in RECMODE == 2.
            r_.interval = iv_plus(p.interval,intval2iv(I_ry)); % p.interval + I_ry, see [E], p.31.
        end  % RECMODE ~= 2
        
        r(i,j) = r_;
    end % j
end %i

end % function pow_
