function r = stdfun(a,fun,funname)
%STDFUN  evaluation of unary standard functions, like exp, log, sin, cos,...,   
%        for Taylor models 
%
%   r = stdfun(a,fun,funname) 
%
% input arguments: 
%   
%         a: Taylor model to which the function fun shall be applied
%       fun: function handle of a standard function, for example fun = @exp
%   funname: name string of the standard function, for example funname = 'exp'
%
% output argument: 
%
%        r = fun(a), r has, like "a", type taylormodel.  
% 
% Description of the algorithm: 
%
%   For simplicity let us assume that the input argument "a" is a single 
%   Taylor model (not a vector or matrix of Taylor models). 
%   Then, the enclosure described by "a" is simply its range R which is 
%   given by R := a.image + a.interval. Recall that a.image is an interval 
%   enclosure of the image of the polynomial part P of "a", and that the 
%   error interval "a.interval" covers all rounding and truncation errors.
%   Let a0 := P(0,...,0) be the constant term of P. Then, since a.interval 
%   is, or at least should be, a tight interval around zero(!), we have 
%   that a0 is an element of the range R. 
%   It is quite reasonable to consider a0 somehow as the "center" of R. 
%   Therefore, Taylor expansion of the standard function "fun" with expansion 
%   point a0 makes most sense for computing fun(a). Thus, if we define 
%  
%       h := a - a0, 
%
%   which is a Taylor model whose polynomial part now has zero constant 
%   term and which is simply a constant shift of "a", then the Taylor
%   series of order n becomes
%
%       r := fun(a0) + fun'(a0)*h + fun''(a0)/2 * h^2 + ... + fun^(n)(a0)/n! * h^n . 
%
%   The right-hand side is evaluated by Horner's scheme and uses 
%   Taylor model arithmetic for summation and multiplication only.  
%   This is the expensive part of the algorithm.
%
%   The remainder term is enclosed by interval evaluation(!) of
%
%       I := fun^(n+1)(xi)/(n+1)! * H^(n+1),  xi := hull(a0,R), 
%                                              H := R-a0 = Range(h) 
%                                                 = h.image + h.interval
%  
%   where xi is the interval hull of a0 and R which encloses
%   the unknown intermediate value. Note that almost always 
%   
%       hull(a0,R) = R
%  
%   holds true, since almost always a0 is an element of R, but we do not 
%   want to exclude pathological, exceptional cases.  
%
%   Updating the error interval of r by 
%
%       r.interval := r.interval + I 
%
%   gives the final verified Taylor model result r = fun(a). 
% 
%   Here, "verified" means the following: For each real vector x in the 
%   domain of "a" which is also the domain of r we have that
%   
%       fun(P(x-x0) + a.interval)   is contained in   Q(x-x0) + r.interval   
%
%   where Q is the polynomial part of r and x0 := a.center = r.center 
%   is the common center point of the common domain a.domain = r.domain.

% written  06/08/18     F. Buenger

global INTLAB_ODE_VARS
global INTLAB_ODE_TCOEFF

ODEMODE = INTLAB_ODE_VARS.ODEMODE;
RECMODE = INTLAB_ODE_VARS.RECMODE;

HORNER_SCHEME = 1; 
USE_H_SQUARES = 2; 
POLY_EVAL_METHOD = HORNER_SCHEME;   % The Taylor polynomial is evaluated by Horner's scheme  
%POLY_EVAL_METHOD = USE_H_SQUARES;  % The Taylor polynomial p(x) = c0 + c1 * h + c2*h^2 + ... + cn * h^n
                                    % is evaluated by computing the powers h,h^2,...,h^n first where
                                    % even powers h^(2k) are obtained by squaring h^k. This diminishes 
                                    % the error interval a little bit compared to Horner's scheme,
                                    % where dependencies of h powers are not taken into account.
                                    % The evaluation in mode USE_H_SQUARES is a bit more expensive than Horner's scheme, 
                                    % since extra multiplications of the constants ci with the h powers must be performed.
if RECMODE ~= 2
    if strcmp(funname,'exp') || strcmp(funname,'sin') || strcmp(funname,'cos') || ...
       strcmp(funname,'sinh') || strcmp(funname,'cosh')
        N = [a.order];
        n_max = max(N(:));
        if length(INTLAB_ODE_TCOEFF.INVFACT_FLP) < (n_max + 2)
            set_inverse_factorial(n_max+1);
        end
    end
    if strcmp(funname,'sqrt')
        N = [a.order];
        n_max = max(N(:));
        if length(INTLAB_ODE_TCOEFF.SQRT_FLP) < (n_max + 2)
            set_sqrt_factor(n_max);
        end
    end
end

S_a = size(a);
r = a; % Initialize result r with a. (just storage preallocation)
C_.inf = 0;
C_.sup = 0;

for i = 1:S_a(1)
    for j = 1:S_a(2)
        a_ = a(i,j);                
        if RECMODE ~= 2
            n = a_.order;
            [a0,i0] = get_constant_term(a_);            % a0 := P(0,...,0) is the constant term of the polynomial part P of a_.
            h = subtract_constant_term(a_,i0); % h := a_ - a0                    
            C = taylor_coeff(a0,fun,funname,n,ODEMODE); % Compute Taylor coefficients C of function fun upto order n with expansion point a0.
            a_image = a_.image;
            h_image = h.image;
            if RECMODE == 1 % record write mode. Note that RECMODE == 1 implies ODEMODE == 1 so that this must not be checked.
                % Write n, a0, i0, C, a_.image to record list.
                % This record mode stuff is only done for performance reasons!
                push_reclist('n', n);
                push_reclist('a0', a0);
                push_reclist('i0', i0);
                push_reclist('C', C);
                push_reclist('a_.image', a_.image);
                push_reclist('h.image', h.image);                
            end
        else % RECMODE == 2, record read mode, only the error interval component is computed!
            % Read n, a0, i0, C, a_.image from the record list.
            n = pull_reclist('n');
            a0 = pull_reclist('a0');
            i0 = pull_reclist('i0');
            C = pull_reclist('C');
            a_image = pull_reclist('a_.image');
            h_image = pull_reclist('h.image');
            h = subtract_constant_term(a_,i0); % h := a_ - a0                    
        end
                
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
        
        %------------------------------------------------------------------------------
        % Horner scheme evaluation of the Taylor series of the
        % standard function fun with expansion point a0 upto order n:
        %
        %   r_ := fun(a0) + fun'(a0)*h + fun''(a0)/2*h^2 + ... + fun^(n)(a0)/n!*h^n .
        %       = C(1)    + C(2)*h     + C(3)*h^2        + ....+ C(n+1)*h^n
        %
        %   (*) C(k) := fun^(k)(a0)/k!
        %
        % In ODEMODE == 0 (non-verified mode) the Taylor coefficients C(i)
        % are floats. In ODEMODE == 1 (verified mode) they are intervals.
        %------------------------------------------------------------------------------
                
        if ODEMODE == 1                               % verified mode.            
            C_.inf = C.inf(n+1);
            C_.sup = C.sup(n+1);
            if POLY_EVAL_METHOD == HORNER_SCHEME
                rs = C_;                              % The Horner scheme starts with the highest Taylor coefficient C(n+1).
                for k = n:-1:1
                    C_.inf = C.inf(k);
                    C_.sup = C.sup(k);
                    rs = C_ + h .* rs;                % Since h_ is a Taylor model, summation and multiplication are done in Taylor model arithmetic.
                end
            else                                      % POLY_EVAL_METHOD == USE_H_POWERS
                rs = C_.*h_(n);                       % Start with rs := C(n+1)*h^n.
                for k = n-1:-1:1
                    C_.inf = C.inf(k+1);
                    C_.sup = C.sup(k+1);
                    rs = rs + C_.*h_(k);              % rs := rs + C(k+1)*h^k , k = n-1,n-2,...,1. Numerically, backward summation might be a little bit more stable than forward summation. 
                end   
                C_.inf = C.inf(1);
                C_.sup = C.sup(1);
                rs = rs + C_;                         % Finally, add constant Taylor coefficient (h power is h^0 = 1). 
                
            end
            %----------------------------------------------------------------
            % Compute enclosure I of the remainder of the Taylor expansion
            %
            %   I := fun^(n+1)(xi)/(n+1)! * H^(n+1),
            %
            % where xi := hull(a0,R), H := R-a0, R := a.image + a.interval .
            %----------------------------------------------------------------
            R = iv_plus(a_image,a_.interval);         % R := a_.image + a_.interval    , range of a_
            %H = iv_minus(R,a0);                      % H := R - a0                    , enclosure of s-a0 for all s in R
            H = iv_plus(h_image,h.interval);          % H := R - a0 = Range(h) = h.image + h.interval. This is numerically better than computing R-a0.
            xi = R;
            xi.inf = min(R.inf,a0); 
            xi.sup = max(R.sup,a0);                   % xi := convex hull of R and a0   , enclosure of the "intermediate" vector for estimating the Taylor series remainder
            H_pow = iv_intpower(H,n+1);               % H_pow := H.^(n+1)
            c = taylor_coeff(xi,fun,funname,n+1,1,1); % c := fun^(n+1)(xi)/(n+1)!          , n+1 - th Taylor coefficient of fun at xi
            I = iv_times(c,H_pow);                    % I := fun^(n+1)(xi)/(n+1)! * H^(n+1), enclosure of the remainder term of the Taylor expansion
            
            if RECMODE ~= 2
                r_ = rs;
                r_.interval = iv_plus(r_.interval,I); % r_.interval := r_.interval + I , update error interval of result r_
                r_.image = image(r_);
            else
                r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL; % Unnecessarily, just for clearity, the component result r_ is cleared in RECMODE == 2.
                r_.interval = iv_plus(rs.interval,I); % Only error intervals are of interest in RECMODE == 2.
            end
        else                                          % non-verified mode.
            if POLY_EVAL_METHOD == HORNER_SCHEME
                r_ = C(n+1);                          % The Horner scheme starts with the highest Taylor coefficient C(n+1).
                for k = n:-1:1
                    r_ = C(k) + h .* r_;              % Since h_ is a Taylor model, summation and multiplication are done in Taylor model arithmetic.
                end
            else                                      % POLY_EVAL_METHOD == USE_H_POWERS
                r_ = C(n+1).*h_(n);                   % Start with r_ := C(n+1)*h^n.
                for k = n-1:-1:1
                    r_ = r_ + C(k+1).*h_(k);          % Since h_ is a Taylor model, summation and multiplication are done in Taylor model arithmetic.
                end                
                r_ = r_ + C(1);                       % Finally, add constant Taylor coefficient (h power is h^0 = 1). 
            end
            r_.interval = INTLAB_ODE_VARS.EMPTYIV;    % Error intervals and images are not of interest in non-verified mode.
            r_.image = INTLAB_ODE_VARS.EMPTYIV;       % For clarity both components are cleared.
        end                
        r(i,j) = r_;
    end
end

end % function stdfun

function C = taylor_coeff(x,fun,funname,n,odemode,remainder)
% Compute Taylor coefficients fun^(i)(x)/i! of function fun at x
% upto order n or just for order n, if "remainder" is specified.
%
% For some functions, like exp, sin, cos those Taylor coefficients 
% are easy to compute explicitly. For many other standard functions, 
% like asin, explicit formulas are complicated and recursion formulas are
% better suited. For that the INTLAB taylor toolbox (data type taylor) or 
% The INTLAB AWA toolbox containing the data type tcoeff can be used
% For non-verified mode (odemode == 0), always the taylor toolbox 
% is used which works with floats where computation is fast. 
% In verified mode (odemode == 1), the data type tcoeff from the 
% AWA toolbox is the default which is 3-4 times faster than
% the Taylor toolbox. Nevertheless, also in verified mode, you can freely 
% switch between both toolboxes by setting the constant TAYLOR_TOOLBOX:
% 
%  TAYLOR_TOOLBOX = 1;  % ---> Use data type tcoeff from the AWA toolbox. 
%  TAYLOR_TOOLBOX = 2;  % ---> Use data type taylor from the taylor toolbox.

global INTLAB_AWA_VARS
global INTLAB_ODE_TCOEFF

% Choose toolbox for computing Taylor coefficients
TCOEFF = 1; 
TAYLOR = 2;   
TAYLOR_TOOLBOX = TCOEFF;  % ---> Use data type tcoeff from the AWA toolbox.    
%TAYLOR_TOOLBOX = TAYLOR; % ---> Use data type taylor from the taylor toolbox.  

remainder = (nargin == 6); % If remainder == true , then only the n-th Taylor coefficient is returned.
                           % If remainder == false, then all Taylor coefficients up to order n are returned 
if remainder
    n_ = n;        
    C.inf = 0;             % Initialize result C with zero. C is a single interval for enclosing the n-th Taylor coefficient f^(n)(x)/k!
    C.sup = 0;             % (Remark: The function is called by stdfun with "n = n'+1" so that in fact the remainder of a Taylor expansion of order n' will be returned.)
else
    n_ = (0:n)'; 
    if odemode == 0
        C = zeros(n+1,1);
    else
        C.inf = zeros(n+1,1);  % Initialize result C with zero. C is an interval vector of length n+1  
        C.sup = C.inf;         % for enclosing the Taylor coefficients f^(k)(x)/k!, k = 0,...,n
    end
end                            
n_idx = n_ + 1;

if odemode == 0 % non-verified mode
    x_ = x;
else % verified mode
    if isfloat(x)
        x_ = intval(x,x,'infsup');
    else
        x_ = iv2intval(x);
    end
end

switch funname
    case 'inv'   % inv(x) := 1/x,  inv^(k)(x)/k! = (-1)^k * (1/x)^(k+1), k = 0,1,...          
        if odemode == 0
            C = -(-1./x_).^(n_+1);                
        else
            C = iv_uminus( iv_intpower(iv_rdivide(-1,x),n_+1) ); % C := -(-1./x).^(n_+1);
        end
    case 'sqrt'  % sqrt^(k)(x)/k! = (-1)^(k-1)*(2k-3)!!/(k!*2^k) * x^(1/2-k)
        if odemode == 0
            C = INTLAB_ODE_TCOEFF.SQRT_FLP;
            C = C(n_idx); 
            C = C .* x_.^(0.5-n_);
        else
            C = INTLAB_ODE_TCOEFF.SQRT_IV;               
            C.inf = C.inf(n_idx);
            C.sup = C.sup(n_idx);
            C = iv_times(C,intval2iv(x_.^(0.5-n_))); % C := C .* x_.^(0.5-n_);           
        end
    case 'exp'  % exp^(k)(x)/k! = exp(x)/k! 
        if odemode == 0
            C = INTLAB_ODE_TCOEFF.INVFACT_FLP;
            C = C(n_idx).*exp(x_); 
        else
            C = INTLAB_ODE_TCOEFF.INVFACT_IV;               
            C.inf = C.inf(n_idx);
            C.sup = C.sup(n_idx);
            C = iv_times(C,intval2iv(exp(x_))); % C := C(n_idx).*exp(x_);             
        end
    case 'log'   % log^(k)(x)/k! = (-1)^(k+1)/(k*x^k) for k >= 1
        if odemode == 0
            n_ = (1:n)';
            C = [log(x_);-(-x_).^(-n_)./n_];                
        else
            if remainder
                C = iv_rdivide(-1,iv_times(iv_intpower(iv_uminus(x),n_),n_));  %  C := -1 / (n_.*(-x).^n_)
            else
                log_x = intval2iv(log(x_));
                C.inf(1) = log_x.inf;
                C.sup(1) = log_x.sup;
                n_ = (1:n)';
                C_ = iv_rdivide(-1,iv_times(iv_intpower(iv_uminus(x),n_),n_)); %  C := -1 / (n_.*(-x).^n_)
                C.inf(2:n+1) = C_.inf;
                C.sup(2:n+1) = C_.sup;
            end
        end                
    case 'sin'   % sin^(k)(x)/k! = (-1)^j * f(x)/k!, where j = 1 if (n-1) mod 4 in {1,2} and f(x) = cos(x) if n is even and f(x) = sin(x) if n is odd. 
        if odemode == 0
            F = INTLAB_ODE_TCOEFF.INVFACT_FLP;
            F = F(n_idx);
            sin_x = sin(x_);
            cos_x = cos(x_);
            C(1:4:n+1) = sin_x;
            C(2:4:n+1) = cos_x;
            C(3:4:n+1) = -sin_x;
            C(4:4:n+1) = -cos_x;  
            C = C.*F;
        else
            F = INTLAB_ODE_TCOEFF.INVFACT_IV;               
            F.inf = F.inf(n_idx);
            F.sup = F.sup(n_idx);
            sin_x = intval2iv(sin(x_));
            cos_x = intval2iv(cos(x_));
            if remainder
                mod_4 = mod(n-1,4);
                if mod(n-1,2) == 0  % n is even
                    if mod_4 == 2
                        C = iv_uminus(cos_x); % C := -cos_x
                    else
                        C = cos_x;                            
                    end
                else              % n is odd
                    if mod_4 == 1
                        C = iv_uminus(sin_x); % C := -sin_x;   
                    else
                        C = sin_x;                            
                    end
                end
            else
                C.inf(1:4:n+1) = sin_x.inf;  
                C.sup(1:4:n+1) = sin_x.sup;   
                C.inf(2:4:n+1) = cos_x.inf;  
                C.sup(2:4:n+1) = cos_x.sup;
                C.inf(3:4:n+1) = -sin_x.sup; 
                C.sup(3:4:n+1) = -sin_x.inf;
                C.inf(4:4:n+1) = -cos_x.sup; 
                C.sup(4:4:n+1) = -cos_x.inf;
            end
            C = iv_times(C,F);
        end
    case 'cos'   % cos^(k)(x)/k! = (-1)^j * f(x)/k!, where j = 1 if (n-1) mod 4 in {0,1} and f(x) = sin(x) if n is even and f(x) = cos(x) if n is odd.           
        if odemode == 0
            F = INTLAB_ODE_TCOEFF.INVFACT_FLP;
            F = F(n_idx);
            sin_x = sin(x_);
            cos_x = cos(x_);
            C(1:4:n+1) = cos_x;
            C(2:4:n+1) = -sin_x;
            C(3:4:n+1) = -cos_x;
            C(4:4:n+1) = sin_x;  
            C = C.*F;
        else
            F = INTLAB_ODE_TCOEFF.INVFACT_IV;               
            F.inf = F.inf(n_idx);
            F.sup = F.sup(n_idx);
            sin_x = intval2iv(sin(x_));
            cos_x = intval2iv(cos(x_));
            if remainder
                mod_4 = mod(n-1,4);
                if mod(n-1,2) == 0  % n-1 is even
                    if mod_4 == 0
                        C = iv_uminus(sin_x); % C := -sin_x   
                    else
                        C = sin_x;                            
                    end
                else              % n-1 is odd
                    if mod_4 == 1
                        C = iv_uminus(cos_x); % C := -cos_x;    
                    else
                        C = cos_x;                            
                    end
                end
            else
                C.inf(1:4:n+1) = cos_x.inf;
                C.sup(1:4:n+1) = cos_x.sup;
                C.inf(2:4:n+1) = -sin_x.sup;
                C.sup(2:4:n+1) = -sin_x.inf;
                C.inf(3:4:n+1) = -cos_x.sup;
                C.sup(3:4:n+1) = -cos_x.inf;
                C.inf(4:4:n+1) = sin_x.inf;
                C.sup(4:4:n+1) = sin_x.sup;
            end
            C = iv_times(C,F);
        end        
    case 'sinh'  % sinh^(k)(x)/k! = cosh(x)/k! if k is odd and sinh^(k)(x)/k! = sinh(x)/k! if k is even
        if odemode == 0
            F = INTLAB_ODE_TCOEFF.INVFACT_FLP;
            F = F(n_idx);
            C(1:2:n+1) = sinh(x_);
            C(2:2:n+1) = cosh(x_);
            C = C.*F;            
        else
            F = INTLAB_ODE_TCOEFF.INVFACT_IV;               
            F.inf = F.inf(n_idx);
            F.sup = F.sup(n_idx);            
            sinh_x = intval2iv(sinh(x_));
            cosh_x = intval2iv(cosh(x_));
            if remainder
                if mod(n-1,2) == 0  % n is even
                    C = cosh_x;
                else
                    C = sinh_x;
                end
            else
                C.inf(1:2:n+1) = sinh_x.inf;
                C.sup(1:2:n+1) = sinh_x.sup;
                C.inf(2:2:n+1) = cosh_x.inf;                
                C.sup(2:2:n+1) = cosh_x.sup;                
            end
            C = iv_times(C,F);
        end        
    case 'cosh'  % cosh^(k)(x)/k! = sinh(x)/k! if k is odd and cosh^(k)(x)/k! = cosh(x)/k! if k is even
        if odemode == 0
            F = INTLAB_ODE_TCOEFF.INVFACT_FLP;
            F = F(n_idx);            
            C(1:2:n+1) = cosh(x_);
            C(2:2:n+1) = sinh(x_);
            C = C.*F;                        
        else
            F = INTLAB_ODE_TCOEFF.INVFACT_IV;               
            F.inf = F.inf(n_idx);
            F.sup = F.sup(n_idx);            
            sinh_x = intval2iv(sinh(x_));
            cosh_x = intval2iv(cosh(x_));
            if remainder
                if mod(n-1,2) == 0  % n is even
                   C = sinh_x;
                else
                   C = cosh_x; 
                end
            else
                C.inf(1:2:n+1) = cosh_x.inf;
                C.sup(1:2:n+1) = cosh_x.sup;
                C.inf(2:2:n+1) = sinh_x.inf;                
                C.sup(2:2:n+1) = sinh_x.sup;                
            end
            C = iv_times(C,F);
        end                      
    otherwise % Compute Taylor coefficients by chosen Taylor toolbox, that is, by automatic differentiation. 
        if odemode == 0 || TAYLOR_TOOLBOX == TAYLOR
            xt = taylorinit(x_,n);         
            C = fun(xt);
            C = struct(C).t;          % If odemode == 0, then Taylor coefficients have type float otherwise (odemode == 1) they have type intval.
            if odemode == 1
                C = intval2iv(C);     % conversion of type intval to interval-like structures with components .inf and .sup.
            end
        elseif TAYLOR_TOOLBOX == TCOEFF
            if isfloat(x)
                xt.inf = [x;1;zeros(n-1,1)];
                xt.sup = xt.inf;
            else
                z = zeros(n-1,1);
                xt.inf = [x.inf;1;z];
                xt.sup = [x.sup;1;z];                
            end
            xt = taylorcoeffinit(xt); % Taylor coefficients are interval-like structures with components .inf and .sup.
            INTLAB_AWA_VARS.STATUS = [];
            C = fun(xt);
            C = struct(C);
        end
        if remainder % Only the last, n+1-th Taylor coefficient is returned.
            C.inf = C.inf(n+1); 
            C.sup = C.sup(n+1); 
        end
end

end % function taylor_coeff