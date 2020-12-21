function r = log(a)
%EXP  Taylor model logarithm log(a)
%
%   r = log(a)

% written  09/09/15     F. Buenger
% modified 11/23/15     F. Buenger  matrix input (componentwise evaluation)
% modified 12/18/15     F. Buenger  switch between verified/non-verified Taylor model arithmetic
% modified 01/20/16     F. Buenger  record feature
% modified 06/01/18     F. Buenger  "intval"-components --> intval-like structures
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

global INTLAB_ODE_VARS

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@log,'log');
else
    %global INTLAB_ODE_VARS
    
    ODEMODE = INTLAB_ODE_VARS.ODEMODE;
    RECMODE = INTLAB_ODE_VARS.RECMODE;
    
    if ODEMODE == 1 || RECMODE == 2
        theta = intval(0,1,'infsup');
    end
    
    S_a = size(a);
    r = a; % Initialize result r with a_. (just storage preallocation)
    
    for i = 1:S_a(1)
        for j = 1:S_a(2)
            a_ = a(i,j);
            if RECMODE ~= 2
                n = a_.order;
                rng_ = iv_plus(a_.image,a_.interval);
                if ~(rng_.inf > 0)
                    error('a.image + a.interval must be contained in (0,+Inf).')
                end
                
                [c,i0] = get_constant_term(a_);
                nn = (1:n+1);
                if ODEMODE == 1
                    c = intval(c,c,'infsup');
                end
                
                % Taylor polynomial p(x) of degree n of the function phi(x):=log(x) at x_0:=c
                % evaluated for x:=a_ . See [E], p.23 and p.29,30,34.
                %
                % p(x) := log(c) + h/c -1/2*h^2/c^2+...+(-1)^(n+1) * 1/n * h^n/c^n
                
                h = subtract_constant_term(a_,i0); % fast computation of a-c
                C = [log(c),-((-1/c).^nn)./nn];
                if RECMODE == 1 % record write mode
                    % Note that RECMODE == 1 implies ODEMODE == 1 so that this has not to be checked.
                    % Store n,c,C,i0,a_.image in record list.
                    push_reclist('n', n);
                    push_reclist('c', c);
                    push_reclist('C', C);
                    push_reclist('i0', i0);
                    push_reclist('a_.image', a_.image);
                end
                % Remark: The previous push_reclist statements in RECMODE == 1 are located in front of the subsequent for loop
                %         because in it Taylor model arithmetic is done which also writes to the record list.
                %         In a later read mode that linear order of recording must be kept.
                p = C(n+1);
                for k = n:-1:1 % Horner scheme evaluation
                    p = C(k)+h.*p;
                end
                r_ = p; % Initialize the result component with p.
                if ODEMODE == 1 % verified mode.
                    % Computation of the error interval r_.interval, see [E], p.32,33, and p. 30,31.
                    I_y = iv2intval(iv_plus(a_.image,a_.interval)); % See [E], p.30.
                    
                    % See [E], p.30,34, error term caused by the remainder of the Taylor expansion
                    % r(x) := (-1)^(n+2)*1/(n+1)*h^(n+1)/c^(n+1)*1/(1+theta*h/c)^(n+1), theta in (0,1).
                    h_ = I_y-c;
                    I_ry = C(n+2) .* h_^(n+1) ./ (1+theta*h_./c).^(n+1); % See [E], p.30,34.
                    r_.interval =  iv_plus(p.interval,intval2iv(I_ry));  % See [E], p.31.
                    r_.image = image(r_);
                else % non-verified mode
                    r_.interval = INTLAB_ODE_VARS.EMPTYIV;
                    r_.image = INTLAB_ODE_VARS.EMPTYIV;
                end
            else % RECMODE == 2, record read mode, only the error interval component is computed !!!
                % Read n,c,C,i0,a_.image from record list.
                n = pull_reclist('n');
                c = pull_reclist('c');
                C = pull_reclist('C');
                i0 = pull_reclist('i0');
                a_image = pull_reclist('a_.image');
                rng_ = iv_plus(a_image,a_.interval);
                if ~(rng_.inf > 0)
                    error('a.image + a.interval must be contained in (0,+Inf).')
                end
                h = subtract_constant_term(a_,i0); % fast computation of a-c
                % Remark: In RECMODE == 2 actually h = a_ would also be o.k. since only interval components
                % are of interest/considered and for h := subtract_constant_term(a_,i0) we have h.interval = a_.interval.
                p = C(n+1);
                for k = n:-1:1 % Horner scheme evaluation
                    p = C(k)+h.*p;
                end
                I_y = iv2intval(iv_plus(a_image,a_.interval)); % See [E], p.30.
                % See [E], p.30,34, error term caused by the remainder of the Taylor expansion
                % r(x) := (-1)^(n+2)*1/(n+1)*h^(n+1)/c^(n+1)*1/(1+theta*h/c)^(n+1), theta in (0,1).
                h_ = I_y-c;
                I_ry = C(n+2) .* h_^(n+1) ./ (1+theta*h_./c).^(n+1); % See [E], p.30,34.
                r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;                % Unnecessarily, just for clarity the component result r_ is cleared in RECMODE == 2.
                r_.interval = iv_plus(p.interval,intval2iv(I_ry));   % See [E], p.31.
            end
            r(i,j) = r_;
        end
    end
end

end % function log


