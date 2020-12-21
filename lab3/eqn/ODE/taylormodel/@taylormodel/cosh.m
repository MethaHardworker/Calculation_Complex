function r = cosh(a)
%SINH  Taylor model hyperbolic cosine  cosh(a)
%
%   r = cosh(a)

% written  09/09/15     F. Buenger
% modified 11/23/15     F. Buenger  matrix input (componentwise evaluation)
% modified 12/18/15     F. Buenger  switch between verified/non-verified Taylor model arithmetic
% modified 01/20/16     F. Buenger  record feature
% modified 06/01/18     F. Buenger  "intval"-components --> intval-like structures
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

global INTLAB_ODE_INVERSE_FACTORIAL
global INTLAB_ODE_VARS

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@cosh,'cosh');
else
    % global INTLAB_ODE_INVERSE_FACTORIAL
    % global INTLAB_ODE_VARS
    
    ODEMODE = INTLAB_ODE_VARS.ODEMODE;
    RECMODE = INTLAB_ODE_VARS.RECMODE;
    
    if ODEMODE == 1 || RECMODE == 2
        theta = intval(0,1,'infsup');
    end
    
    if RECMODE ~= 2
        N = [a.order];
        n_max = max(N(:));
        if isempty(INTLAB_ODE_INVERSE_FACTORIAL) || length(INTLAB_ODE_INVERSE_FACTORIAL.FLP) < (n_max + 2)
            set_inverse_factorial_(n_max+1);
        end
    end
    
    S_a = size(a);
    r = a; % Initialize result r with a. (just storage preallocation)
    
    for i = 1:S_a(1)
        for j = 1:S_a(2)
            a_ = a(i,j);
            if RECMODE ~= 2
                n = a_.order;
                [c,i0] = get_constant_term(a_);
                C = zeros(1,n+2);
                if ODEMODE == 1 % verified mode.
                    c = intval(c,c,'infsup');
                    C = intval(C,C,'infsup');
                    F = INTLAB_ODE_INVERSE_FACTORIAL.INTVAL(1:n+2);
                else
                    F = INTLAB_ODE_INVERSE_FACTORIAL.FLP(1:n+2);
                end
                % Taylor polynomial p(x) of degree n of the function phi(x):=cosh(x) at x_0:=c
                % evaluated for x:=a_ . See [E], p.23 and p.29,30,35.
                %
                % p(x) := cosh(c) + sinh(c)*h + 1/2!*cosh(c)*h^2 + 1/3!*sinh(c)*h^3 ...
                %
                n = a_.order;
                h = subtract_constant_term(a_,i0); % fast computation of a_ - c
                C(1:2:n+2) = cosh(c);
                C(2:2:n+2) = sinh(c);
                C = C.*F;
                if RECMODE == 1 % record write mode
                    % Note that RECMODE == 1 implies ODEMODE == 1 so that this has not to be checked.
                    % Store n,c,i0,C,a_.image in record list.
                    push_reclist('n', n);
                    push_reclist('c', c);
                    push_reclist('i0', i0);
                    push_reclist('C', C);
                    push_reclist('a_.image', a_.image);
                end
                % Remark: The previous push_reclist statements in RECMODE == 1 are located in front of the subsequent for loop
                %         because in it Taylor model arithmetic is done which also writes to the record list.
                %         In a later read mode that linear order of recording must be kept.
                p = C(n+1);
                for k = n:-1:1 % Horner scheme evaluation
                    p = C(k)+h.*p;
                end
                r_ = p; %initialize the component result with p
                if ODEMODE == 1 % verified mode.
                    % Computation of the error interval r.interval, see [E], p.32,33, and p. 30,31.
                    I_y = iv2intval(iv_plus(a_.image,a_.interval)); % See [E], p.30.
                    
                    % see [E], p.30,35 error term caused by the remainder of the Taylor expansion
                    % r_(x) := 1/(n+1)!*h^(n+1)*s1
                    % s1 = -s2 if n mod 4 = 1,2
                    % s1 = -s2 else
                    % s2:= cos(c+theta*h) if n is even,
                    % s2:= sin(c+theta*h) if n is odd,
                    %      theta in (0,1)
                    
                    h_ = I_y-c;
                    if even(n)
                        I_ry = F(n+2)*h_^(n+1)*sinh(c+theta*h_);       % See [E], p.30,35.
                    else
                        I_ry = F(n+2)*h_^(n+1)*cosh(c+theta*h_);       % See [E], p.30,35.
                    end
                    r_.interval = iv_plus(p.interval,intval2iv(I_ry)); % See [E], p.31.
                    r_.image = image(r_);
                else % non-verified mode.
                    r_.interval = INTLAB_ODE_VARS.EMPTYIV;
                    r_.image = INTLAB_ODE_VARS.EMPTYIV;
                end
            else % RECMODE == 2, record read mode, only the error interval component is computed !!!
                % Read n,c,i0,C,a_.image from record list.
                n = pull_reclist('n');
                c = pull_reclist('c');
                i0 = pull_reclist('i0');
                C = pull_reclist('C');
                a_image = pull_reclist('a_.image');
                
                h = subtract_constant_term(a_,i0); % fast computation of a_ - c
                % Remark: In RECMODE == 2 actually h = a_ would also be o.k.since only interval components
                % are of interest/considered and for h := subtract_constant_term(a_,i0) we have h.interval = a_.interval.
                p = C(n+1);
                for k = n:-1:1 % Horner scheme evaluation
                    p = C(k)+h.*p;
                end
                I_y = iv2intval(iv_plus(a_image,a_.interval)); % See [E], p.30.
                h_ = I_y-c;
                F = INTLAB_ODE_INVERSE_FACTORIAL.INTVAL(n+2); % From the record write call in the probably far past, still
                % length(INTLAB_ODE_INVERSE_FACTORIAL.INTVAL) >= (n + 2) holds true.
                if even(n)
                    I_ry = F * h_^(n+1)*sinh(c+theta*h_); % See [E], p.30,35.
                else
                    I_ry = F * h_^(n+1)*cosh(c+theta*h_); % See [E], p.30,35.
                end
                r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;              % Unnecessarily, just for clarity the component result r_ is cleared in RECMODE == 2.
                r_.interval = iv_plus(p.interval,intval2iv(I_ry)); % See [E], p.31.
            end
            r(i,j) = r_;
        end
    end
end

end % function cosh
