function r = inv(a)
%INV  Taylor model elementwise inverse 1/a  (same as a.^-1)
%
%   r = inv(a)

% written  09/08/15     F. Buenger
% modified 11/23/15     F. Buenger  matrix input (componentwise evaluation)
% modified 12/18/15     F. Buenger  switch between verified/non-verified Taylor model arithmetic
% modified 01/20/16     F. Buenger  record feature
% modified 06/13/18     F. Buenger  encapsulation in general function stdfun

global INTLAB_ODE_VARS

TRIGGER = true;
%TRIGGER = false;

if TRIGGER
    r = stdfun(a,@inv,'inv');
else    
    % global INTLAB_ODE_VARS
    
    ODEMODE = INTLAB_ODE_VARS.ODEMODE;
    RECMODE = INTLAB_ODE_VARS.RECMODE;
    
    if ODEMODE == 1 || RECMODE == 2
        theta = intval(0,1,'infsup');
    end
    
    S_a = size(a);
    r = a; % Initialize result r with a. (just storage preallocation)
    
    for i = 1:S_a(1)
        for j = 1:S_a(2)
            a_ = a(i,j);
            if RECMODE ~= 2
                n = a_.order;
                [c,i0] = get_constant_term(a_);
                if ODEMODE == 1
                    c = intval(c,c,'infsup');
                end
                
                % Taylor polynomial p(x) of degree n of the function phi(x) := 1/x at x_0 := c
                % evaluated for x := a_. See [E], p.23 and p.30.
                %
                % p(x) := 1/c * { 1 - (x-c)/c + (x-c)^2/c^2 - ... + (-1)^n (x-c)^n/c^n}
                %
                h = subtract_constant_term(a_,i0); % fast computation of a_-c
                C = - (-1/c).^(1:n+1);
                if RECMODE == 1 % record write mode
                    % Note that RECMODE == 1 implies ODEMODE == 1 so that this has not to be checked.
                    % Store n,c,i0,a_.image in record list
                    push_reclist('n', n);
                    push_reclist('c', c);
                    push_reclist('i0', i0);
                    push_reclist('C', C);
                    push_reclist('a_.image', a_.image);
                end
                % Remark: The previous push_reclist statements in RECMODE == 1 are located in front of the subsequent for loop
                %         because in it Taylor model arithmetic is done which also writes to the record list
                %         in a later read mode that linear order of recording must be kept
                p = C(n+1);
                for k = n:-1:1 % Horner scheme evaluation
                    p = C(k)+h.*p;
                end
                r_ = p; %Initialize the result component with p.
                if ODEMODE == 1 % verified mode.
                    % Computation of the error interval r.interval, see [E], p.32,33, and p. 30,31
                    I_y = iv2intval(iv_plus(a_.image,a_.interval)); % See [E], p.30.
                    
                    % See [E], p.30,32, error term caused by the remainder of the Taylor expansion
                    % r_(x) := (-1)^(n+1)/c * (x-c)^(n+1)/c^(n+1) * 1/(1+theta*(x-c)/c)^(n+2)
                    %        = h^(n+1) / ( c * (1-theta*h)^(n+2) ) , theta in (0,1), h:= 1-x/c
                    h_ = 1-I_y/c;
                    I_ry = h_.^(n+1) ./ ( c.*(1-theta.*h_).^(n+2) );
                    
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
                a_image = pull_reclist('a_.image');
                
                h = subtract_constant_term(a_,i0); % fast computation of a-c
                % Remark: In RECMODE == 2 actually h = a_ would also be o.k.since only interval components
                % are of interest/considered and for h := subtract_constant_term(a_,i0) we have h.interval = a_.interval.
                p = C(n+1);
                for k = n:-1:1 % Horner scheme evaluation
                    p = C(k)+h.*p;
                end
                I_y = iv2intval(iv_plus(a_image,a_.interval));     % See [E], p.30.
                h_ = 1-I_y/c;
                I_ry = h_.^(n+1) ./ ( c.*(1-theta.*h_).^(n+2) );
                r_ = INTLAB_ODE_VARS.ZEROTAYLORMODEL;              % Unnecessarily, just for clarity the component result r_ is cleared in RECMODE == 2.
                r_.interval = iv_plus(p.interval,intval2iv(I_ry)); % p.interval + I_ry, see [E], p.31.
            end
            
            r(i,j) = r_;
        end
    end
end

end % function inv
