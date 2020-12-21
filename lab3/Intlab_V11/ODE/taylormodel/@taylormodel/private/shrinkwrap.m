function r = shrinkwrap(a,base_tm)
%SHRINKWRAP   shrink wrapping of Taylor model vector a.
%
%   res = shrinkwrap(a,base_tm)
%
% The implementation is based on 
% [Bue] F. Buenger, "Shrink wrapping for Taylor models revisited", Numerical Algorithms 78(4), pp. 1001-1017, 2018

% written  06/29/16     F. Buenger

global INTLAB_ODE_OPTIONS

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

n = a(1).dim-1; % order of the underlying ODE system
% Split the Taylor model in constant, linear, and remainder parts.
A = get_linear_terms(a); % real nx(n+1)-matrix of the linear part of "a"
A = A(:,1:n); % The (n+1)-th column of A corresponds to the linear part of the time variable
              % which stands at the last position (n+1) for each component of "a".
              % This column should be zero since for shrink wrapping, i.e., "a" should not depend
              % on the time variable. This final column is cut of so that A becomes an nxn matrix.
              
if INTLAB_ODE_OPTIONS.blunting
    A = blunt(A);
end
                          
[c,idx] = get_constant_term(a); % c := constant terms of a, idx stores the corresponding monomial indices.
a0 = subtract_constant_term(a,idx); % a0 := a - (constant terms of "a")
                       
cond_max = 1E16; % heuristic bound for the condition number of A. Feel free to change this constant!
%cond_max = 1E3; % heuristic bound for the condition number of A. Feel free to change this constant!
                 % The value stems from [E], p.139, for Van-der-Pol equation.
cond_A = cond(A);
if cond_A > cond_max
    err_code = 2; % Shrink wrapping seems not promising for such a bad conditioned A. Try other strategies.
else 
    [q,err_code] = shrinkwrap_factor(a0,A);
    if err_code == 0 % shrink wrap factors were successfully determined.
        r = a; % Initialize result r(x):=a(q*x) with a.
               % Note that a already contains the constant term c again in contrast to a0. 
        for i = 1:n
            r_ = r(i);            
            M = r_.monomial(:,1:n); % Neglect last (n+1-th) column for time exponents which are all zero here.  
            
            % Compute a verified enclosure of the coefficients of r(x) := a(q*x) 
            % Note that all shrink wrap factors q(i) >= 1 so that sparsity of coefficients cannot occur.
            % Recall also that rounding is upwards!

            P = prod((q').^M,2);
            c_lower =  -(r_.coefficient .* (-P));
            c_upper = r_.coefficient .* P;
            c_mid = 0.5 * (c_lower + c_upper); % estimate for interval midpoints 
            r_.coefficient = c_mid ;           % Store the corresponding coefficients as those of the result r = a(q*b).
                        
            % Determine error interval 
            rest = r_;
            coeff.inf = -(c_mid-c_lower); % centered lower bound for coefficients
            coeff.sup = c_upper-c_mid;    % centered upper bound for coefficients 
            rest.coefficient = coeff;     % interval coefficients
            
            r_.interval = image(rest);    % error interval                 
            r_.image = image(r_);
            r(i) = r_;
        end
    else
        err_code = 1; % Shrink wrap factors could not be determined.
    end
end  

% If shrink wrapping was not possible, then calculate either a box or a parallelotope enclosure.

if err_code ~= 0
    % First check if we really have to use a box or parallelotope enclosure or if we can use the
    % given Taylor model for the next integration step without doing shrink wrapping.
    
    % Calculate an approximate of the volume of the remainder box.
    fac = 0.01^n; % needed for the condition later (see below)
    vol = 2.0^n;  % 2.0 is the length of the interval [-1,1].      

    a0_interval = get_interval(a0); % Note that the error interval of a0 equals the error interval of a.
    a0_image = get_image(a0); % image of polynomial part of a0    
    
    vol_interval = prod(a0_interval.sup-a0_interval.inf);  % non-verified computation of the volume of the error interval box    
    vol_parallelotope = det(A) .* vol;                     % non-verified calculation of the volume of the parallelotope A*[-1,1]^n
    
    % Now we check if the volume of the remainder box is in relation to the parallelotope volume too big.
    % In this case a box or parallelotope enclosure is used. Otherwise nothing is done and a is returned. 
    
    if vol_interval >= (fac * vol_parallelotope)
        % The remainder term has, relative to the parallelotope volume, much bigger size.
        % So we calculate a verified parallelotope and a verified box enclosure and take the one with the smaller volume.
        
        box = iv_plus(a0_image,a0_interval); % verified computation of the box enclosure:  "box = a0_image + a0_interval"
        vol_box = prod(box.sup-box.inf);     % non-verified calculation of the volume of the box enclosure
        
        % Calculate the parallelotope enclosure and the corresponding volume if possible.
        
        if err_code ~= 2 % Matrix inversion didn't fail. The calculation of a parallelotope enclosure seems possible.
            
            R = intval2iv( inv(intval(A,A,'infsup')) ); % R is a verified inverse of A computed with INTLAB.
            
            % Compute the nonlinear polynomial part h = a0(x)-Ax of a0 (without remainder, just polynomials) simply by deleting the linear terms in a0.
            h = a0; % initialization               
            for i = 1:n
                idx = find(sum(h(i).monomial,2)==1); % Find linear terms of Taylor model a(i).
                if ~isempty(idx)
                    h(i).monomial(idx,:) = [];       % delete linear monomials
                    h(i).coefficient(idx) = [];      % ... and corresponding coefficients
                    h(i).interval.inf = 0;           % error interval is zero, i.e., [0,0]
                    h(i).interval.sup = 0;                    
                    % Empty Taylor models shall be avoided.
                    if isempty(h(i).monomial)
                        h(i).monomial = zeros(1,n+1);
                        h(i).coefficient = 0;        % zero-polynomial 
                        h(i).image.inf = 0;          % The image is directly set to the zero-interval [0,0].
                        h(i).image.sup = 0;
                    else
                        h(i).image = image(h(i)); % Compute image of h(i).
                    end
                end
            end
            
            Rh = R*h;            % verified enclosure of A^{-1}*h
            E = a0_interval;     % E = error interval of a0 =  error interval of a
            RE = iv_mtimes(R,E); % RE := R*E verified matrix-vector-product 
            
            % Compute verified enclosure [-r,r] of Rh.image + Rh.interval + RE.
            Rh_image = get_image(Rh); 
            Rh_interval = get_interval(Rh); 
            hlp = iv_plus(iv_plus(Rh_image,Rh_interval),RE);  % Rh_image + Rh_interval + RE
            r = max(abs(hlp.inf),abs(hlp.sup));

            q = 1+r; % Simplified shrink wrap factor. Recall that rounding is upwards.
            vol_parallelotope = vol_parallelotope * prod(q); % non-verified volume of parallelotope enclosure using det(A*diag(q)) = det(A)*prod(q)
        end        
        if err_code == 2 || vol_box < vol_parallelotope % Note that the non-verified computed quantities vol_box and vol_parallelotope
                                                        % are only used for deciding which one of the two verified enclosures is taken.
                                                        % Thus in any case the chosen enclosure is verified!
            % Take the box enclosure. 
            % Calculate the box representation with Taylor models each defined on [-1,1].
            [box_mid,box_rad] = iv_getmidrad(box);
            r = box_mid + box_rad .* base_tm + c; % verified "box enclosure" of a
        else % Take the parallelotope enclosure.
            r = A * (diag(q) * base_tm) + c; % r := A * diag(q) + c computed in (verified) Taylor model arithmetic. 
        end 
    else % Do nothing, i.e., return the original input Taylor model "a" without any changes.
        r = a; 
    end 
end

if rndold ~= 1 
    setround(rndold)
end

end  % function shrinkwrap
