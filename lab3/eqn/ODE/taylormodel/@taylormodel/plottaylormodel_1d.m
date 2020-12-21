function [x,y] = plottaylormodel_1d(a,idx,kmax1,kmax2,kmax3,mode)
%PLOTTAYLORMODEL  verified 1d-plots of Taylor model ranges
%
%  [x,y] = plottaylormodel_1d(a,idx,kmax1,kmax2,kmax3,mode)
%
% The plots are overall verified! This means that the displayed enclosures
% are also valid for all points which are not discretization points.
%
% Input arguments:  
%      a : Taylor model vector of length m. A row vector is converted to a column vector.
%          All Taylor models a(i) must also have the same domain except for component idx.
%          The domains for the component idx must be consecutive w.r.t. i, i.e., 
%          a(i).domain.sup(idx) = a(i+1).domain.inf(idx) for i = 1... m-1.
%          Typically, these domains are the time intervals 
%          [t_0 t_1],[t_1 t_2],...,[t_m-1 t_m] discretizing the integration 
%          period from t_0 to t_e =: t_m of an ODE. 
%          The Taylor models a(i) are plotted with respect to component idx (x-axis).
%          Lower and upper bounds of a(i) and a(i+1) are continuously connected   
%          at the boundary points a(i).domain.sup(idx) = a(i+1).domain.inf(idx).
%  idx   : domain index of Taylor models, default: last domain index n:=dim(a(i)).
%  kmax1 : optional splitting number. Each domain interval 
%          [a(i,j).domain.inf(k),a(i,j).domain.sup(k)] 
%          k <> idx, is subdivided into kmax1 subintervals. 
%          If kmax1 is not specified, then the default is kmax1 := 1.  
%  kmax2 : optional splitting number. Each domain interval 
%          [a(i).domain.inf(idx),a(i).domain.sup(idx)] 
%          is subdivided into kmax2 subintervals. 
%          If kmax2 is not specified, then the default is kmax2 := 10.          
%          The boundary points of these subintervals are the endpoints 
%          of the piecewise linear upper and lower bound curves.
%          kmax2 is irrelevant for "mode" = 0 (pointwise evaluation).
%  kmax3 : optional splitting number. Each of the kmax2 subintervals of 
%          [a(i).domain.inf(idx),a(i).domain.sup(idx)] 
%          is again subdivided into kmax3 subintervals in order to compute 
%          better values for first and second derivatives. These derivatives 
%          are needed to plot the verified, piecewise linear lower and upper 
%          bound curves. The default is kmax3 = 10.
%          kmax3 is irrelevant for "mode" = 0 (pointwise evaluation).
%  mode:  plotting mode, 
%         "mode" = 0 means "pointwise" evaluation. Upper and lower bounds
%            are plotted as points only at the boundaries
%            a(i).domain.inf(idx), a(i).domain.sup(idx), i=1,...,m.
%         "mode" = 1 (default) means "interval" evaluation of variable "idx". 
%            upper and lower bounds are plotted as piecewise linear lines. 
%
% Output arguments: 
%         x,y are the plotted x- and y-values, where the lower bound function values 
%         are given by y.inf and the upper bound function values by y.sup. 
%
%         If mode = 0, then the lower bound points (x(i),y.inf(i)) are plotted as red stars '*' 
%         and the upper bound points (x(i),y.sup(i)) are plotted as green stars '*'.
%
%         If mode = 1, then the lower bound points (x(i),y.inf(i)) are linearly connected by
%         red lines and the upper bound points (x(i),y.sup(i)) are linearly connected by
%         green lines. We emphasize that these piecewise linear upper and lower bound curves
%         are rigorously verified, i.e., the enclosure is also verified for function values 
%         of intermediate x-values which are not grid points x(i)!
%
%--------
%
% WARNING: Plotting verified inclusions of Taylor model images is quite 
%          critical to performance. Please use this function only for 
%          moderate matrix sizes, moderate Taylor model orders, small
%          domain dimensions, and small splitting numbers kmax1, kmax2, kmax3! 

% written  11/07/16     F. Buenger

e = 1e-30;

if 1+e > 1 % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards
end
if nargin < 1 
    error('wrong call.')
end
if isempty(a)
    error('Taylor model is empty.')
end
S_a = size(a); 
m = numel(a);
if length(S_a) > 2
  error('maximally two dimensions for type taylormodel')  
end
if S_a(1) > 1 &&  S_a(2) > 1
  error('Taylor model must be a vector.')
end
if S_a(2) > 1 
  a = a'; % Transform row vector to column vector.  
end
if nargin < 2 || isempty(idx)
    idx = a(1,1).dim; % The default index is the last index of the domain.
end
if nargin < 3 || isempty(kmax1)
    kmax1 = 1; % The default of number of domain subdivisions for component idx.
end
if nargin < 4 || isempty(kmax2)       
    kmax2 = 100; % The default of number of domain subdivisions for component idx.
end
if nargin < 5 || isempty(kmax3)       
    kmax3 = 10; 
end
if nargin < 6 || isempty(mode)
    mode = 1; % default mode
end

dim_a = get_dim(a);
n = max(dim_a(:)); 
if min(min(dim_a)) < n || idx > n || idx < 1 
    error('inconsistent input')
end

IDX = [1:idx-1,idx+1:n]; % complementary indices

% Initialization of idx-domains. 
% Later: x_idx.inf(i) := a(i).domain.inf(idx) and x_idx.sup(i) := a(i).domain.sup(idx), i=1,...,m
x_idx.inf = zeros(m,1); 
x_idx.sup = zeros(m,1);

for i = 1:m
    if any(a(1).domain.inf(IDX) ~= a(i).domain.inf(IDX)) || ...
            any(a(1).domain.sup(IDX) ~= a(i).domain.sup(IDX))
        error('Taylor model domains except for component idx must be equal.');
    end
    if i < S_a(1) && a(i).domain.sup(idx) ~= a(i+1).domain.inf(idx)
        error('Taylor model domains for component idx must be consecutive.');
    end
    x_idx.inf(i) = a(i).domain.inf(idx);
    x_idx.sup(i) = a(i).domain.sup(idx);
end

% Subdivide all domains for components distinct from idx into kmax1 subintervals.
M = kmax1^(n-1);
X0.inf = zeros(M,n); 
X0.sup = X0.inf;
for j = IDX
    x =  linspace(a(1).domain.inf(j),a(1).domain.sup(j),kmax1+1)'; 
    w_inf = x(1:kmax1);
    w_sup = x(2:kmax1+1);  
    subs1 = kron((1:kmax1)',ones(kmax1^(n-1-j),1));
    Mj = kmax1^(j-1);
    X0.inf(:,j) = repmat(w_inf(subs1),Mj,1); 
    X0.sup(:,j) = repmat(w_sup(subs1),Mj,1); 
end

X.inf = repmat(X0.inf,kmax2+1,1);
X.sup = repmat(X0.sup,kmax2+1,1);
subs2 = kron((1:kmax2+1)',ones(M,1));

K = kmax2*kmax3;
XX.inf = repmat(X0.inf,K,1);
XX.sup = repmat(X0.sup,K,1);
subs3 = kron((1:K)',ones(M,1));
subs4 = kron((1:kmax2)',ones(kmax3*M,1));


y_low = zeros(kmax2,2); % Initialize array for intervals of piecewise linear, lower bound curve.
y_up  = y_low;          % Initialize array for intervals of piecewise linear, upper bound curve.

da = deriv(a,1,idx);    % Compute first derivatives.
dda = deriv(a,2,idx);   % Compute second derivatives.

% Prepare linear dummy-Taylor model dummy(x) = x_idx-center_idx .
dummy = a(1);
dummy.interval.inf = 0;
dummy.interval.sup = 0;
dummy.monomial = zeros(1,n);
dummy.monomial(1,idx) = 1;   
dummy.coefficient = 1;

x = zeros(m*kmax2+1,1); % Initialize x-values for the whole interval D := [a(1).domain.inf(idx),a(m).domain.sup(idx)].
y.inf = x;              % Initialize corresponding function value inclusions.
y.sup = x;

for i = 1:m
    a_ = a(i);
    dummy.domain = a_.domain;      
    dummy.center = a_.center;  
    dummy.image = image(dummy);
    da_ = da(i);
    dda_ = dda(i);
    c_idx = a_.center(idx);
    switch mode
        case 0 % point plot 
            x_ =  linspace( x_idx.inf(i),x_idx.sup(i),kmax2+1)'; % Create kmax2+1 equispaced grid points of the [domain a_.domain.inf(idx),domain a_.domain.sup(idx)]. 
            X.inf(:,idx) = x_(subs2);                            % Insert x in column idx of X for each partition of the remaining domains with index distinct from idx.
            X.sup(:,idx) = x_(subs2);  
            y_ = evaltaylormodel(a_,X);                          % Evaluate a_ on all subdomains contained in the rows of X.  
            y_.inf = accumarray(subs2,y_.inf,[],@min);           % Accumulate minima for all coarse grid points. 
            y_.sup = accumarray(subs2,y_.sup,[],@max);           % Accumulate maxima for all coarse grid points. 
            I = (i-1)*kmax2+1:i*kmax2+1;                         % Index set for grid points of Taylor model a(i).
            x(I) = x_;
            if i == 1
              y.inf(I) = y_.inf;
              y.sup(I) = y_.sup;
            else % For subsequent Taylor models the enclosure at the left boundary of of the idx-domain the enclosure can only become wider. 
                % Therefore, the enclosure of the previous Taylor model at its right boundary of the idx-domain is taken / not overwritten.
              y.inf(I(2:end)) = y_.inf(2:end);
              y.sup(I((2:end))) = y_.sup(2:end);            
            end
        case 1 % piecewise linear plot
            % Subdivide domain for component idx into kmax1 subintervals.
            xx = linspace(x_idx.inf(i),x_idx.sup(i),K+1)';     % fine grid
            x_  = xx(1:kmax3:K+1);                             % coarse grid            
            X.inf(:,idx) = x_(subs2);                          % Insert coarse grid points (point intervals) in column idx of X for each partition of the remaining domains.
            X.sup(:,idx) = x_(subs2);  
            
            y_ = evaltaylormodel(a_,X);                        % Evaluate a_ on all subdomains contained in the rows of X.  
            y_.inf = accumarray(subs2,y_.inf,[],@min);         % Accumulate minima for all coarse grid points. 
            y_.sup = accumarray(subs2,y_.sup,[],@max);         % Accumulate maxima for all coarse grid points. 
            
            dy_ = evaltaylormodel(da_,X);                      % Evaluate derivatives da_ on all subdomains contained in the rows of X.  
            dy_.inf = accumarray(subs2,dy_.inf,[],@min);       % Accumulate minima for all coarse grid points. 
            dy_.sup = accumarray(subs2,dy_.sup,[],@max);       % Accumulate maxima for all coarse grid points. 

            diff_x = diff(x_);
            s_lower = (y_.inf(2:end)-y_.inf(1:end-1))./diff_x; % non-verified computation of the slopes of piecewise linear lower bound curve
            s_upper = (y_.sup(2:end)-y_.sup(1:end-1))./diff_x; % non-verified computation of the slopes of piecewise linear upper bound curve                                                
            
            xx_inf = xx(1:end-1);                              % left boundaries of fine grid intervals
            xx_sup = xx(2:end);                                % right boundaries of fine grid intervals
            XX.inf(:,idx) = xx_inf(subs3);                     % Insert fine grid intervals in column idx of X for each partition of the remaining domains.
            XX.sup(:,idx) = xx_sup(subs3);  
            ddy = evaltaylormodel(dda_,XX);  
            ddy.inf = accumarray(subs4,ddy.inf,[],@min);       % Accumulate minima of second derivatives on all coarse(!) subintervals. 
            ddy.sup = accumarray(subs4,ddy.sup,[],@max);       % Accumulate maxima of second derivatives on all coarse(!) subintervals.  
            
            idx_convex = (ddy.inf >= 0);                       % Indices of all coarse subintervals on which a_ with respect to variable idx is convex.
            idx_concave = (ddy.inf <= 0);                      % Indices of all coarse subintervals on which a_ with respect to variable idx is concave.
            
            % On convex subintervals the straight line between the upper values at the grid points build an upper bound for all intermediate function values.
            % The lower values at the left interval boundaries are taken as starting points of straight lines having the lower derivatives at these points as their slopes. 
            % Convexity implies that these straight lines are lower bounds.
            I1 = [idx_convex;false];
            I2 = [false;idx_convex];
            y_up(idx_convex,:) = [y_.sup(I1),y_.sup(I2)];            
            y_low(idx_convex,1) = y_.inf(I1);
            v = iv_plus(y_.inf(I1),iv_times(dy_.inf(I1),iv_minus(x_(I2),x_(I1)))); % interval evaluation of y_.inf(I1) + dy_.inf(I1)*(x_(I2)-x_(I1));
            y_low(idx_convex,2) = v.inf; 
            
            % On concave subintervals the straight line between the lower values at the grid points build a lower bound for all intermediate function values.
            % The upper values at the left interval boundaries are taken as starting points of straight lines having the upper derivatives at these points as their slopes. 
            % Concavity implies that these straight lines are upper bounds.
            I1 = [idx_concave;false];
            I2 = [false;idx_concave];
            y_low(idx_concave,:) = [y_.inf(I1),y_.inf(I2)];
            y_up(idx_concave,1) = y_.sup(I1);
            v = iv_plus(y_.sup(I1),iv_times(dy_.sup(I1),iv_minus(x_(I2),x_(I1)))); % interval evaluation of y_.sup(I1) + dy_.sup(I1)*(x_(I2)-x_(I1));
            y_up(idx_concave,2) = v.sup;
            
            % Now the (hopefully) few remaining intervals on which a_ is neither convex nor concave w.r.t. x_idx are considered.
            % (Sadly a vectorization seems not possible so that this is done in a loop.)
            jj = 1:kmax2; 
            jj = jj(~or(idx_convex,idx_concave));
            for j = jj
                J = ((j-1)*M*kmax3+1:j*M*kmax3);
                z.inf = XX.inf(J,:); % z is a partition of subdomains corresponding to the j-th coarse grid subinterval.
                z.sup = XX.sup(J,:); 
                
                s = s_upper(j);
                d = a_-s*dummy;
                b = evaltaylormodel(d,z);
                b = max(b.inf);                                      % This implies a_(x-c_idx) + [err] <= s*(x-c_idx)+b for all x in [x(j),x(j+1)].
                y1 = iv_plus(iv_times(s,iv_minus(x_(j),c_idx)),b);   % s*(x(j)-c_idx) + b, left upper bound
                y2 = iv_plus(iv_times(s,iv_minus(x_(j+1),c_idx)),b); % s*(x(j+1)-c_idx) + b, right upper bound
                y_up(j,:) = [y1.sup,y2.sup];
                
                s = s_lower(j);
                d = a_-s*dummy;
                b = evaltaylormodel(d,z);
                b = min(b.inf);                                      % This implies a_(x-c_idx) + [err] >= s*(x-c_idx)+b for all x in [x(j),x(j+1)].
                y1 = iv_plus(iv_times(s,iv_minus(x_(j),c_idx)),b);   % s*(x(j)-c_idx) + b, left lower bound
                y2 = iv_plus(iv_times(s,iv_minus(x_(j+1),c_idx)),b); % s*(x(j+1)-c_idx) + b, right lower bound
                y_low(j,:) = [y1.inf,y2.inf];
            end
            
            % At the gridpoints the piecewise linear lower bound shall be continuous.
            y_.inf = [y_low(1,1);
                min(y_low(2:end,1),y_low(1:end-1,2));
                y_low(end,2)];
            
            % At the gridpoints the piecewise linear upper bound shall be continuous.
            y_.sup = [y_up(1,1);
                max(y_up(2:end,1),y_up(1:end-1,2));
                y_up(end,2)];
            
            % The connections between two consecutive Taylor models a(i) and a(i+1) shall be continuous.
            I = (i-1)*kmax2+1:i*kmax2+1; % Index set for grid points of Taylor model a(i).
            x(I) = x_;
            if i == 1
                y.inf(I) = y_.inf;
                y.sup(I) = y_.sup;
            else
                % y.inf(I(1)) is the lower bound of the previous Taylor model a(i-1) at the right boundary b of the domain for the variable with index idx.
                % y_.inf(1) is the lower bound of the actual Taylor model a(i) at the left boundary b'=b of the domain for the variable with index idx.
                % In order to make the connection continuous, the minimum of both values  y.inf(I(1)) and y_.inf(1).
                y.inf(I(1)) = min(y.inf(I(1)),y_.inf(1));
                y.inf(I(2:end)) = y_.inf(2:end);
                % same procedure as before for upper bound
                y.sup(I(1)) = max(y.sup(I(1)),y_.sup(1));
                y.sup(I(2:end)) = y_.sup(2:end);
            end
    end
end
switch mode
    case 0 % pointwise plot
        plot(x,y.inf,'*r',x,y.sup,'*g');
    case 1 % piecewise linear plot
        plot(x,y.inf,'r',x,y.sup,'g');
end

if rndold ~= 1
    setround(rndold)
end

end % function plottaylormodel