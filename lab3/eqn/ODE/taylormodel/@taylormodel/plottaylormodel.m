function plottaylormodel(a,T,idx,parts,color,fillcolor,mode)
%PLOTTAYLORMODEL  plotting for Taylor models,
%                 two- and three-dimensional phase plots are possible. 
%
%   plottaylormodel(a,T,idx,m,color,fillcolor,mode)
%
% Input arguments:  
%  a : mxn-Taylor model matrix. 
%      All Taylor models in one row a(i,:) must have the same domain.   
%      All Taylor models a(:,j) must also have the same domain except for 
%      j = idx.
%
%  T : vector of floating point values at which the variable with index idx 
%      shall be evaluated.
%      If "mode" = 0, then for each T(k) the first Taylor model line a(i,:)
%      which contains T(k) in its common domains a(i,:).domain(idx) 
%      is evaluated and displayed for T(k).
%      If T is empty or not specified, then by default T is taken to be
%      the set of boundary points of a(i,:).domain(idx).
%
%      If "mode" = 1, then for each T(k) the first Taylor model line a(i,:)
%      which contains T(k) in a(i,:).domain(idx) is evaluated and displayed 
%      on the  whole interval [T(k), min(T(k+1),sup(a(i,:).domain(idx)))].
%      If T is initial or not specified, then by default 
%      T := linspace(t0,tf,500), where t0 is the smallest left boundary of 
%      some a(i,:).domain(idx) and tf is the largest right boundary of some
%      a(i',:).domain(idx).
% 
%  idx   : domain index of Taylor models, default: last domain index n.
%  parts : optional splitting number. Each domain interval a(i,j).domain(k), 
%          k <> idx, is subdivided into "parts" subintervals. If "parts" is  
%          not specified, then the default is parts := 8.  
%  
%  color, fillcolor : Same meaning as in function 'plotintval'.  
%
%  mode: plotting mode, 
%        "mode" = 0 means "pointwise" evaluation and 
%        "mode" = 1 (default) means "interval" evaluation of variable "idx", 
%          see description of parameter "T".
%          In mode "1" upper and lower bounds are plotted as piecewise
%          constant step functions! This mode makes most sense if 
%          sup(a(i,:).domain(idx)) = inf(a(i+1,:).domain(idx))
%          for all i. 
%
%--------
%
% WARNING: Plotting verified inclusions of Taylor model images is quite tedious.
%          Please use this function only for moderate matrix sizes, 
%          moderate Taylor model orders, small domain dimensions, 
%          and small splitting numbers "parts"!

% written  27/11/15     F. Buenger


S_a = size(a); 
if length(S_a) > 2
  error('maximally two dimensions for taylormodel')  
end
if S_a(2) > 3 
  error('Only 2d or 3d-plots are possible.')
end
if nargin < 1 
    error('wrong call.')
end
if nargin < 2 
   T = []; % default T = [].
else
    T = T(:); % T must be a column vector.
end
if nargin < 3 || isempty(idx)
    idx = a(1,1).dim;           % The default index is the last index of the domain.
end
if nargin < 4 || isempty(parts)
    parts = 8;                  % The default of number of domain subdivisions.
end
if nargin < 5 || isempty(color) % Use plot with blue boundary.
    color = 'b';
end
if nargin < 6 || isempty(fillcolor)
    fillcolor = 1;
end
if nargin < 7 || isempty(mode)
    mode = 1; 
end

dim_a = get_dim(a);
n = max(dim_a(:)); 
if min(dim_a(:)) < n || idx > n || idx < 1 
    error('inconsistent input')
end

IDX = [1:idx-1,idx+1:n];
n = a(1,1).dim;

% A priori check that for all Taylor models a(i,j) all domains
% a(i,j).domain are the same with the only possible exception
% a(i,j).domain(idx) <> a(i',j').domain(idx) for i<>i'.
% In particular, this means that all Taylor models in one row a(i,:) 
% must have the same domain.
if isempty(T)
    switch mode
        case 0 % pointwise evaluation
            TT = zeros(2*S_a(1),1); % Initialize evaluation point vector with zeros.
        case 1 % interval evaluation
            t0 = a(1,1).domain.inf(idx); % Initialize minimal value for variable idx.
            tf = a(1,1).domain.sup(idx); % Initialize maximal value for variable idx.
    end
end
for i = 1:S_a(1)
    for j = 1:S_a(2)
      if any(a(1,1).domain.inf(IDX) ~= a(i,j).domain.inf(IDX)) || ...
         any(a(1,1).domain.sup(IDX) ~= a(i,j).domain.sup(IDX)) || ...
             a(i,1).domain.inf(idx) ~= a(i,j).domain.inf(idx)  || ...
             a(i,1).domain.sup(idx) ~= a(i,j).domain.sup(idx)
          error('Taylor model domains may only differ for the specified index idx');
      end
    end
    if isempty(T)
        switch mode
            case 0 % pointwise evaluation
                TT(2*i-1) = a(i,1).domain.inf(idx);  % left boundary
                TT(2*i) = a(i,1).domain.sup(idx);    % right boundary
            case 1 % interval evaluation
                t0 = min(t0,a(i,1).domain.inf(idx)); % actual minimal value for variable idx
                tf = max(tf,a(i,1).domain.sup(idx)); % actual maximal value for variable idx
        end
    end
end

if isempty(T)
    switch mode
        case 0                     % pointwise evaluation
            T = TT;                % all boundary points
        case 1                     % interval evaluation
            N = 500;               % quite arbitrary number. Feel free to change that!
            T = linspace(t0,tf,N); % fixed number of equal size subintervals
            T = T';                % work with column vectors
    end
end

% Subdivide all domains V(j):=a(1,1).domain(k), j=1,...,n, k<>idx into m subdomains V(i,j) of equal diameter. 
% Then, generate all combinations of subdomains: 
% V(i_1,1) x ... x V(i_{idx-1},idx-1) x {0} x V(i_{idx+1},idx+1) x ... x V(i_n,n) 
% and store the result in a big matrix X. The number of such combinations is m^(n-1). 
% Therefore X has size (m^(n-1),n).

X = intval(zeros(parts^(n-1),n)); % Initialize X with zeros.
for j = IDX                                                           % Subdivide only domain intervals with index distinct from idx. 
    x =  linspace(a(1,1).domain.inf(j),a(1,1).domain.sup(j),parts+1); % Build m subdomains of domain j.
    x = x';
    w = infsup(x(1:parts),x(2:parts+1));                              % Transfer the m subdomains to type 'intval'. 
    w = w( kron((1:parts)',ones(parts^(n-1-j),1)) );                  % Build multiple copies of w according to its index j.
    X(:,j) = repmat(w,parts^(j-1),1);                                 % Insert these multiple copies in column j of X.
end

Y = intval( zeros(S_a(2),size(X,1)) ); % initialization for Taylor model images + error intervals

if S_a(2) == 1  % 2d - plot of variable idx (x-direction) and Taylor model values (y-direction)
    Y_min = []; % Initialize plot values in y-direction, lower bound.
    Y_max = []; % Initialize plot values in y-direction, upper bound.
    T_plt = []; % Initialize plot values in x-direction.
else
    YY = [];    % Initialize phase plot values.
end

for i = 1:S_a(1)
    D.inf = a(i,1).domain.inf(idx); % D is the common idx-domain-component of all Taylor models in row a(i,:).
    D.sup = a(i,1).domain.sup(idx);
    t_idx = and(T >= D.inf,T <= D.sup);
    t = T( t_idx ); % t := (T \cap D ), intersection of T and D. In typical applications t has small too moderate length.           
    if mode == 1 % interval evaluation
        if isempty(t)
            t = [a(i,1).domain.inf(idx);a(i,1).domain.sup(idx)]; % evaluate the whole domain.
        else
            if D.inf < t(1)
                t = [D.inf;t]; % Add left domain boundary.
            end
            if D.sup > t(end)
                t = [t;D.sup]; % Add right domain boundary.
            end
        end
        n_t = length(t)-1;
    else % pointwise evaluation
        n_t = length(t);  
        T(t_idx) = []; % Delete the elements stored in t from T to prevent double evaluations in pointwise plot mode. 
    end
    for k = 1:n_t % Loop over all values t(k) in t.
        switch mode
            case 0 % pointwise evaluation
                X(:,idx) = intval(t(k),t(k),'infsup'); % Set column idx of X to constant point interval t(k).
            case 1 % interval evaluation
                X(:,idx) = intval(t(k),t(k+1),'infsup');
        end
        for j = 1:S_a(2)  % Compute images of all polynomial sets a(i,j), 1<=j<=S_a(2)<=3.            
            Y(j,:) = iv2intval(evaltaylormodel(a(i,j),intval2iv(X))); % Y(j) := (image + error) of Taylor model a(i,j) on subdomain X.
        end % j
        if S_a(2) == 1 % 2d-plot (t,a(i,1))
            y_min = min(inf_(Y(1,:)));
            y_max = max(sup(Y(1,:)));
            switch mode
                case 0  
                    T_plt = [T_plt;t(k)];        % Append t to plot values in x-direction.
                    Y_min = [Y_min,y_min];
                    Y_max = [Y_max,y_max];
                case 1   
                    T_plt = [T_plt;t(k);t(k+1)]; % Append t to plot values in x-direction.
                    Y_min = [Y_min,y_min,y_min];
                    Y_max = [Y_max,y_max,y_max];       
            end                        
        else  % 2d-phase plot (a(i,1),a(i,2)) or 3d-phase plot (a(i,1),a(i,2),a(i,3))
            if length(YY) > 1E6
                % plotintval(YY,[],color,fillcolor);     
                plotintval(YY,color,[],fillcolor);     
                YY = [];
            else
                YY = [YY,Y];
            end
        end % S_a
    end % k
end % i

if S_a(2) == 1
    switch mode
        case 0
            plot(T_plt,Y_min,'*r',T_plt,Y_max,'*g'); % pointwise evaluation => point plot
        case 1
            plot(T_plt,Y_min,'r',T_plt,Y_max,'g');   % interval evaluation => line plot
    end
elseif ~isempty(YY)
   % plotintval(YY,[],color,fillcolor);  
   plotintval(YY,color,[],fillcolor);  
end

end % function plottaylormodel