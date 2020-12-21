function [t,y] = awa_disp(T,Y,Z,t,print,fmt)
%AWA_DISP  computes and displays the inclusions of an ODE at specified time grid points.             
%
%   [t,y] = awa_disp(T,Y,t,print,fmt)
%                       
% input parameters: 
%
%   T,Y,Z: output of a function call [T,Y,Z] = awa(...)
%   T_ : optional parameter. If T_ is not empty, then inclusions and 
%        diameters at all time points T_(k) are displayed. Only entries of 
%        T_ are considered that are contained in the integration interval 
%        [T(1),T(end)], respectively [T(end),T(1)] if integration is done 
%        in negative direction. All these time points are returned in t
%        sorted in integration direction order. 
%
%        If T_ is empty, then simply the inclusions Y and their diameters 
%        at all, by awa automatically determined time grid points T are 
%        displayed. In this case t := T. 
% 
%   print: optional parameter. Inclusion intervals and their diameters are 
%          displayed. Default = true;
%
%   fmt: optional parameter. Display format of the enclosure bounds.
%        Example: fmt = 'longe' 
%   
% output parameters:
% 
%   t,y :  [y.inf(j,i),y.sup(j,i)] is an enclosure of u_i(t(j)) where   
%          u_i is the i-th component of the exact solution u of the ode.
%          

% written  05/24/17     F. Buenger
% modified 04/09/18     F. Buenger, additional output parameter Z for continuous solution inclusion
% modified 08/29/18     F. Buenger, INTLAB display formats

global INTLAB_CONST

fmt_old = getformat;
if nargin >= 6 && ~isempty(fmt)
    eval(['format ' fmt]);
end

if nargin  < 5 || isempty(print) 
    print = true;
end

n = length(Y.inf(1,:)); % length of initial conditions = dimension of ODE system
direction = sign(T(end)-T(1));

if nargin  < 4 || isempty(t)
    t = T;
    y = Y;
    if print
      y_ = Y; 
    end
else  
    Z = intval2iv(Z);
    t = t(:); % Convert t to column vector.
    if isempty(Z)
        error('empty Taylor coefficient array')
    else        
        if direction > 0
            t = t( and(t >= T(1),t <= T(end)) ) ;
            t = sort(t);
        else
            t = t( and(t >= T(end),t <= T(1)) ) ;
            t = sort(t,'descend');            
        end
        y_ = intval(zeros(length(t),n)); % just preallocation
        for k = 1:length(t)
            tk = t(k);
            if direction > 0
                j = find(tk >= T,1,'last');
            else
                j = find(tk <= T,1,'last');                
            end
            Tj = T(j);
            if tk == Tj % t(k) is a grid point at which the enclosure is already given in Y.
                y_(k,:) = Y(j,:);
            else        % t(k) is an intermediate point contained in (T(j),T(j+1)) if direction > 0, respectively in (T(j+1),T(j)) if direction < 0.
                Tj1 = T(j+1);
                x = (intval(tk)-Tj) / (intval(Tj1)-Tj);  
                for i = 1:n
                    p_.inf = flipud(squeeze(Z.inf(j,i,:))); 
                    p_.sup = flipud(squeeze(Z.sup(j,i,:))); 
                    p = polynom(iv2intval(p_)); % generalized Taylor polynomial. Its last coefficient encloses the remainder term on the whole interval [T(j),T(j+1)].  
                    y_(k,i) = polyval(p,x);     % y(t(k)) is contained in p(x).    
                end
            end
        end
        y = y_;        
    end
end

if print
    y_ = y_';
    switch INTLAB_CONST.INTVAL_DISPLAY
        case 'DisplayMidrad'
            y_str = midrad(y_(:),[],[]);
        otherwise
            y_str = infsup(y_(:),[],[]);
    end
    y_str_exp = strtrim(y_str.exp);    
    y_diam = diam(y_(:));
    K = length(t);
    i_ = []; 
    
    for k = 1:K
        t_str = strtrim(evalc('disp(t(k))'));
        % Note that t_str might be rounded!
        % Thus, the "display version" t_str of t(k) may not
        % exactly represent the floating-point number t(k).
        % We explicitly did not(!) want to display t(k)
        % artificially as an interval in such a case.
        
        disp(['t = ',t_str]);        % Display time at which an inclusion shall be printed.
        % Display inclusions and their diameters of inclusions.        
        for i = 1:n
            idx = (k-1)*n + i;
            if n > 1 
                i_ = ['_',num2str(i)];
            end
            y_string = strtrim(y_str.str(idx,:));
            if ~isempty(y_str_exp)
                y_string = [y_str_exp,' ', y_string];
            end
            switch INTLAB_CONST.INTVAL_DISPLAY
                case 'DisplayMidrad'
                    y_str_diam = []; % In midrad notation the diameter must not bedisplayed separately. 
                otherwise
                    y_str_diam = ['    d([y',i_,']) = ',num2str(y_diam(idx),'%.2e')];                    
            end
            disp(['    [y',i_,'] = ',y_string,y_str_diam]);            
        end
        disp(' ')
    end
end

eval(['format ' fmt_old]);
end  % function awa_disp