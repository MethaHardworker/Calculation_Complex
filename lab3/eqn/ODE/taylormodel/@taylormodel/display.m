function display(a,name,newline)
%DISPLAY  command window display of Taylor models
%    
%   display(a,name,newline)
%
%  Single Taylor models are displayed with the following notation:
% 
%      dim   :  a.dim, number of variables of the polynomial part
%      order :  a.order, maximum degree of the polynomial part
%      type  :  a.type  0 : general Taylor model (default), 
%                       1 : ODE - Taylor model 
%
%      <iv_mid,iv_rad> : error interval a.interval in mid-/rad-notation
%      [im_inf,im_sup] : image interval a.image in infsup-notation 
%      [min,max]       : interval vector a.domain in infsup-notation 
%      center          : center point a.center of the polynomial part
%      coeff           : coefficient vector a.coefficent of the polynomial part.
%                        The corresponding monomials (exponents of variables)
%                        are displayed to the left.  
%
%  Matrices of Taylor models are displayed in a MATLAB-generic, 
%  more technical way.
%  

% written  08/27/15     F. Buenger
% modified 11/19/15     F. Buenger  display arrays of Taylor models
% modified 02/10/16     F. Buenger  "intval"-components --> intval-like structures 

if nargin < 2 || isempty(name)
    name = inputname(1);
    if isempty(name)
        name = 'ans';
    end
end
if nargin < 3 || isempty(newline)
    newline = true;
end

space = '  ';
disp(['taylormodel ' name ' = ' ])

if size(a,1) == 1 && size(a,2) == 1 % special display for single Taylor model
    
    % Display dimension, order, and type.
    n = a.dim;   
%     T = array2table([n a.order a.type mid(a.interval) rad(a.interval) inf(a.image) sup(a.image)],...
%         'VariableNames',{'dim' 'order' 'type' 'iv_mid' 'iv_rad' 'im_inf' 'im_sup'});
    [a_iv_mid,a_iv_rad] = iv_getmidrad(a.interval);   
    M = cell(1,7);
    M{1} = uint32(n);
    M{2} = uint32(a.order);
    M{3} = uint32(a.type);
    M{4} = a_iv_mid;
    M{5} = a_iv_rad;
    M{6} = a.image.inf;
    M{7} = a.image.sup;
%     M = [num2cell(uint32(n)),...
%          num2cell(uint32(a.order)),...
%          num2cell(uint32(a.type)),...
%          num2cell(a_iv_mid),...
%          num2cell(a_iv_rad),...
%          num2cell(a.image.inf),...
%          num2cell(a.image.sup) ];
    X = {'dim' 'order' 'type' 'iv_mid' 'iv_rad' 'im_inf' 'im_sup'}; 
    T = cell2table(M,'VariableNames',X); 
    disp(T);
    
    % Display domains and evaluation points.
    if newline 
        disp(' '); 
    end; % newline, but not for disp(a)
    
    M = [a.domain.inf, a.domain.sup,a.center];
    if a.type == 1 % Taylor models for ODEs get variable names y1,...,ym,t, m:=n-1.
        if n == 2
            X = {'x' 't'};
        else
            X = cell(1,n);
            for i = 1:n-1
                X{i} = ['x',num2str(i)];
            end;
            X{n} = 't'; % name for time variable
        end
    else % General Taylor models get variable names x1,...,xn.        
       fmt = lower(matlab.internal.display.format);
       % Formats short and long are transformed by disp() for type table
       % to shortg and longg which does not look very nice since the decimal 
       % points are not aligned. Therefore we switch to shorte and longe, 
       % respectively, in these cases for displaying the domains and centers. 
        switch fmt
            case {'short'}
                format 'shorte';
            case {'long'}
                format 'longe';
        end
        if n == 1
            X = {'x'};
        else
            X = cell(1,n);
            for i = 1:n
                X{i} = ['x',num2str(i)];
            end
        end
    end
    T = array2table(M,'VariableNames',{'min' 'max' 'center'},...
                      'RowNames',X');
    disp(T);
    
    % Display monomials and coefficients.
    if newline 
        disp(' ');
    end; % newline, but not for disp(a)
    fmt = lower(matlab.internal.display.format);
    % Formats short and long are transformed by disp() for type table
    % to shortg and longg which does not look very nice since the decimal 
    % points are not alligned. Therefore we switch to shorte and longe, 
    % respectively, in these cases for displaying coefficients. 
    switch fmt
        case {'short'}
            format 'shorte';
        case {'long'}
            format 'longe';
    end
    if isfield(a.coefficient,'inf') && isfield(a.coefficient,'sup')
        M = [num2cell(uint32(a.monomial)) , num2cell(a.coefficient.inf) , num2cell(a.coefficient.sup)];
        T = cell2table(M,'VariableNames',[X,{'coeff_inf'},{'coeff_sup'}]);
    else
        M = [num2cell(uint32(a.monomial)) , num2cell(a.coefficient)];
        T = cell2table(M,'VariableNames',[X,{'coeff'}]);
    end
    disp(T);
    format (fmt);
       
else % generic, more technical display for non-trivial Taylor model vector/matrix
    s.type = '.';
    F = fields(a);
    for i = 1:length(F)
        s.subs = F{i};
        disp([space,F{i},' : ']);
        disp(subsref_(a,s)); % a.(fname) is not wanted here.
        % 'subsref_' already performs certain type conversions
        % which are more suitable for displaying.
    end
end

end  % function display

  
  