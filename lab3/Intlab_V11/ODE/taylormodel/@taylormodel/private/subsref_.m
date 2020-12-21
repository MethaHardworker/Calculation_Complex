function r = subsref_(a,s)
%SUBSREF  implements subscripted references for Taylor models
%
%   r = subsref_(a,s)

% written  08/27/15     F. Buenger
% modified 11/19/15     F. Buenger  subsref for matrix input
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures

istaylormodel_a = true;
len_s = length(s);

while len_s
    if istaylormodel_a
        if strcmp(s(1).type,'()')     % index reference a(i)
            r = a(s(1).subs{:});
        elseif strcmp(s(1).type,'{}') % reference a{i}
            r = a{s(1).subs{:}};
            %   error('reference a{idx} not implemented for taylormodel')
        elseif strcmp(s(1).type,'.')  % index reference a.t
            fname =  s(1).subs;
            if strcmp(fname,'center') || strcmp(fname,'domain') || strcmp(fname,'coefficient') || strcmp(fname,'monomial')
                S = size(a);
                if length(S) > 2
                    error('maximally two dimensions for type taylormodel.');
                end
                if S(1) == 1 && S(2) == 1 % In the one-dimensional case the field value is directly returned.
                    r = a.(fname);
                else                      % In the multidimensional case each field value is stored in a cell of a cell array.
                    r = cell(S);
                    for i = 1:S(1)        % The double loop seems not to be avoidable for the present implementation of the type 'taylormodel'.
                        for j = 1:S(2)
                            r{i,j} = a(i,j).(fname); % Note that a(i,j).(fname) is (in general) an array and for distinct pairs (i,j)
                                                     % and (i',j') the array sizes of a(i,j).(fname) and a(i',j').(fname)
                                                     % might be different.
                        end
                    end
                end
            elseif strcmp(fname,'order') || strcmp(fname,'dim') || strcmp(fname,'type')
                S = size(a);
                if length(S) > 2
                    error('maximally two dimensions for type taylormodel.');
                end
                r = zeros(S);
                for i = 1:S(1) % The double loop seems not to be avoidable for the present implementation of the type 'taylormodel'.
                    for j = 1:S(2)
                        x =  a(i,j).(fname);
                        if ~isempty(x)
                            r(i,j) = x; % Scalar field values are stored in an array.
                        else
                            r = x;
                            return;
                        end
                    end
                end
            elseif strcmp(fname,'interval') || strcmp(fname,'image')
                S = size(a);
                if length(S) > 2
                    error('maximally two dimensions for type taylormodel.');
                end
                r_inf = zeros(S);
                r_sup = r_inf;
                for i = 1:S(1) % The double loop seems not to be avoidable for the present implementation of the type 'taylormodel'.
                    for j = 1:S(2)
                        x = a(i,j).(fname);
                        %[x_inf,x_sup] = getinfsup(x);
                        % faster implementation than the previous line in comments
                        
                        if ~(isempty(x.inf) || isempty(x.sup)); % faster than isempty(x)
                            r_inf(i,j) = x.inf; % Scalar field values are stored in an array.
                            r_sup(i,j) = x.sup; % Scalar field values are stored in an array.
                        else
                            r = x; % x is an empty interval.
                            return;
                        end
                    end
                end
                r.inf = r_inf;
                r.sup = r_sup;
            else
                error('invalid subscript reference for taylormodel')
            end
        else
            error('invalid index reference for taylormodel')
        end
        
    else  % for a.inf(i) etc.
        r = subsref(a,s(1));        
    end
    len_s = len_s-1;
    
    if len_s
        s = s(2:end);
        a = r;
        clar r;
        istaylormodel_a = isa(a,'taylormodel');
    end
end % function subsref_
