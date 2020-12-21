function res = sqr_interval(rest_im,a_iv,a_im)
%SQR_INTERVAL fast computation of the error interval of the square of a Taylor model.
%
%   res = sqr_interval(rest_im,a_iv,a_im)

% written  01/04/16     F. Buenger
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures 

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

rest_im_inf = rest_im.inf;
rest_im_sup = rest_im.sup;

a_iv_inf = a_iv.inf;
a_iv_sup = a_iv.sup;

a_im_inf = a_im.inf;
a_im_sup = a_im.sup;

% The following is a fast implementation of
% I := rest_im + (a_iv).^2 + 2*a_im.*a_iv;

% Recall that rounding is upwards. 

% Compute [s1_inf,s1_sup] := (a_iv).^2 
if (a_iv_inf <= 0) && (a_iv_sup >= 0) % a_iv contains 0 
    s1_inf = 0;
else
    s1_inf = min( -((-a_iv_inf)*a_iv_inf) , -((-a_iv_sup)*a_iv_sup) ) ;
end
s1_sup = max(a_iv_inf^2,a_iv_sup^2);

% Compute [s2_inf,s2_sup] := 2*a_im.*a_iv  
s2_inf = 2 * min( [-((-a_im_inf).*a_iv_inf), ... % Note that no rounding error is caused by multiplication with 2.
                   -((-a_im_sup).*a_iv_sup), ... 
                   -((-a_im_inf).*a_iv_sup), ...
                   -((-a_im_sup).*a_iv_inf)] );
           
s2_sup = 2 * max( [a_im_inf.*a_iv_inf,...        % Note that no rounding error is caused by multiplication with 2.
                   a_im_sup.*a_iv_sup,... 
                   a_im_inf.*a_iv_sup,... 
                   a_im_sup.*a_iv_inf] );

res.inf = -( (-rest_im_inf) - s1_inf - s2_inf ) ;  
res.sup = rest_im_sup + s1_sup + s2_sup;

if rndold ~= 1
    setround(rndold)
end

end % function sqr_interval