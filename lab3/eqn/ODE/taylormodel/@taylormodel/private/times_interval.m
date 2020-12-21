function res = times_interval(a_iv,a_im,b_iv,b_im,rest_im)
%TIMES_INTERVAL fast computation of the error interval of the product of two Taylor models.
%
%   res = times_interval(a_iv,a_im,b_iv,b_im,rest_im)

% written  01/18/16     F. Buenger
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures 

e = 1e-30;
if 1+e > 1 % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards in verified mode
end

a_iv_inf = a_iv.inf;
a_iv_sup = a_iv.sup;

a_im_inf = a_im.inf;
a_im_sup = a_im.sup;

b_iv_inf = b_iv.inf;
b_iv_sup = b_iv.sup;

b_im_inf = b_im.inf;
b_im_sup = b_im.sup;

rest_im_inf = rest_im.inf;
rest_im_sup = rest_im.sup;

% I_ab = rest.image + a_interval .* (b_.image + b_.interval) + a_image .* b_.interval;
% I_ba = rest.image + b_.interval .* (a_image + a_interval) + b_.image .* a_interval;
% r_.interval = intersect(I_ab,I_ba); % error interval r_.interval of a_ * b_ according to p.107,108 of [E]

% Compute interval bounds for I_ab := rest.image + a_interval .* (b_.image + b_.interval) + a_image .* b_.interval;
% Recall that rounding is upwards.
s_inf = -((-b_im_inf) - b_iv_inf);  % Compute bounds of b_im + b_iv.
s_sup = b_im_sup + b_iv_sup; 

p1_inf = min( [-((-a_iv_inf)*s_inf),-((-a_iv_inf)*s_sup),-((-a_iv_sup)*s_inf),-((-a_iv_sup)*s_sup)] );
p1_sup = max([a_iv_inf*s_inf,a_iv_inf*s_sup,a_iv_sup*s_inf,a_iv_sup*s_sup]);

p2_inf = min([-((-a_im_inf)*b_iv_inf),-((-a_im_inf)*b_iv_sup),-((-a_im_sup)*b_iv_inf),-((-a_im_sup)*b_iv_sup)]);
p2_sup = max([a_im_inf*b_iv_inf,a_im_inf*b_iv_sup,a_im_sup*b_iv_inf,a_im_sup*b_iv_sup]);

I_ab_inf = -( (-rest_im_inf) - p1_inf - p2_inf ); 
I_ab_sup = rest_im_sup + p1_sup + p2_sup;

% Analogously compute interval bounds for I_ba := rest.image + b_.interval .* (a_image + a_interval) + b_.image .* a_interval;

s_inf = -((-a_im_inf) - a_iv_inf); % Compute bounds of a_im + a_iv.
s_sup = a_im_sup + a_iv_sup; 

p1_inf = min( [-((-b_iv_inf)*s_inf),-((-b_iv_inf)*s_sup),-((-b_iv_sup)*s_inf),-((-b_iv_sup)*s_sup)] );
p1_sup = max([b_iv_inf*s_inf,b_iv_inf*s_sup,b_iv_sup*s_inf,b_iv_sup*s_sup]);

p2_inf = min([-((-b_im_inf)*a_iv_inf),-((-b_im_inf)*a_iv_sup),-((-b_im_sup)*a_iv_inf),-((-b_im_sup)*a_iv_sup)]);
p2_sup = max([b_im_inf*a_iv_inf,b_im_inf*a_iv_sup,b_im_sup*a_iv_inf,b_im_sup*a_iv_sup]);

I_ba_inf = -( (-rest_im_inf) - p1_inf - p2_inf ); 
I_ba_sup = rest_im_sup + p1_sup + p2_sup;

% res := r_.interval := intersect(I_ab,I_ba) , the error interval of a_ * b_ according to p.107,108 of [E].
res.inf = max(I_ab_inf,I_ba_inf);
res.sup = min(I_ab_sup,I_ba_sup);

if rndold ~= 1
    setround(rndold)
end
end % function times_interval