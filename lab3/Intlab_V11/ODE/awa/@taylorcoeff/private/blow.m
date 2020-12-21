function r = blow(a,eps_rel,eps_abs)
%BLOW  epsilon inflation of the interval a. 
% 
%   r = blow(a,eps_rel,eps_abs) 
%
% Corresponds to (but slightly different from) Pascal-XSC function BLOW, 
% see folder pxsc-3.6.2, file ari/i_ari.p:
% 
%       GLOBAL FUNCTION BLOW ( VAR A : INTERVAL; EPS : REAL ) : INTERVAL;
% 
%       VAR B : INTERVAL;
% 
%       BEGIN
%           B:= (1.0+EPS)*A - EPS*A;
%           BLOW.INF:= PRED( B.INF );
%           BLOW.SUP:= SUCC( B.SUP );
%       END
% 
%  AWA uses EPS := VFAKT = 0.1.

% written  07/17/17     F. Buenger

% AWA-style epsilon inflation
r = iv_plus( iv_minus(iv_times(1.0+eps_rel,a) , iv_times(eps_rel,a)) , iv_midrad(0,eps_abs)); % r = (1+eps_rel)*a-eps_rel*a + [-eps_abs,eps_abs];

% Rump-style epsilon inflation
%r = iv_plus( iv_times(iv_midrad(1,eps_rel),a) , iv_midrad(0,eps_abs) ); % r = [1-eps_rel,1+eps_rel]*a + [-eps_abs,eps_abs];

end % function blow




