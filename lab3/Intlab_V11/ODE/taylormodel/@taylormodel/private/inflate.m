function res = inflate(a,inclusion_)
%INFLATE        epsilon inflation of the interval a 
% 
%   res = inflate(a) 

% written  11/11/15     F. Buenger
% modified 02/11/16     F. Buenger  "intval"-components --> intval-like structures
% modified 04/25/16     F. Buenger  new parameter inclusion_


% The values of eps_rel and eps_abs are taken from Riot, function DGLSolver::inflate, file dglalg.cpp.
% See also [E], p.61, line 2.

%  eps_rel = 1E-2;   % constant for relative epsilon inflation
%  eps_abs = 1E-6;   % constant for absolute epsilon inflation

% eps_rel = 1.2E-1;    
% eps_rel_incl = 1E-1; % Those components where inclusion was already established, are inflated with a smaller relative factor.
% 
% eps_abs = 4.2E-16;
% eps_abs_incl = 2E-16;


%eps_rel = 1E-2; % Lotka-Volterra optimal for shrink wrapping
eps_rel = 1E-1;   % In AWA file awa.p we have: 
eps_abs = 2E-16;  %
                  %   457|  VFAKT:= 0.1;
                  %   462|  write('Enter error tolerances E_a, E_r : '); read(infile,E_a,E_r);
                  %   501|  UTAU:= BLOW( UTAU, VFAKT );
                  %   502|  FOR I:= 1 TO n DO  utau[i]:= utau[i] + intval( -(E_a+E_r), E_a+E_r );
                  % 
                  % The absolute and relative error tolerances E_a and E_r 
                  % are user specified and read from the awa-input file "infile".
                  % Typical values are E_a = E_r = 1e-16, see the examples in [E], Chap.5, p.135 ff.
                  % Thus E_a+E_r = 2E-16 so that
                  %                  
                  %   eps_rel := VFAKT = 0.1 
                  %   eps_abs := 2E-16  
                  %
                  % are comman AWA-like settings for the epsilon-inflation.                   
                  
% AWA/RIOT-style-epsilon-inflation, 
%   Riot: function DGLSolver::inflate, file dglalg.cpp
%   AWA : function BLOW, pxsc-3.6.2, file ari/i_ari.p

% res = (1+eps_rel) * a - eps_rel * a + intval(-eps_abs,eps_abs,'infsup'); 

% Rump-style-epsilon-inflation, see S.M. Rump. A Note on Epsilon-Inflation. 
% Reliable Computing, Volume 4, Issue 4, pp. 371-375, 1998.

% res = midrad(1,eps_rel).*a + midrad(0,eps_abs);

f_rel.inf = -((-1)+eps_rel); % f_rel.inf := 1-eps_rel; would also be o.k., since eps_rel >> eps('double') and inflation is heuristic anyway.
f_rel.sup = 1+eps_rel;

f_abs.inf = -eps_abs;
f_abs.sup = eps_abs;

res = iv_plus(iv_times(f_rel,a),f_abs); % interval evaluation of f_rel * a + f_abs

% % Alternative implementation using input argument inclusion_. Only for testing.
%
% f_rel.sup = ~inclusion_ + eps_rel;
% f_rel.sup(inclusion_) = 1 + eps_rel_incl;
% 
% f_rel.inf = ~inclusion_ - eps_abs;
% f_rel.inf(inclusion_) = 1 - eps_abs_incl;
% 
% f_abs.sup = ~inclusion_ + eps_abs;
% f_abs.sup(inclusion_) = 1 + eps_abs_incl;
% 
% f_abs.inf = ~inclusion_ - eps_abs;
% f_abs.inf(inclusion_) = 1 - eps_abs_incl;
%
% res = iv_plus(iv_times(f_rel,a),f_abs); % interval evaluation of f_rel * a + f_abs


end % function inflate