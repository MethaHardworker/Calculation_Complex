function a = uminus(a)
%UMINUS  Taylor coefficient unary minus  - a
%
%   r = uminus(a)
%
% Corresponds to Lohner's function TKOEFF, case OPER = UMINUS, file awa.p.

% written  08/01/17     F. Buenger
% modified 01/18/18     F. Buenger  vectorized version added

numel_a = numel(a);
numel_bound = 30;   % numel(a) < numel_bound, then the non-vectorized uminus is faster than vectorized uminus. Feel free to change this heuristic value!
trigger = (numel_a < numel_bound);
%trigger = true;    % only for testing !!!
%trigger = false;   % only for testing !!!

if trigger    
    for i = 1:numel_a
        x = a(i).inf;
        a(i).inf = -a(i).sup;
        a(i).sup = -x;
    end
else 
    a_ = tcoeff2tc(a);
    x = a_.inf;
    a_.inf = -a_.sup;
    a_.sup = -x;
    a = tc2tcoeff(a_,a(1));
end
    
    
end % function uminus