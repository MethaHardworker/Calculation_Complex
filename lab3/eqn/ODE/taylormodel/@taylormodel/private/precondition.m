function [yl,yr] = precondition(a,b,base_tm)
%PRECONDITION   preconditioning of Taylor models
%
%   res = precondition(a,base_tm,t,type)
%
%   input: 
%       a:       Integrated current left Taylor model.
%                (Its time dependence is already eliminated by evaluating at t = t_next.)
%       b:       Current right Taylor model
%       base_tm: "identity" Taylor model matrix       
%       type:    type of preconditioning. At the moment two types, namely QR preconditioning and PE (parallel epiped) preconditioning are implemented:
%                   1: QR preconditioning. The matrix A for the linear part of "a" is decomposed by a reordered QR factorization. 
%                      The reordering is simply the ordering of the columns of A in descending order of Euclidean length.
%                      The left Taylor model gets linear part Q, while R sticks at the right Taylor model.     
%                   2: PE preconditioning. The left Taylor model simply gets linear part A of "a".     
% 
%   output: 
%       yl: left Taylor model for the next time step
%       yr: right Taylor model for the next time step
% 
%  The implementation is based on: 
%
%  [NJN] M. Neher, K.R. Jackson, N.S. Nedialkov, "On Taylor model based integration of ODEs", 
%          SIAM J. Numer. Anal. 45(1), pp. 236-262, 2007
%  [MB]  Kyoko Makino and Martin Berz, "Suppression of the wrapping effect by Taylor model - based validated integrators",
%          MSU HEP Report 40910, 2003 

% written  05/15/17     F. Buenger

e = 1e-30;
if 1+e > 1      % fast check for rounding upwards
    rndold = 1;
else
    rndold = getround;
    setround(1) % rounding upwards 
end

global INTLAB_ODE_OPTIONS 

% constants for preconditioning types
PREC_QR  = 1;  % QR preconditioning
PREC_PE  = 2;  % parallelepiped preconditioning
           
n = a(1).dim-1; % order of the underlying ODE system
A = get_linear_terms(a); % real nx(n+1)-matrix of the linear part of "a"
A = A(:,1:n); % The (n+1)-th column of A corresponds to the linear part of the time variable
              % which stands at the last position (n+1) for each component of "a".
              % This column should be zero since "a" should not depend on the time variable.
              % This final column is cut of so that A becomes an nxn-matrix.

switch INTLAB_ODE_OPTIONS.precondition 
    case PREC_QR %  QR preconditioning      
        % Compute a QR decomposition of As := A*P, where P is a permutation matrix such that
        % the columns of As are sorted by Euclidean length in descending order.
        % This is step (i) of Algorithm 6.1, p. 257, [NJN], see also [MB], Definition 12 (QR preconditioning), p.28.
        
        normsA = sqrt(sum(A.*A));  % non-verified computation of the Euclidean norms of the column vectors of A
        [normsA,iA] = sort(normsA,'descend');  % The variable "normsA" is not used later on. Only the permutation iA is of interest.
    
        As = A(:,iA);          % The columns of As are sorted by Euclidean length in descending order.
        [Q,R] = qr(As);        % Compute a non-verified QR decomposition As = Q*R 
        dummy = sign(diag(R));
        dummy(dummy == 0) = 1;
        Q = Q.*dummy';         % The columns of Q are multiplied with the non-zero signs of the diagonal of R. 
                               % In other words, Q is chosen such that R has non-negative diagonal.
                               % Only the orthogonal part Q is used later on, the triangular part R is not explicitly used.                                   
    case PREC_PE % parallelepiped preconditioning        
        if INTLAB_ODE_OPTIONS.blunting
            Q = blunt(A);
        else
            cond_bound = 1e6; % heuristic condition number bound at which blunting is automatically used, even though it is not explicitly switched on by INTLAB_ODE_OPTIONS.blunting = true.
                              % Feel free to change that value!
            cond_A = cond(A);
            if cond_A > cond_bound
                Q = blunt(A);
            else
                Q = A;                
            end
        end
end  
                
Q_inv = intval2iv( inv(intval(Q)) ); % Compute verified inverse of Q.
                
% All but the constant term of the integrated left Taylor model "a" is shifted to the right Taylor model. 
% Q will become the linear part of the left Taylor model and Q^{-1} is applied to the right Taylor model.
% This is step (ii) of Algorithm 6.1, [NJN].

% Additional explanations to this step: 
%   Let a0 := a - a(0) be "all but the constant term of the integrated left Taylor model", 
%   where a(0) is the constant term of the integrated left Taylor model "a". 
%   Then the intermediate setting of the right Taylor model "yr" in this step is done as follows:  
%   
%       (*)  yr = Q^{-1} * a0 o b =  (Q^{-1} * a0) o b  ,  note that b has range in [-1,1] componentwise.
% 
%   It is numerically important(!) that the evaluation is done in the following order
%  
%       1. z  = Q^{-1} * a0
%       2. yr = z o b
%
%   The mathematically equivalent evaluation of (*) yr := Q^{-1} * (a0 o b), i.e.:
%
%       1. z  = a0 o b
%       2. yr = Q^{-1}*z
%
%   is numerically instable.

[c,idx] = get_constant_term(a);     % c := constant terms of a, idx stores the corresponding monomial indices.
a0 = subtract_constant_term(a,idx); % a0 := a - (constant terms of "a")

if INTLAB_ODE_OPTIONS.shrinkwrap
    b = shrinkwrap(b,base_tm);      % Try shrinkwrapping in order to get rid off the remainder interval.
end

yr = concatenate(Q_inv*a0,b);       % Compute yr = Q^{-1} * a0 o b =  (Q^{-1} * a0) o b , note that b has range in [-1,1] componentwise.
%yr = Q_inv*concatenate(a0,b);      %            = Q^{-1} * (a0 o b), mathematically equivalent expression but numerically instable

% Bound the range of the new right Taylor model.
% This is step (iii) of Algorithm 6.1, [NJN].
I = iv_plus(get_image(yr),get_interval(yr));
[m,s] = iv_getmidrad(I);
s_min = 1E-16; % minimum scaling constant. The value is chosen quite arbitrarily small (near zero). Feel free to change it!
               % For stability of rescaling with 1./s, each component of a scaling s is bounded from below by s_min.                                    
s = max(s_min,s,'includenan'); % For stability reasons of scaling with 1./s, each component of s is bounded from below. 

% Apply a scaling matrix S^-1 on the right Taylor model such that afterwards 
% each component is contained in [-1,1] and spans [-1,1] approximately. 
% This is step (iv) of Algorithm 6.1 in [NJN].
yr = (yr-m)./s; 

% Apply the inverse scaling S to the new left Taylor model (so that scaling of left and right Taylor models cancels over all: "S*S^-1 = I").
% This is step (v) of Algorithm 6.1, [NJN].
yl = Q*(s.*base_tm + m) + c;    

if rndold ~= 1 
    setround(rndold)
end

end  % function precondition
