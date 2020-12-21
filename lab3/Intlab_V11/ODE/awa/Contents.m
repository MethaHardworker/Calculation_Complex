%INTLAB-AWA
%
%ODE-solver 
%  awa          - main function, ODE-solver
%  awaTC        - kernel of main function awa (for performance reasons 
%                 placed in local area of type taylorcoeff)  
%  awaset       - set global parameters used by awa,  
%                 analogous to MATLAB function odeset
%
%Taylor coefficient constructors (type taylorcoeff) 
%  taylorcoeffinit   - Initialization of dependent variables
%  taylorcoeff       - constructor for type taylorcoeff
%
%Display of Taylor coefficient (rigorous)
%  display      - Command window display of Taylor coefficients
%  disp         - Display function for pop-up windows in debugger
%  awa_disp     - Display of awa results
%
%Taylor coefficient arithmetic operations
%  plus         - Plus                              +
%  uplus        - Unary plus                        +
%  minus        - Minus                             -
%  uminus       - Unary minus                       -
%  mtimes       - Matrix multiply                   *
%  times        - Elementwise multiply              .*
%  mrdivide     - Slash or right division           /
%  mldivide     - Backslash or left division        \
%  rdivide      - Elementwise right division        ./
%  ldivide      - Elementwise left division         .\
%  mpower       - Matrix power                      ^
%  power        - Elementwise power                 .^
%
%Other Taylor coefficient operations
%  sum          - Sum
%  prod         - Product
%  inv          - Inversion
%
%Taylor coefficient trigonometric functions
%  sin          - Sine
%  cos          - Cosine
%  tan          - Tangent
%  cot          - Cotangent
%  sec          - Secant
%  csc          - Cosecant
%  asin         - Inverse sine
%  acos         - Inverse cosine
%  atan         - Inverse tangent
%  acot         - Inverse cotangent
%  asec         - Inverse secant
%  acsc         - Inverse cosecant
%  sinh         - Hyperbolic sine
%  cosh         - Hyperbolic cosine
%  tanh         - Hyperbolic tangent
%  coth         - Hyperbolic cotangent
%  asinh        - Inverse hyperbolic sine
%  acosh        - Inverse hyperbolic cosine
%  atanh        - Inverse hyperbolic tangent
%  acoth        - Inverse hyperbolic cotangent
%
%Taylor coefficient exponential functions (rigorous, real and complex)
%  exp          - Exponential
%  log          - Natural logarithm
%  log10        - Logarithm to base 10
%  sqr          - Square
%  sqrt         - Square root
%
%Verification routines and auxiliary
%  awaget       - get awa-parameters (mainly for internal use)   
%  typeadjust   - Type adjustment
%  
%Demonstration, samples
%  dawa   - Some examples for using INTLAB-AWA
%
% INTLAB-AWA (see main function awa.m) is a MATLAB/INTLAB-implementation 
% of the well-known verified ODE-solver AWA written in Pascal by 
% Rudolf Lohner, Institute for Applied Mathematics, Univ. of Karlsruhe. 
%
% For details see:
%
% [L]   R. Lohner, Einschliessung der Loesung gewoehnlicher Anfangs- und
%       Randwertaufgaben und Anwendungen, Diss. Univ. Karlsruhe, 1988
%
% [AWA] R. Lohner, Pascal implementation according to [L], available from 
%       http://www2.math.uni-wuppertal.de/~xsc/xsc/pxsc_software.html#awa 
%
% INTLAB-AWA was written by Florian Buenger,
% Institute for Reliable Computing, Hamburg University of Technology.

% written  03/29/18     F. Buenger
