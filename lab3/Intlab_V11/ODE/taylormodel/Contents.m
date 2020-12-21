%INTLAB Taylor model toolbox
%
%ODE-solver 
%  verifyode    - main function, ODE-solver
%  verifyodeTM  - kernel of main function verifyode (for performance reasons 
%                 placed in local area of type taylormodel)  
%  verifyodeset - set global parameters used by verifyode,  
%                 analogous to MATLAB function odeset

%Taylor model constructors (type taylormodel)
%  taylormodelinit   - Initialization 
%  taylormodel       - Constructor for type taylormodel
%
%Display of Taylor models (rigorous)
%  display            - Command window display of Taylor models
%  disp               - Display function for pop-up windows in debugger
%  verifyode_disp     - Display function for verifyode results
%  plottaylormodel    - visualization of Taylor models 
%                       2D- and 3D-phase plots are possible
%  plottaylormodel_1d - 1D-visualization of Taylor models
%
%Taylor model arithmetic operations
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
%  inv          - Elementwise inverse               .^-1
%
%Other Taylor model operations
%  sum          - Sum
%  prod         - Product
%  inv          - Elementwise inversion
%
%Taylor model trigonometric functions (rigorous, real and complex)
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
%Taylor model exponential functions 
%  exp          - Exponential
%  log          - Natural logarithm
%  log10        - Logarithm to base 10
%  sqr          - Square
%  sqrt         - Square root
%
%Taylor model comparison 
%  eq           - Equal                             ==
%  ne           - Not equal                         ~=
%  
%ODE utility routines
%  evaltaylormodel  - Verified evaluation of Taylor models on subdomains
%  evaltaylormodel2 - Verified evaluation  of Taylor models on subdomains 
%                     (alternative implementation to evaltaylormodel)
%  image            - Verified enclosure of the image of the polynomial
%                     part of a Taylor model
%  tie              - Concatenation of left and right Taylor models
%  typeadjust       - Type adjustment
%  verifyodeget     - Get verifyode options (mainly for internal use) 
%
%Demonstration, samples
%  dtaylormodel    - Some examples for using the INTLAB Taylor model toolbox
%
%The Taylor model toolbox is mainly designed for solving ODEs in a 
% verified manner. The main function for that purpose is verifyode.
% 
%The implementation is based on: 
%  
% [E]  I. Eble, "Über Taylor-Modelle", Dissertation at Karlsruhe Institute of Technology, 2007 (written in German),
%        Riot, C++-implementation, http://www.math.kit.edu/ianm1/~ingo.eble/de
% [M]  K. Makino, "Rigorous analysis of nonlinear motion in particle accelerators", 
%        Dissertation at Michigan State University, 1998
% [MB] K. Makino and M. Berz, "Suppression of the wrapping effect by Taylor model - based validated integrators",
%        MSU HEP Report 40910, 2003 
% [NJN] M. Neher, K.R. Jackson, N.S. Nedialkov, "On Taylor model based integration of ODEs", 
%        SIAM J. Numer. Anal. 45(1), pp. 236-262, 2007
% [Bue] F. Buenger, "Shrink wrapping for Taylor models revisited", 
%        Numerical Algorithms 78(4), pp. 1001-1017, 2018
%
% The INTLAB Taylor model toolbox was written by Florian Buenger,
% Institute for Reliable Computing, Hamburg University of Technology.

% written  06/04/18     F. Buenger
