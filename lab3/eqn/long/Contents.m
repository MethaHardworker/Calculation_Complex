%INTLAB long toolbox (rudimentary, slow but correct)
%
%Long constructors
%  long         - Long constructor
%  long2dble    - Long to double
%  long2intval  - Long to intval
%  vertcat      - Vertical concatenation            [ ; ]
%  subsasgn     - Subscripted assignment A(2:4) = 1
%  subsref      - Subscripted reference r = A(3)
%
%Display of long numbers
%  display      - Command window display for long
%  disp         - Display function for pop-up windows in debugger
%
%Long arithmetic operations
%  plus         - Plus                              +
%  uplus        - Unary plus                        +
%  minus        - Minus                             -
%  uminus       - Unary minus                       -
%  mtimes       - Elementwise multiply              *
%  times        - Elementwise multiply              .*
%  mrdivide     - Elementwise right division        /
%  mldivide     - Elementwise left division         \
%  rdivide      - Elementwise right division        ./
%  ldivide      - Elementwise left division         .\
%  mpower       - Elementwise power                 ^
%  power        - Elementwise power                 .^
%  longshift    - Shift by r bits
%
%Other long operations
%  min          - Minimum
%  max          - Maximum
%  abs          - Absolute value
%  mid          - Midpoint
%  rad          - Radius
%  diam         - Diameter
%  inf          - Infimum
%  sup          - Supremum
%
%Long comparisons
%  eq           - Equal                             ==
%  ne           - Not equal                         ~=
%  le           - Less equal                        <=
%  lt           - Less than                         <
%  ge           - Greater equal                     >=
%  gt           - Greater than                      >
%
%Utility routines
%
%  isempty      - Long is empty in Matlab sense, i.e. []
%  isnan        - True for Not a Number
%  numel        - Number of elements in array (default 1)
%  numels       - Number of elements in array
%  end          - Last index
%
%Structural operations
%  length       - Length
%  size         - Size
%
%Other long operations
%  longprecision- Sets/gets current precision
%  addlongerror - Add error to long number
%
%Some sample functions
%  exp          - Exponential
%  longpi       - Long computation of Pi
%
%Initialization of INTLAB long package and system variables
%  longinit     - Initialization and definition of defaults
%
%Demonstration, samples
%  demolong     - Some examples for using INTLAB long package
%

% written  12/30/98     S.M. Rump
%
% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
