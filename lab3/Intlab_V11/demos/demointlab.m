function demointlab
%DEMOINTLAB   Wrapper routine to call INTLAB demos
%
%A selection of INTLAB demos, call
%
%  demointlab
%

% written  10/13/12     S.M. Rump  
% modified 11/07/12     S.M. Rump  INTLAB_larger added
% modified 04/04/14     S.M. Rump  end function
% modified 05/13/14     S.M. Rump  new demos
% modified 01/21/15     S.M. Rump  error if no web available
% modified 09/03/18     S.M. Rump  new order and update
%

  if ~exist('web','file')
    error(['Command "web" not supported. ', ...
           'To see INTLAB demo-files please display the corresponding ' ...
           'html-files in directory demos/html directly by some browser.     '])
  end    

  d = which('demointlab');
  wng = warning;
  warning off
  addpath([ d(1:end-13) '\html' ])
  warning(wng)
  
  clc
  disp('Welcome to INTLAB, the Matlab toolbox for reliable computing.')
  disp(' ')
  disp('The current version consists of more than 1700 .m-functions with more ')
  disp('  than 57 thousand lines of Matlab-code (more than 99 KLOC with comments). ')
  disp('The test suite for INTLAB consists of another 102 KLOC. ')
  disp(' ')
  while 1
    name = displaycomments;
    str = input('select demo ','s');
    if str=='0'
      break
    end
    try
      web([name{index(upper(str))} '.html'])
    catch
      disp('invalid input')
    end
  end
  
  disp(' ')
  disp('Enjoy INTLAB. Comments and suggestions always welcome to rump (at) tuhh.de .')
  disp(' ')
  
end  % function demointlab
  
function name = displaycomments
  disp(' ')
  disp('This is a wrapper routine to call several INTLAB demos, selected by numbers. ')
  disp(' ')
  num = '0';
  name = {};
  num = succ(num); disp([num '  A general demo of some features of INTLAB']); name{index(num)} = 'dintlab';
  num = succ(num); disp([num '  Details about interval arithmetic']); name{index(num)} = 'darithmetic';
  num = succ(num); disp([num '  Some examples of interval computations']); name{index(num)} = 'dintval';
  num = succ(num); disp([num '  Some larger examples with INTLAB']); name{index(num)} = 'dintlab_larger';
  num = succ(num); disp([num '  Accurate real and complex interval standard functions in INTLAB']); name{index(num)} = 'dstdfcts';
  num = succ(num); disp([num '  ODEs: The AWA toolbox']); name{index(num)} = 'dawa';
  num = succ(num); disp([num '  ODEs: The Taylor model toolbox']); name{index(num)} = 'dtaylormodel';
  num = succ(num); disp([num '  Global (un-)constrained optimization and all roots of nonlinear functions']); name{index(num)} = 'dglobal';
  num = succ(num); disp([num '  The polynomial toolbox (univariate and multivariate polynomials']); name{index(num)} = 'dpolynom';
  num = succ(num); disp([num '  Gradients: automatic differentiation of multivariate functions']); name{index(num)} = 'dgradient';
  num = succ(num); disp([num '  Hessians: automatic differentiation with second derivative']); name{index(num)} = 'dhessian';
  num = succ(num); disp([num '  Taylor series: automatic Taylor coefficients']); name{index(num)} = 'dtaylor';  
  num = succ(num); disp([num '  Affine interval arithmetic']); name{index(num)} = 'daffari';
  num = succ(num); disp([num '  Slopes: automatic slope generation']); name{index(num)} = 'dslope';
  num = succ(num); disp([num '  Utility routines']); name{index(num)} = 'dutility';
  num = succ(num); disp([num '  fl-numbers: k-bit point and interval arithmetic']); name{index(num)} = 'dfl';
  num = succ(num); disp([num '  Long numbers: a non-optimal multiple precision package']); name{index(num)} = 'dlong';
  num = succ(num); disp([num '  Accurate summation and dot products']); name{index(num)} = 'daccsumdot';
  disp(' ')
  disp('0  exit this wrapper')
  disp(' ')
end  % function displaycomments

function num = succ(num)
% next in numbering
  if isequal(num,'9')
    num = 'A';
  else
    num = char(num+1);
  end
end  % function succ

function i = index(num)
% convert num into index
  if ( '0' <= num ) & ( num <= '9' )
    i = num - '0';    
  else
    i = num - 'A' + 10;
  end
end  % function index

function ch = upper(ch)
% switch to upper case if letter
  if ( 'a' <= ch ) & ( ch <= 'z' )
    ch = ch - 'a' + 'A';
  end
end  % function upper
