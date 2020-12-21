function coeff = taylorcoeffinit(coeff)
%TAYLORCOEFFININT  Initialization of generalized Taylor coefficients
%
%  r = taylorcoeffinit(coeff)
%
% It is assumed that coeff is an at most 2-dimensional structure array with 
% components "inf" and "sup" such that coeff(i,j).inf and coeff(i,j).sup are 
% column vectors of equal size k, independent of i, j.
% The interval [coeff(i,j).inf(s),coeff(i,j).sup(s)] represents an 
% inclusion of the (s-1)-th generalized Taylor coefficient of component (i,j).

% written  12/11/17     F. Buenger

global INTLAB_AWA_VARS % global parameters 

if nargin == 0    % taylorcoeffinit should only be called by startintlab.m without parameters  
                  % in order to initialize the subsequent global variables once.                  
                  
    INTLAB_AWA_VARS.STATUS = [];
    INTLAB_AWA_VARS.TREENR = [];
    INTLAB_AWA_VARS.VERTEXNR = [];
    INTLAB_AWA_VARS.TREE = {[];[];[];[]};
    return;
end

coeff = taylorcoeff(coeff ,'taylorcoeffinit');

end % function taylorcoeffinit