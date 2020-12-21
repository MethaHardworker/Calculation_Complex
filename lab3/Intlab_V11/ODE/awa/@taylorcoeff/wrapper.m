function [varargout] = wrapper(dummy,fname,vargin)
%WRAPPER  This wrapper is only used by test programs for calling 
%         private awa functions from outside. 
%
%   [varargout] =  wrapper(dummy,fname,vargin)   
% 
%   dummy:   Taylor coefficient for accessing the wrapper
%   fname:   name of the private function to be called
%   vargin:  cell array containing the input for function "fname"
%            (For some unknown reasons in Matlab 2016b varargin{:} returns a cell array
%            and not a comma-separated list. Therefore, we use our own cell array "vargin".) 

% written  05/11/17     F. Buenger

    [varargout{1:nargout}] = feval(fname,vargin{:});
    
end % function wrapper

