function cleardata(param)
%CLEAR        store/retrieve global data
%   param  0  store data
%          1  retrieve data
%

% written  05/19/14     S.M. Rump
% modified 09/12/18     S.M. Rump  ODE toolbox added
%

  % store INTLAB_CONST at a safe place - out of reach of clear all
  global INTLAB_CONST
  global INTLAB_AWA_VARS
  global INTLAB_ODE_OPTIONS
  global INTLAB_ODE_VARS
  global INTLAB_ODE_TCOEFF
  var = { 'CONST','AWA_VARS','ODE_OPTIONS','ODE_VARS','ODE_TCOEFF' };
  if param
    for i=1:length(var)
      v = var{i};
      eval(['INTLAB_' v ' = getappdata(0,''INTLAB_' v ''');'])
    end
    %     INTLAB_CONST = getappdata(0,'INTLAB_CONST');
  else
    %     setappdata(0,'INTLAB_CONST',INTLAB_CONST);
    for i=1:length(var)
      v = var{i};
      eval(['setappdata(0,''INTLAB_' v ''',INTLAB_' v ');'])
    end

  end
  