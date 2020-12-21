function progress(str,width_window)
% display progress
%
% Call:
%   progress(input,width_window)
%
% input  string        string to be displayed after some time
%        real >=0      start display after mintime [sec] (inf means no display at all)
%        real <0       close popup window
%
%        width_window  (optional) default 600 pixels
%

% written  07/15/18  S.M. Rump
%

  global INTLAB_CONST
  persistent old_text
  
  if ~isstr(str)
    if str<0
      try
        close(INTLAB_CONST.POPUPWINDOW)
      end
      INTLAB_CONST.POPUPWINDOW = [];
      INTLAB_CONST.SHOWPOPUPWINDOW = true;
    else
      INTLAB_CONST.MINTIME = str;
    end
    return
  end
  
  if isempty(INTLAB_CONST.POPUPWINDOW)
    INTLAB_CONST.SHOWPOPUPWINDOW = true;
    screen = get(0,'ScreenSize');
    width = screen(3);
    height = screen(4);
    if nargin==1
      width_window = 600;
    end
    INTLAB_CONST.POPUPWINDOW = figure('CloseRequestFcn',@closereq);
    set(INTLAB_CONST.POPUPWINDOW, ...
      'Name','busy ...   [to avoid popup, call progress(inf)]', ...
      'MenuBar','none', ...
      'ToolBar', 'none', ...
      'NumberTitle','off', ...
      'Interruptible', 'off', ...
      'Position',[width-width_window-50 height-100 width_window 30])
    t = text;
    old_text = t;
    t.String = str;
    t.FontSize = 12;
    t.Position = [0.02 0.6];
    axis off
  else
    if INTLAB_CONST.SHOWPOPUPWINDOW
      figure(INTLAB_CONST.POPUPWINDOW);
      old_text.String = str;
%       drawnow
    end
  end
end  % progress


function closereq(src,callbackdata)
% do not stop program when popup window is closed
  global INTLAB_CONST
  INTLAB_CONST.SHOWPOPUPWINDOW = false;
  delete(INTLAB_CONST.POPUPWINDOW);
end  % function closereq
