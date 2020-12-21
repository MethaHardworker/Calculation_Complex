function restore
  % clean up function for sudden death or error
  global INTLAB_CONST
  INTLAB_CONST.RealStdFctsExcptnIgnore = 0;
  progress(-1);
end  % function restore
