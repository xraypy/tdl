@echo off
REM ******************* INSTALL PATH  *****************************REM below is the install directory, point to dir with tdl script
set TDL_DIR=C:\\user\\trainor\\work\\code\\tdl_src\\tdl
if "%1" =="" goto :default
python tdl %* 
goto :end

:default
python tdl 
goto :end

:end
set TDL_DIR=
