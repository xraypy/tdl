@echo off

REM ******************* INSTALL PATH  *****************************
REM below is the install directory 

set TDL_DIR=C:\\temp\\tdl
set TDL=tdl

REM ***************************************************************


REM ********************* PY PATHS  ********************************
REM below are the paths to search for modules
set TDL_PATH=%TDL_DIR%
set TDL_PATH=%TDL_PATH%;%TDL_DIR%\\%TDL%
set TDL_PATH=%TDL_PATH%;%TDL_DIR%\\modules
set TDL_PATH=%TDL_PATH%;%TDL_DIR%\\modules\\IO
set TDL_PATH=%TDL_PATH%;%TDL_DIR%\\modules\\GUI

set TDL_INI_PATH=%TDL_DIR%

set WX_RSRC=%TDL_DIR%\\modules\\GUI\\tdl_wxGUI.rsrc.py

REM ******************** RUN IT  **********************************

set PYTHONPATH=%TDL_PATH%

if "%1" =="" goto :default

python tdl.py %* 
goto :end

:default
python tdl.py 
goto :end


REM ***************************************************************

:end
set TDL_DIR=
set TDL_PATH=
set TDL_INI_PATH=

