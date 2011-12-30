@echo off
REM ##############################################################
REM # NOTE for this to work you must have the following 
REM # Environment Variables defined:
REM #   - MSVSHOME : path to MS Visual Studio 
REM #   - TDLDIR = /path/to/tdl
REM ##############################################################

REM ##############################################################
REM # Add vc stuff to path
REM # note clear path first since if rerun this file path gets big!
path = ''
path = %path%;%MSVSHOME%\\Vc7
path = %path%;%MSVSHOME%\\Vc7\\bin
path = %path%;%MSVSHOME%\\Common7\\IDE
REM ##############################################################

REM ##############################################################
REM # call vc bat to setup compiler variables 
REM # MSVSHOME\\Vc7\\bin\\vcvars32.bat
call vcvars32.bat
REM ##############################################################

REM ##############################################################
REM # Set up environment variables used in makefiles
REM # The important ones are 'COMPILER', 'COMPFLAGS', 
REM # 'INCLUDE' and 'LIB'
set TOOLS= %MSVSHOME%\\Vc7
set TOOLS2= %MSVSHOME%\\Vc7\\PlatformSDK
set COMPILER= %TOOLS%\\bin\\cl.exe
set COMPFLAGS= /Z7 /Od /c /nologo
set INCLUDE32= -I%TOOLS%\\include -I%TOOLS2%\\include
REM 
set LIB1=%TOOLS%\\lib
set LIB2=%TOOLS2%\\lib
set LIB= %LIB1%; %LIB2%
REM 
set GSLINCLUDE= -IC:\\user\\trainor\\work\\code\\gsl\\GSL\\gsl_1.3\\include
set TDLINCLUDE= -I%TDLDIR%\\lib\\src %INCLUDE32% %GSLINCLUDE% 
REM ##############################################################

REM ##############################################################
REM # run make

NMAKE /F makefile.msvs %1

REM ##############################################################
