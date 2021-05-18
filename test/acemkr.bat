@echo off
REM batch file to run acemaker
REM
REM usage: runacemk filename
REM        filename: full name of acemaker input file
REM        Default=acemaker.inp
REM 
REM acemkd : full path to acemaker executables
REM
set acemkd=\ACEMAKER\exe\
if %1.==. goto run
copy %1 acemaker.inp
:run
%acemkd%acemaker
