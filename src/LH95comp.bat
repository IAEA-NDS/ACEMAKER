@ echo off
if exist sixlin.exe   del sixlin.exe
if exist doace.exe    del doace.exe
if exist acemaker.exe del acemaker.exe
lf95 sixlin.f
lf95 doace.f
lf95 acemaker.f
if exist sixlin.exe   copy sixlin.exe ..\exe\
if exist doace.exe    copy doace.exe  ..\exe\
if exist acemaker.exe copy acemaker.exe ..\exe\
if exist *.obj del *.obj
if exist *.map del *.map