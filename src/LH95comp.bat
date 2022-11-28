@ echo off
if exist sixlin.exe   del sixlin.exe
if exist gamlin.exe   del gamlin.exe
if exist doace.exe    del doace.exe
if exist dotsl.exe    del dotsl.exe
if exist dodos.exe    del dodos.exe
if exist dophn.exe    del dophn.exe
if exist acemaker.exe del acemaker.exe
lf95 sixlin.f
lf95 gamlin.f
lf95 doace.f
lf95 dotsl.f
lf95 dodos.f
lf95 dophn.f
lf95 acemaker.f
if exist sixlin.exe   copy sixlin.exe   ..\exe\
if exist gamlin.exe   copy gamlin.exe   ..\exe\
if exist doace.exe    copy doace.exe    ..\exe\
if exist dotsl.exe    copy dotsl.exe    ..\exe\
if exist dodos.exe    copy dodos.exe    ..\exe\
if exist dophn.exe    copy dophn.exe    ..\exe\
if exist acemaker.exe copy acemaker.exe ..\exe\
if exist *.obj del *.obj
if exist *.map del *.map