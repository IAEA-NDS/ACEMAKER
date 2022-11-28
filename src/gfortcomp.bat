@ echo off
if exist sixlin.exe   del sixlin.exe
if exist gamlin.exe   del gamlin.exe
if exist doace.exe    del doace.exe
if exist dotsl.exe    del dotsl.exe
if exist dodos.exe    del dodos.exe
if exist dophn.exe    del dophn.exe
if exist acemaker.exe del acemaker.exe
gfortran -o sixlin.exe   sixlin.f
gfortran -o gamlin.exe   gamlin.f
gfortran -o doace.exe    doace.f
gfortran -o dotsl.exe    dotsl.f
gfortran -o dodos.exe    dodos.f
gfortran -o dophn.exe    dophn.f
gfortran -o acemaker.exe acemaker.f
if exist sixlin.exe   copy sixlin.exe   ..\exe\
if exist gamlin.exe   copy gamlin.exe   ..\exe\
if exist doace.exe    copy doace.exe    ..\exe\
if exist dotsl.exe    copy dotsl.exe    ..\exe\
if exist dodos.exe    copy dodos.exe    ..\exe\
if exist dophn.exe    copy dophn.exe    ..\exe\
if exist acemaker.exe copy acemaker.exe ..\exe\
if exist *.obj del *.obj
if exist *.map del *.map
if exist *.o   del *.o
if exist a     del a
 