@echo off
if exist x.f del x.f
if exist groupie.exe del groupie.exe
copy groupie.f + endfio.f + scratchb.f + timer.f + seconds.f   x.f
lf95 -out groupie.exe x.f
if exist groupie.exe copy groupie.exe ..\exe\ 
del x.*
if exist groupie.map del groupie.map