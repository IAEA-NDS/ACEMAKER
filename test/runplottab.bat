@echo off
rem
rem usage plot.inp plot.cur plot.pnt
rem
set plottabd=d:\Plottab\
rem usage plot.inp plot.cur
if exist PLOTTAB.INP del PLOTTAB.INP
if exist PLOTTAB.CUR del PLOTTAB.CUR
if exist PLOTTAB.PNT del PLOTTAB.PNT
if exist PLOT*.ps    del PLOT*.ps
if exist PLOTTAB.LST del PLOTTAB.LST
copy %plottabd%PLOT.CHR
copy %plottabd%PAGE.DAT
copy %plottabd%PLOT.SYM
copy %1 PLOTTAB.INP
copy %2 PLOTTAB.CUR
if not %3.==. copy %3 PLOTTAB.PNT
%plottabd%plottab
%plottabd%plotsave
if exist PLOT0001.ps del PLOT0001.ps
if exist PLOT*.ps copy PLOT*.ps %2.ps
if exist PLOTTAB.INP del PLOTTAB.INP
if exist PLOTTAB.CUR del PLOTTAB.CUR
if exist PLOTTAB.PNT del PLOTTAB.PNT
if exist PLOT*.ps    del PLOT*.ps
if exist PLOTTAB.LST del PLOTTAB.LST
if exist PLOT.CHR del PLOT.CHR
if exist PAGE.DAT del PAGE.DAT
if exist PLOT.SYM del PLOT.SYM
ps2pdf %2.ps %2.pdf
set plottabd=