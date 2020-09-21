      program acemaker
c     version 1.0
c
c     Driver program to produce ACE-formatted files for Monte Carlo
c     calculations.
c
c     The code uses the modules LINEAR, RECENT, SIGMA1, LEGEND, SPECTRA,
c     FIXUP, GROUPIE, MERGER and DICTIN from the PREPRO-2019 code
c     package. The module SIXLIN is used to linearize and process MF6
c     data. The code DOACE produces the ACE-formatted file.
c
c     In general the driver invokes the following calculation sequence:
c     LINEAR + RECENT + LEGEND + SPECTRA + SIXLIN +
c     + SIGMA1 + FIXUP + GROUPIE + MERGER + DICTIN + DOACE
c
c     Configuration file: acemaker.cfg
c
c     By default the fullpath for the PREPRO and ACEMAKER executables
c     is \ACEMAKER\exe\ saved on the internal variables preprod and
c     acemakerd.
c
c     The configuration file  acemaker.cfg can be created on the
c     working directory to redefine the fullpaths for PREPRO and
c     ACEMAKER. Two lines are needed: the first line specifies the
c     PREPRO directory and the second one is for the ACEMAKER
c     executables.
c
c     Example of acemaker.cfg:
c
c     D:\PREPRO19\
c     D:\ACEMAKER\exe\
c
c     Note: The use of  '\' or '/' as ending character is compulsory.
c
c     It is also possible to redefine fullpaths changing the values
c     of preprod and acemakerd in the corresponding DATA statement on
c     the main program and recompiling ACEMAKER.
c
c     Input file: acemaker.inp
c
c     The input options for running ACEMAKER should be entered on the
c     acemaker.inp file. It is fixed inside the code. acemaker.inp is a
c     text file containing the input options given by keywords.The first
c     line should be the title of the case. It is followed by a set of
c     keywords, one by line.  The keywords should be entered in the
c     first 6 character positions. The data fields start in column 7,
c     if require. Some keywords could require additional lines. The
c     keywords could be written in lower- or uppercase, but not in a
c     combination. The last keyword should be the card END or end.
c
c     The keyword set of ACEMAKER is described below:
c
c     Keywrd    Data-field(7-72)        Comments
c     =================================================================
c
c     endf      full-filename       ENDF input data filename (A66)
c     ENDF
c
c     mat       NMAT                Number of material
c     MAT                           Should be followed by NMAT lines
c                                   containing material numbers (MAT),
c                                   one by line
c
c     temp      NTEMP               Number of temperatures
c     TEMP                          Should be followed by NTEMP lines
c                                   containing the temperatures in
c                                   Kelvin, one by line
c
c     tol       EPS                 Maximum fraction tolerance for
c     TOL                           linearization/reconstruction of data
c
c     ymin      FYMIN               Minimum allowable value of
c     YMIN                          cross section or distribution
c
c     dxmu      FMUMAX              Maximum cosine interval(mu) for
c     DXMU                          angular distribution reconstruction
c
c     suff      .XX                 ZAID suffix for the ACE formatted
c     SUFF                          file (ZAID.XXc). Should start with
c                                   '.', followed by two digits (XX)
c
c     ptab      IPTAB               Probability table trigger in the URR
c     PTAB                          IPTAB = 0/1 = no/yes
c
c     pndf      IPNDF               trigger to keep the PENDF tape
c     PNDF                          IPNDF = 0/1 = no/yes
c
c     mon       IMON                Monitor printing trigger
c     MON                           IMON = 0/1 = min./max. print out
c
c     keep      IKEEP               Trigger to keep intermediate files
c     KEEP                          IKEEP = 0/1 = no/keep
c
c     end                           End of case
c     END
c
c     The default input options are equivalent to the following
c     acemaker.inp file:
c
c     ACEMAKER default title (default acemaker.inp)
c     ENDF   ENDFB.IN
c     MAT    0
c     TEMP   1
c      293.6
c     TOL    0.001
c     YMIN   1.0D-30
c     DXMU   0.001
c     SUFF   .00
c     PTAB   1
c     PNDF   0
c     MON    0
c     KEEP   0
c     END
c
c     The keyword:
c     MAT    0
c     in the third input line above means the first material on the
c     input ENDF tape (ENDFB.IN in this case)
c
c     The keywords can be entered in any order. The END keyword
c     plays the role of case terminator. If a keyword is missing, the
c     default value is applied. So, it is only necessary to include
c     the keywords that change defaults. Several cases can be stacked in
c     one acemaker.inp file by appending after the current END a new
c     title line and new input options finishing with a new END keyword.
c     ACEMAKER will stop after the last END. An example is given below:
c
c     Case 1 title
c     ENDF   FENDLEN.DAT
c     MAT    3
c      5728
c      725
c      825
c     SUFF   .32
c     END
c     Case 2 title
c     ENDF   \ENDF-B-VIII.0\n-001_H_001.endf
c     TEMP   2
c      293.6
c      600.0
c     SUFF   .80
c     END
c     Case 3 title
c     ENDF   \ENDF\ENDFBVIII.endf
c     MAT    2
c      9228
c      9237
c     TEMP   2
c      293.6
c      600.0
c     SUFF   .80
c     END
c
c     In case 1 above, evaluated data will be retrieved from FENDLEN.DAT
c     for 3 materials with MAT number 5728(La-139), 725(N-14) and 825
c     (O-16). The ZAID on the ACE formatted files will have the suffix
c     .32c for each material. Default values apply to the rest of input
c     options.
c     In case 2, evaluated nuclear data are read for H-1 from
c     \ENDF-B-VIII.0\n-001_H_001.endf (NMAT=0, wich means first material
c     on the ENDF tape). Two temperatures are requested, then two ace-
c     formatted files will be generated: one at 293.6 K and other at
c     600.0 K. The ZAID number will be 1001.80c and 1001.81c. Default
c     values are used for the rest of input options.
c     In Case 3 evaluated nuclear data are retrieved from tape
c     \ENDF\ENDFBVIII.endf for U-235 (9228) and U-238 (9237). Two
c     temperatures are requested for each material. Four ACE-formatted
c     files are generated with ZAID numbers:
c      92235.80c at 293.6K and 92235.81c at 600.0K for U-235 (9228)
c      92238.80c at 293.6K and 92238.81c at 600.0K for U-238 (9237)
c     Default values are taken for the rest of input options.
c
c     Each ACE-formatted file generated is named according to the
c     pattern ZAzzzaaa.XXc.acef where zzzaaa is the ZA number printed
c     in a field of six positions and filled with zeroes in place of
c     spaces. It follows the GROUPIE-2019 conventions for 2-band data
c     filenames. The symbol .XX is the ZAID suffix entered as input for
c     the first temperature and it is eventually increased in .01 for
c     each additional temperature(NTEMP>1).
c
c     Filenames for the cases previously presented:
c
c     Case 1
c     ZA057139.32c.acef  for  La-139
c     ZA007014.32c.acef  for  N-14
c     ZA008016.32c.acef  for  O-16
c     Case 2
c     ZA001001.80c.acef  and  ZA001001.81c.acef  for H-1
c     Case 3
c     ZA092235.80c.acef  and  ZA092235.81c.acef  for U-235
c     ZA092238.80c.acef  and  ZA092238.81c.acef  for U-238
c
c     The output files ZAzzzaaa.XXc.xsd and ZAzzzaaa.XXc.lst are also
c     prepared. The first file contains the directory information for
c     the XSDIR file of MCNP and the second one is a compilation of the
c     listing files produced by all the modules called by ACEMAKER. The
c     files acemaker.log summaries ACEMAKER processing.
c
      implicit real*8 (a-h,o-z)
      parameter (nmatx=300, ntempx=10)
      dimension mat(nmatx), temp(ntempx)
      character*80 cmd
      character*72 fendf
      character*66 title
      character*62 preprod,acemakerd
      character*11 ctime
      character*10 cdate
      character*8  fza
      character*5  str5
      character*3  suff
      real*4 xsuff
      data nin/20/,nou/21/,nlst/22/
      data preprod/'\ACEMAKER\exe\'/,acemakerd/'\ACEMAKER\exe\'/
      data emax/1.0d+09/,fangmin/1.0d-10/,flegmin/1.0d-02/
c
c      Open acemaker log file
c
      call getdtime(cdate,ctime)
      open(nou,file='acemaker.log')
      write(nou,'(a)')'========================================='
      write(nou,'(a,1x,a,1x,a)')'ACEMAKER log file ',cdate,ctime
      write(nou,'(a)')'========================================='
      write(*,'(1x,a)')'========================================='
      write(*,'(1x,a,1x,a,1x,a)')'ACEMAKER log file ',cdate,ctime
      write(*,'(1x,a)')'========================================='
c
c      Get fullpath for PREPRO and ACEMAKER
c
      open(nin,file='acemaker.cfg',iostat=icfg,err=5,status='old')
    5 if (icfg.eq.0) then
        read(nin,'(a)')preprod
        read(nin,'(a)')acemakerd
        write(nou,'(a)')'Read filepaths will be used: '
        write(*,*)'Read filepaths will be used: '
      else
        write(nou,'(a)')'Default filepaths will be used: '
        write(*,*)'Default filepaths will be used: '
      endif
      close(nin)
      write(nou,'(a,a)')'PREPRO   filepath: ',trim(preprod)
      write(nou,'(a,a)')'ACEMAKER filepath: ',trim(acemakerd)
      write(*,'(1x,a,a)')'PREPRO   filepath: ',trim(preprod)
      write(*,'(1x,a,a)')'ACEMAKER filepath: ',trim(acemakerd)
c
c      Get input options for first case
c
      ic=1
      jc=0
      call getinp(ic,title,fendf,nmat,mat,ntemp,temp,tol,ymin,dxmu,
     &            suff,iptab,ipndf,imon,ikeep)
c
c      Input case cycle
c
      do while (ic.ge.0)
c
c      Write input options for case jc
c
      jc=jc+1
      write(nou,*)
      write(nou,'(a6,i4)')'Case: ',jc
      write(nou,'(a10)')'=========='
      write(nou,'(a7,a)')'Title: ',trim(title)
      write(nou,*)
      write(nou,'(a17,a)')'ENDF input tape: ',trim(fendf)
      write(nou,*)
      write(nou,'(a33,1pe12.5)')'Minimum allowable cross section: ',ymin
      write(nou,'(a25,8x,1pe12.5)')'Linearization tolerance: ',tol
      write(nou,'(a25,8x,1pe12.5)')'Maximum cosine interval: ',dxmu
      write(nou,'(a33,i2)')'Probability tables 0/1 = no/yes: ',iptab
      write(nou,'(a28,5x,i2)')'Monitor option 0/1 = no/yes:',imon
      write(nou,'(a29,4x,i2)')'Keep PENDF tape 0/1 = no/yes:',ipndf
      write(nou,'(a28,5x,i2)')'Keep all tapes 0/1 = no/yes:',ikeep
      write(nou,'(a30,4x,a3)')'ACE-formatted file ID suffix: ',suff
      write(nou,*)
      write(nou,'(a21,i6)')'Number of materials: ',nmat
      write(nou,'(a20)')'===================='
      write(nou,'(a5,a12)')'  No.','    MATERIAL'
      do i=1,nmat
        write(nou,'(i4,i11)')i,mat(i)
      enddo
      write(nou,*)
      write(nou,'(a24,i3)')'Number of temperatures: ',ntemp
      write(nou,'(a23)')'======================='
      write(nou,'(a5,a16)')'  No.',' TEMPERATURE [K]'
      do i=1,ntemp
        write(nou,'(i4,1pe13.5)')i,temp(i)
      enddo
c
c      Deleting intermediate files
c
      call delfile('merger.inp')
      call delfile('merger.lst')
      call delfile('linear.inp')
      call delfile('linear.lst')
      call delfile('recent.inp')
      call delfile('recent.lst')
      call delfile('legend.inp')
      call delfile('legend.lst')
      call delfile('legend.tmp')
      call delfile('spectra.inp')
      call delfile('spectra.lst')
      call delfile('sixlin.inp')
      call delfile('sixlin.lst')
      call delfile('sixlin.tmp')
      call delfile('sigma1.inp')
      call delfile('sigma1.lst')
      call delfile('fixup.inp')
      call delfile('fixup.lst')
      call delfile('groupie.inp')
      call delfile('groupie.lst')
      call delfile('dictin.inp')
      call delfile('dictin.lst')
      call delfile('doace.inp')
      call delfile('doace.lst')
      call delfile('endf6.endf')
      call delfile('linear.pendf')
      call delfile('recent.pendf')
      call delfile('legend.pendf')
      call delfile('spectra.pendf')
      call delfile('sixlin.pendf')
      call delfile('sigma1.pendf')
      call delfile('fixup.pendf')
      call delfile('groupie.endf')
      call delfile('groupie.pendf')
      call delfile('endf6.pendf')
c
c      Material cycle
c
      do im=1,nmat
        mati=mat(im)
        write(nou,*)
        write(nou,'(a5,i4)')'MAT= ',mati
        write(nou,'(a)')'========='
        write(*,*)
        write(*,'(1x,a5,i4)')'MAT= ',mati
        write(*,'(1x,a)')'========='
c
c       MERGER: get endf-6 formatted data from endf input tape
c
        open(nin,file='merger.inp')
        write(nin,'(a)')'endf6.endf'
        write(nin,'(a19,i4)')'ENDF tape for MAT= ',mati
        write(nin,'(a)')trim(fendf)
        write(nin,'(a3)')'END'
        write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
        write(nin,*)
        close(nin)
        call setcmd(nou,preprod,'merger',cmd)
        call system(cmd)
        call delfile('merger.inp')
        call delfile('merger.lst')
c
c       Get ZA
c
        matj=0
        za=0.0d0
        open(nin,file='endf6.endf',err=10)
        read(nin,*,err=10,end=10)
        read(nin,'(d11.0,55x,i4)',err=10,end=10)za,matj
   10   if (matj.ne.mati) then
          write(nou,'(a,i4,a)')'*** ERROR: MAT= ',mati,' not found'
          write(*,'(1x,a,i4,a)')'*** ERROR: MAT= ',mati,' not found'
          close (nin,err=15)
   15     call delfile('endf6.endf')
          read(nin,'(d11.0,55x,i4)',err=40,end=40)za,matj
        endif
        close(nin)
        nza=za
        fza=''
        write(fza,'(a2,i6)')'ZA',nza
        do ii=3,8
          if (fza(ii:ii).eq.' ') fza(ii:ii)='0'
        enddo
c
c       LINEAR: linearize cross sections at 0.0 K
c
        open(nin,file='linear.inp')
        write(nin,'(2i11,1pd11.4,i11)')0,imon,ymin,1
        write(nin,'(a)')'endf6.endf'
        write(nin,'(a)')'linear.pendf'
        write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
        write(nin,*)
        write(nin,'(1p2d11.4)')0.0,0.5*tol
        write(nin,'(1p2d11.4)')1.0,0.5*tol
        write(nin,'(1p2d11.4)')2.0,tol
        write(nin,'(1p2d11.4)')emax,tol
        write(nin,*)
        close(nin)
        call setcmd(nou,preprod,'linear',cmd)
        call system(cmd)
        if (ikeep.ne.1) then
          call delfile('linear.inp')
          call delfile('endf6.endf')
        endif
c
c       RECENT: reconstruct resonance cross sections at 0.0 K
c
        open(nin,file='recent.inp')
        write(nin,'(i11,1pd11.4,4i11)')0,ymin,1,1,1,imon
        write(nin,'(a)')'linear.pendf'
        write(nin,'(a)')'recent.pendf'
        write(nin,'(2i11)')mati,mati
        write(nin,*)
        write(nin,'(1p2d11.4)')0.0,0.5*tol
        write(nin,'(1p2d11.4)')1.0,0.5*tol
        write(nin,'(1p2d11.4)')2.0,tol
        write(nin,'(1p2d11.4)')emax,tol
        write(nin,*)
        close(nin)
        call setcmd(nou,preprod,'recent',cmd)
        call system(cmd)
        if (ikeep.ne.1) then
          call delfile('recent.inp')
          call delfile('linear.pendf')
        endif
c
c       LEGEND: Process/linearize angular distributions on MF4
c
        open(nin,file='legend.inp')
        write(nin,'(1pd11.4,5i11)')tol,0,2,1,2,0
        write(nin,'(a)')'recent.pendf'
        write(nin,'(a)')'legend.pendf'
        write(nin,'(i6,i2,i3,i6,i2,i3,1p4d11.4)')mati,1,1,mati,99,999,
     &              0.0d0,emax,fangmin,flegmin
        write(nin,*)
        close(nin)
        call setcmd(nou,preprod,'legend',cmd)
        call system(cmd)
        call delfile('legend.tmp')
        if (ikeep.ne.1) then
          call delfile('legend.inp')
          call delfile('recent.pendf')
        endif
c
c       SPECTRA: Process/linearize energy distributions on MF5
c
        open(nin,file='spectra.inp')
        write(nin,'(2i11,1pd11.4,i11)')0,imon,ymin,1
        write(nin,'(a)')'legend.pendf'
        write(nin,'(a)')'spectra.pendf'
        write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
        write(nin,*)
        write(nin,'(1p2d11.4)')0.0,tol
        write(nin,*)
        close(nin)
        call setcmd(nou,preprod,'spectra',cmd)
        call system(cmd)
        if (ikeep.ne.1) then
          call delfile('spectra.inp')
          call delfile('legend.pendf')
        endif
c
c       SIXLIN: Process/linearize angular-energy distributions on MF6
c
        open(nin,file='sixlin.inp')
        write(nin,'(2i11,1pd11.4)')0,imon,ymin
        write(nin,'(a)')'spectra.pendf'
        write(nin,'(a)')'sixlin.pendf'
        write(nin,'(2i11)')mati,mati
        write(nin,'(1p2d11.4)')tol,dxmu
        write(nin,*)
        close(nin)
        call setcmd(nou,acemakerd,'sixlin',cmd)
        call system(cmd)
        if (ikeep.ne.1) then
          call delfile('sixlin.inp')
          call delfile('spectra.pendf')
        endif
c
c       Merge temperature independent listing files
c
        call delfile('temp1.lst')
        open(nlst,file='temp1.lst')
        write(nlst,'(a)')'ACEMAKER listing file'
        write(nlst,'(a)')'====================='
        call cpfile(nlst,'linear.lst')
        call delfile('linear.lst')
        call cpfile(nlst,'recent.lst')
        call delfile('recent.lst')
        call cpfile(nlst,'legend.lst')
        call delfile('legend.lst')
        call cpfile(nlst,'spectra.lst')
        call delfile('spectra.lst')
        call cpfile(nlst,'sixlin.lst')
        call delfile('sixlin.lst')
        close(nlst)
c
c        Temperature cycle
c
        do it=1,ntemp
          tempi=temp(it)
          write(nou,'(a5,i4,a7,1pd11.4)')'MAT= ',mati,' TEMP= ',tempi
          write(*,'(1x,a5,i4,a7,1pd11.4)')'MAT= ',mati,' TEMP= ',tempi
c
c         SIGMA1: Doppler broadening cross sections
c
          open(nin,file='sigma1.inp')
          write(nin,'(2i11,1p2d11.4,2i11)')0,imon,tempi,ymin,1,0
          write(nin,'(a)')'sixlin.pendf'
          write(nin,'(a)')'sigma1.pendf'
          write(nin,'(2i11)')mati,mati
          write(nin,*)
          write(nin,'(1p2d11.4)')0.0,0.5*tol
          write(nin,'(1p2d11.4)')1.0,0.5*tol
          write(nin,'(1p2d11.4)')2.0,tol
          write(nin,'(1p2d11.4)')emax,tol
          write(nin,*)
          close(nin)
          call setcmd(nou,preprod,'sigma1',cmd)
          call system(cmd)
          if (ikeep.ne.1) call delfile('sigma1.inp')
c
c         FIXUP: fix formats and XS. Prepare unified energy grid
c
          open(nin,file='fixup.inp')
          write(nin,'(a14)')'10002111100011'
          write(nin,'(a)')'sigma1.pendf'
          write(nin,'(a)')'fixup.pendf'
          close(nin)
          call setcmd(nou,preprod,'fixup',cmd)
          call system(cmd)
          if (ikeep.ne.1) then
            call delfile('fixup.inp')
            call delfile('sigma1.pendf')
          endif
c
c         Get probability tables
c
          if (iptab.gt.0) then
            write(cmd,'(a8,a)')fza,'.SHIELD.LST'
            call delfile(cmd)
            write(cmd,'(a8,a)')fza,'.UNSHIELD.LST'
            call delfile(cmd)
            write(cmd,'(a8,a)')fza,'.MULTBAND.LST'
            call delfile(cmd)
            write(cmd,'(a8,a)')fza,'.MULTBAND.TAB'
            call delfile(cmd)
            write(cmd,'(a8,a)')fza,'.PLOT.CUR'
            call delfile(cmd)
            write(cmd,'(a8,a)')fza,'.URR.ENDF'
            call delfile(cmd)
c
c           GROUPIE: Generate 2-bands(bins) probability tables in URR
c
            open(nin,file='groupie.inp')
            write(nin,'(4i11,1pd11.4,i11)')0,-11,2,0,1.0d-3,0
            write(nin,'(a)')'fixup.pendf'
            write(nin,'(a)')'groupie.endf'
            write(nin,'(5i11)')1,1,1,1,1
            write(nin,'(a20,i5,a6,1pd11.4)')'GROUPIE 2-BANDS MAT=',mati,
     &                ' TEMP=',tempi
            write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
            write(nin,*)
            close(nin)
            call setcmd(nou,preprod,'groupie',cmd)
            call system(cmd)
            if (ikeep.ne.1) then
              call delfile('groupie.inp')
              call delfile('groupie.endf')
              write(cmd,'(a8,a)')fza,'.SHIELD.LST'
              call delfile(cmd)
              write(cmd,'(a8,a)')fza,'.UNSHIELD.LST'
              call delfile(cmd)
              write(cmd,'(a8,a)')fza,'.MULTBAND.LST'
              call delfile(cmd)
              write(cmd,'(a8,a)')fza,'.MULTBAND.TAB'
              call delfile(cmd)
              write(cmd,'(a8,a)')fza,'.PLOT.CUR'
              call delfile(cmd)
            endif
c
c           Checking fza=ZAzzzaaa.URR.ENDF
c
            write(cmd,'(a8,a)')fza,'.URR.ENDF'
            title=trim(cmd)
            ichk=1
            open(nin,file=title,status='old',err=20)
            read(nin,'(a)',err=20,end=20)cmd
            read(nin,'(a)',err=20,end=20)cmd
            ichk=0
            close(nin)
   20       if (ichk.eq.0) then
c
c             MERGER: Merge sections 2/152 & 2/153
c
              open(nin,file='merger.inp')
              write(nin,'(a)')'groupie.pendf'
              write(nin,'(a16,i4,a7)')'PENDF tape MAT= ',mati,' PTAB-2'
              write(nin,'(a)')trim(title)
              write(nin,'(a)')'fixup.pendf'
              write(nin,'(a3)')'END'
              write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
              write(nin,*)
              close(nin)
              call setcmd(nou,preprod,'merger',cmd)
              call system(cmd)
              if (ikeep.ne.1) then
                call delfile('merger.inp')
                call delfile(title)
              endif
            else
              close(nin)
              write(nou,*)'** Warning: No unresolved probability tables'
              write(*,*)' ** Warning: No unresolved probability tables'
            endif
          endif
c
c         DICTIN: Update dictionary section of pendf tape
c
          open(nin,file='dictin.inp')
          if (iptab.eq.1.and.ichk.eq.0) then
            write(nin,'(a)')'groupie.pendf'
          else
            write(nin,'(a)')'fixup.pendf'
          endif
          write(nin,'(a)')'endf6.pendf'
          close(nin)
          call setcmd(nou,preprod,'dictin',cmd)
          call system(cmd)
          if (ikeep.ne.1) then
            call delfile('dictin.inp')
            call delfile('fixup.pendf')
            if (ichk.eq.0) call delfile('groupie.pendf')
          endif
c
c         DOACE: Prepare ACE-formatted file for MC (ZAaaazzz.xxc.acef)
c
          read(suff,'(f3.0)')xsuff
          xsuff=xsuff+0.01*(it-1)
          write(str5,'(f5.2)')xsuff
          write(cmd,'(a8,a3,a6)')fza,str5(3:5),'c.acef'
          title=trim(cmd)
          call delfile(title)
          open(nin,file='doace.inp')
          write(nin,'(2i11)')0,imon
          write(nin,'(a)')'endf6.pendf'
          write(nin,'(a)')trim(title)
          write(nin,'(i11)')mati
          write(nin,'(a3)')str5(3:5)
          write(nin,*)
          close(nin)
          call setcmd(nou,acemakerd,'doace',cmd)
          call system(cmd)
          if (ikeep.ne.1) call delfile('doace.inp')
c
c         Simple file access checking
c
          ichk=1
          open(nin,file=title,status='old',err=30)
          read(nin,'(a)',err=30,end=30)cmd
          read(nin,'(a)',err=30,end=30)cmd
          ichk=0
          close(nin)
   30     if(ichk.eq.0) then
            write(nou,'(a11,a,a16,i4,a7,1pd11.4)')' ACE-file: ',
     &            trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
            write(*,'(1x,a11,a,a16,i4,a7,1pd11.4)')' ACE-file: ',
     &            trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
          else
            write(nou,'(a24,a,a16,i4,a7,1pd11.4)')
     &            ' Access error ACE-file: ',
     &            trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
            write(*,'(1x,a24,a,a16,i4,a7,1pd11.4)')
     &            ' Access error ACE-file: ',
     &            trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
          endif
c
c         Save pendf tape onto ZAaaazzz.xxc.pendf, if ipndf=1
c
          if (ipndf.gt.0.or.ikeep.eq.1) then
            write(cmd,'(a8,a3,a7)')fza,str5(3:5),'c.pendf'
            title=trim(cmd)
            call delfile(title)
            open(nlst,file=title)
            call cpfile(nlst,'endf6.pendf')
            close(nlst)
          endif
          call delfile('endf6.pendf')
c
c         Save *.lst files onto ZAaaazzz.xxc.lst
c
          write(cmd,'(a8,a3,a5)')fza,str5(3:5),'c.lst'
          title=trim(cmd)
          call delfile(title)
          open(nlst,file=title)
          call cpfile(nlst,'temp1.lst')
          call cpfile(nlst,'sigma1.lst')
          call delfile('sigma1.lst')
          call cpfile(nlst,'fixup.lst')
          call delfile('fixup.lst')
          if (iptab.eq.1) then
            call cpfile(nlst,'groupie.lst')
            call delfile('groupie.lst')
            call cpfile(nlst,'merger.lst')
            call delfile('merger.lst')
          endif
          call cpfile(nlst,'dictin.lst')
          call delfile('dictin.lst')
          call cpfile(nlst,'doace.lst')
          call delfile('doace.lst')
          write(nlst,'(a)')'End of ACEMAKER listing file'
          write(nlst,'(a)')'============================'
          close(nlst)
          title=''
c
c         End temperature cycle
c
        enddo
        call delfile('temp1.lst')
        if (ikeep.ne.1) then
          call delfile('sixlin.pendf')
          call delfile('a.tmp')
        endif
c
c         End material cycle
c
   40 enddo
c
c      Check for ending case cycle
c
      if (ic.eq.0) then
c
c       EOF was found the case before, force to end case cycle
c
        ic=-1
      else
c
c       Get input options for a possible next case
c
        call getinp(ic,title,fendf,nmat,mat,ntemp,temp,tol,ymin,dxmu,
     &              suff,iptab,ipndf,imon,ikeep)
      endif
c
c      End of case cycle
c
      enddo
      call getdtime(cdate,ctime)
      write(nou,'(a)')'==============================================='
      write(nou,'(a,1x,a,1x,a)')'End of ACEMAKER log file',cdate,ctime
      write(nou,'(a)')'==============================================='
      write(*,'(1x,a)')'==============================================='
      write(*,'(1x,a,1x,a,1x,a)')'End of ACEMAKER log file',cdate,ctime
      write(*,'(1x,a)')'==============================================='
      close(nou)
      stop
      end
c======================================================================
      subroutine getinp(ic,title,fendf,nmat,mat,ntemp,temp,tol,ymin,
     &                  dxmu,suff,iptab,ipndf,imon,ikeep)
c
c      Read ACEMAKER input options
c
      implicit real*8 (a-h,o-z)
      dimension mat(*),temp(*)
      character*80 line
      character*72 fendf
      character*66 title
      character*6 keyw
      character*3 suff
      data inp/1/,ndat/2/
      title='ACEMAKER default title (default acemaker.inp)'
      fendf='endfb.in'
      nmat=0
      ntemp=1
      temp(1)=293.6d0
      tol=1.0d-3
      ymin=1.0d-30
      dxmu=1.0d-3
      suff='.00'
      iptab=1
      ipndf=0
      imon=0
      ikeep=0
      if (ic.eq.1) then
        open(inp,file='acemaker.inp', iostat=ista, status='old', err=5)
    5   if (ista.ne.0) then
          write(*,*)' Default values used'
          open(ndat,file=fendf,status='old',err=20)
          read(ndat,*,err=20,end=20)
          read(ndat,'(a66,i4)',err=20,end=20)line(1:66),mat(1)
          close(ndat)
          nmat=1
          ic=0
          return
        endif
      endif
      read(inp,'(a66)',end=15,err=15)title
      keyw=title(1:6)
      istop=0
      if (index(keyw,'end').gt.0.or.index(keyw,'END').gt.0) istop=2
      do while (istop.eq.0)
        read(inp,'(a80)',end=10,err=20)line
        keyw=line(1:6)
        if (index(keyw,'endf').gt.0.or.index(keyw,'ENDF').gt.0) then
          i=7
          do while (line(i:i).eq.' '.and.i.le.80)
            i=i+1
          enddo
          if (i.le.80) then
            fendf=''
            k=1
            do while (line(i:i).ne.' '.and.i.le.80)
              fendf(k:k)=line(i:i)
              k=k+1
              i=i+1
            enddo
          endif
        elseif (index(keyw,'mat').gt.0.or.index(keyw,'MAT').gt.0) then
          read(line(7:80),*,err=20)nmat
          if (nmat.le.0) then
            nmat=0
          else
            do i=1,nmat
              read(inp,*,err=20)mat(i)
            enddo
          endif
        elseif (index(keyw,'temp').gt.0.or.index(keyw,'TEMP').gt.0) then
          read(line(7:80),*,err=20)ntemp
          if (ntemp.le.0) then
            ntemp=1
            temp(1)=293.6d0
          else
            do i=1,ntemp
              read(inp,*,err=20)temp(i)
            enddo
          endif
        elseif (index(keyw,'tol').gt.0.or.index(keyw,'TOL').gt.0) then
          read(line(7:80),*,err=20)tol
          if (tol.le.0.0d0) tol=1.0d-3
        elseif (index(keyw,'ymin').gt.0.or.index(keyw,'YMIN').gt.0) then
          read(line(7:80),*,err=20)ymin
          if (ymin.le.0.0d0) ymin=1.0d-30
        elseif (index(keyw,'dxmu').gt.0.or.index(keyw,'DXMU').gt.0) then
          read(line(7:80),*,err=20)dxmu
          if (dxmu.le.0.0d0) dxmu=1.0d-3
        elseif (index(keyw,'suff').gt.0.or.index(keyw,'SUFF').gt.0) then
          i=index(line,'.')
          if (i.gt.6) then
            suff=line(i:i+2)
            if (suff(2:2).lt.'0'.or.suff(2:2).gt.'9') suff(2:2)='0'
            if (suff(3:3).lt.'0'.or.suff(3:3).gt.'9') suff(3:3)='0'
          else
            suff='.00'
          endif
        elseif (index(keyw,'ptab').gt.0.or.index(keyw,'PTAB').gt.0) then
          read(line(7:80),*,err=20)iptab
          if (iptab.ge.1) then
            iptab=1
          else
            iptab=0
          endif
        elseif (index(keyw,'pndf').gt.0.or.index(keyw,'PNDF').gt.0) then
          read(line(7:80),*,err=20)ipndf
          if (ipndf.ge.1) then
            ipndf=1
          else
            ipndf=0
          endif
        elseif (index(keyw,'mon').gt.0.or.index(keyw,'MON').gt.0) then
          read(line(7:80),*,err=20)imon
          if (imon.ge.1) then
            imon=1
          else
            imon=0
          endif
        elseif (index(keyw,'keep').gt.0.or.index(keyw,'KEEP').gt.0) then
          read(line(7:80),*,err=20)ikeep
          if (ikeep.ge.1) then
            ikeep=1
          else
            ikeep=0
          endif
        elseif (index(keyw,'end').gt.0.or.index(keyw,'END').gt.0) then
          istop=1
        endif
      enddo
   10 if (istop.lt.2) then
c
c       Checking input options
c
        open (ndat,file=fendf,status='old',err=20)
        if (nmat.eq.0) then
          read(ndat,*,err=20,end=20)
          read(ndat,'(a66,i4)',err=20,end=20)line(1:66),mat(1)
          nmat=1
        endif
        close(ndat)
        if (istop.eq.1) then
          ic=2
        else
          ic=0
        endif
      else
        ic=-1
      endif
      if (ikeep.eq.1) ipndf=1
      return
c
c      EOF found
c
   15 ic=-2
      return
   20 write(*,*)
      write(*,*)' *** Fatal error reading acemaker input'
      write(*,*)' *** Check input options and input endf tape'
      close(inp)
      stop
      end
c======================================================================
      subroutine delfile(fname)
c
c      Delete a file
c
      character*(*) fname
      data ndel/91/
      open(ndel,file=fname,err=10)
   10 close(ndel,status='DELETE',err=20)
   20 return
      end
c======================================================================
      subroutine cpfile(ncpy,fname)
c
c      Copy a file into unit ncpy
c
      character*(*) fname
      character*120 line
      data lst/90/
      open(lst,file=fname,status='old',err=10)
      do while (.true.)
        read(lst,'(a)',err=10,end=10)line
        write(ncpy,'(a)')trim(line)
      enddo
   10 write(ncpy,*)
      close(lst,err=20)
   20 return
      end
c======================================================================
      subroutine setcmd(nou,fpath,codenm,cmd)
      character*(*) cmd
      character*(*) fpath
      character*(*) codenm
      character*11  ctime
      character*10  cdate
      call getdtime(cdate,ctime)
      write(cmd,'(a,a)')trim(fpath),trim(codenm)
      if (nou.gt.0) then
        write(nou,'(1x,a,1x,a,1x,a)')ctime,cdate,trim(cmd)
      endif
      write(*,'(1x,1x,a,1x,a,1x,a)')ctime,cdate,trim(cmd)
      write(cmd,'(a,a,a)')trim(fpath),trim(codenm),' 1>a.tmp 2>&1'
      return
      end
c======================================================================
      subroutine getdtime(cdate,ctime)
      character*11  ctime
      character*10  cdate,time
      character*8   date
      call date_and_time(date,time)
      ctime(1:2)=time(1:2)
      ctime(3:3)=':'
      ctime(4:5)=time(3:4)
      ctime(6:6)=':'
      ctime(7:8)=time(5:6)
      ctime(9:9)=':'
      ctime(10:11)=time(8:9)
      cdate(1:4)=date(1:4)
      cdate(5:5)='/'
      cdate(6:7)=date(5:6)
      cdate(8:8)='/'
      cdate(9:10)=date(7:8)
      return
      end
