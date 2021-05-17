      program acemaker
c     version 2.0
c
c     Driver program to produce ACE-formatted files for Monte Carlo
c     calculations.
c
c     The code uses the modules LINEAR, RECENT, SIGMA1, LEGEND, SPECTRA,
c     FIXUP, GROUPIE, MERGER and DICTIN from the PREPRO-2019 code
c     package.
c
c     The package addtionally includes:
c       1. Driver program ACEMAKER
c       2. Module SIXLIN for linearizing and processing MF6 data
c       3. Module DOACE for generating fast ACE-formatted file.
c       4. Module DOTSL for preparing ACE-formatted file describing the
c          thermal scattering law (TSL)
c
c     For preparing fast ACE-formatted files, ACEMAKER invokes the
c     calculation sequence:
c     LINEAR + RECENT + LEGEND + SPECTRA + SIXLIN + + SIGMA1 + FIXUP +
c     GROUPIE + MERGER + DICTIN + DOACE
c
c     For generating a thermal ACE-formatted file from TSL data just
c     the module DOTSL is required.
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
c     Note the use of  '\' (windows) or '/' (linux) at the end.
c
c     It is also possible to redefine fullpaths changing the variables
c     preprod and acemakerd in the corresponding DATA statement on the
c     main program and recompiling ACEMAKER.
c
c     Input file: acemaker.inp
c
c     The input options for running ACEMAKER should be entered on the
c     acemaker.inp file. The name is fixed inside the code and it is a
c     text file containing the input options given by keywords.The first
c     line should be the title of the case. It is followed by a set of
c     keywords, one by line.  The keywords should be entered in the
c     first 6 character positions. The data fields start in column 7,
c     if require. Some keywords could require additional lines. The
c     keywords could be written in lower- or uppercase, but not in a
c     combination. The last keyword should be the card END or end.
c
c     The keyword set of ACEMAKER is described below for each type of
c     calculation:
c
c==============   For all types of calculations   =====================
c
c     Keyword  Data-field(7-72)  Description
c     -----------------------------------------------------------------
c
c     ace       IACE             Type of ACE-formatted file to be
c     ACE                        generated:
c                                  IACE=0 - fast    (Default)
c                                  IACE=1 - thermal
c
c     endf      full-filename    ENDF-6 formatted data filename
c     ENDF
c
c     mat       NMAT             Number of material
c     MAT                        It should be followed by NMAT lines
c                                containing material numbers (MAT),
c                                one by line (NMAT=0 means first MAT
c                                in the input ENDF tape, if IACE=0)
c
c     temp      NTEMP            Number of temperatures
c     TEMP                       It should be followed by NTEMP lines
c                                containing the temperatures in
c                                Kelvin, one by line
c                                (Default: NTEMP=1, 293.6 K)
c
c     tol       TOL              Maximum fraction tolerance for
c     TOL                        linearization/reconstruction of data
c                                (Default: TOL=0.001)
c
c     mcnp      MCNPX            Trigger to use extended ZAID for MCNPX
c     MCNP                       (MCNPX = 0/1 = MCNP/MCNPX, Default = 0)
c
c     suff      .XX              input ID suffix for ZAID on the ACE
c     SUFF                       file (ZAID.XXc). It should start with
c                                '.', followed by two digits (MCNPX=0)
c                                or three digits (MCNPX=1)
c                                (Default: .XX=.00)
c
c     mon       IMON             Monitor printing trigger
c     MON                        IMON = 0/1/2 = min./max./max.+plot
c                                IMON=2 is only allowed for thermal
c                                calculation (IACE=1)
c                                (Default: IMON=0)
c
c     keep      IKEEP            Trigger to keep intermediate files
c     KEEP                       IKEEP = 0/1 = no/yes (Default: IKEEP=0)
c
c     end                        End of case
c     END                        It should be the last keyword
c
c
c========   For fast ACE-formatted file generation (ACE 0)   ==========
c
c     Keyword  Data-field(7-72)  Description
c     -----------------------------------------------------------------
c
c     ymin      FYMIN            Minimum allowable value of
c     YMIN                       cross section or distribution
c                                (Default: FYMIN=1.0d-30)
c
c     dxmu      FMUMAX           Maximum cosine interval(mu) for
c     DXMU                       angular distribution reconstruction
c                                (Default: FMUMAX=0.001)
c
c     ptab      IPTAB            Probability table trigger in the URR
c     PTAB                       IPTAB = 0/1 = no/yes (Default: IPTAB=1)
c
c     pndf      IPNDF            trigger to keep the PENDF tape
c     PNDF                       IPNDF = 0/1 = no/yes (Default: IPNDF=0)
c
c     urr       IURR             trigger to keep the PTAB intermediate
c     URR                        file ZAzzzaaa.URR.ENDF from groupie
c                                IURR = 0/1 = no/yes (Default: IURR=0)
c
c========   For thermal ACE-formatted file generation (ACE 1)   =======
c
c     Keyword  Data-field(7-72)  Description
c     -----------------------------------------------------------------
c
c     nbin      NBIN             Number of equi-probable cosines for
c     NBIN                       incoherent elastic and inelastic
c                                scattering (Default: NBIN=16)
c
c     ethm      ETHMAX           Maximun energy for thermal treatment
c     ETHM                       (Default: ETHMAX=4 eV)
c
c     eps       EPS              Fractional tolerance for preparing
c     EPS                        incident energy grid assuming "1/E"
c                                behaviour (Default: 0.003 = 0.3%)
c
c     thdat     NMATH            Thermal data description by material.
c     THDAT                      NMATH must be equal to NMAT on the MAT
c                                keyword (NMATH=NMAT).
c                                It should be followed by NMATH sets 
c                                of NZAM(i)+1 lines containing the data
c                                presented below:
c
c                                         THDAT NMATH
c     line         1:                       THZAID(1) NMIX(1) NZAM(1)
c     line         2:                         IZAM(1,1)
c     line         3:                         IZAM(2,1)
c     ..............                             ...
c     ..............                          IZAM(j,1)
c     ..............                             ...
c     line NZAM(1)+1:                         IZAM(NZAM(1),1)
c     line NZAM(1)+2:                       THZAID(2) NMIX(2) NZAM(2)
c     line NZAM(1)+3:                         IZAM(1,2)
c     line NZAM(1)+4:                         IZAM(2,2)
c     ..............                             ...
c     ..............                          IZAM(j,2)
c     ..............                             ...
c     line NZAM(1)+NZAM(2)+2:                 IZAM(NZAM(2),2)
c     .............                              ...
c     .............                              ...
c     line NZAM(1)+...+NZAM(i-1)+i:         THZAID(i) NMIX(i) NZAM(i)
c     line NZAM(1)+...+NZAM(i-1)+i+1:         IZAM(1,i)
c     line NZAM(1)+...+NZAM(i-1)+i+2:         IZAM(2,i)
c     ..............                             ...
c     ..............                          IZAM(j,i)
c     ..............                             ...
c     line NZAM(1)+...+NZAM(i-1)+NZAM(i)+i:   IZAM(NZAM(i),i)
c     ..............                             ...
c     ..............                             ...
c
c     where i=1,2,3, ... NMATH  and
c
c     THZAID(i): Thermal ZAID for material MAT(i), up to a maximum of
c                6 characters (Examples: lwtr, hwtr, H_ZrH, grph0)
c     NMIX(i):   Number of atom types mixed in the TSL of material
c                MAT(i) (Default: 1)
c     NZAM(i):   Number of moderator materials for which the TSL of
c                material MAT(i) is applicable (up to a maximum of 16)
c     IZAM(j,i): ZA=1000*Z+A value of moderator material j for which
c                the TSL of material MAT(i) is applicable.
c
c     The first line of each NMATH set contains the values of THZAID(i),
c     NMIX(i) and NZAM(i) for material MAT(i). The following NZAM(i)
c     lines contain the ZA value of each moderator material(one by line)
c
c======================================================================
c
c     The default input options for ACEMAKER are equivalent to the
c     following acemaker.inp file:
c
c     123456789012345678901234567890123456789012345678901234.. (columns)
c
c     ACEMAKER default title (default options)
c     ACE    0
c     ENDF   ENDFB.IN
c     MAT    0
c     TEMP   1
c      293.6
c     TOL    0.001
c     YMIN   1.0e-30
c     DXMU   0.001
c     SUFF   .00
c     PTAB   1
c     PNDF   0
c     MON    0
c     KEEP   0
c     MCNP   0
c     END
c
c     The default is a fast ACE-formatted file generation.The keyword
c     MAT    0 in the 4th input line above, means the first material on
c     the input ENDF tape (ENDFB.IN in this case)
c
c     The keywords can be entered in any order. The END keyword
c     plays the role of case terminator. If a keyword is missing, the
c     default value is applied. So, it is only necessary to include
c     the keywords that change defaults. Several cases can be stacked in
c     one acemaker.inp file by appending after the keyword END a new
c     title line and new input options finishing with a new END keyword.
c     ACEMAKER will stop after the last END.
c
c     For the thermal option ACE 1 (IACE=1) it is compulsory to enter
c     the keywords: ACE, MAT and THDAT. If the rest of input options
c     are omitted, then the following default values are applied:
c
c     ENDF   ENDFB.IN
c     TEMP   1
c      293.6
c     NBIN   16
c     ETHM   4.0
c     TOL    0.001
c     EPS    0.003
c     SUFF   .00
c     MON    0
c     KEEP   0
c     MCNP   0
c
c     Input examples:
c     ===============
c
c     a) Input examples for option ACE 0 (fast)
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
c     In the case 1 above, evaluated data will be retrieved from
c     FENDLEN.DAT for 3 materials with MAT numbers 5728(La-139),
c     725(N-14)and 825(O-16). The ZAID on the ACE formatted files will
c     have the suffix .32c for each material. Default values apply to
c     the rest of input options.
c     In the case 2, evaluated nuclear data are read for H-1 from
c     \ENDF-B-VIII.0\n-001_H_001.endf (NMAT=0, wich means first material
c     on the ENDF tape). Two temperatures are requested, then two ace-
c     formatted files will be generated: one at 293.6 K and other at
c     600.0 K. The ZAID number will be 1001.80c and 1001.81c. Default
c     values are used for the rest of input options.
c     In the case 3 evaluated nuclear data are retrieved from tape
c     \ENDF\ENDFBVIII.endf for U-235 (9228) and U-238 (9237). Two
c     temperatures are requested for each material. Four ACE-formatted
c     files are generated with ZAID numbers:
c      92235.80c at 293.6K and 92235.81c at 600.0K for U-235 (9228)
c      92238.80c at 293.6K and 92238.81c at 600.0K for U-238 (9237)
c     Default values are taken for the rest of input options.
c
c     For ACE 0 (fast option), each produced ACE-formatted file  is
c     named according to the pattern ZAzzzaaa.XXc.acef where zzzaaa is
c     the ZA number printed in a field of six positions and filled with
c     zeroes in place of spaces. It follows the GROUPIE-2019 conventions
c     for 2-band data filenames. The symbol .XX is the ZAID suffix
c     entered as input for the first temperature and it is eventually
c     increased in .01 or .001 for each additional temperature(NTEMP>1).
c     For the cases previously presented the filenames are:
c
c     Case 1:
c     ZA057139.32c.acef  for  La-139
c     ZA007014.32c.acef  for  N-14
c     ZA008016.32c.acef  for  O-16
c     Case 2:
c     ZA001001.80c.acef  and  ZA001001.81c.acef  for H-1
c     Case 3:
c     ZA092235.80c.acef  and  ZA092235.81c.acef  for U-235
c     ZA092238.80c.acef  and  ZA092238.81c.acef  for U-238
c
c     The output files ZAzzzaaa.XXc.xsd and ZAzzzaaa.XXc.lst are also
c     prepared. The first file contains the directory information for
c     the XSDIR file of MCNP and the second one is a compilation of the
c     listing files produced by all the modules called by ACEMAKER. The
c     files acemaker.log summaries ACEMAKER processing.
c
c     b) Input examples for the option ACE 1 (thermal):
c
c     123456789012345678901234567890123456789012345678901234.. (columns)
c
c     H bound in H2O
c     ACE    1
c     ENDF   \ENDF-B-VIII.0\tsl\tsl-HinH2O.endf
c     MAT    1
c      1
c     TEMP   1
c      293.6
c     MON    2
c     NBIN   64
c     ETHM   4.0
c     TOL    0.001
c     EPS    0.003
c     THDAT  1
c      H_H2O 1 1
c      1001
c     SUFF   .00
c     MCNP   0
c     END
c
c     Multiple cases are also allowed:
c
c     H bound in H2O
c     ACE    1
c     ENDF   \ENDF-B-VIII.0\tsl\tsl-HinH2O.endf
c     MAT    1
c      1
c     TEMP   1
c      293.6
c     MON    2
c     NBIN   64
c     ETHM   4.0
c     TOL    0.001
c     EPS    0.003
c     THDAT  1
c      H_H2O 1 1
c      1001
c     SUFF   .00
c     MCNP   0
c     END
c     D bound in D2O
c     ACE    1
c     ENDF   \ENDF-B-VIII.0\tsl\tsl-DinD2O.endf
c     MAT    1
c      11
c     TEMP   1
c      293.6
c     MON    2
c     NBIN   64
c     ETHM   4.0
c     TOL    0.001
c     EPS    0.003
c     THDAT   1
c      D_D2O  1 1
c      1002
c     SUFF   .00
c     MCNP   0
c     END
c     Graphite 0% porosity
c     ACE    1
c     ENDF   \ENDF-B-VIII.0\tsl\tsl-crystalline-graphite.endf
c     MAT    1
c      30
c     TEMP   2
c      296.0
c      400.0
c     MON    2
c     NBIN   64
c     ETHM   4.0
c     TOL    0.001
c     EPS    0.003
c     THDAT  1
c      grph0 1 3
c      6000
c      6012
c      6013
c     SUFF   .00
c     MCNP   0
c     END
c     ZrH(H in ZrH and Zr in ZrH are contained on tsl-ZrH.endf)
c     ACE    1
c     ENDF   \ENDF-B-VIII.0\tsl\tsl-ZrH.endf
c     MAT    2
c      7
c      58
c     TEMP   1
c      296.0
c     MON    2
c     NBIN   64
c     ETHM   4.0
c     TOL    0.001
c     EPS    0.003
c     THDAT  2
c      H_ZrH  1 1
c      1001
c      Zr_ZrH 1 3
c      40000
c      40090
c      40091
c     SUFF   .00
c     MCNP   0
c     END
c
c     For ACE 1 (thermal option), each ACE-formatted file generated is
c     named according to the pattern .acef where  THZAID is
c     the six character ID entered as input on the THDAT keyword. The
c     symbol .XX is the ZAID suffix entered as input for the first
c     temperature and eventually increased in .01 or .001 for each
c     additional temperature, if NTEMP>1.
c
c     In the previous examples the following files are generated:
c      H_H2O.00t.acef  at 293.6 K
c      D_D2O.00t.acef  at 293.6 K
c      grph0.00t.acef  at 296.0 K
c      grph0.01t.acef  at 400.0 K
c      H_ZrH.00t.acef  at 296.0 K
c     Zr_ZrH.00t.acef  at 296.0 K
c
c     with the corresponding THZAID.XXt.xsd files for XSDIR.
c     Additionaly, if IMON=2, the files THZAID.XXt.plt and
c     THZAID.XXt.cur are saved for PLOTTAB plottings.
c
      implicit real*8 (a-h,o-z)
      parameter (nmatx=300, ntempx=25)
      dimension temp(ntempx)
      dimension mat(nmatx),nmix(nmatx),nzam(nmatx),izam(16,nmatx)
      character*80 cmd
      character*72 fendf
      character*66 title
      character*62 preprod,acemakerd
      character*11 ctime
      character*10 cdate
      character*8  fza
      character*6  thzaid(nmatx)
      character*4  suff
      data nin/20/,nou/21/,nlst/22/,nplt/30/,ncur/31/
      data preprod/'\ACEMAKER\exe\'/,acemakerd/'\ACEMAKER\exe\'/
      data emax/1.0d+09/,fangmin/1.0d-10/,flegmin/1.0d-02/
c
c      Open acemaker log file
c
      call getdtime(cdate,ctime)
      call delfile('acemaker.log')
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
      call getinp(ic,title,iace,fendf,nmat,mat,ntemp,temp,tol,mcnpx,
     &  ymin,dxmu,xsuff,iptab,ipndf,iurr,imon,ikeep,nbin,ethmax,eps,
     &  thzaid,nmix,nzam,izam)
c
c      case cycle
c
      do while (ic.ge.0)
c
c      Print input options for case jc
c
       jc=jc+1
       write(*,*)
       write(*,'(1x,a6,i4)')'Case: ',jc
       write(*,'(1x,a10)')'=========='
       write(*,'(1x,a7,a)')'Title: ',trim(title)
       if (iace.eq.0) then
         write(*,'(1x,a)')'Fast ACE-formatted file'
       elseif (iace.eq.1) then
         write(*,'(1x,a)')'Thermal ACE-formatted file'
       endif
       write(nou,*)
       write(nou,'(a6,i4)')'Case: ',jc
       write(nou,'(a10)')'=========='
       write(nou,'(a7,a)')'Title: ',trim(title)
       write(nou,*)
       write(nou,'(a17,a)')'ENDF input tape: ',trim(fendf)
       write(nou,*)
       if (iace.eq.0) then
         write(nou,'(a)')'Fast ACE-formatted file'
       elseif (iace.eq.1) then
         write(nou,'(a)')'Thermal ACE-formatted file'
       endif
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
       write(nou,*)
       write(nou,'(a33,i2)')'MCNP trigger: 0/1 = MCNP/MCNPX : ',mcnpx
       write(nou,'(a,i2)')'Monitor option 0/1/2=min/max/plt:',imon
       write(nou,'(a25,8x,1pe12.5)')'Linearization tolerance: ',tol
       nsuff=nint(1000.0d0*xsuff+1.0d-5)
       if (nsuff.ge.1000) nsuff=nsuff-(nsuff/1000)*1000
       write(suff,'(a1,i3)')'.',nsuff
       if (suff(2:2).eq.' ') suff(2:2)='0'
       if (suff(3:3).eq.' ') suff(3:3)='0'
       if (mcnpx.eq.1) then
         write(nou,'(a,7x,a4)')'Entered ID-suffix for ZAID:',suff(1:4)
       else
         write(nou,'(a,7x,a3)')'Entered ID-suffix for ZAID:',suff(1:3)
       endif
       suff='    '
       xsuff0=xsuff
       if (iace.eq.0) then
         write(nou,'(a,4x,1pe12.5)')'Min. allowable cross section:',ymin
         write(nou,'(a,12x,1pe12.5)')'Max. cosine interval:',dxmu
         write(nou,'(a,1x,i2)')'Probability tables 0/1 = no/yes:',iptab
         write(nou,'(a,4x,i2)')'Keep PENDF tape 0/1 = no/yes:',ipndf
         write(nou,'(a,5x,i2)')'Keep all tapes 0/1 = no/yes:',ikeep
       elseif (iace.eq.1) then
         write(nou,'(a,i6)')'Number of equiprobable cosines:',nbin
         write(nou,'(a,1pe12.5)')'Maximum thermal energy:',ethmax
         write(nou,'(a,a,1pe12.5)')'Fractional tolerance for',
     &     ' incident energy grid:',eps
         write(nou,*)
         write(nou,'(a)')'Thermal data by material:'
         write(nou,'(a)')'========================='
         write(nou,'(a5,a10,a8,a6,a5,a10)')'  No.',' MATERIAL ',
     &     ' THZAID ',' NMIX ',' NZA ',' ZA-values'
         do i=1,nmat
           write(nou,'(i4,i8,4x,a6,i5,i5,i10)')i,mat(i),thzaid(i),
     &       nmix(i),nzam(i),izam(1,i)
           if (nzam(i).gt.1) then
             do j=2,nzam(i)
               write(nou,'(32x,i10)')izam(j,i)
             enddo
           endif
         enddo
       endif
c
c      Deleting internal files
c
       call delfile('MERGER.INP')
       call delfile('MERGER.LST')
       call delfile('LINEAR.INP')
       call delfile('LINEAR.LST')
       call delfile('RECENT.INP')
       call delfile('RECENT.LST')
       call delfile('LEGEND.INP')
       call delfile('LEGEND.LST')
       call delfile('LEGEND.TMP')
       call delfile('SPECTRA.INP')
       call delfile('SPECTRA.LST')
       call delfile('SIXLIN.INP')
       call delfile('SIXLIN.LST')
       call delfile('SIXLIN.TMP')
       call delfile('SIGMA1.INP')
       call delfile('SIGMA1.LST')
       call delfile('FIXUP.INP')
       call delfile('FIXUP.LST')
       call delfile('GROUPIE.INP')
       call delfile('GROUPIE.LST')
       call delfile('DICTIN.INP')
       call delfile('DICTIN.LST')
       call delfile('DOACE.INP')
       call delfile('DOACE.LST')
       call delfile('DOTSL.INP')
       call delfile('DOTSL.LST')
       call delfile('DOTSL.PLT')
       call delfile('DOTSL.CUR')
       call delfile('ENDF6.ENDF')
       call delfile('LINEAR.PENDF')
       call delfile('RECENT.PENDF')
       call delfile('LEGEND.PENDF')
       call delfile('SPECTRA.PENDF')
       call delfile('SIXLIN.PENDF')
       call delfile('SIGMA1.PENDF')
       call delfile('FIXUP.PENDF')
       call delfile('GROUPIE.ENDF')
       call delfile('GROUPIE.PENDF')
       call delfile('ENDF6.PENDF')
       call delfile('TEMP1.LST')
       call delfile('a.tmp')
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
c        MERGER: get endf-6 formatted data from endf input tape
c
         open(nin,file='MERGER.INP')
         write(nin,'(a)')'ENDF6.ENDF'
         write(nin,'(a19,i4)')'ENDF tape for MAT= ',mati
         write(nin,'(a)')trim(fendf)
         write(nin,'(a3)')'END'
         write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
         write(nin,*)
         close(nin)
         call setcmd(nou,preprod,'merger',cmd)
         call system(cmd)
         call delfile('MERGER.INP')
         call delfile('MERGER.LST')
c
c        Get ZA
c
         matj=0
         za=0.0d0
         open(nin,file='ENDF6.ENDF',err=10)
         read(nin,*,err=10,end=10)
         read(nin,'(d11.0,55x,i4)',err=10,end=10)za,matj
   10    if (matj.ne.mati) then
           write(nou,'(a,i4,a)')'*** ERROR: MAT= ',mati,' not found'
           write(*,'(1x,a,i4,a)')'*** ERROR: MAT= ',mati,' not found'
           close (nin,err=15)
   15      call delfile('ENDF6.ENDF')
           read(nin,'(d11.0,55x,i4)',err=40,end=40)za,matj
         endif
         close(nin)
         nza=za
         fza=' '
         write(fza,'(a2,i6)')'ZA',nza
         do ii=3,8
           if (fza(ii:ii).eq.' ') fza(ii:ii)='0'
         enddo
         xsuff=xsuff0
c
c        Fast ACE-file
c
         if (iace.eq.0) then
c
c          LINEAR: linearize cross sections at 0.0 K
c
           open(nin,file='LINEAR.INP')
           write(nin,'(2i11,1pe11.4,i11)')0,imon,ymin,1
           write(nin,'(a)')'ENDF6.ENDF'
           write(nin,'(a)')'LINEAR.PENDF'
           write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
           write(nin,*)
           write(nin,'(1p2e11.4)')0.0,0.5d0*tol
           write(nin,'(1p2e11.4)')1.0,0.5d0*tol
           write(nin,'(1p2e11.4)')2.0,tol
           write(nin,'(1p2e11.4)')emax,tol
           write(nin,*)
           close(nin)
           call setcmd(nou,preprod,'linear',cmd)
           call system(cmd)
           if (ikeep.ne.1) then
             call delfile('LINEAR.INP')
             call delfile('ENDF6.ENDF')
           endif
c
c          RECENT: reconstruct resonance cross sections at 0.0 K
c
           open(nin,file='RECENT.INP')
           write(nin,'(i11,1pe11.4,4i11)')0,ymin,1,1,1,imon
           write(nin,'(a)')'LINEAR.PENDF'
           write(nin,'(a)')'RECENT.PENDF'
           write(nin,'(2i11)')mati,mati
           write(nin,*)
           write(nin,'(1p2e11.4)')0.0,0.5d0*tol
           write(nin,'(1p2e11.4)')1.0,0.5d0*tol
           write(nin,'(1p2e11.4)')2.0,tol
           write(nin,'(1p2e11.4)')emax,tol
           write(nin,*)
           close(nin)
           call setcmd(nou,preprod,'recent',cmd)
           call system(cmd)
           if (ikeep.ne.1) then
             call delfile('RECENT.INP')
             call delfile('LINEAR.PENDF')
           endif
c
c          Merge temperature independent listing files
c
           open(nlst,file='TEMP1.LST')
           write(nlst,'(a)')'ACEMAKER listing file'
           write(nlst,'(a)')'====================='
           call cpfile(nlst,'LINEAR.LST')
           call delfile('LINEAR.LST')
           call cpfile(nlst,'RECENT.LST')
           call delfile('RECENT.LST')
         endif
c
c        Fast ACE-file
c
         if (iace.eq.0) then
c
c          LEGEND: Process/linearize angular distributions on MF4
c
           open(nin,file='LEGEND.INP')
           write(nin,'(1pe11.4,5i11)')tol,0,2,1,2,0
           write(nin,'(a)')'RECENT.PENDF'
           write(nin,'(a)')'LEGEND.PENDF'
           write(nin,'(i6,i2,i3,i6,i2,i3,1p4e11.4)')mati,1,1,
     &                 mati,99,999,0.0d0,emax,fangmin,flegmin
           write(nin,*)
           close(nin)
           call setcmd(nou,preprod,'legend',cmd)
           call system(cmd)
           call delfile('LEGEND.TMP')
           if (ikeep.ne.1) then
             call delfile('LEGEND.INP')
             call delfile('RECENT.PENDF')
           endif
c
c          SPECTRA: Process/linearize energy distributions on MF5
c
           open(nin,file='SPECTRA.INP')
           write(nin,'(2i11,1pe11.4,i11)')0,imon,ymin,1
           write(nin,'(a)')'LEGEND.PENDF'
           write(nin,'(a)')'SPECTRA.PENDF'
           write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
           write(nin,*)
           write(nin,'(1p2e11.4)')0.0,tol
           write(nin,*)
           close(nin)
           call setcmd(nou,preprod,'spectra',cmd)
           call system(cmd)
           if (ikeep.ne.1) then
             call delfile('SPECTRA.INP')
             call delfile('LEGEND.PENDF')
           endif
c
c          SIXLIN: Process/linearize angular-energy distributions on MF6
c
           open(nin,file='SIXLIN.INP')
           write(nin,'(2i11,1pe11.4)')0,imon,ymin
           write(nin,'(a)')'SPECTRA.PENDF'
           write(nin,'(a)')'SIXLIN.PENDF'
           write(nin,'(2i11)')mati,mati
           write(nin,'(1p2e11.4)')tol,dxmu
           write(nin,*)
           close(nin)
           call setcmd(nou,acemakerd,'sixlin',cmd)
           call system(cmd)
           call delfile('SIXLIN.TMP')
           if (ikeep.ne.1) then
             call delfile('SIXLIN.INP')
             call delfile('SPECTRA.PENDF')
           endif
c
c          Merge temperature independent listing files
c
           call cpfile(nlst,'LEGEND.LST')
           call delfile('LEGEND.LST')
           call cpfile(nlst,'SPECTRA.LST')
           call delfile('SPECTRA.LST')
           call cpfile(nlst,'SIXLIN.LST')
           call delfile('SIXLIN.LST')
         endif
c
c         Temperature cycle
c
         do it=1,ntemp
           tempi=temp(it)
           write(nou,'(a5,i4,a7,1pe11.4)')'MAT= ',mati,' TEMP= ',tempi
           write(*,'(1x,a5,i4,a7,1pe11.4)')'MAT= ',mati,' TEMP= ',tempi
c
c          Fast ACE-file
c
           if (iace.eq.0) then
c
c            SIGMA1: Doppler broadening cross sections
c
             open(nin,file='SIGMA1.INP')
             write(nin,'(2i11,1p2e11.4,2i11)')0,imon,tempi,ymin,1,0
             if (iace.eq.0) then
               write(nin,'(a)')'SIXLIN.PENDF'
             else
               write(nin,'(a)')'RECENT.PENDF'
             endif
             write(nin,'(a)')'SIGMA1.PENDF'
             write(nin,'(2i11)')mati,mati
             write(nin,*)
             write(nin,'(1p2e11.4)')0.0,0.5d0*tol
             write(nin,'(1p2e11.4)')1.0,0.5d0*tol
             write(nin,'(1p2e11.4)')2.0,tol
             write(nin,'(1p2e11.4)')emax,tol
             write(nin,*)
             close(nin)
             call setcmd(nou,preprod,'sigma1',cmd)
             call system(cmd)
             if (ikeep.ne.1) call delfile('SIGMA1.INP')
           endif
c
c          Fast ACE-file
c
           if (iace.eq.0) then
c
c            FIXUP: fix formats and XS. Prepare unified energy grid
c
             open(nin,file='FIXUP.INP')
             write(nin,'(a14)')'10002111100011'
             write(nin,'(a)')'SIGMA1.PENDF'
             write(nin,'(a)')'FIXUP.PENDF'
             close(nin)
             call setcmd(nou,preprod,'fixup',cmd)
             call system(cmd)
             if (ikeep.ne.1) then
               call delfile('FIXUP.INP')
               call delfile('SIGMA1.PENDF')
             endif
c
c            Get probability tables
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
c              GROUPIE: Generate 2-bands(bins) probability tables in URR
c
               open(nin,file='GROUPIE.INP')
               write(nin,'(4i11,1pe11.4,i11)')0,-11,2,0,1.0d-3,0
               write(nin,'(a)')'FIXUP.PENDF'
               write(nin,'(a)')'GROUPIE.ENDF'
               write(nin,'(5i11)')1,1,1,1,1
               write(nin,'(a20,i5,a6,1pe11.4)')'GROUPIE 2-BANDS MAT=',
     &           mati,' TEMP=',tempi
               write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
               write(nin,*)
               close(nin)
               call setcmd(nou,preprod,'groupie',cmd)
               call system(cmd)
               if (ikeep.ne.1) then
                 call delfile('GROUPIE.INP')
                 call delfile('GROUPIE.ENDF')
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
c              Checking fza=ZAzzzaaa.URR.ENDF
c
               write(cmd,'(a8,a)')fza,'.URR.ENDF'
               title=trim(cmd)
               ichk=1
               open(nin,file=title,status='old',err=20)
               read(nin,'(a)',err=20,end=20)cmd
               read(nin,'(a)',err=20,end=20)cmd
               ichk=0
   20          close(nin)
               if (ichk.eq.0) then
c
c                MERGER: Merge sections 2/152 & 2/153
c
                 open(nin,file='MERGER.INP')
                 write(nin,'(a)')'GROUPIE.PENDF'
                 write(nin,'(a,i4,a)')'PENDF tape MAT= ',mati,' PTAB-2'
                 write(nin,'(a)')title
                 write(nin,'(a)')'FIXUP.PENDF'
                 write(nin,'(a3)')'END'
                 write(nin,'(i6,i2,i3,i6,i2,i3)')mati,1,1,mati,99,999
                 write(nin,*)
                 close(nin)
                 call setcmd(nou,preprod,'merger',cmd)
                 call system(cmd)
                 if (ikeep.ne.1) then
                   call delfile('MERGER.INP')
                   if (iurr.ne.1) call delfile(title)
                 endif
               else
                 write(nou,'(a,a)')' ** Warning:',
     &             ' No unresolved probability tables'
                 write(*,'(1x,a,a)')' ** Warning:',
     &             ' No unresolved probability tables'
               endif
             endif
c
c            DICTIN: Update dictionary section of pendf tape
c
             call delfile('DICTIN.INP')
             call delfile('DICTIN.LST')
             open(nin,file='DICTIN.INP')
             if (iptab.eq.1.and.ichk.eq.0) then
               write(nin,'(a)')'GROUPIE.PENDF'
             else
               write(nin,'(a)')'FIXUP.PENDF'
             endif
             write(nin,'(a)')'ENDF6.PENDF'
             close(nin)
             call setcmd(nou,preprod,'dictin',cmd)
             call system(cmd)
             if (ikeep.ne.1) then
               call delfile('DICTIN.INP')
               call delfile('FIXUP.PENDF')
               if (ichk.eq.0) call delfile('GROUPIE.PENDF')
             endif
c
c            DOACE: Prepare fast ACE-formatted file for MC
c                   (Filename: ZAaaazzz.xxc.acef)
c
             if (mcnpx.eq.1) then
               xsuff=xsuff+0.001d0*dble(it-1)
             else
               xsuff=xsuff+0.01d0*dble(it-1)
             endif
             nsuff=nint(1000.0d0*xsuff+1.0d-5)
             if (nsuff.ge.1000) nsuff=nsuff-(nsuff/1000)*1000
             write(suff,'(a1,i3)')'.',nsuff
             if (suff(2:2).eq.' ') suff(2:2)='0'
             if (suff(3:3).eq.' ') suff(3:3)='0'
             if (mcnpx.eq.1) then
               write(cmd,'(a8,a4,a7)')fza,suff(1:4),'nc.acef'
             else
               write(cmd,'(a8,a3,a6)')fza,suff(1:3),'c.acef'
             endif
             title=trim(cmd)
             call delfile(title)
             open(nin,file='DOACE.INP')
             write(nin,'(3i11)')0,imon,mcnpx
             write(nin,'(a)')'ENDF6.PENDF'
             write(nin,'(a)')trim(title)
             write(nin,'(i11)')mati
             if (mcnpx.eq.1) then
               write(nin,'(a4)')suff(1:4)
             else
               write(nin,'(a3,a1)')suff(1:3),' '
             endif
             close(nin)
             call setcmd(nou,acemakerd,'doace',cmd)
             call system(cmd)
             if (ikeep.ne.1) call delfile('DOACE.INP')
c
c            Simple file access checking
c
             ichk=1
             open(nin,file=title,status='old',err=30)
             read(nin,'(a)',err=30,end=30)cmd
             read(nin,'(a)',err=30,end=30)cmd
             ichk=0
   30        close(nin)
             if(ichk.eq.0) then
               write(nou,'(a11,a,a16,i4,a7,1pe11.4)')' ACE-file: ',
     &           trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
               write(*,'(1x,a11,a,a16,i4,a7,1pe11.4)')' ACE-file: ',
     &           trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
             else
               write(nou,'(a24,a,a16,i4,a7,1pe11.4)')
     &           ' Access error ACE-file: ',
     &           trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
               write(*,'(1x,a24,a,a16,i4,a7,1pe11.4)')
     &           ' Access error ACE-file: ',
     &           trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
             endif
c
c            Save pendf tape as ZAaaazzz.xxc.pendf, if ipndf=1
c
             if (ipndf.gt.0.or.ikeep.eq.1) then
               if (mcnpx.eq.1) then
                 write(cmd,'(a8,a4,a8)')fza,suff(1:4),'nc.pendf'
               else
                 write(cmd,'(a8,a3,a7)')fza,suff(1:3),'c.pendf'
               endif
               title=trim(cmd)
               call delfile(title)
               open(nlst,file=title)
               call cpfile(nlst,'ENDF6.PENDF')
               close(nlst)
             endif
             call delfile('ENDF6.PENDF')
c
c            Save *.lst files as ZAaaazzz.xxc.lst
c
             close(nlst)
             if (mcnpx.eq.1) then
               write(cmd,'(a8,a4,a6)')fza,suff(1:4),'nc.lst'
             else
               write(cmd,'(a8,a3,a5)')fza,suff(1:3),'c.lst'
             endif
             title=trim(cmd)
             call delfile(title)
             open(nlst,file=title)
             call cpfile(nlst,'TEMP1.LST')
             call cpfile(nlst,'SIGMA1.LST')
             call delfile('SIGMA1.LST')
             call cpfile(nlst,'FIXUP.LST')
             call delfile('FIXUP.LST')
             if (iptab.eq.1) then
               call cpfile(nlst,'GROUPIE.LST')
               call delfile('GROUPIE.LST')
               call cpfile(nlst,'MERGER.LST')
               call delfile('MERGER.LST')
             endif
             call cpfile(nlst,'DICTIN.LST')
             call delfile('DICTIN.LST')
             call cpfile(nlst,'DOACE.LST')
             call delfile('DOACE.LST')
           elseif (iace.eq.1) then
c
c            DOTSL: Prepare thermal ACE-formatted file for MC
c                   (Filename: THZAID.xxt.acef)
c
             if (mcnpx.eq.1) then
               xsuff=xsuff+0.001d0*dble(it-1)
             else
               xsuff=xsuff+0.01d0*dble(it-1)
             endif
             nsuff=nint(1000.0d0*xsuff+1.0d-5)
             if (nsuff.ge.1000) nsuff=nsuff-(nsuff/1000)*1000
             write(suff,'(a1,i3)')'.',nsuff
             if (suff(2:2).eq.' ') suff(2:2)='0'
             if (suff(3:3).eq.' ') suff(3:3)='0'
             i0=6
             do i=1,6
               if (thzaid(im)(i:i).ne.' ') then
                 i0=i
                 exit
               endif
             enddo
             fza='        '
             fza=thzaid(im)(i0:6)
             if (mcnpx.eq.1) then
               write(cmd,'(a,a4,a7)')trim(fza),suff(1:4),
     &           'nt.acef'
             else
               write(cmd,'(a,a3,a6)')trim(fza),suff(1:3),
     &           't.acef'
             endif
             title=trim(cmd)
             call delfile(title)
             open(nin,file='DOTSL.INP')
             write(nin,'(a)')'ENDF6.ENDF'
             write(nin,'(a)')trim(title)
             write(nin,'(i11,1pe11.4,2i11)')mati,tempi,nmix(im),imon
             write(nin,'(i11,1p3e11.4)')nbin,ethmax,tol,eps
             if (mcnpx.eq.1) then
               write(nin,'(5x,a6,7x,a4,2i11)')trim(fza),suff(1:4),
     &           mcnpx,nzam(im)
             else
               write(nin,'(5x,a6,8x,a3,2i11)')trim(fza),suff(1:3),
     &           mcnpx,nzam(im)
             endif
             write(nin,'(2(8i7))')(izam(i,im),i=1,nzam(im))
             close(nin)
             call setcmd(nou,acemakerd,'dotsl',cmd)
             call system(cmd)
             write(nou,'(a,a,a,i4,a,1pe11.4)')' ACE-file: ',
     &         trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
             write(*,'(1x,a,a,a,i4,a,1pe11.4)')' ACE-file: ',
     &         trim(title),' saved for MAT= ',mati,' TEMP= ',tempi
             if (ikeep.ne.1) call delfile('DOTSL.INP')
c
c            Save *.lst file as THZAID.xxt.lst
c
             if (mcnpx.eq.1) then
               write(cmd,'(a,a4,a6)')trim(fza),suff(1:4),'nt.lst'
             else
               write(cmd,'(a,a3,a5)')trim(fza),suff(1:3),'t.lst'
             endif
             title=trim(cmd)
             call delfile(title)
             open(nlst,file=title)
             write(nlst,'(a)')'ACEMAKER listing file'
             write(nlst,'(a)')'====================='
             call cpfile(nlst,'DOTSL.LST')
             call delfile('DOTSL.LST')
c
c            Save *.log file as THZAID.xxt.log
c             
             if (mcnpx.eq.1) then
               write(cmd,'(a,a4,a6)')trim(fza),suff(1:4),'nt.log'
             else
               write(cmd,'(a,a3,a5)')trim(fza),suff(1:3),'t.log'
             endif
             title=trim(cmd)
             call delfile(title)
             open(nin,file=title)
             write(nin,'(a)')'DOTSL log file'
             write(nin,'(a)')'=============='
             call cpfile(nin,'a.tmp')
             close(nin)
c
c            Save *.plt and *.cur files
c                          
             if (imon.eq.2) then
               if (mcnpx.eq.1) then
                 write(cmd,'(a,a4,a6)')trim(fza),suff(1:4),'nt.plt'
               else
                 write(cmd,'(a,a3,a5)')trim(fza),suff(1:3),'t.plt'
               endif
               title=trim(cmd)
               call delfile(title)
               open (nplt,file=title)
               call cpfile(nplt,'DOTSL.PLT')
               close(nplt)
               call delfile('DOTSL.PLT')
               write(nou,'(a,a,a,i4,a,1pe11.4)')' PLT-file: ',
     &           trim(title),'  saved for MAT= ',mati,' TEMP= ',tempi
               write(*,'(1x,a,a,a,i4,a,1pe11.4)')' PLT-file: ',
     &           trim(title),'  saved for MAT= ',mati,' TEMP= ',tempi
               if (mcnpx.eq.1) then
                 write(cmd,'(a,a4,a6)')trim(fza),suff(1:4),'nt.cur'
               else
                 write(cmd,'(a,a3,a5)')trim(fza),suff(1:3),'t.cur'
               endif
               title=trim(cmd)
               call delfile(title)
               open (ncur,file=title)
               call cpfile(ncur,'DOTSL.CUR')
               close(ncur)
               call delfile('DOTSL.CUR')
               write(nou,'(a,a,a,i4,a,1pe11.4)')' CUR-file: ',
     &           trim(title),'  saved for MAT= ',mati,' TEMP= ',tempi
               write(*,'(1x,a,a,a,i4,a,1pe11.4)')' CUR-file: ',
     &           trim(title),'  saved for MAT= ',mati,' TEMP= ',tempi
             endif
           endif
           call delfile('a.tmp')
           write(nlst,'(a)')'End of ACEMAKER listing file'
           write(nlst,'(a)')'============================'
           close(nlst)
           title=' '
c
c          End temperature cycle
c
         enddo
         if (ikeep.ne.1) then
           if (iace.eq.0) then
             call delfile('TEMP1.LST')
             call delfile('SIXLIN.PENDF')
           elseif (iace.eq.1) then
             call delfile('ENDF6.ENDF')
           else
             call delfile('TEMP1.LST')
             call delfile('RECENT.PENDF')
           endif
         endif
c
c        End material cycle
c
   40    continue
       enddo
c
c      Check for ending case cycle
c
       if (ic.eq.0) then
c
c       EOF was found before, then case cycle trigger set to -1
c
        ic=-1
       else
c
c       Get input options for next case
c
        call getinp(ic,title,iace,fendf,nmat,mat,ntemp,temp,tol,mcnpx,
     &    ymin,dxmu,xsuff,iptab,ipndf,iurr,imon,ikeep,nbin,ethmax,eps,
     &    thzaid,nmix,nzam,izam)
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
      subroutine getinp(ic,title,iace,fendf,nmat,mat,ntemp,temp,tol,
     &  mcnpx,ymin,dxmu,xsuff,iptab,ipndf,iurr,imon,ikeep,nbin,ethmax,
     &  eps,thzaid,nmix,nzam,izam)
c
c      Read ACEMAKER input options
c
      implicit real*8 (a-h,o-z)
      character*80 line
      character*72 fendf
      character*66 title
      character*6 keyw,thzaid(*)
      dimension temp(*),mat(*),nmix(*),nzam(*),izam(16,*)
      data inp/1/,ndat/2/
      title='ACEMAKER (default input options)'
      iace=0
      fendf='endfb.in'
      nmat=0
      ntemp=1
      temp(1)=293.6d0
      tol=1.0d-3
      ymin=1.0d-30
      dxmu=1.0d-3
      xsuff=0.0d0
      iptab=1
      ipndf=0
      iurr=0
      imon=0
      ikeep=0
      nbin=16
      ethmax=4.0d0
      eps=3.0d-3
      mcnpx=0
      nmath=-1
      keyw='      '
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
        if (index(keyw,'ace').gt.0.or.index(keyw,'ACE').gt.0) then
          read(line(7:80),*,err=20)iace
          if (iace.lt.0.or.iace.gt.2) iace=0
        elseif (index(keyw,'endf').gt.0.or.index(keyw,'ENDF').gt.0) then
          i=7
          do while (line(i:i).eq.' '.and.i.le.80)
            i=i+1
          enddo
          if (i.le.80) then
            fendf=' '
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
          if (tol.le.0.0d0) then
            tol=1.0d-3
          elseif (tol.lt.1.0d-9) then
            tol=1.0d-9
          endif
        elseif (index(keyw,'mcnp').gt.0.or.index(keyw,'MCNP').gt.0) then
          read(line(7:80),*,err=20)mcnpx
          if (mcnpx.le.0.or.mcnpx.gt.1) mcnpx=0
        elseif (index(keyw,'suff').gt.0.or.index(keyw,'SUFF').gt.0) then
          read(line(7:80),*,err=20)xsuff
          if (xsuff.lt.0.0d0.or.xsuff.ge.1.0d0) xsuff=0.0d0
        elseif (index(keyw,'ymin').gt.0.or.index(keyw,'YMIN').gt.0) then
          read(line(7:80),*,err=20)ymin
          if (ymin.lt.1.0d-30) ymin=1.0d-30
        elseif (index(keyw,'dxmu').gt.0.or.index(keyw,'DXMU').gt.0) then
          read(line(7:80),*,err=20)dxmu
          if (dxmu.le.0.0d0) then
            dxmu=1.0d-3
          elseif (dxmu.lt.1.0d-9) then
            dxmu=1.0d-9
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
        elseif (index(keyw,'urr').gt.0.or.index(keyw,'URR').gt.0) then
          read(line(7:80),*,err=20)iurr
          if (iurr.ge.1) then
            iurr=1
          else
            iurr=0
          endif
        elseif (index(keyw,'mon').gt.0.or.index(keyw,'MON').gt.0) then
          read(line(7:80),*,err=20)imon
          if (imon.lt.0.or.imon.gt.2) imon=0
        elseif (index(keyw,'keep').gt.0.or.index(keyw,'KEEP').gt.0) then
          read(line(7:80),*,err=20)ikeep
          if (ikeep.ge.1) then
            ikeep=1
          else
            ikeep=0
          endif
        elseif (index(keyw,'nbin').gt.0.or.index(keyw,'NBIN').gt.0) then
          read(line(7:80),*,err=20)nbin
          if (nbin.lt.4) nbin=16
        elseif (index(keyw,'ethm').gt.0.or.index(keyw,'ETHM').gt.0) then
          read(line(7:80),*,err=20)ethmax
          if (ethmax.le.1.0d-5) ethmax=4.0d0
        elseif (index(keyw,'eps').gt.0.or.index(keyw,'EPS').gt.0) then
          read(line(7:80),*,err=20)eps
          if (eps.le.0.0d0) then
            eps=3.0d-3
          elseif (eps.lt.1.0d-9) then
            eps=1.0d-9
          endif
        elseif (index(keyw,'thdat').gt.0.or.
     &          index(keyw,'THDAT').gt.0) then
          read(line(7:80),*,err=20)nmath
          do i=1,nmath
            read(inp,*,err=20)thzaid(i),nmix(i),nzam(i)
            if (nmix(i).le.0) nmix(i)=1
            if (nzam(i).le.0.or.nzam(i).gt.16) then
              write(*,*)' *** Fatal error reading acemaker input'
              write(*,*)' *** NZA value out of range ',nzam(i)
              write(*,*)' *** Check input options'
              stop
            endif
            do j=1,nzam(i)
              read(inp,*,err=20)izam(j,i)
            enddo
          enddo
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
      if (ikeep.eq.1) then
        ipndf=1
        if (iptab.eq.1) then
          iurr=1
        else
          iurr=0
        endif
      endif
      if (iace.ne.1) then
        if (imon.eq.2) imon=1
      else
        if (nmat.ne.nmath.or.nmath.le.0) then
          write(*,*)' *** Fatal error reading acemaker input'
          write(*,*)' *** NMAT and NMATH are not consistent',nmat,nmath
          write(*,*)' *** Check input options'
          stop
        endif
      endif
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
c======================================================================
