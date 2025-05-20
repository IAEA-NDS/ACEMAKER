      program dotsl
c     version 1.0
c
c     Process thermal scattering data for neutrons to prepare ACE files
c     for Monte Carlo simulations
c
c     DOTSL reads the scattering law data from ENDF-6 formatted file and
c     prepares coherent elastic, incoherent elastic and inelastic
c     scattering data in ACE-format to be used in MCNP and others
c     Monte Carlo codes
c
c     The code belongs to ACEMAKER code package devoloped by NDS/IAEA
c
c     INPUT data (DOTSL.INP):
c
c     Input data options should be entered on the DOTSL.INP text file.
c
c     line  1: Input ENDF-6 formatted filename containing TSL(MF7) (a72)
c     line  2: Output ACE-formatted filename                       (a72)
c     line  3: MAT     TEMP    NMIX   IMON              (i11,e11.0,2i11)
c     line  4: NBIN    ETHMAX  TOL    TOLE                  (i11,3e11.0)
c     line  5: THZAID  SUFF    MCNPX   NZA            (5x,a6,f11.0,2I11)
c     line 6a: (IZA(i),i=1, min(8,NZA))                            (8i7)
c     line 6a: (IZA(i),i=9, NZA)   required if NZA>8               (8i7)
c
c     where
c          MAT: Material number in the ENDF-formatted file containing
c               the thermal scattering law (TSL) on MF7 (matsl)
c         TEMP: Temperature [K]
c               (Default = 296 K)
c         NMIX: Number of atom types in mixed moderator
c               (Default = 1)
c         IMON: Monitor printing trigger (0/1/2) = min/max/max+plot
c               (Default = 0)
c         NBIN: Number of equi-probable cosines for incoherent elastic
c               and inelastic scattering (Default = 16)
c       ETHMAX: Maximun energy for thermal treatment
c               (Default = 4 eV)
c          TOL: Fractional tolerance for linearization/reconstruction
c               (Default = 0.001 = 0.1%)
c         TOLE: Fractional tolerance for preparing incident energy grid
c               (Default = 0.001 = 0.1%)
c       THZAID: Thermal ZAID name on the ACE-formatted file
c               (Default = ZA number)
c         SUFF: Thermal suffix for THID (thsuff)
c               (Default = .00)
c        MCNPX: Trigger to use extended ZAID for MCNPX
c               (0/1 = MCNP/MCNPX, Default = 0)
c          NZA: Number of ZA values of moderator components
c               (Default = 1, Max. = 16)
c       IZA(i): ZA values of the moderator components, given in two
c               lines, with a maximum of 8 values in each
c               (Max. = 16 values of ZA)
c               if NZA < 9, then the second line is not required.
c
c     Example input (next text line shows column numbers on DOTSL.INP):
c
c     12345678901234567890123456789012345678901234567891012345 (columns)
c     -----------------------------------------------------------------
c     \TSL\tsl-HinH2O.endf
c     \ACE\HinH2O.acef
c               1      293.6          1          2
c              64        4.0      0.001      0.003
c            lwtr        .00          0          1
c        1001
c     -----------------------------------------------------------------
c
c     Output files:
c       1. Output ACE-formatted file (full filename: input line 2)
c       2. DOTSL.LST listing file (fix name)
c      if imon=2
c       3. DOTSL.PLT PLOTTAB option file
c       4. DOTSL.CUR PLOTTAB curve file
c
      implicit real*8(a-h, o-z)
      parameter (bk=8.6173303d-5,ev2mev=1.0d-6)
      parameter (nfixx=150, nei0=nfixx+1 )
      parameter (ethmin=1.0d-5, ethmaxd=4.0d0, tempd=296.0d0)
      parameter (nbind=16, told=0.001d0, toled=0.001d0)
      character*1  lin130(130)
      character*3  tht
      character*4  suff
      character*6  thzaid
      character*10 hd,hm,cdate
      character*11 zsymam,str11,ctime
      character*13 hz
      character*66 line
      character*70 hk
      character*72 fin2,fout
      common/acetxt/hz,hd,hm,hk
      common/acecte/awrth,tmev,awm(16),izam(16)
      common/acepnt/nxs(16),jxs(32)
      common/acexss/xss(100000000),nxss
      allocatable ei(:)
      data nin/1/,in2/2/,lst/20/
      data thzaid/'      '/,thsuff/.00/,suff/'    '/
      data lin130/130*'='/
c
c     Initialization (nxss should be set to xss dimension)
c
      nxss=100000000
      do i=1,16
        nxs(i)=0
        izam(i)=0
        awm(i)=0.0d0
      enddo
      do i=1,32
        jxs(i)=0
      enddo
      do i=1,70
        hk(i:i)=' '
      enddo
c
c     Open main input/listing/plotting units with fix filenames
c       DOTSL.INP: DOTSL input option
c       DOTSL.LST: DOTSL listing file
c
      open (nin, file='DOTSL.INP', status='OLD')
      open (lst, file='DOTSL.LST')
      write(lst,*)
      write(lst,'(a)')' PROGRAM DOTSL: Process Thermal Scattering Law'
      write(lst,'(a)')' ============================================='
      write(lst,*)
      call getdtime(cdate,ctime)
      write(lst,'(a,a11,a,a10)')' Started at ',ctime,' on ',cdate
      write(lst,*)
      write(*,'(a)')' PROGRAM DOTSL: Process Thermal Scattering Law'
      write(*,'(a,a11,a,a10)')' Started at ',ctime,' on ',cdate
c
c     Read input data from DOTSL.INP
c
      read(nin,'(a)')fin2
      read(nin,'(a)')fout
      read(nin,'(i11,e11.0,2i11)')matsl,temp,nmix,imon
      read(nin,'(i11,3e11.0)')nbin,ethmax,tol,tole
      read(nin,'(5x,a6,f11.0,2i11)')thzaid,thsuff,mcnpx,nza
      if (nza.ge.1.and.nza.le.8) then
        read(nin,'(8i7)')(izam(i),i=1,nza)
      elseif (nza.ge.9.and.nza.le.16) then
        read(nin,'(8i7)')(izam(i),i=1,8)
        read(nin,'(8i7)')(izam(i),i=9,nza)
      elseif (nza.gt.16) then
        write(lst,'(a,i4)')' Input fatal error: NZA greater than 16',nza
        write(lst,'(a)')' Only 16 moderator components are allowable'
        write(*,'(a,i4)')' Input fatal error: NZA greater than 16',nza
        write(*,'(a)')' Only 16 moderator components are allowable'
        close(nin)
        close(lst)
        stop
      else
        write(lst,'(a,i4)')' Input fatal error: NZA less than 1',nza
        write(lst,'(a)')' At least one moderator component is needed'
        write(*,'(a,i4)')' Input fatal error: NZA less than 1',nza
        write(*,'(a)')' At least one moderator component is needed'
        close(nin)
        close(lst)
        stop
      endif
      close(nin)
c
c     Open input data file and assigning defaults, if required
c
      open (in2,file=fin2, status='OLD')
      call findmat(in2,matsl,icod)
      if (icod.ne.0) then
        write(lst,'(a,i4,a)' )' MAT=',matsl,' not found on ',fin2
        write(*,'(a,i4,a)' )' MAT=',matsl,' not found on ',fin2
        close(in2)
        close(lst)
        stop
      else
        call readcont(in2,zath,awrth,l1,l2,n1,n2,mat,mf,mt,nsi)
        izath=nint(zath)
      endif
      if (temp.le.0.0d0) temp=tempd
      if (imon.lt.0.or.imon.gt.2) imon=0
      if (nbin.lt.4) nbin=nbind
      if (ethmax.le.0.0d0) ethmax=ethmaxd
      if (tol.le.0.0d0) tol=told
      if (tole.le.0.0d0) tole=toled
      if (thzaid.eq.''.or.thzaid.eq.'      ') write(thzaid,'(i6)')izath
      if (thsuff.lt.0.0d0.or.thsuff.ge.1.0d0) thsuff=0.0d0
      if (mcnpx.ne.1) mcnpx=0
      if (nmix.le.0) nmix=1
c
c     Printing input data
c
      write(lst,'(a)')' Input options'
      write(lst,'(a)')' ============='
      write(lst,'(a,a)')' TSL-ENDF file =',fin2
      write(lst,'(a,a)')' TSL-ACE  file =',fout
      write(lst,'(a,i4)')' TSL-ENDF MAT =',matsl
      write(lst,'(a,1pe12.5)')' Temperature [K] =',temp
      write(lst,'(a,i3)')' Number of atom types in mixed moderator =',
     &  nmix
      write(lst,'(a,i2)')' Printing option =',imon
      write(lst,'(a,i6)')' Number of equiprobable cosines =',nbin
      write(lst,'(a,1pe12.5)')' Maximum thermal energy =',ethmax
      write(lst,'(a,a,1pe12.5)')' Fractional tolerance for',
     &  ' linearization =',tol
      write(lst,'(a,a,1pe12.5)')' Fractional tolerance for',
     &  ' incident energy grid =',tole
      write(lst,'(a,a)')' Thermal ZAID name = ',thzaid
      write(suff,'(a1,i3)')'.',nint(1000*thsuff)
      if (suff(2:2).eq.' ') suff(2:2)='0'
      if (suff(3:3).eq.' ') suff(3:3)='0'
      if (mcnpx.eq.1) then
        tht='nt '
        write(lst,'(a,a4)')' Thermal suffix = ',suff(1:4)
        write(lst,'(a,a3)')' MCNPX class of data identifier = ',tht
      else
        tht='t  '
        write(lst,'(a,a3)')' Thermal suffix = ',suff(1:3)
        write(lst,'(a,a3)')' MCNP class of data identifier = ',tht
      endif
      write(lst,'(a,a,i3)')' Number of ZA values of moderator',
     &  ' components =',nza
      write(lst,'(a)')' List of ZA values of moderator components:'
      write(lst,'(2(8i7))')(izam(i),i=1,nza)
c
c     Explore input tape MF1/MT451
c
      call readcont(in2,elis,sta,lis,liso,n1,nfor,mat,mf,mt,nsi)
      call readcont(in2,awi,emax,lrel,l2,nsub,nver,mat,mf,mt,nsi)
      call readcont(in2,temp1,c2,ldrv,l2,nwd,nxc,mat,mf,mt,nsi)
      if (nsub.ne.12) then
        write(lst,*)
        write(lst,*)' === Error: NSUB=12 and IPART=1 required'
        write(lst,*)'            NSUB=',nsub,' IPART=',(nsub/10)
        write(*,*)' === Error: NSUB=12 and IPART=1 required'
        write(*,*)'            NSUB=',nsub,' IPART=',(nsub/10)
        close(in2)
        close(lst)
        stop
      endif
      if (nfor.lt.6) then
        write(lst,*)
        write(lst,*)' === Error: ENDF-6 format required'
        write(lst,*)'            nfor=',nfor
        write(*,*)' === Error: ENDF-6 format required'
        write(*,*)'            nfor=',nfor
        close(in2)
        close(lst)
        stop
      endif
c
c     Preparing heading lines for ACE-formatted file
c
      call readtext(in2,line,mat,mf,mt,nsi)
      zsymam=line(1:11)
      hk(1:11)=zsymam
      if (line(11:11).ne.' '.and.line(11:11).ne.'m'.and.
     &    line(11:11).ne.'M'.and.line(12:12).ne.' ') then
        hk(12:12)=line(12:12)
      else
        hk(12:12)=' '
      endif
      hk(13:15)=' T='
      if (temp.lt.1.0d6) then
        write(str11,'(f11.2)')temp
      else
        write(str11,'(1pe11.4)')temp
      endif
      i=index(str11,' ',.true.)
      if (i.eq.0) i=1
      k=11-i
      hk(16:16+k)=str11(i:11)
      hk(17+k:24+k)=' K from '
      call readtext(in2,line,mat,mf,mt,nsi)
      call readtext(in2,line,mat,mf,mt,nsi)
      hk(25+k:42+k)=line(5:22)
      hk(43+k:53+k)=' (ACEMAKER)'
      hk(67:69)='TSL'
      str11=' '
      call getdtime(hd,str11)
      if (mcnpx.eq.1) then
        write(hz,'(a6,a4,a3)')thzaid,suff(1:4),tht
      else
        write(hz,'(a6,a3,a3,a1)')thzaid,suff(1:3),tht,' '
      endif
      write(hm,'(a6,i4)')'   mat',matsl
c
c     Prepare incident energy grid for incoherent scattering
c
      allocate(ei(nei0))
      call egrid(ei,nei,ethmax,ethmin,temp)
c
c     Thermal inelastic scattering (always present)
c
      call sigine(in2,lst,matsl,temp,nmix,nbin,ethmax,tol,tole,ei,nei,
     &  xnatom,imon)
c
c     Thermal elastic scattering
c
      call sigela(in2,lst,matsl,temp,nmix,nbin,tole,ei,nei,xnatom,imon)
c
c     Writing ACE-formatted file and its xsdir information
c
      tmev=bk*temp*ev2mev
      write(lst,*)
      write(lst,'(a)')' ACE-formatted file heading lines:'
      if (mcnpx.eq.1) then
        write(lst,'(a13,f12.6,1x,1pe11.4,1x,a10)')hz(1:13),awrth,tmev,hd
      else
        write(lst,'(a10,f12.6,1x,1pe11.4,1x,a10)')hz(1:10),awrth,tmev,hd
      endif
      write(lst,'(a70,a10)')hk,hm
      write(lst,*)
      call throut(fout,mcnpx)
      write(lst,'(a,a)')' Output ACE-formatted file: ',fout
c
c     Preparing some plots for PLOTTAB
c
      if (imon.eq.2) then
        call thrplot
      endif
c
c     Finishing
c
      write(lst,*)
      call getdtime(cdate,ctime)
      write(lst,'(a,a11,a,a10)')' DOTSL ended at ',ctime,' on ',cdate
      write(lst,*)
      write(*,'(a,a11,a,a10)')' DOTSL ended at ',ctime,' on ',cdate
      close (in2)
      close (lst)
      deallocate(ei)
      stop
      end
C======================================================================
      subroutine egrid(e,n,ethmax,ethmin,temp)
c
c      Generate the initial incident energy grid
c
      implicit real*8(a-h, o-z)
      parameter (nfixx=150,temp0=293.6d0,thigh=3000.0d0)
      dimension e(*)
      dimension efix(nfixx)
      data efix/
     &   1.00d-5,   1.78d-5,   2.50d-5,   3.50d-5,   5.00d-5,   7.00d-5,
     &   1.00d-4,   1.26d-4,   1.60d-4,   2.00d-4, .000253d0, .000297d0,
     & .000350d0, .000420d0, .000506d0, .000615d0, .000750d0, .000870d0,
     & .001000d0, .001012d0, .001230d0, .001500d0, .001800d0, .002030d0,
     & .002277d0, .002600d0, .003000d0, .003500d0, .004048d0, .004500d0,
     & .005000d0, .005600d0, .006325d0, .007200d0, .008100d0, .009108d0,
     & .010000d0, .010630d0, .011500d0, .012397d0, .013300d0, .014170d0,
     & .015000d0, .016192d0, .018200d0, .019900d0, .020493d0, .021500d0,
     & .022800d0, .025300d0, .028000d0, .030613d0, .033800d0, .036500d0,
     & .039500d0, .042757d0, .046500d0, .050000d0, .056925d0, .062500d0,
     & .069000d0, .075000d0, .081972d0, .090000d0, .096000d0, .100000d0,
     & .103500d0, .107000d0, .111573d0, .120000d0, .128000d0, .135500d0,
     & .145728d0, .152000d0, .160000d0, .172000d0, .184437d0, .190000d0,
     &.2000000d0,.2277000d0,.2400000d0,.2510392d0,.2600000d0,.2705304d0,
     &.2907501d0,.3011332d0,.3206421d0,.3576813d0,.3750000d0,.3900000d0,
     &.4170351d0,.4350000d0,.4500000d0,.4750000d0,.5032575d0,.5350000d0,
     &.5600000d0,.5800000d0,.6000000d0,.6250000d0,.7000000d0,.7400000d0,
     &.7800000d0,.8200000d0,.8600000d0,.9000000d0,.9500000d0,1.000000d0,
     &1.050000d0,1.100000d0,1.160000d0,1.220000d0,1.280000d0,1.350000d0,
     &1.420000d0,1.480000d0,1.550000d0,1.620000d0,1.700000d0,1.800000d0,
     &1.855000d0,1.900000d0,1.950000d0,2.020000d0,2.100000d0,2.180000d0,
     &2.250000d0,2.360000d0,2.450000d0,2.590000d0,2.700000d0,2.855000d0,
     &2.920000d0,3.000000d0,3.120000d0,3.270000d0,3.420000d0,3.750000d0,
     &4.070000d0,4.460000d0,4.900000d0,5.350000d0,5.850000d0,6.400000d0,
     &7.000000d0,7.650000d0,8.400000d0,9.150000d0,9.850000d0,10.00000d0/
      if (temp.gt.thigh) then
        e1=efix(1)
        re=efix(nfixx)/e1
        r=log((temp/temp0)*re)/log(re)
        do k=2,nfixx
          efix(k)=e1*exp(r*log(efix(k)/e1))
        enddo
      endif
      if (ethmin.gt.0.0d0.and.ethmin.lt.efix(1)) efix(1)=ethmin
      ethmax=fix8dig(ethmax,0)
      n=0
      do k=1,nfixx
        ek=fix8dig(efix(k),0)
        if (ek.gt.ethmax) exit
        n=n+1
        e(n)=ek
      enddo
      if (e(n).lt.ethmax) then
        n=n+1
        e(n)=ethmax
      endif
      return
      end
C======================================================================
      subroutine sigine(in2,lst,matsl,temp,nmix,nbin,ethmax,tol,tole,
     &  ei,nei,xnatom,imon)
c
c      processing of inelastic thermal scattering
c
      implicit real*8(a-h, o-z)
      parameter (nemax=2000,nmax=7000, nkks=20)
      parameter (slgmin=-228.0d0)
      parameter (epmin=0.0d0, tolde=1.0d-6)
      parameter (tzref=0.0253d0, bk=8.6173303d-5)
      dimension ei(*)
      common/acecte/awrth,tmev,awm(16),izam(16)
      common/acepnt/nxs(16),jxs(32)
      common/acexss/xss(100000000),nxss
      dimension nbt(20),ibt(20)
      character*1 lin130(130)
      allocatable sigb(:),xat(:),aws(:),isl(:),teff(:)
      allocatable alpha(:),beta(:),sab(:,:)
      allocatable x(:),y(:),ubar(:,:)
      allocatable xs1(:),ys1(:),ubar1(:,:)
      allocatable xs(:,:),ys(:,:),ubars(:,:,:)
      allocatable es(:),sigss(:),uus(:),neps(:)
      allocatable w0(:),w1(:)
      data lin130/130*'='/
c
c     Initio (the order of smin is 1.0e-99)
c
      smin=exp(slgmin)
      allocate(x(nmax),y(nmax),ubar(nbin,nmax))
      allocate(xs1(nmax),ys1(nmax),ubar1(nbin,nmax))
      allocate(xs(nmax,nkks),ys(nmax,nkks),ubars(nbin,nmax,nkks))
      allocate(es(nkks),sigss(nkks),uus(nkks),neps(nkks))
      allocate(w0(nbin),w1(nbin))
c
c     Read thermal inelastic scattering law from ENDF-6 file(MF7/MT4)
c
      call findmt(in2,matsl,7,4,icod)
      if (icod.ne.0) then
        write(lst,*)
        write(lst,'(a,a,i4)')' === Fatal error: Inelastic scattering',
     &    ' MT=4 not found on TSL unit',in2
        write(*,'(a,a,i4)')' === Fatal error: Inelastic scattering',
     &    ' MT=4 not found on TSL unit',in2
        close(in2)
        close(lst)
        stop
      endif
      write(lst,*)
      write(lst,'(a)')' Thermal inelastic scattering'
      write(lst,*)
      write(*,*)' Thermal inelastic scattering'
c
c     Get general information on moderator atoms and TSL
c
      call readcont(in2,za,awr,l1,lat,lasym,n2,mat,mf,mt,nsi)
      tz=bk*temp
      if (lat.eq.1) then
        tz0=tzref
      else
        tz0=tz
      endif
      call readlist(in2,c1,c2,lln,l2,ni,ns,y)
      nsa=ns+1
      if (ni.ne.6*nsa.or.nsa.lt.1) then
        write(lst,*)' Incorrect ni value ni=',ni,' ns=',ns
        write(*,*)' Incorrect ni value ni=',ni,' ns=',ns
        close(in2)
        close(lst)
        stop
      endif
      allocate(sigb(nsa),xat(nsa),aws(nsa),isl(nsa),teff(nsa))
      if (y(1).gt.0.0d0) then
        isl(1)=-1
        ehigh=y(2)
        aws(1)=y(3)
        elim=y(4)
        xat(1)=y(6)
        awrth=aws(1)
        xnatom=xat(1)
        fb=(awrth+1.0d0)/awrth
        sigb(1)=y(1)*fb*fb/xnatom
        i0=1
      else
        isl(1)=0
        aws(1)=y(3)
        xat(1)=y(6)
        xnatom=1.0d0
        fb=(aws(1)+1.0d0)/aws(1)
        sigb(1)=y(2)*fb*fb/xat(1)
      endif
      if (ns.gt.0) then
        do k=1,ns
          j0=6*k
          i=k+1
          isl(i)=nint(y(j0+1)+1.0d-5)
          aws(i)=y(j0+3)
          xat(i)=y(j0+6)
          fb=(aws(i)+1.0d0)/aws(i)
          sigb(i)=y(j0+2)*fb*fb
          if (isl(1).lt.0.and.isl(i).eq.0) then
            sigb(i)=sigb(i)/xnatom
          else
            sigb(i)=sigb(i)/xat(i)
          endif          
        enddo
      endif
      if (isl(1).lt.0) then
c
c       Tabulated Sab given for principal atom
c
        ib0=0
        call readtab2(in2,c1,c2,l1,l2,nr,nb,nbt,ibt)
c
c       Reading first beta at requested temperature
c
        call readtab1(in2,temp1,b0,lt,l2,nr,na,nbt,ibt,x,y)
        if (na.gt.nmax) then
          write(lst,'(a)')' === Error: Buffer memory is not enough'
          write(lst,'(a,i6)')'     Increase nmax at least to ',na
          write(*,'(a)')' === Error: Buffer memory is not enough'
          write(*,'(a,i6)')'     Increase nmax at least to ',na
          close(in2)
          close(lst)
          stop
        endif
        allocate (alpha(na),beta(nb),sab(na,nb))
        deltmax=dtemp(temp,0)
        itemp=0
        nt=lt+1
        do k=1,nt
          if (k.ne.1) then
            call readlist(in2,temp1,bb,li,l2,np,n2,y)
            if (np.ne.na) then
              write(lst,'(a)')' === Error: TSL format error NP <> NA'
              write(lst,'(a,i6,a,i6)')'            NP=',np, '  NA=',na
              write(*,'(a)')' === Error: TSL format error NP <> NA'
              write(*,'(a,i6,a,i6)')'            NP=',np, '  NA=',na
              close(in2)
              close(lst)
              stop
            endif
          endif
          delt=abs(temp1-temp)
          if (delt.le.deltmax) then
            do j=1,na
              sab(j,1)=y(j)
            enddo
            deltmax=delt
            itemp=k
            temp0=temp1
          endif
        enddo
        if (itemp.gt.0) then
          beta(1)=b0
          if (b0.eq.0.0d0) ib0=1
          do j=1,na
            alpha(j)=x(j)
          enddo
          if (abs(temp-temp0).gt.0.001d0) then
            write(lst,'(a,a,1pe12.5,a)')' === Warning:',
     &        ' Requested temperature = ',temp,' K not found'
            write(lst,'(a,a,1pe12.5,a)')'             ',
     &        ' Closest temperature found = ',temp0,
     &        ' K will be used'
          else
            write(lst,'(a,1pe12.5,a)')' Requested temperature = ',
     &        temp,' K (found)'
          endif
        else
          write(lst,'(a,a,1pe12.5,a)')' === Error:',
     &      ' Requested temperature = ',temp,' K not found'
          write(*,'(a,a,1pe12.5,a)')' === Error:',
     &      ' Requested temperature = ',temp,' K not found'
          close(in2)
          close(lst)
          stop
        endif
c
c       Reading rest of beta values at requested temperature
c
        do i=2,nb
          call readtab1(in2,temp1,beta(i),lt,l2,nr,np,nbt,ibt,x,y)
          if (ib0.eq.0.and.beta(i).eq.0.0d0) ib0=i
          if ((lt+1).ne.nt.or.na.ne.np) then
            if (np.ne.na) then
              write(lst,'(a)')' === Error: TSL format error NP <> NA'
              write(lst,'(a,i6,a,i6)')'            NP=',np, '  NA=',na
              write(*,'(a)')' === Error: TSL format error NP <> NA'
              write(*,'(a,i6,a,i6)')'            NP=',np, '  NA=',na
            endif
            if ((lt+1).ne.nt) then
              write(lst,'(a,a,2i5)')' === Error: Different number of',
     &          ' temperatures found for beta ',lt+1,nt
              write(*,'(a,a,2i5)')' === Error: Different number of',
     &          ' temperatures found for beta ',lt+1,nt
            endif
            close(in2)
            close(lst)
            stop
          endif
          do k=1,nt
            if (k.ne.1) then
              call readlist(in2,temp1,bb,li,l2,np,n2,y)
              if (np.ne.na) then
                write(lst,'(a)')' === Error: TSL format error NP <> NA'
                write(lst,'(a,i6,a,i6)')'            NP=',np, '  NA=',na
                write(*,'(a)')' === Error: TSL format error NP <> NA'
                write(*,'(a,i6,a,i6)')'            NP=',np, '  NA=',na
                close(in2)
                close(lst)
                stop
              endif
            endif
            if (k.eq.itemp) then
              if (abs(temp1-temp0).le.delt) then
                do j=1,na
                  sab(j,i)=y(j)
                enddo
              else
                write(lst,'(a,a,1p2e13.5)')' === Error: Different',
     &            ' temperature set found for beta',temp0,temp1
                write(*,'(a,a,1p2e13.5)')' === Error: Different',
     &            ' temperature set found for beta',temp0,temp1
                close(in2)
                close(lst)
                stop
              endif
            endif
          enddo
        enddo
c
c       Check and modify Sab to deal with potential numerical problems
c
        write(lst,*)
        write(lst,'(2(a,i5))')' alpha values =',na,'  beta values =',nb
        nsmall=0
        if (lln.eq.1) then
          do i=1,nb
            do j=1,na
              if (sab(j,i).eq.0.0d0.or.sab(j,i).le.slgmin) then
                if (sab(j,i).ne.0.0d0) nsmall=nsmall+1
                sab(j,i)=slgmin
              endif
            enddo
          enddo
        else
          do i=1,nb
            do j=1,na
              if (sab(j,i).eq.0.0d0.or.sab(j,i).le.smin) then
                if (sab(j,i).ne.0.0d0) nsmall=nsmall+1
                sab(j,i)=slgmin
              else
                sab(j,i)=log(sab(j,i))
              endif
            enddo
          enddo
        endif
        if (nsmall.gt.0) then
          write(lst,*)
          write(lst,'(a,i9,a,i10,a)')' === Warning: ',nsmall,
     &      ' small values found out of ',na*nb,' Sab values.'
          write(lst,'(14x,a,a,1pe12.5,a)')
     &      ' It could affect energy transfers greater',
     &      ' than ',2.0d0*(-(slgmin+4.0d0))*tz0,' eV'
        endif
        etmax=beta(nb)*tz0
        if (etmax.lt.ethmax) then
          write(lst,*)
          write(lst,'(a,a,1pe12.5,a,a)')
     &      ' Sab tabulated data will be used for energy transfers',
     &      ' below ',etmax,' eV. SCT approximation will be applied',
     &      ' otherwise'
        endif
c
c       Trigger special issue for interpolating in liquid materials
c
        if (sab(1,ib0).gt.sab(2,ib0)) then
          liq=ib0
          write(lst,*)
          write(lst,'(a,a,1pe15.8,a,i5)')' liquid behavior detected,',
     &      ' beta=',beta(liq),' for ib=',liq
        else
          liq=0
        endif
      else
c
c       No principal atom. Use analytic expressions. a-priori beta grid
c       is generated for adaptively reconstruction of Sab
c
        temp0=temp
        liq=0
        na=1
        nb=200
        allocate(beta(nb),alpha(na),sab(na,nb))
        beta(1)=0.0d0
        beta(2)=max(1.0d-4*tz0/tzref,1.0d-5)
        beta(nb)=ethmax/tz0
        rb=exp(log(beta(nb)/beta(2))/dble(nb-2))
        do i=3,nb-1
          beta(i)=beta(i-1)*rb
        enddo
      endif
c
c     Reading effective temperatures for SCT
c
      delt=dtemp(temp,1)
      do i=1,nsa
        if (isl(i).le.0) then
          call readtab1(in2,c1,c2,l1,l2,nr,np,nbt,ibt,x,y)
          if ((x(1)-delt).gt.temp.or.(x(np)+delt).lt.temp) then
            write(lst,*)
            write(lst,'(a,a)')' === Error: Requested temperature is',
     &        ' out of effective temperature range'
            write(*,'(a,a)')' === Error: Requested temperature is',
     &        ' out of effective temperature range'
            close(in2)
            close(lst)
            stop
          elseif (x(1).ge.temp) then
            teff(i)=y(1)
          elseif (x(np).le.temp) then
            teff(i)=y(np)
          else
            teff(i)=fvalue(nr,nbt,ibt,np,x,y,temp)
          endif
        else
          teff(i)=temp
        endif
      enddo
      write(lst,*)
      write(lst,'(1x,a,7x,a,6x,a,6x,a,7x,a,4x,a)')
     &  'ATOM','SIGB','NATOM','AWR','TSL','Teff(T)'
      do i=1,nsa
        if (isl(i).lt.1) then
          write(lst,'(i3,1pe18.8,i4,e16.6,i4,e15.6)')
     &      i,sigb(i),nint(xat(i)+0.01d0),aws(i),isl(i),teff(i)
        else
          write(lst,'(i3,1pe18.8,i4,e16.6,i4)')
     &      i,sigb(i),nint(xat(i)+0.01d0),aws(i),isl(i)
        endif
      enddo
c
c      initial ACE pointers for inelastic scattering
c
      itie=1
      itix0=1+itie+nemax
      itxe0=itix0+nemax
      itnep0=itxe0+nemax
      ixss=itnep0+nemax-1
c
c      Compute initial number of secondary energies nodos
c      Initialize incident energies cycle
c
      if (lasym.eq.1) then
        nbb=nb
      else
        nbb=2*nb-1
      endif
      je=0
      tolu=2.0d0*tol
      do ie=1,nei
c
c       Initial incident energies cycle
c
        e=ei(ie)
        ke=1
        do while (ke.gt.0)
c
c       Incident energy(e) grid adaptively reconstruction
c
          x1m=fix8dig(e,-1)
          x1p=fix8dig(e,1)
c
c         Secondary-energy(ep) grid adaptively reconstruction
c
          x0=epmin
          call flinmu(y0,w0,nbin,e,x0,temp,tz,tz0,lasym,liq,
     &      na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,slgmin,tol)
          j=1
          x(1)=x0
          y(1)=y0
          do l=1,nbin
            ubar(l,1)=w0(l)
          enddo
          do ib=1,nbb
            if(lasym.eq.0) then
              ibb=ib-nb
              isg=sign(1,ibb)
              b=dble(isg)*beta(isg*ibb+1)
            else
              b=beta(ib)
            endif
            if (b.eq.0.0d0) then
              ep=e
              x1=x1m
            else
              ep=e+b*tz0
              x1=fix8dig(ep,0)
            endif
            if (x1.gt.x0.and.((x1.lt.x1m.and.b.lt.0.0d0).or.
     &         (x1.gt.x1p.and.b.gt.0.0d0).or.(b.eq.0.0d0))) then
              call flinmu(y1,w1,nbin,e,x1,temp,tz,tz0,lasym,liq,
     &          na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,slgmin,tol)
              call flinep(x0,y0,w0,x1,y1,w1,nbin,j,x,y,ubar,nmax,e,temp,
     &          tz,tz0,lasym,liq,na,alpha,nb,beta,sab,nsa,isl,sigb,aws,
     &          teff,slgmin,tol)
              if (b.eq.0.0d0) then
                x1=x1p
                call flinmu(y1,w1,nbin,e,x1,temp,tz,tz0,lasym,liq,
     &            na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,slgmin,tol)
                call inc(j,nmax)
                x(j)=x1
                y(j)=y1
                do l=1,nbin
                  ubar(l,j)=w1(l)
                enddo
              endif
              x0=x1
              y0=y1
              do l=1,nbin
                w0(l)=w1(l)
              enddo
            endif
          enddo
          nep=j
c
c         Removing leading zeros, if any
c
          iep0=nep
          do i=2,nep
            if (y(i).gt.0.0d0) then
              iep0=i-1
              exit
            endif
          enddo
          if (iep0.lt.nep) then
c
c           Removing trailing zeros, if any
c
            knep=nep-1
            do i=knep,iep0,-1
              if (y(i).gt.0.0d0) then
                nep=i+1
                exit
              endif
            enddo
            iep0=iep0-1
            nep=nep-iep0
            do i=1,nep
              x(i)=x(iep0+i)
              y(i)=y(iep0+i)
              do j=1,nbin
                ubar(j,i)=ubar(j,iep0+i)
              enddo
            enddo
c
c           Calculating cross section and average cosine at e
c
            sigs=0.0d0
            uu=0.0d0
            x0=x(1)
            y0=y(1)
            u0=avecos(ubar(1,1),nbin)
            do i=2,nep
              x1=x(i)
              y1=y(i)
              u1=avecos(ubar(1,i),nbin)
              sigsi=0.5d0*(x1-x0)*(y1+y0)
              sigs=sigsi+sigs
              uu=0.5d0*(u0+u1)*sigsi+uu
              x0=x1
              y0=y1
              u0=u1
            enddo
            uu=uu/sigs
          else
c
c           Null cross section at incident energy e
c
            nep=3
            x(1)=fix8dig(e,-1)
            y(1)=0.0d0
            x(2)=fix8dig(e,0)
            y(2)=2.0d0/(x(3)-x(1))
            x(3)=fix8dig(e,1)
            y(3)=0.0d0
            do i=1,nep
              do j=1,nbin
                ubar(j,i)=1.0d0
              enddo
            enddo
            sigs=0.0d0
            uu=0.0d0
          endif
c
c         Printing control
c
c          write(lst,'(130a1)')lin130
c          do i=iep0,nep
c            call printine(lst,i,i,x(i),y(i),y(i),ubar,nbin)
c          enddo
c          write(lst,'(130a1)')lin130
c
c         Checking incident energy grid adaptively reconstruction
c
          if (ie.eq.1) then
            e0=e
            sigs0=sigs
            uu0=uu
            call thrload(lst,je,e,sigs,nep,nbin,x,y,ubar,nmix,
     &        imon,itie,itix0,itxe0,itnep0,ixss,nemax)
            ke=0
          else
            if (ke.eq.1) then
              e1=e
              sigs1=sigs
              uu1=uu
              nep1=nep
              do i=1,nep
                xs1(i)=x(i)
                ys1(i)=y(i)
                do j=1,nbin
                  ubar1(j,i)=ubar(j,i)
                enddo
              enddo
              e=0.5d0*(e0+e1)
              e=fix8dig(e,0)
              ke=2
              kk=0
            else
              if (e.le.e0.or.e.ge.e1) then
                dsige=0.0d0
                duu=0.0d0
                slope0=1.0d0
                slope1=1.0d0
                slopu0=1.0d0
                slopu1=1.0d0
                tests=1.0d0
                testu=1.0d0
              else
                sigsm=0.5d0*(sigs0+sigs1)
                uum=0.5d0*(uu0+uu1)
                dsige=abs(sigsm-sigs)
                duu=abs(uum-uu)
                de0=e-e0
                de1=e1-e
                slope0=(sigs-sigs0)/de0
                slope1=(sigs1-sigs)/de1
                slopu0=(uu-uu0)/de0
                slopu1=(uu1-uu)/de1
                tests=tole*abs(sigs)+1.0d-12
                testu=tolu*abs(uu)+1.0d-6
              endif
              if (((dsige.lt.tests.and.duu.lt.testu).and.
     &          (slope0*slope1.gt.0.0d0.and.slopu0*slopu1.gt.0.0d0)).or.
     &          ((e1-e0).le.tolde*e1).or.(kk.ge.nkks)) then
                e0=e1
                sigs0=sigs1
                uu0=uu1
                call thrload(lst,je,e1,sigs1,nep1,nbin,xs1,ys1,ubar1,
     &            nmix,imon,itie,itix0,itxe0,itnep0,ixss,nemax)
                if (kk.eq.0) then
                  ke=0
                else
                  e1=es(kk)
                  sigs1=sigss(kk)
                  uu1=uus(kk)
                  nep1=neps(kk)
                  do i=1,nep1
                    xs1(i)=xs(i,kk)
                    ys1(i)=ys(i,kk)
                    do j=1,nbin
                      ubar1(j,i)=ubars(j,i,kk)
                    enddo
                  enddo
                  kk=kk-1
                  e=0.5d0*(e0+e1)
                  e=fix8dig(e,0)
                endif
              else
                kk=kk+1
                es(kk)=e1
                sigss(kk)=sigs1
                uus(kk)=uu1
                neps(kk)=nep1
                do i=1,nep1
                  xs(i,kk)=xs1(i)
                  ys(i,kk)=ys1(i)
                  do j=1,nbin
                    ubars(j,i,kk)=ubar1(j,i)
                  enddo
                enddo
                e1=e
                sigs1=sigs
                uu1=uu
                nep1=nep
                do i=1,nep
                  xs1(i)=x(i)
                  ys1(i)=y(i)
                  do j=1,nbin
                    ubar1(j,i)=ubar(j,i)
                  enddo
                enddo
                e=0.5d0*(e0+e1)
                e=fix8dig(e,0)
              endif
            endif
          endif
c
c         loop over incident energy grid adaptively reconstruction
c
        enddo
c
c       loop over initial incident energies
c
      enddo
c
c       Check XSS dimension
c
      if (ixss.gt.nxss) then
        write(lst,'(a,a,i10,a,i10)')' === Fatal error: Not',
     &    ' enough memory. Increase XSS array size in more than',
     &    ixss-nxss,' units. Current size ',nxss
        write(*,'(a,a,i10,a,i10)')' === Fatal error: Not',
     &    ' enough memory. Increase XSS array size in more than',
     &    ixss-nxss,' units. Current size ',nxss
        close(in2)
        close(lst)
        stop
      endif
c
c      repacking XSS array and correct pointers
c
      if (je.ne.nemax) then
        itix=1+itie+je
        itxe=itix+je
        itnep=itxe+je
        itxss=itnep+je-1
        itxss0=itnep0+nemax-1
        nde=itxss0-itxss
        nxss0=ixss-itxss0
        do i=1,je
          xss(itix+i-1)=xss(itix0+i-1)
        enddo
        do i=1,je
          xss(itxe+i-1)=nint(xss(itxe0+i-1)+1.0d-6)-nde
        enddo
        do i=1,je
          xss(itnep+i-1)=xss(itnep0+i-1)
        enddo
        do i=1,nxss0
          xss(itxss+i)=xss(itxss0+i)
        enddo
        ixss=itxss+nxss0
      else
        itix=itix0
        itxe=itxe0
      endif
c
c     Assigning triggers and pointers for inelastic scattering
c
      nxs(1)=ixss
      nxs(2)=3
      nxs(3)=nbin+1
      nxs(4)=nbin
      nxs(7)=2
      jxs(1)=itie
      jxs(2)=itix
      jxs(3)=itxe
      write(lst,*)
      write(lst,'(a,2i10)')' Length of inelastic data and XSS array: ',
     &  ixss,ixss
      deallocate (beta,alpha,sab)
      deallocate (sigb,xat,aws,isl,teff)
      deallocate(x,y,ubar)
      deallocate(xs1,ys1,ubar1)
      deallocate(xs,ys,ubars)
      deallocate(es,sigss,uus,neps)
      deallocate (w0,w1)
      return
      end
C======================================================================
      subroutine flinep(x0,y0,w0,x1,y1,w1,nbin,j,x,y,ubar,nmax,e,temp,
     &  tz,tz0,lasym,liq,na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,
     &  arglim,tol)
      implicit real*8 (a-h, o-z)
      parameter (ns=22, azero=1.0d-10, xszero=1.0d-12, uzero=1.0d-6)
      dimension alpha(*),beta(*),sab(na,*),isl(*),sigb(*),aws(*),teff(*)
      dimension w0(*),w1(*),x(*),y(*),ubar(nbin,*)
      dimension wm(nbin)
      dimension xs(ns),ys(ns),zs(nbin,ns)
      tolw=2.0d0*tol
      k=0
      nostop=1
      do while (nostop.eq.1)
        yl=0.5d0*(y0+y1)
        area0=(x1-x0)*yl
        xm=0.5d0*(x0+x1)
        xm=fix8dig(xm,0)
        iflag=0
        if (area0.lt.azero.or.xm.le.x0.or.xm.ge.x1.or.k.eq.ns) then
          iflag=1
        else
          call flinmu(ym,wm,nbin,e,xm,temp,tz,tz0,lasym,liq,
     &      na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim,tol)
          if (abs(yl-ym).lt.(tol*abs(ym)+xszero)) then
            mflag=0
            sumu=0.0d0
            sumw=0.0d0
            do imu=1,nbin
              um=0.5d0*(w0(imu)+w1(imu))
              wmm=wm(imu)
              if (abs(um-wmm).gt.tol) then
                mflag=imu
                exit
              endif
              sumu=sumu+um
              sumw=sumw+wmm
            enddo
            if (mflag.eq.0) then
              if (abs(sumu-sumw).lt.(tolw*abs(sumw)+uzero)) iflag=2
            endif
          endif
        endif
        if (iflag.gt.0) then
          call inc(j,nmax)
          x(j)=x1
          y(j)=y1
          do l=1,nbin
            ubar(l,j)=w1(l)
          enddo
          if (k.eq.0) then
            nostop=0
          else
            x0=x1
            y0=y1
            do l=1,nbin
              w0(l)=w1(l)
            enddo
            x1=xs(k)
            y1=ys(k)
            do l=1,nbin
              w1(l)=zs(l,k)
            enddo
            k=k-1
          endif
        else
          k=k+1
          xs(k)=x1
          ys(k)=y1
          do l=1,nbin
            zs(l,k)=w1(l)
          enddo
          x1=xm
          y1=ym
          do l=1,nbin
            w1(l)=wm(l)
          enddo
        endif
      enddo
      return
      end
C======================================================================
      subroutine flinmu(f0,ubar,nmu,ein,eou,temp,tz,tz0,lasym,liq,
     &  na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim,tol)
      implicit real*8 (a-h, o-z)
      real*16 fbin,fbinlo,sf0,sf1,sf0i,sf1i,u1,v1,u2,v2,v0,du,dv,slope
      real*16 dsf0,p1,pm,sd,d2,one,one3
      parameter(nmumax=7000,nmu00=4,one=1.0d0,one3=1.0d0/3.0d0)
      parameter(tolmu=1.0d-5, rtolmu=1.0d-8, d2min=1.0d-9)
      parameter(vtol=1.0d-6, vtol2=vtol/(0.5d0+vtol))
      parameter(f0min=1.0d-30, tollow=0.999999999d0)
      dimension ubar(*)
      dimension alpha(*),beta(*),sab(na,*),isl(*),sigb(*),aws(*),teff(*)
      dimension xmu00(nmu00)
      dimension xmu(nmumax),fmu(nmumax)
      allocatable xmu0(:)
      data xmu00/-1.00d0,0.0d0,0.99d0,1.00d0/
c
c     a-priori cosine grid preparation
c
      sqei=sqrt(ein)
      sqeo=sqrt(eou)
      sqme=sqei*sqeo
      allocate(xmu0(4+nsa+na))
      ymax=tol
      do imu0=1,nmu00
        x0=xmu00(imu0)
        x0=fix8dig(x0,0)
        xmu0(imu0)=x0
        y0=sigtsl(ein,eou,x0,temp,tz,tz0,lasym,liq,
     &     na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim)
        ymax=max(y0,ymax)
      enddo
      imu0=nmu00
      if (na.gt.1.and.sqme.gt.0.0d0) then
c
c       Calculate possible cosines from the alpha grid
c
        c0=aws(1)*tz0
        esm=sqei-sqeo
        alfmin=esm*esm/c0
        esm=sqei+sqeo
        alfmax=esm*esm/c0
        esum=eou+ein
        do i=na,1,-1
          alfa=alpha(i)
          if (alfa.gt.alfmin.and.alfa.lt.alfmax) then
            x0=0.5d0*(esum-c0*alfa)/sqme
            x0=fix8dig(x0,0)
            if (x0.gt.-1.0d0.and.x0.lt.1.0d0) then
              imu0=imu0+1
              xmu0(imu0)=x0
              y0=sigtsl(ein,eou,x0,temp,tz,tz0,lasym,liq,
     &          na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim)
              ymax=max(y0,ymax)
            endif
          endif
        enddo
      endif
c
c     Add the epithermal static scattering cosine for each atom
c
      if (sqme.gt.0.0d0) then
        c1=0.5d0*(eou-ein)/sqme
        c2=0.5d0*(ein+eou)/sqme
        do i=1,nsa
          alf=(aws(i)-1)/(aws(i)+1)
          alf=alf*alf
          if (eou.gt.alf*ein.and.eou.lt.ein) then
            x0=aws(i)*c1+c2
            x0=fix8dig(x0,0)
            if (x0.gt.-1.0d0.and.x0.lt.1.0d0) then
              imu0=imu0+1
              xmu0(imu0)=x0
              y0=sigtsl(ein,eou,x0,temp,tz,tz0,lasym,liq,
     &           na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim)
              ymax=max(y0,ymax)
            endif
          endif
        enddo
      endif
      nmu0=imu0
c
c     Ordering cosines in ascending order and removing repeated values
c
      call orderx(xmu0,nmu0,rtolmu,1)
      if (xmu0(1).ne.-1.0d0) xmu0(1)=-1.0d0
      if (xmu0(nmu0).ne.1.0d0) xmu0(nmu0)=1.0d0
c
c     Adaptively cosine grid preparation for lin-lin interpolation
c
      dymax=0.01d0*ymax
      toly=0.5d0*tol
      x0=xmu0(1)
      y0=sigtsl(ein,eou,x0,temp,tz,tz0,lasym,liq,
     &   na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim)
      j=1
      xmu(1)=x0
      fmu(1)=y0
      do i=2,nmu0
        x1=xmu0(i)
        y1=sigtsl(ein,eou,x1,temp,tz,tz0,lasym,liq,
     &     na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim)
        call flin(x0,y0,x1,y1,j,xmu,fmu,nmumax,
     &    ein,eou,temp,tz,tz0,lasym,liq,
     &    na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim,
     &    toly,tolmu,dymax)
        x0=x1
        y0=y1
      enddo
      mmu=j
c
c     Removing leading and trailing zeros, if any
c
      kmu=mmu-1
      im0=mmu
      do i=1,kmu
        if (fmu(i+1).gt.0.0d0) then
          im0=i
          exit
        endif
      enddo
      if (im0.lt.mmu) then
        do i=kmu,im0,-1
          if (fmu(i).gt.0.0d0) then
            mmu=i+1
            exit
          endif
        enddo
      endif
c
c     Calculation of f0=f(Ein,Eou)
c
      f0=0.0d0
      i0=im0+1
      if (i0.le.mmu) then
        do i=i0,mmu
          i1=i-1
          f0=0.5d0*(xmu(i)-xmu(i1))*(fmu(i)+fmu(i1))+f0
        enddo
      endif
c
c     Calculation of equiprobable cosines
c       mu(imu)=cosine(Ein,Eou,imu) , imu=1,nmu.
c       pdf(mu(imu))=1/nmu
c
      if (f0.gt.f0min) then
        fbin=f0/dble(nmu)
        fbinlo=tollow*fbin
        nmu1=nmu-1
        imu=0
        sf0=0.0d0
        sf1=0.0d0
        u1=xmu(im0)
        v1=fmu(im0)
        i=im0+1
        do while(i.le.mmu)
          u2=xmu(i)
          v2=fmu(i)
          du=u2-u1
          dv=v2-v1
          slope=dv/du
          v0=v1-slope*u1
          sf0i=0.5d0*du*(v2+v1)
          sf0=sf0i+sf0
          sf1i=du*(slope*one3*(u2*u2+u1*u2+u1*u1)+0.5d0*v0*(u2+u1))
          sf1=sf1i+sf1
          if (imu.lt.nmu1.and.sf0.ge.fbinlo) then
            if (sf0.gt.fbin) then
              dsf0=fbin+sf0i-sf0
              sf1=sf1-sf1i
              if (abs(slope).lt.abs(vtol2*v1).and.v1.gt.f0min) then
                d2=dsf0/v1
                icod=1
              else
                p1=v1/slope
                pm=p1-u1
                sd=pm*pm+(u1*(v1+v0)+2.0d0*dsf0)/slope
                if (sd.lt.0.0d0) then
                  sdm=abs(sd)
                  if (sdm.le.tolmu) then
                    write(*,'(a,a,1pe18.10,a,i3)')
     &                '   Warning: Small negative discriminant.',
     &                ' Set to its absolute value ',sdm,
     &                ' for cosine interval i=',imu+1
                    sd=sdm
                  else
                    write(*,'(a,i5,a)')' === Fatal error calculating',
     &                nmu,' equiprobable cosines. No real roots found'
                    stop
                  endif
                endif
                sd=sqrt(sd)
                d2=sign(one,slope)*sd-p1
                icod=2
              endif
              if (d2.gt.du) then
                if ((d2-du).lt.tolmu) then
                  write(*,'(a,a,i4,a,i3)')'   Warning: Solution',
     &              ' too close to u2, i=',imu+1,' icod=',icod
                  d2=tollow*du
                else
                  write(*,'(a,i5,a)')' === Fatal error calculating',
     &              nmu,' equiprobable cosines. Root out of range'
                  write(*,*)' u1,u2,d2,du,cod,imu = ',u1,u2,d2,du,
     &              icod,imu+1,i,mmu
                  stop
                endif
              elseif (d2.le.0.0d0) then
                if (abs(d2).lt.tolmu) then
                  write(*,'(a,a,i4,a,i3)')'   Warning: Solution',
     &              ' too close to u1, i=',imu+1,' icod=',icod
                  d2=d2min*u1
                else
                  write(*,'(a,i5,a)')' === Fatal error calculating',
     &              nmu,' equiprobable cosines. Root out of range'
                  write(*,*)' u1,u2,d2,du,cod,imu = ',u1,u2,d2,du,
     &              icod,imu+1,i,mmu
                  stop
                endif
              endif
              u2=u1+d2
              v2=slope*d2+v1
              sf1i=d2*(slope*one3*(u2*u2+u1*u2+u1*u1)+0.5d0*v0*(u2+u1))
              sf1=sf1i+sf1
              i=i-1
            endif
            imu=imu+1
            ubar(imu)=sf1/fbin
            sf0=0.0d0
            sf1=0.0d0
          endif
          u1=u2
          v1=v2
          i=i+1
        enddo
        ubar(nmu)=sf1/fbin
c
c       Checking equi-probable cosines
c
        call chkcos(ubar,nmu,ineg,ipos,uneg,upos)
        tolwrt=1.0d0+0.5*tol
        if (ineg.gt.1.and.uneg.le.-tolwrt) then
          write(*,'(a,1pe14.7,a,e14.7,a,i4,a,0pf7.4,a)')
     &      '   Warning: ein=',ein,' eou=',eou,' nn=',
     &      ineg,' umin=',uneg,' corrected'
        endif
        if (ipos.gt.1.and.upos.ge.tolwrt) then
          write(*,'(a,1pe14.7,a,e14.7,a,i4,a,0pf7.4,a)')
     &      '   Warning: ein=',ein,' eou=',eou,' np=',
     &      ipos,' umax=',upos,' corrected'
        endif
      else
        f0=0.0d0
        do imu=1,nmu
          ubar(imu)=0.0d0
        enddo
      endif
      deallocate(xmu0)
      return
      end
C======================================================================
      subroutine flin(x0,y0,x1,y1,j,x,y,npmax,
     &  ein,eou,temp,tz,tz0,lasym,liq,
     &  na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim,
     &  tol,tolmu,dymax)
      implicit real*8 (a-h, o-z)
      parameter (ns=32, dumax=0.5d0)
      dimension alpha(*),beta(*),sab(na,*),isl(*),sigb(*),aws(*),teff(*)
      dimension x(*),y(*)
      dimension xs(ns),ys(ns)
      dymmax=2.0d0*dymax
      k=0
      nostop=1
      do while (nostop.eq.1)
        h=x1-x0
        yl=0.5d0*(y0+y1)
        xm=0.5d0*(x0+x1)
        xm=fix8dig(xm,0)
        iflag=0
        if (h.lt.tolmu.or.k.eq.ns.or.xm.le.x0.or.xm.ge.x1) then
          iflag=1
        else
          ym=sigtsl(ein,eou,xm,temp,tz,tz0,lasym,liq,
     &       na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim)
          if (abs(yl-ym).lt.tol*abs(ym+dymmax).and.h.lt.dumax.and.
     &        abs(y1-y0).lt.abs(yl+dymax)) iflag=2
        endif
        if (iflag.gt.0) then
          call inc(j,npmax)
          x(j)=x1
          y(j)=y1
          if (k.eq.0) then
            nostop=0
          else
            x0=x1
            y0=y1
            x1=xs(k)
            y1=ys(k)
            k=k-1
          endif
        else
          k=k+1
          xs(k)=x1
          ys(k)=y1
          x1=xm
          y1=ym
        endif
      enddo
      return
      end
C======================================================================
      real*8 function sigtsl(ein,eou,u,temp,tz,tz0,lasym,liq,
     &       na,alpha,nb,beta,sab,nsa,isl,sigb,aws,teff,arglim)
      implicit real*8 (a-h, o-z)
      parameter (amin0=1.0d-6, smin0=1.0d-16)
      dimension alpha(*),beta(*),sab(na,*),isl(*),sigb(*),aws(*),teff(*)
      aws0=aws(1)
      b=(eou-ein)/tz0
      bz=(eou-ein)/tz
      a=max((ein+eou-2.0d0*u*sqrt(ein*eou))/(aws0*tz0),amin0)
      c=0.5d0*sqrt(eou/ein)/tz
      if (isl(1).lt.0) then
c
c       Principal atom has Sab tabulation
c
        s=stab(aws0,a,b,bz,lasym,liq,na,alpha,nb,beta,sab,arglim)
        if (s.lt.0.0d0) then
          s=0.0d0
          aa=a*aws0*tz0/tz
          do i=1,nsa
            if (isl(i).lt.1) then
              az=aa/aws(i)
              sn=sct(az,bz,temp,teff(i),arglim)
              s=sigb(i)*sn+s
            endif
          enddo
          s=c*s
        else
          s=sigb(1)*c*s
        endif
      else
c
c       No principal atom (analytic functions are used)
c
        s=0.0d0
        aa=a*aws0*tz0/tz
        do i=1,nsa
          az=aa/aws(i)
          if (isl(i).lt.1) then
            sn=sct(az,bz,temp,teff(i),arglim)
          else
            sn=sfree(az,bz,arglim)
          endif
          s=sigb(i)*sn+s
        enddo
        s=c*s
      endif
      if (s.lt.smin0) then
        sigtsl=0.0d0
      else
        sigtsl=s
      endif
      return
      end
C======================================================================
      real*8 function stab(aw0,a,b,bz,lasym,liq,na,alpha,nb,beta,sab,
     &  arglim)
      implicit real*8 (a-h, o-z)
      parameter (bmin=0.2d0, bmax=30.0d0)
      dimension alpha(*),beta(*),sab(na,*)
      bm=abs(b)
      if (a.gt.alpha(na).or.
     &   (lasym.eq.0.and.bm.gt.beta(nb)).or.
     &   (lasym.eq.1.and.(b.lt.beta(1).or.b.gt.beta(nb)))) then
c
c       Out of tabulated range (SCT will be applied)
c
        stab=-1.0d0
        return
      else
c
c       Liquid case (b**2)/a interpolation for small values of a and b
c
        if (liq.ne.0.and.a.lt.alpha(1).and.bm.le.bmin) then
          if (lasym.eq.1.and.b.lt.0.0d0) then
            ib1=liq-1
          else
            ib1=liq+1
          endif
          ar=alpha(1)/a
          br=b/beta(ib1)
          br=br*br
          sabx=sab(1,liq)+0.5d0*log(ar)-(sab(1,liq)-sab(1,ib1))*ar*br
          if (sabx.lt.arglim) sabx=arglim
        else
c
c         Interpolation on the Sab table will be tried
c
          bb=b
          if (lasym.eq.0.and.b.lt.0.0d0) bb=bm
          nbm=nb-1
          ib=nbm
          do i=1,nbm
            if (bb.lt.beta(i+1)) then
              ib=i
              exit
            endif
          enddo
          nam=na-1
          ia=nam
          do i=1,nam
            if (a.lt.alpha(i+1)) then
              ia=i
              exit
            endif
          enddo
c
c          Check for big transfer of momentum or energy
c
          if (a*aw0.ge.bmax.or.bm.ge.bmax) then
            ia1=ia+1
            ib1=ib+1
            if (sab(ia,ib).le.arglim.or.sab(ia1,ib).le.arglim.or.
     &          sab(ia,ib1).le.arglim.or.sab(ia1,ib1).le.arglim) then
c
c             SCT will be used
c
              stab=-2.0d0
              return
            endif
          endif
c
c         Interpolation on the Sab table will be used
c
          if (ia.eq.nam) ia=ia-1
          if (ib.eq.nbm) ib=ib-1
          ia1=ia+1
          ia2=ia+2
          ib1=ib+1
          ib2=ib+2
          sab0=sinter(a,alpha(ia),sab(ia,ib),alpha(ia1),sab(ia1,ib),
     &                alpha(ia2),sab(ia2,ib),arglim)
          sab1=sinter(a,alpha(ia),sab(ia,ib1),alpha(ia1),sab(ia1,ib1),
     &                alpha(ia2),sab(ia2,ib1),arglim)
          sab2=sinter(a,alpha(ia),sab(ia,ib2),alpha(ia1),sab(ia1,ib2),
     &                alpha(ia2),sab(ia2,ib2),arglim)
          sabx=sinter(bb,beta(ib),sab0,beta(ib1),sab1,
     &                beta(ib2),sab2,arglim)
        endif
        x=sabx-0.5d0*bz
        if (x.gt.arglim) then
          stab=exp(x)
        else
          stab=0.0d0
        endif
      endif
      return
      end
C======================================================================
      real*8 function sinter(x,x1,y1,x2,y2,x3,y3,ylim)
      implicit real*8 (a-h, o-z)
      parameter (delty=2.0d0)
      if (x.lt.x1) then
        if (y1.gt.y2) then
          y=y1
        else
          call terp1m(x1,y1,x2,y2,x,y,3)
        endif
      elseif (x.gt.x3) then
        if (y3.gt.y2) then
          y=y3
        else
          call terp1m(x2,y2,x3,y3,x,y,2)
        endif
      elseif (abs(y1-y2).gt.delty.or.abs(y2-y3).gt.delty) then
        if (x.lt.x2) then
          call terp1m(x1,y1,x2,y2,x,y,2)
        else
          call terp1m(x2,y2,x3,y3,x,y,2)
        endif
      else
        dx1=x-x1
        dx2=x-x2
        dx3=x-x3
        d12=x1-x2
        d13=x1-x3
        d23=x2-x3
        z1=dx2*dx3/(d12*d13)
        z2=-dx1*dx3/(d12*d23)
        z3=dx1*dx2/(d13*d23)
        y=z1*y1+z2*y2+z3*y3
      endif
      if (y.lt.ylim) then
        sinter=ylim
      else
        sinter=y
      endif
      return
      end
C======================================================================
      real*8 function sct(a,b,t,tef,arglim)
      implicit real*8 (a-h, o-z)
      parameter (c=3.54490770181103d0)
      bm=abs(b)
      tr=tef/t
      d=a-bm
      x=-(0.25d0*d*d/(a*tr)+0.5d0*(b+bm))
      if (x.gt.arglim) then
        sct=exp(x)/(c*sqrt(a*tr))
      else
        sct=0.0d0
      endif
      return
      end
C======================================================================
      real*8 function sfree(a,b,arglim)
      implicit real*8 (a-h, o-z)
      parameter (c=3.54490770181103d0)
      x=a+b
      x=-0.25d0*x*x/a
      if (x.gt.arglim) then
        sfree=exp(x)/(c*sqrt(a))
      else
        sfree=0.0d0
      endif
      return
      end
C======================================================================
      subroutine thrload(lst,je,es1,sigs1,nep1,nbin,xs1,ys1,us1,nmix,
     &  imon,itie,itix,itxe,itnep,ixss,nemax)
c
c     load thermal XSS data for incident energy Ein=es1
c
      implicit real*8(a-h, o-z)
      parameter (cdfstp=1.0d-6, ustp=0.00500d0, epstp=2.5d-7)
      parameter (ev2mev=1.0d-6)
      character*1 lin130(130)
      allocatable cdf(:),icdf(:)
      dimension xs1(*),ys1(*),us1(nbin,*)
      common/acexss/xss(100000000),nxss
      data lin130/130*'='/
      if (je.ge.nemax) then
        write(lst,*)
        write(lst,*)' Fatal error: Too many incident energy points'
        write(lst,*)' Current nemax for inelastic:',nemax
        stop
      endif
      allocate(cdf(nep1),icdf(nep1))
      if (sigs1.gt.0.0d0) then
c
c       correcting angle distribution at edge energies
c
        if (ys1(1).le.0.0d0) then
          do i=1,nbin
            us1(i,1)=us1(i,2)
          enddo
        endif
        if (ys1(nep1).le.0.0d0) then
          nepm=nep1-1
          do i=1,nbin
            us1(i,nep1)=us1(i,nepm)
          enddo
        endif
c
c       renormalizing pdf(Ein,Eou) and calculating cdf(Ein,Eou)
c
        x0=xs1(1)
        y0=ys1(1)/sigs1
        ys1(1)=y0
        cdf(1)=0.0d0
        do i=2,nep1
          x1=xs1(i)
          y1=ys1(i)/sigs1
          ys1(i)=y1
          cdf(i)=0.5d0*(x1-x0)*(y1+y0)+cdf(i-1)
          x0=x1
          y0=y1
        enddo
        cdf(nep1)=1.0d0
c
c       thinning cdf(Ein,Eou)
c
        i0=3
        do i=3,nep1
          if (cdf(i).gt.cdfstp) then
            i0=i
            exit
          endif
        enddo
        icdf(1)=i0-2
        icdf(2)=i0-1
        icdf(3)=i0
        ep0=xs1(i0)
        cdf0=cdf(i0)
        u0=avecos(us1(1,i0),nbin)
        i0=i0+1
        j=3
        ncdf1=nep1-1
        do i=i0,nep1
          ep1=xs1(i)
          cdf1=cdf(i)
          ys1i=ys1(i)
          im=i-1
          xs1m=xs1(im)
          ys1m=ys1(im)
          ip=i+1
          xs1p=xs1(ip)
          ys1p=ys1(ip)
          slopm=(ys1i-ys1m)/(ep1-xs1m)
          slopp=(ys1p-ys1i)/(xs1p-ep1)
          u1=avecos(us1(1,i),nbin)
          if ((((cdf1-cdf0).ge.cdfstp.or.(abs(u1-u0).ge.ustp.and.
     &       (ep1-ep0).gt.epstp*ep0)).and.i.lt.ncdf1).or.
     &       (slopm*slopp.le.0.0d0.and.(xs1p-xs1m).gt.epstp*xs1m).or.
     &       i.eq.nep1) then
            j=j+1
            icdf(j)=i
            ep0=ep1
            cdf0=cdf1
            u0=u1
          endif
        enddo
        mcdf=j
      else
        mcdf=3
        icdf(1)=0
        icdf(2)=2
        icdf(3)=3
        cdf(1)=0.0d0
        cdf(2)=0.0d0
        cdf(3)=1.0d0
      endif
      nnep=mcdf-1
      je=je+1
      sigs=sigs1/dble(nmix)
c
c     Loading the XSS array for energy Ein=es1
c
      write(*,'(a,i5,a,1p,e15.8,a,e15.8)')' ie=',je,
     &  ' incident energy=',es1,' inelastic scattering=',sigs
      write(lst,*)
      write(lst,'(a,i5,a,1p,e15.8,a,e15.8)')' ie=',je,
     &  ' incident energy=',es1*ev2mev,' inelastic scattering=',sigs
      write(lst,'(a,i5,a,i5,a)')
     &  ' cumulative probability distribution (cdf) given at ',nnep,
     &  ' points from ',nep1,' initial outgoing energies'
      if (imon.gt.0) then
        write(lst,'(5a)')' iep ',' outgoing energy',
     &    '      pdf      ','       cdf     ',
     &    '   equi-probable cosines'
        write(lst,'(130a1)')lin130
      endif
c
c      Inelastic itie block (incident energy/cross section)
c
      xss(itie)=je
      xss(itie+je)=es1*ev2mev
      xss(itix+je-1)=sigs
c
c      Inelastic itxe block (energy/angle distribution)
c
      if ((ixss+(3+nbin)*mcdf).gt.nxss) then
        write(*,*)' Fatal error in thrload: increase xss dimension'
        write(*,*)' xss array dimension: ',nxss
        write(*,*)' min.  required size: ',(ixss+3+nbin)*mcdf
        write(*,*)' je=',je,' e=',es1,' sigs=',sigs
        write(lst,*)' Fatal error in thrload: increase xss dimension'
        write(lst,*)' xss array dimension: ',nxss
        write(lst,*)' min.  required size: ',(ixss+3+nbin)*mcdf
        write(lst,*)' je=',je,' e=',es1,' sigs=',sigs        
        stop
      endif
      xss(itxe+je-1)=ixss
      xss(itnep+je-1)=nnep
      do i=2,mcdf
        j=icdf(i)
        epx=xs1(j)*ev2mev
        pdfx=ys1(j)/ev2mev
        cdfx=cdf(j)
        xss(ixss+1)=epx
        xss(ixss+2)=pdfx
        xss(ixss+3)=cdfx
        do l=1,nbin
          xss(ixss+3+l)=us1(l,j)
        enddo
        ixss=ixss+3+nbin
        if (imon.gt.0) then
          call printine(lst,i-1,j,epx,pdfx,cdfx,us1,nbin)
        endif
      enddo
      deallocate(cdf,icdf)
      return
      end
C======================================================================
      subroutine sigela(in2,lst,matsl,temp,nmix,nbin,tole,ei,nei,
     &  xnatom,imon)
      implicit real*8 (a-h, o-z)
      parameter (nemax=2000, nkks=20)
      parameter (tolbrg=5.0d-7, tolde=1.0d-6, tolwrt=1.0005d0)
      parameter (ev2mev=1.0d-6)
      character*1 line115(115)
      dimension ei(*)
      common/acepnt/nxs(16),jxs(32)
      common/acexss/xss(100000000),nxss
      dimension nbt(20),ibt(20)
      allocatable eb(:),s(:),w(:),t(:)
      allocatable ee(:),xse(:),uus(:),ue(:,:)
      allocatable uei(:),ue2(:),ees(:),xses(:),ues(:,:)
      data line115/115*'='/
c
c     searching for thermal elastic scattering data
c
      write(lst,*)
      call findmt(in2,matsl,7,2,icod)
      if (icod.ne.0) then
        write(*,'(a,i4)')
     &    ' No thermal elastic scattering data for MAT=',matsl
        write(lst,'(a,i4)')
     &    ' No thermal elastic scattering data for MAT=',matsl
        return
      else
        write(lst,'(a)')' Thermal elastic scattering'
      endif
c
c       reading thermal elastic scattering in ENDF-6 format
c
      call readcont(in2,za,awr,lthr,l2,n1,n2,mat,mf,mt,nsi)
      if (lthr.eq.1.or.lthr.eq.3) then
c
c       reading coherent elastic data (lthr=1 or 3)
c
        call readcont(in2,temp1,c2,lt,l2,nr,np,mat,mf,mt,nsi)
        backspace(in2)
        allocate (eb(np),s(np),t(np))
        call readtab1(in2,temp1,c2,lt,l2,nr,np,nbt,ibt,eb,t)
        call checklaw(nr,ibt,icod)
        if (icod.ne.1) then
          write(*,'(a,a)')
     &      ' === Warning: Bragg edges given with interpolation law',
     &      ' different of INT=1 (constant). INT=1 assumed'
          write(lst,*)
          write(lst,'(a,a)')
     &      ' === Warning: Bragg edges given with interpolation law',
     &      ' different of INT=1 (constant). INT=1 assumed'
        endif
        deltmax=dtemp(temp,2)
        itemp=0
        nt=lt+1
        do k=1,nt
          if (k.ne.1) then
            call readlist(in2,temp1,c2,li,l2,np1,n2,t)
            if (np1.ne.np) then
              write(lst,*)
              write(lst,'(a)')' === Error: list/tab1 NP not match'
              write(lst,'(11x,a,i6,a,i6)')' NP1=',np1,'  NP=',np
              write(*,'(a)')' === Error: list/tab1 NP not match'
              write(*,'(11x,a,i6,a,i6)')' NP1=',np1,'  NP=',np
              close(in2)
              close(lst)
              stop
            endif
          endif
          delt=abs(temp1-temp)
          if (delt.le.deltmax) then
            do j=1,np
              s(j)=t(j)
            enddo
            deltmax=delt
            itemp=k
            temp0=temp1
          endif
        enddo
        if (itemp.gt.0) then
          write(lst,*)
          if (abs(temp-temp0).gt.0.001d0) then
            write(lst,'(a,a,1pe12.5,a)')' === Warning:',
     &                ' Requested temperature = ',temp,' K not found'
            write(lst,'(a,a,1pe12.5,a)')'             ',
     &                ' Closest temperature found = ',temp0,
     &                ' K will be used'
          else
            write(lst,'(a,1pe12.5,a)')' Requested temperature = ',
     &                  temp,' K (found)'
          endif
        else
          write(lst,'(a,a,1pe12.5,a)')' === Error:',
     &                ' Requested temperature = ',temp,' K not found'
          write(*,'(a,a,1pe12.5,a)')' === Error:',
     &                ' Requested temperature = ',temp,' K not found'
          close(in2)
          close(lst)
          stop
        endif
        deallocate(t)
      endif
      if (lthr.eq.2.or.lthr.eq.3) then
c
c       reading incoherent elastic data (lthr=2 or 3)
c
        call readcont(in2,sb,c2,l1,l2,nr,nt,mat,mf,mt,nsi)
        backspace(in2)
        allocate (t(nt),w(nt))
        call readtab1(in2,sb,c2,l1,l2,nr,nt,nbt,ibt,t,w)
        delt=dtemp(temp,2)
        if ((t(1)-delt).gt.temp.or.(t(nt)+delt).lt.temp) then
          write(lst,*)
          write(lst,'(a,a,1pe12.5,a,2e12.5)')' === Error: Requested',
     &      ' temperature',temp,' is out range',t(1),t(nt)
          write(*,'(a,a,1pe12.5,a,2e12.5)')' === Error: Requested',
     &      ' temperature',temp,' is out range',t(1),t(nt)
          close(in2)
          close(lst)
          stop
        elseif (t(1).ge.temp) then
          wp=w(1)
        elseif (t(nt).le.temp) then
          wp=w(nt)
        else
          wp=fvalue(nr,nbt,ibt,nt,t,w,temp)
        endif
        deallocate(t,w)
      endif
      dnmix=dble(nmix)
c
c       Processing thermal elastic scattering data
c
      if (lthr.eq.1.or.lthr.eq.3) then
c
c       Processing thermal coherent elastic scattering
c
        write(*,'(a)')' Coherent elastic scattering'
        write(lst,*)
        write(lst,'(a)')' Thermal coherent elastic scattering'
        write(lst,*)
        allocate (t(np), w(np))
        write(lst,'(a,a)')'   ie   Bragg energy     S(E,T)    ',
     &      ' cross section'
        write(lst,'(50a1)')line115(1:50)
        j=0
        b0=0.0d0
        btol=tolbrg
        emax=ei(nei)
        i=1
        do while (i.le.np.and.eb(i).le.emax)
          b1=s(i)
          if ((b1-b0).gt.btol.or.i.eq.np) then
            j=j+1
            t(j)=eb(i)*ev2mev
            w(j)=b1*ev2mev/dnmix
            b0=b1
            btol=(b0+1.0d0)*tolbrg
          endif
          i=i+1
        enddo
        np=j
c
c       Set flags and triggers for coherent elastic scattering
c
        if (lthr.eq.1) then
          nxs(5)=4
          nxs(8)=0
          jxs(7)=0
          jxs(8)=0
          jxs(9)=0
        else
          nxs(5)=5
        endif
        nxs(6)=-1
        itce=nxs(1)+1
        itcx=itce+np+1
        jxs(4)=itce
        jxs(5)=itcx
        jxs(6)=0
c
c       Loading XSS array
c
        j0=itcx-1
        xss(itce)=np
        do j=1,np
          xss(itce+j)=t(j)
          xss(j0+j)=w(j)
          xscoh=w(j)/t(j)
          write(lst,'(i5,1p3e15.8)')j,t(j),w(j),xscoh
          write(*,'(a,i5,1p,2(a,e15.8))')' ie=',j,' incident energy=',
     &      t(j)/ev2mev,' coherent elastic xs=',xscoh
        enddo
        ixss=2*np+1
        nxs(1)=nxs(1)+ixss
        write(lst,*)
        write(lst,'(a,a,2i10)')' Length of coherent elastic data and',
     &    ' XSS array: ',ixss,nxs(1)
        deallocate(eb,s,t,w)
      endif
c
c       Processing thermal incoherent elastic scattering
c
      if (lthr.eq.2.or.lthr.eq.3) then
        write(*,'(a)')' Incoherent elastic scattering'
        write(lst,*)
        write(lst,'(a)')' Thermal incoherent elastic scattering'
        write(lst,*)
        write(lst,'(a,1pe14.7,a,e14.7,a,e14.7,a,i3)')
     &    ' Bound cross section SB=',sb,'  T=',temp,
     &    '  Debye-Waller parameter W''(T)=',wp,'  NATOM=',nint(xnatom)
        write(lst,*)
        if (imon.gt.0) then
          write(lst,'(a,a)')'   ie     energy     cross section ',
     &      '   equiprobable cosines'
          write(lst,'(115a1)')line115
        else
          write(lst,'(a)')'   ie     energy     cross section '
          write(lst,'(35a1)')line115(1:35)
        endif
c
c       computing cross section and angle distribution in a linearly
c       interpolable incident energy grid from the initial incident
c       energy grid
c
        allocate(ee(nemax),xse(nemax),ue(nbin,nemax))
        allocate(ees(nkks),xses(nkks),uus(nkks),ues(nbin,nkks))
        allocate(uei(nbin),ue2(nbin))
        sb2=0.5d0*sb/(xnatom*dnmix)
        w2=2.0d0*wp
        dbin=dble(nbin)
        j=0
        tolu=2.0d0*tole
        do i=1,nei
          e=ei(i)
          ke=1
          do while (ke.gt.0)
            c2=w2*e
            c2inv=1.0d0/c2
            xe=1.0d0-exp(-2.0d0*c2)
            xeinv=1.0d0/xe
            xsei=sb2*xe*c2inv
            u0=-1.0d0
            do k=1,nbin
              xu=exp(-c2*(1.0d0-u0))
              u1=1.0d0+c2inv*log(xe/dbin+xu)
              uei(k)=dbin*c2inv*(exp(-c2*(1.0d0-u1))*(c2*u1-1.0d0)-
     &          xu*(c2*u0-1))*xeinv
              u0=u1
            enddo
            call chkcos(uei,nbin,ineg,ipos,uneg,upos)
            if (ineg.gt.1.and.uneg.le.-tolwrt) then
              write(*,'(a,1pe14.7,a,i4,a,0pf7.4,a)')'   Warning: ein=',
     &          e,' nn=',ineg,' umin=',uneg,' corrected'
            endif
            if (ipos.gt.1.and.upos.ge.tolwrt) then
              write(*,'(a,1pe14.7,a,i4,a,0pf7.4,a)')'   Warning: ein=',
     &          e,' np=',ipos,' umax=',upos,' corrected'
            endif
            uui=avecos(uei,nbin)
            if (i.eq.1) then
              ee1=e
              xse1=xsei
              uu1=uui
              j=j+1
              ee(j)=e
              xse(j)=xsei
              do k=1,nbin
                ue(k,j)=uei(k)
              enddo
              ke=0
            elseif (ke.eq.1) then
              ee2=e
              xse2=xsei
              uu2=uui
              do k=1,nbin
                ue2(k)=uei(k)
              enddo
              e=0.5d0*(ee1+ee2)
              e=fix8dig(e,0)
              ke=2
              kk=0
            else
              if (e.le.ee1.or.e.ge.ee2) then
                dsige=0.0d0
                duu=0.0d0
                tests=1.0d0
                testu=1.0d0
              else
                xsem=0.5d0*(xse1+xse2)
                uum=0.5d0*(uu1+uu2)
                dsige=abs(xsem-xsei)
                duu=abs(uum-uui)
                tests=tole*abs(xsei)+1.0d-12
                testu=tolu*abs(uui)+1.0d-6
              endif
              if ((dsige.lt.tests.and.duu.lt.testu).or.
     &          ((ee2-ee1).le.tolde*ee2).or.(kk.ge.nkks)) then
                ee1=ee2
                xse1=xse2
                uu1=uu2
                j=j+1
                if (j.gt.nemax) then
                  write(lst,*)
                  write(lst,*)' Fatal error: Too many incident energies'
                  write(lst,*)' Current nemax for elastic:',nemax
                  stop
                else
                  ee(j)=ee2
                  xse(j)=xse2
                  do k=1,nbin
                    ue(k,j)=ue2(k)
                  enddo
                endif
                if (kk.eq.0) then
                  ke=0
                else
                  ee2=ees(kk)
                  xse2=xses(kk)
                  uu2=uus(kk)
                  do k=1,nbin
                    ue2(k)=ues(k,kk)
                  enddo
                  kk=kk-1
                  e=0.5d0*(ee1+ee2)
                  e=fix8dig(e,0)
                endif
              else
                kk=kk+1
                ees(kk)=ee2
                xses(kk)=xse2
                uus(kk)=uu2
                do k=1,nbin
                  ues(k,kk)=ue2(k)
                enddo
                ee2=e
                xse2=xsei
                uu2=uui
                do k=1,nbin
                  ue2(k)=uei(k)
                enddo
                e=0.5d0*(ee1+ee2)
                e=fix8dig(e,0)
              endif
            endif
          enddo
        enddo
        nee=j
c
c       Set flags and triggers for incoherent elastic scattering
c
        itce=nxs(1)+1
        itcx=itce+nee+1
        itca=itcx+nee
        if (lthr.eq.2) then
          nxs(5)=3
          nxs(6)=nbin-1
          nxs(8)=0
          jxs(4)=itce
          jxs(5)=itcx
          jxs(6)=itca
          jxs(7)=0
          jxs(8)=0
          jxs(9)=0
        else
          nxs(8)=nbin-1
          jxs(7)=itce
          jxs(8)=itcx
          jxs(9)=itca
        endif
c
c       Prepare itce and itca block for incoherent elastic
c
        xss(itce)=nee
        do i=1,nee
          xss(itce+i)=ee(i)*ev2mev
          xss(itcx+i-1)=xse(i)
          ij0=itca+(i-1)*nbin-1
          do j=1,nbin
            xss(ij0+j)=ue(j,i)
          enddo
          if (imon.gt.0) then
            call printeli(lst,i,xss(itce+i),xss(itcx+i-1),
     &        xss(ij0+1),nbin)
          else
            write(lst,'(i5,1p2e15.8)')i,xss(itce+i),xss(itcx+i-1)
          endif
          write(*,'(a,i5,1p, 2(a,e15.8))')' ie=',i,' incident energy=',
     &      ee(i),' incoherent elastic xs=',xss(itcx+i-1)
        enddo
        ixss=nee*(nbin+2)+1
        nxs(1)=nxs(1)+ixss
        write(lst,*)
        write(lst,'(a,a,2i10)')' Length of incoherent elastic data and',
     &    ' XSS array: ',ixss,nxs(1)
        deallocate (ee,xse,ue,uei,ue2)
        deallocate (ees,xses,uus,ues)
      endif
      return
      end
C======================================================================
      subroutine chkcos(ubar,nmu,ineg,ipos,uneg,upos)
      implicit real*8 (a-h, o-z)
      parameter (rtolu=1.0d-10, test=1.0d0-rtolu)
      dimension ubar(*)
      call orderx(ubar,nmu,rtolu,-1)
      uneg=-1.0d0
      upos=1.0d0
      ineg=0
      ipos=0
      do imu=1,nmu
        u0=ubar(imu)
        if (u0.le.-test) then
          if (u0.lt.uneg) uneg=u0
          ineg=ineg+1
          ubar(imu)=-1.0d0
        elseif (u0.ge.test) then
          if (u0.gt.upos) upos=u0
          ipos=ipos+1
          ubar(imu)=1.0d0
        endif
      enddo
      if (ineg.gt.1) then
        i0=ineg+1
        u0=ubar(i0)
        du=(u0+1.0d0)/dble(ineg)
        do i=1,ineg-1
          ubar(1+i)=-1.0d0+du*i
        enddo
      endif
      if (ipos.gt.1) then
        i0=nmu-ipos
        u0=ubar(i0)
        du=(1.0d0-u0)/dble(ipos)
        do i=1,ipos-1
          ubar(i0+i)=u0+du*i
        enddo
      endif
      return
      end
C======================================================================
      real*8 function avecos(ubar,nbin)
      implicit real*8 (a-h, o-z)
      dimension ubar(*)
      sumu=0.0d0
      do l=1,nbin
        sumu=ubar(l)+sumu
      enddo
      avecos=sumu/dble(nbin)
      return
      end
C======================================================================
      real*8 function dtemp(temp,ktemp)
      implicit real*8 (a-h, o-z)
      if (ktemp.eq.0) then
        dtemp=0.001d0*temp
      elseif (ktemp.eq.1) then
        dtemp=0.05d0*temp
      else
        dtemp=0.001d0*temp+3.0d0
      endif
      return
      end
C======================================================================
c     General routines
C======================================================================
      subroutine checklaw(nr,ibt,icod)
c
c      Check interpolation law
c       icod = 0 if interpolation law is linear everywhere
c       icod = 1 if interpolation law is constant everywhere
c       icod =-1 if one interpolation law, but not linear nor constant
c       icod =-2 otherwise
c
      dimension ibt(*)
      lawmin=9
      lawmax=0
      do i=1,nr
        ilaw=ibt(i)
        law2=ilaw-(ilaw/10)*10
        lawmin=min(lawmin,law2)
        lawmax=max(lawmax,law2)
      enddo
      if (lawmin.eq.lawmax) then
        if (lawmax.eq.2) then
          icod=0
        elseif (lawmax.eq.1) then
          icod=1
        else
          icod=-1
        endif
      else
        icod=-2
      endif
      return
      end
C======================================================================
      subroutine orderx(x,n,tol,irem)
c
c     the n elements of array x are ordered in ascending order
c     if irem is greater than 0, then elements within relative
c     tolerance tol are removed
c
      implicit real*8 (a-h, o-z)
      dimension x(*)
      if (n.gt.1) then
        m=n
        i=0
        do while (i.lt.m-1)
          i=i+1
          j=i
          do while (j.lt.m)
            j=j+1
            if (x(j).lt.x(i)) then
              temp=x(j)
              x(j)=x(i)
              x(i)=temp
            endif
          enddo
          if (i.gt.1) then
            if (irem.gt.0.and.abs(x(i)-x(i-1)).le.abs(tol*x(i))) then
              m=m-1
              if (i.ge.m) then
                 n=m
                 return
              endif
               do k=i,m
                 x(k)=x(k+1)
               enddo
               i=i-1
            endif
          endif
        enddo
        if (irem.gt.0.and.abs(x(m)-x(m-1)).le.abs(tol*x(m))) m=m-1
        n=m
      endif
      return
      end
C======================================================================
      subroutine inc(i,nmax)
c
c     Increase i in one unit and check range
c
      i=i+1
      if (i.gt.nmax) then
        write(*,*)' Error: dimension = ',nmax,' is not enough'
        write(*,*)'        increase array size'
        stop
      endif
      return
      end
C======================================================================
      real*8 function fix8dig(x,ndig)
      implicit real*8 (a-h, o-z)
      character*14 str14
      character*10 str10
      write(str14,'(1pe14.7)')x
      if (ndig.ne.0) then
        str10=str14(1:10)
        read(str14(12:14),'(i3)')nn
        read(str10,'(f10.0)')xmant
        xx=(xmant+dble(ndig)*0.0000001d0)*(10.0d0)**nn
        write(str14,'(1pe14.7)')xx
      endif
      read(str14,'(e14.0)')fix8dig
      return
      end
C======================================================================
C      Printing routines
C======================================================================
      subroutine printine(lst,i,j,ep,pdf,cdf,ubar,nbin)
      implicit real*8 (a-h, o-z)
      dimension ubar(nbin,*)
      nn=min(nbin,8)
      write(lst,'(i5,1p,3e15.8,0p,8f10.6)')
     &  i,ep,pdf,cdf,(ubar(l,j),l=1,nn)
      if (nbin.gt.8) then
        nn=nbin/8
        mm=nbin-nn*8
        nn=nn-1
        if (nn.gt.0) then
          do k=1,nn
            l0=8*k
            write(lst,'(50x,8f10.6)')(ubar(l0+l,j),l=1,8)
          enddo
        endif
        if (mm.gt.0) then
          l0=8*(nn+1)
          write(lst,'(50x,7f10.6)')(ubar(l0+l,j),l=1,mm)
        endif
      endif
      return
      end
C======================================================================
      subroutine printeli(lst,i,e,sig,ubar,nbin)
      implicit real*8 (a-h, o-z)
      dimension ubar(*)
      nn=min(nbin,8)
      write(lst,'(i5,1p2e15.8,0p,8f10.6)')i,e,sig,(ubar(l),l=1,nn)
      if (nbin.gt.8) then
        nn=nbin/8
        mm=nbin-nn*8
        nn=nn-1
        if (nn.gt.0) then
          do k=1,nn
            l0=8*k
            write(lst,'(35x,8f10.6)')(ubar(l0+l),l=1,8)
          enddo
        endif
        if (mm.gt.0) then
          l0=8*(nn+1)
          write(lst,'(35x,7f10.6)')(ubar(l0+l),l=1,mm)
        endif
      endif
      return
      end
C======================================================================
C      Interpolation routines
C======================================================================
      real*8 function cohxs(n,es,s,e)
c
c      Thermal coherent elastic scattering interpolation/calculation
c
      implicit real*8 (a-h,o-z)
      dimension es(*),s(*)
      if (e.lt.es(1)) then
        cohxs=0.0d0
      elseif (e.ge.es(n)) then
        cohxs=s(n)/e
      else
        n1=n-1
        do i=1,n1
          if (e.ge.es(i).and.e.lt.es(i+1)) then
            ie=i
            exit
          endif
        enddo
        cohxs=s(ie)/e
      endif
      return
      end
C=====================================================================
      real*8 function xsvalue(np,x1,y1,x)
c
c       Linear interpolation
c
        implicit real*8 (a-h, o-z)
        dimension x1(*),y1(*)
        if(x.lt.x1(1).or.x.gt.x1(np)) then
          xsvalue=0.0d0
        else
          nn=np-1
          do i=1,nn
            xa=x1(i)
            xb=x1(i+1)
            if (x.ge.xa.and.x.le.xb) then
              i0=i
              exit
            endif
          enddo
          if (x.eq.xa) then
            xsvalue=y1(i0)
          elseif (x.eq.xb) then
            xsvalue=y1(i0+1)
          else
            ii=i0+1
            xsvalue=y1(i0)+(x-xa)/(xb-xa)*(y1(ii)-y1(i0))
          endif
        endif
        return
      end
C=====================================================================
      real*8 function fvalue(nr,nbt,ibt,np,x,y,x0)
c
c      Return f(x0): function value at x0.
c      Function f is given by an ENDF-6/TAB1 record
c
        implicit real*8 (a-h, o-z)
        dimension nbt(*),ibt(*),x(*),y(*)
        character*11 cx0,cx1
        if (x0.lt.x(1).or.x0.gt.x(np)) then
          fvalue=0.0d0
        else
          call chendf(x0,cx0)
          i=2
          do while (i.le.np.and.x(i).lt.x0)
            i=i+1
          enddo
          i1=i-1
          call chendf(x(i1),cx1)
          if (cx0.eq.cx1) then
            fvalue=y(i1)
          elseif (x0.eq.x(i)) then
            fvalue=y(i)
          else
            j=1
            do while (nbt(j).lt.i)
              j=j+1
            enddo
            ilaw=ibt(j)
            call terp1m(x(i1),y(i1),x(i),y(i),x0,y0,ilaw)
            fvalue=y0
          endif
        endif
        return
      end
C======================================================================
C      General routines for reading ENDF-6 formatted files
C======================================================================
      subroutine findmat(nin,mat,icod)
c
c      Find material mat on endf6 formatted tape
c      on return if icod=0, material found
c                if icod=1, material not found
c
      character*66 line
      read(nin,'(a66,i4,i2,i3,i5)',iostat=iosnin)line,mat0,mf,mt,ns
      if (iosnin.lt.0.or.mat0.eq.-1.or.mat0.ge.mat) then
        rewind(nin)
        read(nin,*)
        mat0=0
      elseif (iosnin.gt.0) then
        icod=1
        return
      elseif (mat0.eq.0) then
        backspace nin
        backspace nin
        read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf,mt,ns
        if (mat0.ge.mat) then
          rewind(nin)
          read(nin,*)
          mat0=0
        endif
      endif
      do while (mat0.lt.mat.and.mat0.ne.-1)
        read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf,mt,ns
      enddo
      if (mat0.eq.mat)then
        icod=0
        backspace nin
      else
        icod=1
      endif
      return
   10 icod=1
      return
      end
C======================================================================
      subroutine findmf(nin,mat,mf,icod)
c
c      Find file mf for material mat on endf6 formatted tape
c      on return if icod=0, mat/mf found
c                if icod=1, mat/mf not found
c
      character*66 line
      read(nin,'(a66,i4,i2,i3,i5)',iostat=iosnin)line,mat0,mf0,mt,ns
      if (mat0.eq.0) then
        read(nin,'(a66,i4,i2,i3,i5)',iostat=iosnin)line,mat0,mf0,mt,ns
      endif
      if (iosnin.ne.0.or.mat0.eq.-1.or.mat0.gt.mat.or.
     &  (mat0.eq.mat.and.mf0.gt.mf)) then
        call findmat(nin,mat,icod)
      elseif (mat0.eq.mat) then
        if (mf0.eq.0.or.mf0.eq.mf) then
          backspace nin
          backspace nin
          read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf0,mt,ns
          if (mf0.ge.mf) then
            call findmat(nin,mat,icod)
          else
            icod=0
          endif
        else
          icod=0
        endif
      elseif (mat0.lt.mat) then
        do while (mat0.lt.mat.and.mat0.ne.-1)
          read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf0,mt,ns
        enddo
        if (mat0.eq.mat)then
          icod=0
          backspace nin
        else
          icod=1
        endif
      endif
      if (icod.eq.0) then
        mat0=mat
        mf0=-1
        do while (mat0.eq.mat.and.mf0.lt.mf)
          read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf0,mt,ns
        enddo
        if (mat0.eq.mat.and.mf0.eq.mf) then
          icod=0
          backspace nin
        else
          icod=1
        endif
      endif
      return
   10 icod=1
      return
      end
C======================================================================
      subroutine findmt(nin,mat,mf,mt,icod)
c
c      Find reaction mt on file mf for material mat on endf6 formatted
c      tape, on return if icod=0, mat/mf/mt found
c                      if icod=1, mat/mf/mt not found
c
      character*66 line
      read(nin,'(a66,i4,i2,i3,i5)',iostat=iosnin)line,mat0,mf0,mt0,ns
      if (iosnin.ne.0.or.mat0.ne.mat.or.
     &  (mat0.eq.mat.and.mf0.ne.mf).or.
     &  (mat0.eq.mat.and.mf0.eq.mf.and.mt0.gt.mt)) then
        call findmf(nin,mat,mf,icod)
      elseif (mat0.eq.mat.and.mf0.eq.mf) then
        if (mt0.eq.0.or.mt0.eq.mt) then
         backspace nin
         backspace nin
         read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf0,mt0,ns
         if (mt0.ge.mt) then
           call findmf(nin,mat,mf,icod)
         else
           icod=0
         endif
        else
         icod=0
        endif
      endif
      if (icod.eq.0) then
        mf0=mf
        mt0=-1
        do while (mt0.lt.mt.and.mf0.eq.mf)
         read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf0,mt0,ns
        enddo
        if (mf0.eq.mf.and.mt0.eq.mt) then
          icod=0
          backspace nin
        else
          icod=1
        endif
      endif
      return
   10 icod=1
      return
      end
C======================================================================
      subroutine readtext(nin,line,mat,mf,mt,ns)
c
c     read a TEXT record
c
      character*66 line
      read(nin,10)line,mat,mf,mt,ns
      return
   10 format(a66,i4,i2,i3,i5)
      end
C======================================================================
      subroutine readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
c
c     read a CONT record
c
      implicit real*8 (a-h, o-z)
      read(nin,10)c1,c2,l1,l2,n1,n2,mat,mf,mt,ns
      return
   10 format(2e11.0,4i11,i4,i2,i3,i5)
      end
C======================================================================
      subroutine readlist(nin,c1,c2,l1,l2,npl,n2,b)
c
c     read a LIST record
c
      implicit real*8 (a-h, o-z)
      dimension b(*)
      read(nin,10)c1,c2,l1,l2,npl,n2,mat,mf,mt,ns
      read(nin,20)(b(n),n=1,npl)
      return
   10 format(2e11.0,4i11,i4,i2,i3,i5)
   20 format(6e11.0)
      end
C======================================================================
      subroutine readtab1(nin,c1,c2,l1,l2,nr,np,nbt,intp,x,y)
c
c     read TAB1 record
c
      implicit real*8 (a-h, o-z)
      dimension nbt(*),intp(*),x(*), y(*)
      read(nin,10)c1,c2,l1,l2,nr,np,mat,mf,mt,ns
      read(nin,20)(nbt(n),intp(n),n=1,nr)
      read(nin,30)(x(n),y(n),n=1,np)
      return
   10 format(2e11.0,4i11,i4,i2,i3,i5)
   20 format(6i11)
   30 format(6e11.0)
      end
C======================================================================
      subroutine readtab2(nin,c1,z,l1,l2,nr,nz,nbt,intp)
c
c     read TAB2 record
c
      implicit real*8 (a-h, o-z)
      dimension nbt(*),intp(*)
      read(nin,10)c1,z,l1,l2,nr,nz,mat,mf,mt,ns
      read(nin,20)(nbt(n),intp(n),n=1,nr)
      return
   10 format(2e11.0,4i11,i4,i2,i3,i5)
   20 format(6i11)
      end
C======================================================================
C     Time routine
C======================================================================
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
C======================================================================
C     Routine to write *.plt & *.cur files for PLOTTAB
C======================================================================
      subroutine thrplot
c
c      Write *.plt and *.cur files for plotting with PLOTTAB
c
      implicit real*8 (a-h, o-z)
      character*10 hd,hm
      character*11 x,y
      character*13 hz
      character*70 hk
      character*72 fout,fdir
      parameter (r0=1.02d0, ncmax=15, ncp=10)
      common/acetxt/hz,hd,hm,hk
      common/acecte/awrth,tmev,awm(16),izam(16)
      common/acepnt/nxs(16),jxs(32)
      common/acexss/xss(100000000),nxss
      allocatable uave(:), eu1(:), eu2(:)
      data nplt/40/,ncur/41/
c
c     Open main plotting units with fix filenames
c       DOTSL.PLT: Plottab plotting options
c       DOTSL.CUR: Plottab curve file
c
      open (nplt,file='DOTSL.PLT')
      open (ncur,file='DOTSL.CUR')
c
c       flags and triggers
c
      len2=nxs(1)
      idpni=nxs(2)
      nil=nxs(3)
      nieb=nxs(4)
      idpnc=nxs(5)
      ncl=nxs(6)
      ifeng=nxs(7)
      nclj=nxs(8)
      itie=jxs(1)
      itix=jxs(2)
      itxe=jxs(3)
      itce=jxs(4)
      itcx=jxs(5)
      itca=jxs(6)
      jtce=jxs(7)
      jtcx=jxs(8)
      jtca=jxs(9)
      nei=0
      ne=0
      nec=0
c
c     Cross section plotting
c
c       Inelastic
c
      nei=nint(xss(itie))
      emax=xss(itie+nei)
      write(ncur,'(a)')'Inelastic'
      do i=1,nei
        call chendf(xss(itie+i),x)
        call chendf(xss(itix+i-1),y)
        write(ncur,'(2a11)')x,y
      enddo
      write(ncur,*)
      nc=1
c
c       Coherent elastic
c
      if ((idpnc.eq.4.or.idpnc.eq.5).and.itce.ne.0) then
        write(ncur,'(a)')'Coherent elastic'
        nec=nint(xss(itce))
        e=xss(itce+1)
        xs=0.0d0
        call chendf(e,x)
        call chendf(xs,y)
        write(ncur,'(2a11)')x,y
        j=2
        do while (e.le.emax)
          xs=cohxs(nec,xss(itce+1),xss(itcx),e)
          call chendf(e,x)
          call chendf(xs,y)
          write(ncur,'(2a11)')x,y
          e=e*r0
          if (e.ge.xss(itce+j).and.j.le.nec) then
            e=xss(itce+j)
            e1=e*0.999999999d0
            xs=cohxs(nec,xss(itce+1),xss(itcx),e1)
            call chendf(e1,x)
            call chendf(xs,y)
            write(ncur,'(2a11)')x,y
            j=j+1
          endif
        enddo
        if (e/r0.lt.emax) then
          e=emax
          xs=cohxs(nec,xss(itce+1),xss(itcx),e)
          call chendf(e,x)
          call chendf(xs,y)
          write(ncur,'(2a11)')x,y
        endif
        write(ncur,*)
        nc=nc+1
      endif
c
c       Incoherent elastic
c
      if ((idpnc.eq.3.and.itce.ne.0).or.(idpnc.eq.5.and.jtce.ne.0)) then
        if (idpnc.eq.3) then
          ktce=itce
          ktcx=itcx
        else
          ktce=jtce
          ktcx=jtcx
        endif
        write(ncur,'(a)')'Incoherent elastic'
        ne=nint(xss(ktce))
        do i=1,ne
          call chendf(xss(ktce+i),x)
          call chendf(xss(ktcx+i-1),y)
          write(ncur,'(2a11)')x,y
        enddo
        write(ncur,*)
        nc=nc+1
      endif
c
c      Total
c
      if (nc.gt.1) then
        write(ncur,'(a)')'Total'
        nu1=nei+ne
        numax=nu1
        allocate(eu1(numax))
        if (ne.gt.0) then
          call union(nei,xss(itie+1),ne,xss(ktce+1),nu1,eu1,numax)
        else
          do i=1,nei
            eu1(i)=xss(itie+i)
          enddo
        endif
        nu2=nu1+nec
        numax=nu2
        allocate(eu2(numax))
        if (nec.gt.0) then
          call union(nu1,eu1,nec,xss(itce+1),nu2,eu2,numax)
        else
          do i=1,nu1
            eu2(i)=eu1(i)
          enddo
        endif
        j=1
        do i=1,nu2
          e=eu2(i)
          if (nec.gt.0.and.e.ge.xss(itce+j).and.j.le.nec) then
            e1=xss(itce+j)
            e0=e1*0.999999999d0
            xs=xsvalue(nei,xss(itie+1),xss(itix),e0)
            if (ne.gt.0) xs=xsvalue(ne,xss(ktce+1),xss(ktcx),e0)+xs
            if (nec.gt.0)xs=cohxs(nec,xss(itce+1),xss(itcx),e0)+xs
            call chendf(e0,x)
            call chendf(xs,y)
            write(ncur,'(2a11)')x,y
            do while (e1.lt.e)
              xs=xsvalue(nei,xss(itie+1),xss(itix),e1)
              if (ne.gt.0) xs=xsvalue(ne,xss(ktce+1),xss(ktcx),e1)+xs
              if (nec.gt.0)xs=cohxs(nec,xss(itce+1),xss(itcx),e1)+xs
              call chendf(e1,x)
              call chendf(xs,y)
              write(ncur,'(2a11)')x,y
              e1=e1*r0
            enddo
            j=j+1
          endif
          xs=xsvalue(nei,xss(itie+1),xss(itix),e)
          if (ne.gt.0) xs=xsvalue(ne,xss(ktce+1),xss(ktcx),e)+xs
          if (nec.gt.0)xs=cohxs(nec,xss(itce+1),xss(itcx),e)+xs
          call chendf(e,x)
          call chendf(xs,y)
          write(ncur,'(2a11)')x,y
        enddo
        write(ncur,*)
        nc=nc+1
        deallocate (eu1,eu2)
      endif
      write(nplt,'(1p,4e11.4,2i11,0p,f4.2)')0.0,13.5,0.0,10.0,1,1,1.5
      write(nplt,'(6i11,i4)')nc,0,1,0,0,0,0
      write(nplt,'(a15,25x,a3)')'Incident Energy','MeV'
      write(nplt,'(a13,27x,a4)')'Cross section','barn'
      write(nplt,'(a70)')hk
      write(nplt,'(a)')'Thermal scattering cross section'
      write(nplt,'(22x,4i11)')0,2,0,0
      write(nplt,'(22x,4i11)')0,2,0,0
c
c     Inelastic pdf energy distribution
c     (ncmax incident energies are picked up)
c
      slope=max(dble(nei-1)/dble(ncmax-1),1.0d0)
      itnep=itxe+nei
      nbin=nil-1
      kxss=nbin+3
      ic=0
      kc=1
      do ie=1,nei
        if (ie.eq.kc) then
          call chendf(xss(itie+ie),x)
          write(ncur,'(a4,a11)')'Ein=',x
          ixss=nint(xss(itxe+ie-1)+1.0d-5)
          nep=nint(xss(itnep+ie-1)+1.0d-5)
          do iep=1,nep
            call chendf(xss(ixss+1),x)
            call chendf(xss(ixss+2),y)
            write(ncur,'(2a11)')x,y
            ixss=ixss+kxss
          enddo
          write(ncur,*)
          ic=ic+1
          kc=nint(slope*dble(ic))+1
        endif
      enddo
      write(nplt,*)
      write(nplt,'(1p,4e11.4,2i11,0p,f4.2)')0.0,13.5,0.0,10.0,1,1,1.5
      write(nplt,'(6i11,i4)')-ic,0,1,0,0,0,0
      write(nplt,'(a16,24x,a3)')'Secondary Energy','MeV'
      write(nplt,'(a3,37x,a9)')'PDF','Prob./MeV'
      write(nplt,'(a70)')hk
      write(nplt,'(a)')'Thermal inelastic secondary energy distribution'
      write(nplt,'(22x,4i11)')0,2,0,0
      write(nplt,'(1pe11.4,11x,4i11)')1.0d-8,0,2,0,0
      write(nplt,'(a70)')hk
      write(nplt,'(a)')'Thermal inelastic secondary energy distribution'
      write(nplt,'(1pe11.4,11x,4i11)')1.0d-7,0,2,0,0
      write(nplt,'(1pe11.4,11x,4i11)')1.0d-8,0,2,0,0
c
c     Inelastic equally cosines as a function of (Ein,Eout)
c     ncmax incident energies and ncp outgoing energies are picked up
c     and the cumulative discrete distribution is plotted
c
      dcdf=1.0d0/dble(nbin)
      ic=0
      kc=1
      do ie=1,nei
        if (ie.eq.kc) then
          ei=xss(itie+ie)
          ixss=nint(xss(itxe+ie-1)+1.0d-5)
          nep=nint(xss(itnep+ie-1)+1.0d-5)
          ep1=xss(ixss+1)
          ep2=xss(ixss+(nep-2)*kxss+1)
          rp=exp(log(ep2/ep1)/dble(ncp-1))
          do iep=1,nep-1
            ep=xss(ixss+1)
            if (ep.ge.ep1) then
              call chendf(ep,x)
              write(ncur,'(a,a11)')'Eou=',x
              cdf0=0.0d0
              call chendf(cdf0,x)
              do l=1,nbin
                u1=xss(ixss+3+l)
                call chendf(u1,y)
                write(ncur,'(2a11)')x,y
                cdf0=cdf0+dcdf
                call chendf(cdf0,x)
                write(ncur,'(2a11)')x,y
              enddo
              write(ncur,*)
              ep1=min(ep1*rp,ep2)
            endif
            ixss=ixss+kxss
          enddo
          write(nplt,*)
          write(nplt,'(1p,4e11.4,2i11,0p,f4.2)')
     &      0.0,13.5,0.0,10.0,1,1,1.5
          write(nplt,'(6i11,i4)')ncp,0,1,0,0,0,0
          write(nplt,'(a)')'Unit interval for random sampling'
          write(nplt,'(a)')'Equally probable discrete cosines'
          write(nplt,'(a70)')hk
          call chendf(ei,x)
          write(nplt,'(a,a,a,i4,a,f9.6)')'Incident energy =',x,
     &      ' MeV  NBIN=',nbin,'  PROB.=',dcdf
          write(nplt,'(1p2e11.4,4i11)') 0.0d0,1.0d0,0,1,0,0
          write(nplt,'(1p2e11.4,4i11)')-1.0d0,1.0d0,0,1,0,0
          ic=ic+1
          kc=nint(slope*dble(ic))+1
        endif
      enddo
c
c     Inelastic average cosine as a function of output energy (Eout)
c     (ncmax incident energies are picked up) and calculation of
c     inelastic average cosine as a function of incident energy (Ein)
c
      allocate (uave(nei))
      ic=0
      kc=1
      do ie=1,nei
        ixss=nint(xss(itxe+ie-1)+1.0d-5)
        nep=nint(xss(itnep+ie-1)+1.0d-5)
        if (ie.eq.kc) then
          call chendf(xss(itie+ie),x)
          write(ncur,'(a4,a11)')'Ein=',x
        endif
        uu=0.0d0
        cdf0=0.0d0
        u0=avecos(xss(ixss+4),nbin)
        do iep=1,nep
          cdf1=xss(ixss+3)
          u1=avecos(xss(ixss+4),nbin)
          uu=0.5d0*(u0+u1)*(cdf1-cdf0)+uu
          if (ie.eq.kc) then
            call chendf(xss(ixss+1),x)
            call chendf(u1,y)
            write(ncur,'(2a11)')x,y
          endif
          cdf0=cdf1
          u0=u1
          ixss=ixss+kxss
        enddo
        uave(ie)=uu
        if (ie.eq.kc) then
          write(ncur,*)
          ic=ic+1
          kc=nint(slope*dble(ic))+1
        endif
      enddo
      write(nplt,*)
      write(nplt,'(1p,4e11.4,2i11,0p,f4.2)')0.0,13.5,0.0,10.0,1,1,1.5
      write(nplt,'(6i11,i4)')ic,0,1,0,0,0,0
      write(nplt,'(a16,24x,a3)')'Secondary Energy','MeV'
      write(nplt,'(a14)')'Average cosine'
      write(nplt,'(a70)')hk
      write(nplt,'(a,a)')'Thermal inelastic average cosine distribution'
      write(nplt,'(22x,4i11)')0,2,0,0
      write(nplt,'(22x,4i11)')0,1,0,0
c
c      Thermal inelastic average cosine plotting
c
      nc=1
      write(ncur,'(a)')'Inelastic'
      do ie=1,nei
        call chendf(xss(itie+ie),x)
        call chendf(uave(ie),y)
        write(ncur,'(2a11)')x,y
      enddo
      write(ncur,*)
      deallocate(uave)
c
c     Calculation of the incoherent elastic average cosine as a
c     function of incident energy (Ein)
c
      if ((idpnc.eq.3.and.itca.ne.0).or.(idpnc.eq.5.and.jtca.ne.0)) then
        if (idpnc.eq.3) then
          nbin=ncl+1
          ktce=itce
          ktca=itca
        else
          nbin=nclj+1
          ktce=jtce
          ktca=jtca
        endif
        write(ncur,'(a)')'Incoherent elastic'
        j0=ktca-1
        do i=1,ne
          ij0=j0+(i-1)*nbin+1
          uu=avecos(xss(ij0),nbin)
          call chendf(xss(ktce+i),x)
          call chendf(uu,y)
          write(ncur,'(2a11)')x,y
        enddo
        write(ncur,*)
        nc=nc+1
      endif
c
c     Calculation of coherent elastic average cosine as a
c     function of incident energy (Ein)
c
      if ((idpnc.eq.4.or.idpnc.eq.5).and.itce.ne.0) then
        write(ncur,'(a)')'Coherent elastic'
        e=xss(itce+1)
        j=2
        do while (e.le.emax)
          uu=0.0d0
          sums=0.0d0
          s0=0.0d0
          do i=1,nec
            ei=xss(itce+i)
            if (ei.gt.e) then
              exit
            else
              uui=1.0d0-2.0d0*ei/e
              s1=xss(itcx+i-1)
              ds=(s1-s0)/e
              uu=ds*uui+uu
              sums=ds+sums
            endif
            s0=s1
          enddo
          if (sums.ne.0.0d0) then
            uu=uu/sums
          else
            uu=1.0d0
          endif
          call chendf(e,x)
          call chendf(uu,y)
          write(ncur,'(2a11)')x,y
          e=e*r0
          if (e.ge.xss(itce+j).and.j.le.nec) then
            e=xss(itce+j)
            j=j+1
          endif
        enddo
        write(ncur,*)
        nc=nc+1
      endif
      write(nplt,*)
      write(nplt,'(1p,4e11.4,2i11,0p,f4.2)')0.0,13.5,0.0,10.0,1,1,1.5
      write(nplt,'(6i11,i4)')nc,0,1,0,0,0,0
      write(nplt,'(a15,25x,a3)')'Incident Energy','MeV'
      write(nplt,'(a14)')'Average cosine'
      write(nplt,'(a70)')hk
      write(nplt,'(a)')'Thermal average cosine'
      write(nplt,'(22x,4i11)')0,2,0,0
      write(nplt,'(22x,4i11)')0,1,0,0
      close(nplt)
      close(ncur)
      return
      end
C======================================================================
      subroutine chendf(ffin,str11)
      implicit real*8 (a-h,o-z)
c
c     pack value into 11-character string
c
      character*11 str11
      character*12 str12
      character*13 str13
      character*15 str15
      write(str15,'(1pe15.8)')ffin
      read(str15,'(e15.0)')ff
      aff=dabs(ff)
      if (aff.lt.1.00000d-99) then
        str11=' 0.0       '
      elseif (aff.lt.9.99999d+99) then
        write(str15,'(1pe15.8)')ff
        read(str15,'(12x,i3)')iex
        select case (iex)
          case (-9,-8,-7,-6,-5,-4, 9)
            write(str13,'(1pe13.6)')ff
            str11(1:9)=str13(1:9)
            str11(10:10)=str13(11:11)
            str11(11:11)=str13(13:13)
          case (-3,-2,-1)
            write(str12,'(f12.9)')ff
            str11(1:1)=str12(1:1)
            str11(2:11)=str12(3:12)
          case (0)
            write(str11,'(f11.8)')ff
          case (1)
            write(str11,'(f11.7)')ff
          case (2)
            write(str11,'(f11.6)')ff
          case (3)
            write(str11,'(f11.5)')ff
          case (4)
            write(str11,'(f11.4)')ff
          case (5)
            write(str11,'(f11.3)')ff
          case (6)
            write(str11,'(f11.2)')ff
          case (7)
            write(str11,'(f11.1)')ff
          case (8)
            write(str11,'(f11.0)')ff
          case default
            write(str12,'(1pe12.5)')ff
            str11(1:8)=str12(1:8)
            str11(9:11)=str12(10:12)
        end select
      else
        if (ff.lt.0.0000d0) then
          str11='-9.99999+99'
        else
          str11=' 9.99999+99'
        endif
      endif
      return
      end
C======================================================================
      subroutine union(np1,x1,np2,x2,np3,x3,npmax)
      implicit real*8 (a-h, o-z)
      parameter(eps=1.0d-9)
c
c     Prepare union grid x3 from x1 and x2
c
      dimension x1(*),x2(*),x3(*)
      i=1
      j=1
      k=0
      do while (i.le.np1.and.j.le.np2)
        call inc(k,npmax)
        xx1=x1(i)
        xx2=x2(j)
        if (abs(xx2-xx1).lt.eps*xx2) then
          x3(k)=xx1
          i=i+1
          j=j+1
        elseif (xx1.lt.xx2) then
          x3(k)=xx1
          i=i+1
        else
          x3(k)=xx2
          j=j+1
        endif
      enddo
      if (i.gt.np1) then
        do i=j,np2
          call inc(k,npmax)
          x3(k)=x2(i)
        enddo
      elseif (j.gt.np2) then
        do j=i,np1
          call inc(k,npmax)
          x3(k)=x1(j)
        enddo
      endif
      np3=k
      return
      end
C======================================================================
C      The following subroutines were taken from NJOY2016 and
C      adapted/modified by D. Lopez Aldama for ACEMAKER:
C       1. subroutine terp1 (renamed as terp1m)
C       2. subroutine throut
C       3. subroutine typen
C
C======================================================================
C     Copyright (c) 2016, Los Alamos National Security, LLC
C     All rights reserved.
C
C     Copyright 2016. Los Alamos National Security, LLC. This software
C     was produced under U.S. Government contract DE-AC52-06NA25396 for
C     Los Alamos National Laboratory (LANL), which is operated by Los
C     Alamos National Security, LLC for the U.S. Department of Energy.
C     The U.S. Government has rights to use, reproduce, and distribute
C     this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL
C     SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C     ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
C     modified to produce derivative works, such modified software
C     should be clearly marked, so as not to confuse it with the
C     version available from LANL.
C
C     Additionally, redistribution and use in source and binary forms,
C     with or without modification, are permitted provided that the
C     following conditions are met:
C     1. Redistributions of source code must retain the above copyright
C        notice, this list of conditions and the following disclaimer.
C     2. Redistributions in binary form must reproduce the above
C        copyright notice, this list of conditions and the following
C        disclaimer in the documentation and/or other materials provided
C        with the distribution.
C     3. Neither the name of Los Alamos National Security, LLC,
C        Los Alamos National Laboratory, LANL, the U.S. Government,
C        nor the names of its contributors may be used to endorse or
C        promote products derived from this software without specific
C        prior written permission.
C
C     THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC
C     AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
C     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
C     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
C     DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC
C     OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C     SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
C     USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
C     AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
C     LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
C     IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
C     THE POSSIBILITY OF SUCH DAMAGE.
C
C======================================================================
      subroutine terp1m(x1,y1,x2,y2,x,y,i)
c
c      interpolate one point using ENDF-6 interpolation laws
c      (x1,y1) and (x2,y2) are the end points
c      (x,y) is the interpolated point
c      i is the interpolation law (1-5)
c
c      Adapted by D. Lopez Aldama for ACEMAKER
c
c     *****************************************************************
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0)
c
c     *** x1=x2
      if (x2.eq.x1) then
        y=y1
        return
      endif
c
c     ***y is constant
      if (i.eq.1.or.y2.eq.y1.or.x.eq.x1) then
         y=y1
c
c     ***y is linear in x
      else if (i.eq.2) then
         y=y1+(x-x1)*(y2-y1)/(x2-x1)
c
c     ***y is linear in ln(x)
      else if (i.eq.3) then
         y=y1+log(x/x1)*(y2-y1)/log(x2/x1)
c
c     ***ln(y) is linear in x
      else if (i.eq.4) then
         y=y1*exp((x-x1)*log(y2/y1)/(x2-x1))
c
c     ***ln(y) is linear in ln(x)
      else if (i.eq.5) then
         if (y1.eq.zero) then
            y=y1
         else
            y=y1*exp(log(x/x1)*log(y2/y1)/log(x2/x1))
         endif
c
c     ***coulomb penetrability law (charged particles only)
      else
        write(*,*) ' Interpolation law: ',i,' not allowed'
        stop
      endif
      return
      end
C======================================================================
      subroutine throut(fout,mcnpx)
c
c      Write a type 1 ACE-formatted file for TSL with ifeng=2
c
c      Adapted by D. Lopez Aldama for ACEMAKER
c
      implicit real*8 (a-h, o-z)
      character*10 hd,hm
      character*13 hz
      character*70 hk
      character*72 fout,fdir
      common/acetxt/hz,hd,hm,hk
      common/acecte/awrth,tmev,awm(16),izam(16)
      common/acepnt/nxs(16),jxs(32)
      common/acexss/xss(100000000),nxss
      data nout/30/,ndir/31/
c
c       open output ACE-formatted file
c
      open(nout,file=fout)
c
c       header block
c
      if (mcnpx.eq.1) then
        write(nout,'(a13,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')
     &    hz(1:13),awrth,tmev,hd,hk,hm
      else
        write(nout,'(a10,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')
     &    hz(1:10),awrth,tmev,hd,hk,hm
      endif
      write(nout,'(4(i7,f11.0))') (izam(i),awm(i),i=1,16)
c
c       flags and triggers
c
      len2=nxs(1)
      idpni=nxs(2)
      nil=nxs(3)
      nieb=nxs(4)
      idpnc=nxs(5)
      ncl=nxs(6)
      ifeng=nxs(7)
      nclj=nxs(8)
      itie=jxs(1)
      itix=jxs(2)
      itxe=jxs(3)
      itce=jxs(4)
      itcx=jxs(5)
      itca=jxs(6)
      jtce=jxs(7)
      jtcx=jxs(8)
      jtca=jxs(9)
      write(nout,'(8i9)')(nxs(i),i=1,16)
      write(nout,'(8i9)')(jxs(i),i=1,32)
c
c       itie block (inelastic)
c
      l=itie
      ne=nint(xss(l))
      call typen(l,nout,1)
      l=l+1
      n=2*ne
      do i=1,n
         call typen(l,nout,2)
         l=l+1
      enddo
c
c       itxe block (inelastic energy/angle distribution)
c
      l=itxe
      n=2*ne
      do i=1,ne
        n=n+nint(xss(l+ne+i-1))*(nil+2)
      enddo
      do i=1,n
        call typen(l,nout,2)
        l=l+1
      enddo
      if (idpnc.ne.0) then
c
c       itce block (elastic)
c
        if (itce.ne.0) then
          l=itce
          ne=nint(xss(l))
          nexe=ne
          call typen(l,nout,1)
          l=l+1
          n=2*ne
          do i=1,n
            call typen(l,nout,2)
            l=l+1
          enddo
        endif
c
c       itca block (elastic angular distribution)
c
        if (itca.ne.0.and.ncl.ne.-1) then
          n=nexe*(ncl+1)
          do i=1,n
            call typen(l,nout,2)
            l=l+1
          enddo
        endif
c
c       mixed TSL (incoherent elastic scattering)
c
        if (idpnc.eq.5) then
c
c         jtce=itce2 block
c
          if (jtce.ne.0) then
            l=jtce
            ne=nint(xss(l))
            nexe=ne
            call typen(l,nout,1)
            l=l+1
            n=2*ne
            do i=1,n
              call typen(l,nout,2)
              l=l+1
            enddo
          endif
c
c         jtca=itca2 block
c
          if (jtca.ne.0.and.nclj.ne.-1) then
            n=nexe*(nclj+1)
            do i=1,n
              call typen(l,nout,2)
              l=l+1
            enddo
          endif
        endif
      endif
c
c       finish
c
      call typen(l,nout,3)
      nern=0
      lrec=0
      close(nout)
c
c      write *.xsd file for xsdir
c
      i=index(fout,'.',.true.)
      if (i.le.0) then
        i=len(trim(fout))
      else
        i=i-1
      endif
      write(fdir,'(a,a)')trim(fout(1:i)),'.xsd'
      fdir=trim(fdir)
      open(ndir,file=fdir)
      if (mcnpx.eq.1) then
        write(ndir,'(a13,f12.6,1x,a,1x,a,i9,2i3,1pe11.4)')
     &    hz(1:13),awrth,trim(fout),'0 1 1 ',len2,lrec,nern,tmev
      else
        write(ndir,'(a10,f12.6,1x,a,1x,a,i9,2i3,1pe11.4)')
     &    hz(1:10),awrth,trim(fout),'0 1 1 ',len2,lrec,nern,tmev
      endif
      close(ndir)
      return
      end
C======================================================================
      subroutine typen(l,nout,iflag)
c
c     -----------------------------------------------------------------
c      Write an integer or a real number to a Type-1 ACE file,
c      or (if nout=0) convert real to integer for type-3 output,
c      or (if nout=1) convert integer to real for type-3 input.
c      Use iflag.eq.1 to write an integer (i20).
c      Use iflag.eq.2 to write a real number (1pe20.11).
c      Use iflag.eq.3 to write partial line at end of file.
c
c      Adapted by D. Lopez Aldama for ACEMAKER
c
c     -----------------------------------------------------------------
c
      implicit real*8 (a-h, o-z)
      parameter (epsn=1.0d-12)
      character*20 hl(4)
      common/acexss/xss(100000000),nxss
      save hl,i
      if (iflag.eq.3.and.nout.gt.1.and.i.lt.4) then
        write(nout,'(4a20)') (hl(j),j=1,i)
      else
        i=mod(l-1,4)+1
        if (iflag.eq.1) write(hl(i),'(i20)') nint(xss(l)+epsn)
        if (iflag.eq.2) write(hl(i),'(1p,e20.11)') xss(l)
        if (i.eq.4) write(nout,'(4a20)') (hl(j),j=1,i)
      endif
      return
      end
C======================================================================