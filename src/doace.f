      program doace
c     version 2.0
c
c     prepare a continuous energy ace formatted file for MCNP
c
c     To run DOACE all endf data should be linearized and pre-processed.
c     The code system ACEMAKER based on PREPRO package and SIXLIN could
c     be used to produce the input PENDF
c
c     INPUT data:
c
c     Input data options should be entered on the DOACE.INP text file.
c
c     line 1:       isel       imon     mcnpx      (3i11)
c     line 2:       input PENDF filename           (a72)
c     line 3:       output ACE-formatted filename  (a72)
c     line 4:       imat                           (i11)
c     line 5:       suff                           (a4)
c
c     where
c       isel: selection criterium      (0/1) = ZA/MAT   Default = 0
c       imon: Monitor printing trigger (0/1) = min/max  Default = 0
c      mcnpx: MCNP trigger (0/1) = MCNP/MCNPX  Default=0
c       imat: Selected material number
c       suff: ZAID suffix for ACE-formatted file
c             (Examples: .00, .32, .80, .067, the dot '.' is required)
c
c     Example of DOACE.INP:
c
c               0          0          0
c     \PENDF\U235.PENDF
c     \MC\U235.ACEF
c            9228
c     .00
c
c     Retrieve material 9228 (U-235) from PENDF tape \PENDF\U235.PENDF
c     and generate ACE-formatted file U235.ACEF on \MC\ sub-directory.
c     ZAID for U235 should be 92235.00c and minimum printout will be
c     produced during DOACE processing.
c
      implicit real*8 (a-h, o-z)
      parameter (nnxc=300,npmax=2000000,nbmax=2000000,nxsmax=200000000)
      parameter (tol=0.001d0,ymin=1.0d-30,bk=8.6173303d-11)
      character*1 line1(80),line2(80)
      character*4 suff
      character*10 hd,hm
      character*11 zsymam,str11
      character*13 hz
      character*66 line
      character*70 hk
      character*72 fin1,fout
      common/acetxt/hz,hd,hm,hk
      common/acecte/awr0,tz,awn(16),izn(16)
      common/acepnt/nxs(16),jxs(32)
      common/acedat/xss(200000000),nxss
      dimension nbt(20),ibt(20)
      dimension mf3(nnxc),mf4(nnxc),mf5(nnxc),mf6(nnxc)
      dimension mtr(nnxc),mtyr(nnxc)
      dimension x(npmax),x1(npmax),x2(npmax)
      dimension y(npmax),y1(npmax),y2(npmax)
      dimension b(nbmax)
      data ev2mev/1.0d-6/,s2shak/1.0d-8/
      data line1/80*'='/,line2/80*'-'/,suff/'    '/
      data in1/2/,nin/3/,iou/11/,nou/12/
c
c      open input file DOACE.INP and list file DOACE.LST
c
      open (in1, file='DOACE.INP')
      open (iou, file='DOACE.LST')
      write(*,*)' PROGRAM DOACE: Prepare continuous energy ACE files'
      write(*,'(1x,80a1)')line1
      write(*,*)
      write(iou,*)' PROGRAM DOACE: Prepare continuous energy ACE files'
      write(iou,'(80a1)')line1
      write(iou,*)
c
c      read input data from DOACE.INP
c
      isel=0
      imon=0
      mcnpx=0
      read(in1,'(3i11)')isel,imon,mcnpx
      if (isel.ne.1) isel=0
      if (imon.ne.1) imon=0
      if (mcnpx.ne.1)mcnpx=0
      read(in1,'(a72)')fin1
      read(in1,'(a72)')fout
      read(in1,'(i11)')nsel
      if (nsel.ne.1) nsel=0
      read(in1,'(a)')suff
      close(in1)
c
c      Printing input data
c
      open (nin, file=fin1)
      write(iou,*)' Input parameters'
      write(iou,'(80a1)')line1
      write(iou,*)' MAT/ZA selection  =',isel
      write(iou,*)' Printing option   =',imon
      write(iou,*)' MCNP trigger      =',mcnpx
      write(iou,*)' Input file name   =',fin1
      write(iou,*)' Output file name  =',fout
      if (nsel.eq.0) then
        write(iou,*)' Selected material = first material'
      else
        write(iou,*)' Selected material =',nsel
      endif
      write(iou,'(a21,a)')'  Suffix            =',suff
c
c      Initialize pointers arrays
c
      do i=1,16
        nxs(i)=0
        izn(i)=0
        awn(i)=0.0d0
      enddo
      do i=1,32
        jxs(i)=0
      enddo
      do i=1,70
        hk(i:i)=' '
      enddo
      nxss=nxsmax
c
c      Read tape header and find selected material
c
      call readtext(nin,line,mat,mf,mt,nsi)
      if (nsel.ne.0) then
        call findmat(nin,nsel,icod)
        if (icod.ne.0) then
          write(*,*)' Material =',nsel,' not found on tape ',fin1
          close(nin)
          close(iou)
          stop
        endif
      endif
c
c      Material to be processed. Reading MF1/MT451
c
      call readcont(nin,za0,awr0,lrp,lfi,nlib,n2,mat0,mf,mt,nsi)
      call readcont(nin,elis,sta,lis,liso,n1,nfor,mat,mf,mt,nsi)
      call readcont(nin,awi,emax,lrel,l2,nsub,nver,mat,mf,mt,nsi)
      call readcont(nin,temp,c2,ldrv,l2,nwd,nxc,mat,mf,mt,nsi)
c
c      Check for incident particle type (library type)
c
      if (nsub.eq.5.or.nsub.eq.10) then
        zai=1.0d0
        izai=nint(zai)
      else
        write(iou,*)'  === Error: incident particle is not coded yet'
        write(iou,*)'  === NSUB=',nsub,' IPART=',(nsub/10)
        close(nin)
        close(iou)
        stop
      endif
      matza=nint(za0+1.0d-6)
      if (liso.gt.0) then
        if (matza.eq.95242) then
          izaid=matza
        else
          izaid=matza+300+100*liso
        endif
      else
        if (matza.eq.95242) then
          izaid=matza+400
        else
          izaid=matza
        endif
      endif
      call readtext(nin,line,mat,mf,mt,nsi)
      zsymam=line(1:11)
      call readtext(nin,line,mat,mf,mt,nsi)
      call getdtime(hd,str11)
      str11=' '
      call readtext(nin,line,mat,mf,mt,nsi)
      hk(1:11)=zsymam
      hk(12:15)='  T='
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
      hk(25+k:42+k)=line(5:22)
      hk(43+k:53+k)=' (ACEMAKER)'
      str11=' '
      if (mcnpx.eq.1) then
        write(hz,'(i6,a4,a3)')izaid,suff(1:4),'nc '
      else
        write(hz,'(i6,a3,a4)')izaid,suff(1:3),'c   '
      endif
      write(hm,'(a6,i4)')'   mat',mat0
      tz=bk*temp
      write(*,*)' Material=',mat0,' read'
      write(iou,*)
      write(iou,'(a,i7,a,1pe15.8)')' ZAI=',izai,' AWI=',awi
      write(iou,'(a,i5,a,i7,a,i3,a,i7,a,1pe15.8)')' Material=',mat,
     &  ' ZA=',matza,' LISO=',liso,' ZAID=',izaid,' AWR=',awr0
      write(iou,'(a,a,a,i3,a,1pe15.8)')' SYM=',zsymam,
     &  ' LFI=',lfi,' ELIS=',elis
      write(iou,'(a,1p,e15.8,a,e13.6,a,e13.6,a)')' EMAX=', emax,
     &  ' Temperature=',temp,' K = ',tz,' MeV'
      izt=matza/1000
      iat=mod(matza,1000)
      nxs(2)=matza
      nxs(9)=liso
      nxs(10)=izt
      nxs(11)=iat
c
c      Explore endf-6 formatted input file
c
      call findmt(nin,mat0,1,452,icod)
      if (icod.eq.0) then
        mt452=1
      else
        mt452=0
      endif
      call findmt(nin,mat0,1,455,icod)
      if (icod.eq.0) then
        mt455=1
      else
        mt455=0
      endif
      call findmt(nin,mat0,1,456,icod)
      if (icod.eq.0) then
        mt456=1
      else
        mt456=0
      endif
      if (lfi.le.0.and.(mt452.ne.0.or.mt456.ne.0)) then
        lfi=1
        write(iou,*)' nu-fission data found. lfi reset to 1'
      elseif (lfi.gt.0.and.mt452.eq.0.and.mt456.eq.0) then
        lfi=0
        write(iou,*)' nu-fission data not found. lfi reset to 0'
      endif
      call findmt(nin,mat0,2,153,icod)
      if (icod.eq.0) then
        mt153=1
        write(iou,*)' probability tables are available in the URR'
      else
        mt153=0
        write(iou,*)' probability tables are not available in the URR'
      endif
      nmf3=0
      call findmf(nin,mat0,3,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,3,mt)
          ival=mtvalid(izai,mt)
          if (ival.gt.0) then
            nmf3=nmf3+1
            mf3(nmf3)=mt
          endif
        enddo
        write(iou,*)' ',nmf3,' sections found on MF3'
      else
        write(iou,*)' ERROR: MF3 not found'
        close(nin)
        close(iou)
        stop
      endif
      nmf4=0
      call findmf(nin,mat0,4,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,4,mt)
          ival=mtvalid(izai,mt)
          if (ival.gt.0) then
            nmf4=nmf4+1
            mf4(nmf4)=mt
          endif
        enddo
        write(iou,*)' ',nmf4,' sections found on MF4'
      else
        write(iou,*)' MF4 not found'
      endif
      nmf5=0
      call findmf(nin,mat0,5,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,5,mt)
          ival=mtvalid(izai,mt)
          if (ival.gt.0) then
            nmf5=nmf5+1
            mf5(nmf5)=mt
          endif
        enddo
        write(iou,*)' ',nmf5,' sections found on MF5'
      else
        write(iou,*)' MF5 not found'
      endif
      nmf6=0
      call findmf(nin,mat0,6,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,6,mt)
          ival=mtvalid(izai,mt)
          if (ival.gt.0) then
            nmf6=nmf6+1
            mf6(nmf6)=mt
          endif
        enddo
        write(iou,*)' ',nmf6,' sections found on MF6'
      else
        write(iou,*)' MF6 not found'
      endif
      if (nmf4.eq.0.and.nmf6.eq.0) then
        write(iou,*)' ERROR: No angular distribution found'
        close(nin)
        close(iou)
        stop
      endif
      m5=iposm(nmf5,mf5,455)
      if (m5.gt.0) then
        mt5455=1
      else
        mt5455=0
      endif
      if ((mt455.eq.1.and.mt5455.eq.0).or.
     &    (mt455.eq.0.and.mt5455.gt.0)) then
        write(iou,*)' ERROR: Delayed neutron data incomplete'
        close(nin)
        close(iou)
        stop
      endif
c
c      Find competitive reactions in the URR, if any
c
      if (mt153.gt.0) then
        call findmt(nin,mat0,2,153,icod)
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
        call readcont(nin,c1,c2,l1,lcomp,n1,n2,mat,mf,mt,ns)
        mtabso=lcomp/1000
        mtinel=lcomp-mtabso*1000
      else
        mtabso=0
        mtinel=0
      endif
c
c      Check mt=16 (z,2n) versus partials mt=875-891
c
      mt16=mtchkd(16,875,891,nmf3,mf3,nmf4,mf4,nmf6,mf6)
c
c      Check mt=18 (fission) vs. mt=19,20,21,38 (partials for neutrons)
c
      if (lfi.gt.0) then
        if (izai.eq.1) then
          mt18=mtchkd(18,19,21,nmf3,mf3,nmf4,mf4,nmf6,mf6)
        else
          mt18=1
        endif
      else
        mt18=0
      endif
c
c      Check if mt=5 produces incident particle
c
      mt05=iposm(nmf3,mf3,5)
      if (mt05.gt.0) then
        mt05=0
        m6=iposm(nmf6,mf6,5)
        if (m6.gt.0) then
          call findmt(nin,mat0,6,5,icod)
          call readcont(nin,c1,c2,l1,lct5,nk5,n2,mat,mf,mt,ns)
          do i=1,nk5
            call readtab1(nin,zap,awp,l1,law,n1,n2,nbt,ibt,x,y)
            izap=nint(zap)
            if (izap.eq.izai) mt05=mt05+1
            call nextsub6(nin,law,nbt,ibt,x,b)
          enddo
        endif
      endif
c
c      Check for mt=103-107 vs. partials
c
      mt103=mtchkd(103,600,649,nmf3,mf3,nmf4,mf4,nmf6,mf6)
      mt104=mtchkd(104,650,699,nmf3,mf3,nmf4,mf4,nmf6,mf6)
      mt105=mtchkd(105,700,749,nmf3,mf3,nmf4,mf4,nmf6,mf6)
      mt106=mtchkd(106,750,799,nmf3,mf3,nmf4,mf4,nmf6,mf6)
      mt107=mtchkd(107,800,849,nmf3,mf3,nmf4,mf4,nmf6,mf6)
c
c      Printing XS control parameters
c
      write(iou,*)' lfi   =',lfi
      write(iou,*)' mt452 =',mt452
      write(iou,*)' mt455 =',mt455
      write(iou,*)' mt456 =',mt456
      write(iou,*)' mt5455=',mt5455
      write(iou,*)' mt153 =',mt153
      write(iou,*)' mtinel=',mtinel
      write(iou,*)' mtabso=',mtabso
      write(iou,*)' mt16  =',mt16
      write(iou,*)' mt18  =',mt18
      write(iou,*)' mt05  =',mt05
      write(iou,*)' mt103 =',mt103
      write(iou,*)' mt104 =',mt104
      write(iou,*)' mt105 =',mt105
      write(iou,*)' mt106 =',mt106
      write(iou,*)' mt107 =',mt107
      write(iou,'(80a1)')line1
c
c                  MTR array
c      First reactions that produce  neutrons
c
      nmtr=0
      do i=1,nmf3
        mt=mf3(i)
        iproc=mtprodi(izai,mt)
        if (izai.eq.1) then
          if ((mt.eq.5.and.mt05.eq.0).or.
     &        (mt.eq.16.and.mt16.eq.0).or.
     &        (mt.ge.875.and.mt.le.891.and.mt16.eq.1).or.
     &        (((mt.ge.19.and.mt.le.21).or.mt.eq.38).and.mt18.eq.1).or.
     &        (mt.eq.18.and.mt18.eq.0)) iproc=0
        endif
        if (iproc.gt.0) then
          nmtr=nmtr+1
          mtr(nmtr)=mt
        endif
      enddo
c
c      Rest of reactions
c
      nmtrt=nmtr
      do i=1,nmf3
        mt=mf3(i)
        ipos=iposm(nmtr,mtr,mt)
        if (ipos.eq.0) then
          iproc=mtrval(izai,mt)
          if (iproc.gt.0) then
            nmtrt=nmtrt+1
            mtr(nmtrt)=mt
          endif
        endif
      enddo
      if (mtinel.gt.0) then
        ipos=iposm(nmtrt,mtr,mtinel)
        if (ipos.eq.0) then
          ipos3=iposm(nmf3,mf3,mtinel)
          if (ipos3.gt.0) then
            nmtrt=nmtrt+1
            mtr(nmtrt)=mtinel
          else
            write(iou,*)' ERROR: competitive MT= ',mtinel,
     &      ' not found on tape'
            stop
          endif
        endif
      endif
      if (mtabso.gt.0) then
        ipos=iposm(nmtrt,mtr,mtabso)
        if (ipos.eq.0) then
          ipos3=iposm(nmf3,mf3,mtabso)
          if (ipos3.gt.0) then
            nmtrt=nmtrt+1
            mtr(nmtrt)=mtabso
          else
            write(iou,*)' ERROR: competitive MT= ',mtabso,
     &      ' not found on tape'
            stop
          endif
        endif
      endif
      nxs(4)=nmtrt
      nxs(5)=nmtr
      write(iou,*)' MTR array nmtrt=',nxs(4),' nmtr=',nxs(5)
      write(iou,'(10i8)')(mtr(i),i=1,nmtrt)
      write(iou,'(80a1)')line1
c
c      SET ESZ BLOCK
c
      if (izai.eq.1) then
        ipos=iposm(nmf3,mf3,1)
        if (ipos.eq.0) then
          write(iou,*)' ERROR: MT=1 data not found for neutrons'
          stop
        endif
        ipos=iposm(nmf3,mf3,2)
        if (ipos.eq.0) then
          write(iou,*)' ERROR: MT=2 data not found for neutrons'
          stop
        endif
        call findmt(nin,mat0,3,1,icod)
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
        call readtab1(nin,qm,qi,l1,lr,nr,ne1,nbt,ibt,x1,y1)
        call checklaw(nr,ibt,icod)
        if (icod.ne.0) then
          write(iou,*)' MF3/MT= 1 linear interpolation expected'
          write(iou,*)' LIN-LIN interpolation assumed'
        endif
        call findmt(nin,mat0,3,2,icod)
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
        call readtab1(nin,qm,qi,l1,lr,nr,ne2,nbt,ibt,x,y)
        call checklaw(nr,ibt,icod)
        if (icod.ne.0) then
          write(iou,*)' MF3/MT= 2 linear interpolation expected'
          write(iou,*)' LIN-LIN interpolation assumed'
        endif
        call checkeq(ne1,x1,ne2,x,icod)
        if (icod.ne.0) then
          write(iou,*)' Not unified grid found. Grid of total XS used'
        endif
        call checkdis(iou,ne1,x1,y1,icod)
        if (icod.gt.0) then
          write(iou,*)' ',icod,' discontinuities found and removed'
        endif
        do k=1,nmtrt
          mti=mtr(k)
          call findmt(nin,mat0,3,mti,icod)
          call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
          call readcont(nin,qm,qi,l1,lr,nr,ne,mat,mf,mt,ns)
          if (qi.lt.0.0d0) then
            xx=-qi*(awr0+awi)/awr0
            call ff2chx(xx,str11)
            read(str11,'(e11.0)')xx
            x(1)=edelta(xx,1.0d0)
            call union(1,x,ne1,x1,ne2,x2,npmax)
            ne1=ne2
            do i=1,ne1
              x1(i)=x2(i)
            enddo
          endif
        enddo
        call checkmon(iou,ne1,x1,icod)
        if (icod.ne.0) then
          write(iou,*)' ERROR: ',icod,' non monotonic energies found'
          write(iou,'(1p5e16.9)')(x1(k),k=1,ne1)
          stop
        endif
        call findmt(nin,mat0,3,2,icod)
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
        call readtab1(nin,qm,qi,l1,lr,nr,ne2,nbt,ibt,x,y)
        do i=1,ne1
          y2(i)=fvalue(nr,nbt,ibt,ne2,x,y,x1(i))
        enddo
      endif
c
c      SET LOCATORS & INITIALIZE ESZ BLOCK
c
      write(iou,*)' Processing ESZ block locators'
      kesz=1
      jxs(1)=kesz
      nxs(3)=ne1
      ktot=kesz+ne1
      kcap=ktot+ne1
      kela=kcap+ne1
      khea=kela+ne1
      do i=1,ne1
        im1=i-1
        xss(kesz+im1)=x1(i)*ev2mev
        xss(ktot+im1)=y2(i)
        xss(kcap+im1)=0.0d0
        xss(kela+im1)=y2(i)
        xss(khea+im1)=0.0d0
      enddo
      lxs=khea+ne1
      write(iou,*)' ESZ block locators prepared'
      write(iou,'(80a1)')line1
c
c      NU BLOCK
c
      write(iou,*)' Processing NU block'
      if (lfi.gt.0) then
        knu=lxs
        if (mt452.gt.0.and.mt456.gt.0) lxs=lxs+1
        if (mt456.gt.0) then
          call findmt(nin,mat0,1,456,icod)
          call readcont(nin,c1,c2,l1,lnu,n1,n2,mat,mf,mt,ns)
          if (lnu.ne.2) then
            write(iou,*)' ERROR: LNU not equal 2 for MT456'
            close(nin)
            close(iou)
            stop
          endif
          call readtab1(nin,c1,c2,l1,l2,nr,ne,nbt,ibt,x,y)
          call checklaw(nr,ibt,icod)
          if (icod.ne.0) then
            write(iou,*)' MF1/MT= 456 linear interpolation expected'
            write(iou,*)' LIN-LIN interpolation assumed'
          endif
          xss(lxs)=lnu
          lxs=lxs+1
          xss(lxs)=0
          lxs=lxs+1
          xss(lxs)=ne
          if (imon.gt.0) then
            write(iou,*)'  LNU=',lnu,' ne=',ne
            write(iou,*)'  Energy [MeV]       ',' Prompt NU-BAR '
          endif
          do i=1,ne
            xss(lxs+i)=x(i)*ev2mev
            xss(lxs+ne+i)=y(i)
            if (imon.gt.0) then
              write(iou,'(1p2e20.11)')xss(lxs+i),xss(lxs+ne+i)
            endif
          enddo
          lxs=lxs+2*ne+1
          if (imon.gt.0) write(iou,'(80a1)')line2
        endif
        if (mt452.gt.0.and.mt456.gt.0) xss(knu)=-(lxs-knu-1)
        if (mt452.gt.0) then
          call findmt(nin,mat0,1,452,icod)
          call readcont(nin,c1,c2,l1,lnu,n1,n2,mat,mf,mt,ns)
          if (lnu.ne.2) then
            write(iou,*)' ERROR: LNU not equal 2 for MT452'
            close(nin)
            close(iou)
            stop
          endif
          call readtab1(nin,c1,c2,l1,l2,nr,ne,nbt,ibt,x,y)
          call checklaw(nr,ibt,icod)
          if (icod.ne.0) then
            write(iou,*)' MF1/MT= 452 linear interpolation expected'
            write(iou,*)' LIN-LIN interpolation assumed'
          endif
          xss(lxs)=lnu
          lxs=lxs+1
          xss(lxs)=0
          lxs=lxs+1
          xss(lxs)=ne
          if (imon.gt.0) then
            write(iou,*)'  LNU=',lnu,' ne=',ne
            write(iou,*)'  Energy [MeV]       ',' Total NU-BAR '
          endif
          do i=1,ne
            xss(lxs+i)=x(i)*ev2mev
            xss(lxs+ne+i)=y(i)
            if (imon.gt.0) then
              write(iou,'(1p2e20.11)')xss(lxs+i),xss(lxs+ne+i)
            endif
          enddo
          lxs=lxs+2*ne+1
          if (imon.gt.0) write(iou,'(80a1)')line2
        endif
      else
        knu=0
      endif
      jxs(2)=knu
      write(iou,*)' LNU=',knu
      write(iou,*)' NU block prepared'
      write(iou,'(80a1)')line1
c
c      SET MTR, LQR, TYR, LSIG & SIG
c
c      first reactions that produce  neutrons
c
      write(iou,*)' Processing MTR, LQR, TYR1, LSIG & SIG blocks'
      kmtr=lxs
      jxs(3)=kmtr
      kqr=kmtr+nmtrt
      jxs(4)=kqr
      ktyr=kqr+nmtrt
      jxs(5)=ktyr
      klsig=ktyr+nmtrt
      jxs(6)=klsig
      ksig=klsig+nmtrt
      jxs(7)=ksig
      lxs=ksig
c
c      reactions producing incident particles
c
      write(iou,*)' Processing reactions producing incident particle'
      do i=1,nmtr
        mti=mtr(i)
        write(iou,*)' Processing reaction MF3/MT=',mti
        write(*,*)' Processing reaction MF3/MT=',mti
        write(iou,*)'  mt= ',mti,' contributes to total'
        call typr(nin,izai,mat0,mti,nmf4,mf4,nmf6,mf6,tyr,mtyr(i),
     &            nbt,ibt,x,y,b)
        call findmt(nin,mat0,3,mti,icod)
        call readcont(nin,za,awr,l1,l2,n1,n2,mat,mf,mt,ns)
        call readtab1(nin,qm,qi,l1,lr,nr,ne,nbt,ibt,x,y)
        call checklaw(nr,ibt,icod)
        if (icod.ne.0) then
          write(iou,*)'  MF3/MT=',mti,' linear interpolation expected'
          write(iou,*)'  LIN-LIN interpolation assumed'
        endif
        call thresh(iou,mti,qi,awr0,awi,nr,ne,nbt,ibt,x,y)
        xss(kmtr+i-1)=mti
        xss(kqr+i-1)=qi*ev2mev
        xss(ktyr+i-1)=tyr
        xss(klsig+i-1)=lxs-ksig+1
        if (mti.eq.18) jxs(21)=lxs
        je=iposx(ne1,x1,x(1))
        xss(lxs)=je
        lxs=lxs+1
        xss(lxs)=ne1-je+1
        lxs=lxs+1
        do j=je,ne1
          yy=fvalue(nr,nbt,ibt,ne,x,y,x1(j))
          xss(lxs)=yy
          lxs=lxs+1
          jtot=ktot+j-1
          xss(jtot)=xss(jtot)+yy
        enddo
        write(*,*)'  ne1=',ne1,' ie=',je,' ne=',ne1-je+1,' qi=',qi
        write(*,*)'  tyr1=',nint(tyr),' lsig=',nint(xss(klsig+i-1))
        write(*,*)'  x(ie)=',x1(je),' y(ie)=',xss(lxs-ne1+je-1)
        write(*,*)'  x(ne)=',x1(ne1),' y(ne)=',xss(lxs-1)
        write(iou,*)'  ne1=',ne1,' ie=',je,' ne=',ne1-je+1,' qi=',qi
        write(iou,*)'  tyr1=',nint(tyr),' lsig=',nint(xss(klsig+i-1))
        write(iou,*)'  x(ie)=',x1(je),' y(ie)=',xss(lxs-ne1+je-1)
        write(iou,*)'  x(ne)=',x1(ne1),' y(ne)=',xss(lxs-1)
        write(iou,'(80a1)')line2
      enddo
c
c      Reactions not producing incident particle
c
      write(iou,*)' Reactions not producing incident particle'
      tyr=0.0d0
      i0=nmtr+1
      do i=i0,nmtrt
        mti=mtr(i)
        if (mti.ne.301) then
          write(iou,*)' Processing reaction MF3/MT=',mti
          write(*,*)' Processing reaction MF3/MT=',mti
          idisa=mtdisa(izai,mti)
          if (izai.eq.1.and.mti.ge.600.and.mti.le.849) then
            if ((mti.ge.600.and.mti.le.649.and.mt103.gt.0).or.
     &          (mti.ge.650.and.mti.le.699.and.mt104.gt.0).or.
     &          (mti.ge.700.and.mti.le.749.and.mt105.gt.0).or.
     &          (mti.ge.750.and.mti.le.799.and.mt106.gt.0).or.
     &          (mti.ge.800.and.mti.le.849.and.mt107.gt.0)) idisa=0
          endif
          if (idisa.gt.0) then
            write(iou,*)'  mt= ',mti,' contributes to capture and total'
          endif
          call findmt(nin,mat0,3,mti,icod)
          call readcont(nin,za,awr,l1,l2,n1,n2,mat,mf,mt,ns)
          call readtab1(nin,qm,qi,l1,lr,nr,ne,nbt,ibt,x,y)
          call checklaw(nr,ibt,icod)
          if (icod.ne.0) then
            write(iou,*)'  MF3/MT=',mti,' linear interpolation expected'
            write(iou,*)'  LIN-LIN interpolation assumed'
          endif
          call thresh(iou,mti,qi,awr0,awi,nr,ne,nbt,ibt,x,y)
          xss(kmtr+i-1)=mti
          xss(kqr+i-1)=qi*ev2mev
          xss(ktyr+i-1)=tyr
          xss(klsig+i-1)=lxs-ksig+1
          if (mti.eq.18.and.izai.ne.1) jxs(21)=lxs
          je=iposx(ne1,x1,x(1))
          xss(lxs)=je
          lxs=lxs+1
          xss(lxs)=ne1-je+1
          lxs=lxs+1
          do j=je,ne1
            yy=fvalue(nr,nbt,ibt,ne,x,y,x1(j))
            if (mti.eq.444) yy=yy*ev2mev
            xss(lxs)=yy
            lxs=lxs+1
            if (idisa.gt.0) then
              jm1=j-1
              jtot=ktot+jm1
              xss(jtot)=xss(jtot)+yy
              jcap=kcap+jm1
              xss(jcap)=xss(jcap)+yy
            endif
          enddo
          write(*,*)'  ne1=',ne1,' ie=',je,' ne=',ne1-je+1,' qi=',qi
          write(*,*)'  tyr1=',nint(tyr),' lsig=',nint(xss(klsig+i-1))
          write(*,*)'  x(ie)=',x1(je),' y(ie)=',xss(lxs-ne1+je-1)
          write(*,*)'  x(ne)=',x1(ne1),' y(ne)=',xss(lxs-1)
          write(iou,*)'  ne1=',ne1,' ie=',je,' ne=',ne1-je+1,' qi=',qi
          write(iou,*)'  tyr1=',nint(tyr),' lsig=',nint(xss(klsig+i-1))
          write(iou,*)'  x(ie)=',x1(je),' y(ie)=',xss(lxs-ne1+je-1)
          write(iou,*)'  x(ne)=',x1(ne1),' y(ne)=',xss(lxs-1)
          write(iou,'(80a1)')line2
        endif
      enddo
c
c      Processing heating factors
c
      ipos=iposm(nmtrt,mtr,301)
      if (ipos.gt.0) then
        write(iou,*)' Processing reaction MF3/MT= 301'
        write(*,*)' Processing reaction MF3/MT= 301'
        call findmt(nin,mat0,3,301,icod)
        call readcont(nin,za,awr,l1,l2,n1,n2,mat,mf,mt,ns)
        call readtab1(nin,qm,qi,l1,lr,nr,ne,nbt,ibt,x,y)
        call checklaw(nr,ibt,icod)
        if (icod.ne.0) then
          write(iou,*)'  MF3/MT=',mti,' linear interpolation expected'
          write(iou,*)'  LIN-LIN interpolation assumed'
        endif
        do i=1,ne1
          yy=fvalue(nr,nbt,ibt,ne,x,y,x1(i))
          im1=i-1
          yt=xss(ktot+im1)
          if (yt.eq.0.0d0) then
            xss(khea+im1)=0.0d0
          else
            xss(khea+im1)=yy/yt*ev2mev
          endif
        enddo
      endif
      if (imon.gt.0) then
        write(iou,*)'  Energy [MeV]       ',' TOTAL XS [barn]    ',
     &  ' Capture XS [barn]  ',' Elastic XS [barn]  ',
     &  ' Heating factor  '
        do i=1,ne1
          im1=i-1
          write(iou,'(1p5e20.11)')xss(kesz+im1),xss(ktot+im1),
     &    xss(kcap+im1),xss(kela+im1),xss(khea+im1)
        enddo
        write(iou,'(80a1)')line2
      endif
      write(iou,*)' MTR BLOCK: kmtr=', kmtr
      write(iou,'(5f16.0)')(xss(kmtr+i-1),i=1,nmtrt)
      write(iou,*)' QR BLOCK: kqr=', kqr
      write(iou,'(1p5e16.7)')(xss(kqr+i-1),i=1,nmtrt)
      write(iou,*)' TYR1 BLOCK: ktyr=', ktyr
      write(iou,'(1p5e16.7)')(xss(ktyr+i-1),i=1,nmtrt)
      write(iou,*)' LSIG BLOCK: klsig=', klsig
      write(iou,'(1p5e16.7)')(xss(klsig+i-1),i=1,nmtrt)
      write(iou,*)' ksig=',ksig,' kfis=',jxs(21),' lxs=',lxs
      write(iou,*)' MTR, LQR, TYR1, LSIG & SIG blocks prepared'
      write(iou,'(80a1)')line1
c
c      LAND and AND block for (MT=2 + first nmtr reactions)
c
      write(iou,*)' Processing LAND & AND blocks'
      kland=lxs
      jxs(8)=kland
      nmtr1=nmtr+1
      kand=kland+nmtr1
      jxs(9)=kand
      lxs=kand
      do k=1,nmtr1
        if (k.eq.1) then
          mti=2
          mparti=1
        else
          mti=mtr(k-1)
          mparti=mtyr(k-1)
        endif
        m4=iposm(nmf4,mf4,mti)
        m6=iposm(nmf6,mf6,mti)
        if (m4.gt.0.and.m6.eq.0) then
          call findmt(nin,mat0,4,mti,icod)
          call readcont(nin,c1,c2,lvt,ltt,n1,n2,mat,mf,mt,ns)
          if (lvt.eq.1.and.mti.eq.2) then
            call readlist(nin,c1,c2,li,lct,n1,n2,b)
          else
            call readcont(nin,c1,c2,li,lct,n1,n2,mat,mf,mt,ns)
          endif
          if (ltt.eq.0.and.li.eq.1) then
            write(*,*)' mti=', mti,' MF4 LTT=0 LI=1 isotropic'
            write(iou,*)' mti=', mti,' MF4 LTT=0 LI=1 isotropic'
            xss(kland+k-1)=0
          elseif (ltt.eq.2) then
            write(*,*)' mti=', mti,' MF4 LTT=2 LI=0'
            write(iou,*)' mti=', mti,' MF4 LTT=2 LI=0'
            call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
            call checklaw(nr,ibt,icod)
            if (icod.ne.0.and.icod.ne.1) then
              write(iou,*)' Incident energy grid is not linearly',
     &        ' interpolable. Linear interpolation assumed'
            endif
            xss(kland+k-1)=lxs-kand+1
            le=lxs
            xss(le)=ne
            lc=le+ne
            lxs=lc+ne+1
            do j=1,ne
              call readtab1(nin,c1,e,l1,l2,nr,np,nbt,ibt,x,y)
              call checklaw(nr,ibt,icod)
              if (icod.eq.0) then
                jj=2
              elseif (icod.eq.1) then
                jj=1
              else
                write(iou,*)' Non a linearly interpolable cosine grid'
                write(iou,*)' Run LEGEND'
                stop
              endif
              xss(le+j)=e*ev2mev
              write(iou,*)' Incident energy j=',j,' E=',xss(le+j),' MeV'
c             iso=isochk(np,y)
              iso=0
              if (iso.eq.1) then
                if (imon.gt.0) then
                  write(iou,*)'  Isotropic angular distribution'
                endif
                xss(lc+j)=0
              else
                if((lxs+2+3*np).gt.nxss)then
                  write(iou,*)' ERROR: Increase size xss ',nxss
                  stop
                endif
                call cdfcal(np,x,y,jj,y2)
                xss(lc+j)=-(lxs-kand+1)
                xss(lxs)=jj
                imu=lxs+1
                xss(imu)=np
                ipdf=imu+np
                icdf=ipdf+np
                if (imon.gt.0) then
                  write(iou,*)' Angular distribution intt=',jj,' np=',np
                  write(iou,*)'  i  ',' xmu                ',
     &            ' pdf                ',' cdf                '
                endif
                do i=1,np
                  xss(imu+i)=x(i)
                  xss(ipdf+i)=y(i)
                  xss(icdf+i)=y2(i)
                  if (imon.gt.0) then
                    write(iou,'(i4,1p3e20.11)')i,xss(imu+i),
     &              xss(ipdf+i),xss(icdf+i)
                  endif
                enddo
                lxs=lxs+2+3*np
              endif
            enddo
          else
            write(iou,*)' ERROR: Only ltt=0 and ltt=2 allowed ltt=',ltt
            write(iou,*)' Run LEGEND'
            stop
          endif
        elseif (m6.gt.0.and.m4.eq.0) then
          call findmt(nin,mat0,6,mti,icod)
          call readcont(nin,c1,c2,l1,lct,nk6,n2,mat,mf,mt,ns)
          kpart=1
          kk=1
          do while (kpart.le.mparti.and.kk.le.nk6)
            call readtab1(nin,zap,awp,lip,law,nr,ne,nbt,ibt,x,y)
            izap=nint(zap)
            if (izap.eq.izai) then
              if (law.eq.1.or.law.eq.6) then
                write(*,*)' mti=',mti,' MF6 LAW=',law,' part=',kpart
                write(iou,*)' mti=',mti,' MF6 LAW=',law,
     &          ' part=',kpart
                if (kpart.eq.1) then
                  xss(kland+k-1)=-1
                elseif (xss(kland+k-1).ne.-1) then
                  write(iou,*)' ERROR: Multiple output particles with ',
     &            'incompatible LAW on MF6'
                  stop
                endif
              elseif (law.eq.3) then
                write(*,*)' mti=',mti,' MF6 LAW= 3 isotropic'
                write(iou,*)' mti=',mti,' MF66 LAW= 3 isotropic'
                if (kpart.eq.1.and.mparti.eq.1) then
                  xss(kland+k-1)=0
                  mf6(m6)=0
                else
                  write(*,*)' ERROR: Use of LAW= 3 isotropic two-',
     &            'body distribution is not expected with mpart=',mparti
                  write(iou,*)' ERROR: Use of LAW= 3 isotropic two-',
     &            'body distribution is not expected with mpart=',mparti
                  stop
                endif
              elseif (law.eq.2) then
                write(*,*)' mti=',mti,' MF6 LAW= 2'
                write(iou,*)' mti=',mti,' MF6 LAW= 2'
                if (kpart.eq.1.and.mparti.eq.1) then
                  call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
                  call checklaw(nr,ibt,icod)
                  if (icod.ne.0.and.icod.ne.1) then
                    write(iou,*)' Incident energy grid is not linearly',
     &              ' interpolable. Linear interpolation assumed'
                  endif
                  xss(kland+k-1)=lxs-kand+1
                  le=lxs
                  xss(le)=ne
                  lc=le+ne
                  lxs=lc+ne+1
                  do j=1,ne
                    call readlist(nin,c1,e,lang,l2,nw,np,b)
                    if (lang.eq.0) then
                      nl=np
                      do i=nl,1,-1
                        b(i+1)=b(i)
                      enddo
                      b(1)=1.0d0
                      call leg2lin(nl,b,np,x,y,tol,ymin,npmax)
                    elseif (lang.eq.12.or.lang.eq.14) then
                      i=0
                      do ij=1,np
                        i=i+1
                        x(ij)=b(i)
                        i=i+1
                        y(ij)=b(i)
                      enddo
                      if (lang.ne.12) then
                        nr=1
                        nbt(1)=np
                        ibt(1)=lang-10
                        call linear(nr,nbt,ibt,np,x,y,tol,ymin,npmax)
                        do i=1,np
                          if (y1(i).lt.0.0d0) y1(i)=1.0d-30
                        enddo
                        jj=2
                        c=1.0d0
                        call renorm(np,x,y,jj,c,fn)
                      endif
                    else
                      write(iou,*)' MF6/LAW=2 LANG=',lang,' not allowed'
                      stop
                    endif
                    xss(le+j)=e*ev2mev
                    write(iou,*)' Incident energy j=',j,' E=',xss(le+j),
     &                ' MeV'
c                   iso=isochk(np,y)
                    iso=0
                    if (iso.eq.1) then
                      if (imon.gt.0) then
                        write(iou,*)'  Isotropic angular distribution'
                      endif
                      xss(lc+j)=0
                    else
                      if((lxs+2+3*np).gt.nxss)then
                        write(iou,*)' ERROR: Increase size xss ',nxss
                        stop
                      endif
                      jj=2
                      call cdfcal(np,x,y,jj,y2)
                      xss(lc+j)=-(lxs-kand+1)
                      xss(lxs)=jj
                      imu=lxs+1
                      xss(imu)=np
                      ipdf=imu+np
                      icdf=ipdf+np
                      if (imon.gt.0) then
                        write(iou,*)' Angular distribution intt=',jj,
     &                  ' np=',np
                        write(iou,*)' i   ',' xmu                ',
     &                  ' pdf                ',' cdf                '
                      endif
                      do i=1,np
                        xss(imu+i)=x(i)
                        xss(ipdf+i)=y(i)
                        xss(icdf+i)=y2(i)
                        if (imon.gt.0) then
                          write(iou,'(i4,1p3e20.11)')i,xss(imu+i),
     &                    xss(ipdf+i),xss(icdf+i)
                        endif
                      enddo
                      lxs=lxs+2+3*np
                    endif
                  enddo
                  mf6(m6)=0
                else
                  write(*,*)' ERROR: Use of LAW= 2 two-body reaction',
     &            'distribution is not expected with mpart=',mparti
                  write(iou,*)' ERROR: Use of LAW= 2 two-body reaction',
     &            'distribution is not expected with mpart=',mparti
                  stop
                endif
              elseif (law.eq.7) then
                write(iou,*)' ERROR: MF6/MT=',mti,' LAW 7 not allowed.',
     &          ' Run SIXLIN'
                stop
              endif
              kpart=kpart+1
            else
              call nextsub6(nin,law,nbt,ibt,x,b)
            endif
            kk=kk+1
          enddo
        endif
      enddo
      write(iou,*)' LAND & AND blocks prepared'
      write(iou,'(80a1)')line1
c
c      ldlw and dlw
c
      write(iou,*)' Processing LDLW & DLW blocks'
      kldlw=lxs
      jxs(10)=kldlw
      kdlw=kldlw+nmtr
      jxs(11)=kdlw
      lxs=kdlw
      do k=1,nmtr
        mti=mtr(k)
        m5=iposm(nmf5,mf5,mti)
        m6=iposm(nmf6,mf6,mti)
        qi=xss(kqr+k-1)
        if (m5.gt.0) then
          call findmt(nin,mat0,5,mti,icod)
          call readcont(nin,c1,c2,l1,l2,nk,n2,mat,mf,mt,ns)
          write(*,*)' MF5 for mti=',mti,' nk=',nk,' qi=',qi,
     &    ' lxs=',lxs,' loc=',lxs-kdlw+1
          write(iou,*)' MF5 for mti=',mti,' nk=',nk,' qi=',qi,
     &    ' lxs=',lxs,' loc=',lxs-kdlw+1
          xss(kldlw+k-1)=lxs-kdlw+1
          do ik=1,nk
            call readtab1(nin,c1,c2,l1,lf,nr,npp,nbt,ibt,x,y)
            if (lf.ne.1) then
              write(iou,*)' ERROR: LF=',lf,' not allowed'
              write(iou,*)' Run SPECTRA'
              stop
            endif
            if ((lxs+2*(nr+npp)+3).gt.nxss) then
              write(iou,*)' ERROR: Increase size xss ',nxss
              stop
            endif
            if (ik.gt.1) xss(lxs0)=lxs-kdlw+1
            lxs0=lxs
            xss(lxs)=0
            xss(lxs+1)=4
            if (imon.gt.0) then
              write(iou,*)' law= 4 data for k=',ik
            endif
            call checklaw(nr,ibt,icod)
            if (icod.ne.0) then
              j1=lxs+3
              xss(j1)=nr
              j2=j1+nr
              if (imon.gt.0) then
                write(iou,*)'  Probability interpolation nr=',nr
                write(iou,*)'  == nbt ==  == ibt == '
              endif
              do j=1,nr
                xss(j1+j)=nbt(j)
                xss(j2+j)=ibt(j)
                if (imon.gt.0) then
                  write(iou,'(1x,i10,1x,i10)')nbt(j),ibt(j)
                endif
              enddo
              lxs=j2+nr+1
            else
              if (imon.gt.0) then
                write(iou,*)'  LIN-LIN probability interpolation nr= 0'
              endif
              xss(lxs+3)=0
              lxs=lxs+4
            endif
            xss(lxs)=npp
            if (imon.gt.0) then
              write(iou,*)'  Law probability npp=',npp
              write(iou,*)'  Energy [MeV]              Probability   '
            endif
            do j=1,npp
              xss(lxs+j)=x(j)*ev2mev
              xss(lxs+npp+j)=y(j)
              if (imon.gt.0) then
                write(iou,'(1p2e21.11)')xss(lxs+j),xss(lxs+npp+j)
              endif
            enddo
            lxs=lxs+2*npp+1
            xss(lxs0+2)=lxs-kdlw+1
            call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
            if((lxs+2*(nr+ne)+3).gt.nxss)  then
              write(iou,*)' ERROR: Increase size xss ',nxss
              stop
            endif
            call checklaw(nr,ibt,icod)
            if (icod.eq.0) then
              if (imon.gt.0) then
                write(iou,*)'  LIN-LIN energy interpolation nr= 0'
              endif
              xss(lxs)=0
              lxs=lxs+1
            else
              xss(lxs)=nr
              if (imon.gt.0) then
                write(iou,*)'  Energy interpolation nr=',nr
                write(iou,*)'  == nbt ==  == ibt == '
              endif
              do j=1,nr
                xss(lxs+j)=nbt(j)
                xss(lxs+nr+j)=min(mod(ibt(j),10),2)
                if (imon.gt.0) then
                  write(iou,'(1x,i10,1x,i10)')nbt(j),nint(xss(lxs+nr+j))
                endif
              enddo
              lxs=lxs+2*nr+1
            endif
            xss(lxs)=ne
            lxsn=lxs+ne
            lxsd=lxsn+ne+1
            do j=1,ne
              call readtab1(nin,c1,e,l1,l2,nr,nep,nbt,ibt,x,y)
              call checklaw(nr,ibt,icod)
              if (icod.lt.0) then
                write(iou,*)' ERROR: E'' data are not lin-lin or',
     &          ' constant interpolable for MF5/MT=',mti
                write(iou,*)' Use SPECTRA'
                stop
              endif
              if (icod.eq.0)icod=2
              if (qi.lt.0.0d0) then
                iep=0
                i=1
                do while (i.le.nep.and.iep.eq.0)
                  if (x(i).gt.e)iep=i
                  i=i+1
                enddo
                if (iep.gt.0) then
                  xx=0.9999999d0*e
                  y(iep)=fvalue(nr,nbt,ibt,nep,x,y,xx)
                  x(iep)=xx
                  nep=iep
                endif
              endif
              if (nep.gt.npmax) then
                write(iou,*)' ERROR: Increase TAB1 arrays size ',npmax
                stop
              endif
              if ((lxs+2+3*nep).gt.nxss) then
                write(iou,*)' ERROR: Increase size xss ',nxss
                stop
              endif
              call cdfcal(nep,x,y,icod,y2)
              xss(lxs+j)=e*ev2mev
              xss(lxsn+j)=lxsd-kdlw+1
              xss(lxsd)=icod
              lxsd=lxsd+1
              xss(lxsd)=nep
              if (imon.gt.0) then
                write(iou,*)'  Incident energy E=',xss(lxs+j),' MeV'
                write(iou,*)'  E'' intt=',icod,' nep=',nep
                write(iou,*)'   i  ','   Outgoing energy  ',
     &          ' pdf                ',' cdf                '
              endif
              nep2=nep+nep
              do i=1,nep
                xss(lxsd+i)=x(i)*ev2mev
                xss(lxsd+nep+i)=y(i)/ev2mev
                xss(lxsd+nep2+i)=y2(i)
                if (imon.gt.0) then
                  write(iou,'(i5,1p3e20.11)')i,xss(lxsd+i),
     &            xss(lxsd+nep+i),xss(lxsd+nep2+i)
                endif
              enddo
              lxsd=lxsd+3*nep+1
            enddo
            lxs=lxsd
          enddo
          write(iou,'(80a1)')line2
        elseif (m6.gt.0) then
c
c         ace law 44,61,66
c
          call findmt(nin,mat0,6,mti,icod)
          call readcont(nin,c1,c2,l1,lct,nk6,n2,mat,mf,mt,ns)
c
c         calculate total yield
c
          mparti=mtyr(k)
          kpart=1
          kk=1
          do while (kpart.le.mparti.and.kk.le.nk6)
            call readtab1(nin,zap,awp,lip,law,nr,npy,nbt,ibt,x,y)
            izap=nint(zap)
            if (izap.eq.izai) then
              call checklaw(nr,ibt,icod)
              if (icod.ne.0) then
                write(iou,*)' ERROR: MF6/MT=',mti,' yields are not',
     &          ' linearly interpolable'
                write(iou,*)' Use SIXLIN'
                stop
              endif
              if (kpart.eq.1) then
                npy2=npy
                do i=1,npy2
                  x2(i)=x(i)
                  y2(i)=y(i)
                enddo
              else
                call union(npy,x,npy2,x2,npy3,y1,npmax)
                do i=1,npy3
                  xx=y1(i)
                  yyy=fvalin(npy,x,y,xx)
                  yy3=fvalin(npy2,x2,y2,xx)
                  b(i)=yyy+yy3
                enddo
                npy2=npy3
                do i=1,npy2
                  x2(i)=y1(i)
                  y2(i)=b(i)
                enddo
              endif
              kpart=kpart+1
            endif
            call nextsub6(nin,law,nbt,ibt,x,b)
            kk=kk+1
          enddo
c
c         Update values of TYR BLOCK
c         Energy dependent yields are also prepared
c
          if (lct.gt.2) then
            if (awi.le.4.0d0) then
              lct=2
            else
              lct=1
            endif
          endif
          iconty=isconst(npy2,y2)
          if ((mti.ge.18.and.mti.le.21).or.mti.eq.38) then
            ytot=19
          elseif (iconty.eq.0) then
            ytot=(lxs-kdlw+101)
            xss(lxs)=0
            lxs=lxs+1
            xss(lxs)=npy2
            do iy=1,npy2
              xss(lxs+iy)=x2(iy)*ev2mev
              xss(lxs+npy2+iy)=y2(iy)
            enddo
            lxs=lxs+2*npy2+1
          else
            ytot=0.5d0*(y2(1)+y2(npy2))
          endif
          xss(ktyr+k-1)=(3-2*lct)*ytot
c
c         Prepare dlw array
c
          call findmt(nin,mat0,6,mti,icod)
          call readcont(nin,c1,c2,l1,lct,nk6,n2,mat,mf,mt,ns)
          xss(kldlw+k-1)=lxs-kdlw+1
          kpart=1
          kk=1
          do while (kpart.le.mparti.and.kk.le.nk6)
            call readtab1(nin,zap,awp,lip,law,nr,npy,nbt,ibt,x,y)
            izap=nint(zap)
            if (izap.eq.izai) then
              write(iou,*)' MF6 for k=',k,' mti=',mti,' nk6=',nk6,
     &        ' qi=',qi,' part=',kpart,' lxs=',lxs,' loc=',lxs-kdlw+1
              write(*,*)' MF6 for k=',k,' mti=',mti,' nk6=',nk6,
     &        ' qi=',qi,' part=',kpart,' lxs=',lxs,' loc=',lxs-kdlw+1
              if (kpart.gt.1) xss(lxs0)=lxs-kdlw+1
              lxs0=lxs
              xss(lxs)=0
              lawxs=lxs+1
              ldatxs=lxs+2
              xss(lxs+3)=0
              lxs=lxs+4
              xss(lxs)=npy
              do iy=1,npy
                xx=x(iy)
                xss(lxs+iy)=xx*ev2mev
                ytot=fvalin(npy2,x2,y2,xx)
                if (ytot.gt.0.0d0) then
                  xss(lxs+npy+iy)=y(iy)/ytot
                else
                  xss(lxs+npy+iy)=0.0d0
                endif
              enddo
              lxs=lxs+2*npy+1
              xss(ldatxs)=lxs-kdlw+1
              if (law.eq.6) then
                xss(lawxs)=66
                call readcont(nin,apxs,c2,l1,l2,n1,npsx,mat,mf,mt,ns)
                xss(lxs)=npsx
                xss(lxs+1)=apxs
                lxs=lxs+2
              else
                call readtab2(nin,c1,c2,lang,lep,nr,ne,nbt,ibt)
                if (law.eq.1) then
                  if(lang.eq.2) then
                    lawn=44
                  elseif (lang.eq.1.or.(lang.ge.11.and.lang.le.15)) then
                    lawn=61
                    if (lang.eq.11) then
                      intmu=1
                    else
                      intmu=2
                    endif
                  else
                    write(iou,*)' MF6/LAW=1 LANG=',lang,' not allowed'
                    stop
                  endif
                  xss(lawxs)=lawn
                  write(iou,*)' LAW=',lawn
                  write(*,*)' LAW=',lawn
                else
                  write(iou,*)' ERROR: law=',law,' not coded'
                  write(iou,*)' Run SIXLIN'
                  stop
                endif
                call checklaw(nr,ibt,icod)
                if (icod.eq.0) then
                  xss(lxs)=0
                  lxs=lxs+1
                else
                  xss(lxs)=nr
                  do j=1,nr
                    xss(lxs+j)=nbt(j)
                    xss(lxs+nr+j)=min(mod(ibt(j),10),2)
                  enddo
                  lxs=lxs+2*nr+1
                endif
                xss(lxs)=ne
                lxsn=lxs+ne
                lxsd=lxsn+ne+1
                do j=1,ne
                  call readlist(nin,c1,e,nd,na,nw,nep,b)
                  if (nw.gt.nbmax) then
                    write(iou,*)' ERROR: Increase size of list arrays'
                    write(iou,*)' Set nbmax greater than ',nbmax
                    stop
                  endif
                  ee=e*ev2mev
                  xss(lxs+j)=ee
                  xss(lxsn+j)=lxsd-kdlw+1
                  na2=na+2
                  do ip=1,nep
                    ip0=na2*(ip-1)
                    x(ip)=b(ip0+1)
                    y(ip)=b(ip0+2)
                  enddo
                  if (qi.lt.0.0d0) then
                    iep=0
                    ip=nd+1
                    do while (ip.le.nep.and.iep.eq.0)
                      if (x(ip).gt.e)iep=ip
                      ip=ip+1
                    enddo
                    if (iep.gt.0) then
                      nr=1
                      nbt(1)=nep
                      ibt(1)=lep
                      xx=0.9999999d0*e
                      y(iep)=fvalue(nr,nbt,ibt,nep,x,y,xx)
                      x(iep)=xx
                      nep=iep
                    endif
                  endif
                  if ((lxs+5*nep).gt.nxss) then
                    write(iou,*)' Increase size XSS', nxss
                    stop
                  endif
                  call cdfcal1(nep,nd,x,y,lep,y1)
                  xss(lxsd)=lep+10*nd
                  lxsd=lxsd+1
                  xss(lxsd)=nep
                  nep2=nep+nep
                  nep3=nep2+nep
                  nep4=nep3+nep
                  do ip=1,nep
                    xss(lxsd+ip)=x(ip)*ev2mev
                    if (ip.gt.nd) then
                      xss(lxsd+nep+ip)=y(ip)/ev2mev
                    else
                      xss(lxsd+nep+ip)=y(ip)
                    endif
                    xss(lxsd+nep2+ip)=y1(ip)
                  enddo
                  if (lawn.eq.44) then
                    do ip=1,nep
                      if (na.eq.0) then
                        xss(lxsd+nep3+ip)=0.0d0
                        xss(lxsd+nep4+ip)=1.0d-20
                      else
                        ip0=na2*(ip-1)
                        xss(lxsd+nep3+ip)=b(ip0+3)
                        if (na.eq.2) then
                         xss(lxsd+nep4+ip)=b(ip0+4)
                        else
                         ep=x(ip)*ev2mev
                         xss(lxsd+nep4+ip)=bachaa(izai,izap,matza,ee,ep)
                        endif
                      endif
                    enddo
                    lxsd=lxsd+5*nep+1
                  elseif (lawn.eq.61) then
                    lxscd=lxsd+nep4+1
                    do ip=1,nep
                      xss(lxsd+nep3+ip)=lxscd-kdlw+1
                      ip0=na2*(ip-1)+2
                      if (na.gt.0.and.b(ip0).gt.0.0d0) then
                        if (lang.eq.1) then
                          call leg2lin(na,b(ip0),nmu,x,y,tol,ymin,npmax)
                        else
                          im0=ip0
                          nmu=na/2
                          do im=1,nmu
                            im0=im0+1
                            x(im)=b(im0)
                            im0=im0+1
                            y(im)=2.0d0*b(im0)*b(ip0)
                          enddo
                          if (lang.gt.12) then
                            nr=1
                            nbt(1)=nmu
                            ibt(1)=lang-10
                            call linear(nr,nbt,ibt,nmu,x,y,tol,ymin,
     &                         npmax)
                            do im=1,nmu
                              if (y(im).lt.0.0d0) y(im)=1.0d-30
                            enddo
                            c=b(i0)
                            ilaw=2
                            call renorm(nmu,x,y,ilaw,c,fn)
                          endif
                        endif
                        call cdfcal(nmu,x,y,intmu,y1)
                        xss(lxscd)=intmu
                        write(iou,*)' mti= ',mti,' ip= ',ip,' nmu= ',nmu
                        lxscd=lxscd+1
                        xss(lxscd)=nmu
                        nmu2=nmu+nmu
                        do im=1,nmu
                          xss(lxscd+im)=x(im)
                          xss(lxscd+nmu+im)=y(im)
                          xss(lxscd+nmu2+im)=y1(im)
                        enddo
                      else
                        intmu=2
                        xss(lxscd)=intmu
                        nmu=2
                        lxscd=lxscd+1
                        xss(lxscd)=nmu
                        xss(lxscd+1)=-1.0d0
                        xss(lxscd+2)=1.0d0
                        xss(lxscd+3)=0.5d0
                        xss(lxscd+4)=0.5d0
                        xss(lxscd+5)=0.0d0
                        xss(lxscd+6)=1.0d0
                      endif
                      lxscd=lxscd+3*nmu+1
                    enddo
                    lxsd=lxscd
                  else
                    write(icod,*)' ERROR: MF6/MT=',mti,' not coded'
                    stop
                  endif
                enddo
                kpart=kpart+1
                lxs=lxsd
              endif
            else
              call nextsub6(nin,law,nbt,ibt,x,b)
            endif
            kk=kk+1
          enddo
        else
c
c         law 3 for inelastics on mf4 or mf6/law=2-3
c
          write(*,*)' LAW 3 for k=',k,' mti=',mti,' nk=',nk,' qi=',qi,
     &    ' lxs=',lxs,' loc=',lxs-kdlw+1
          write(iou,*)' LAW 3 for k=',k,' mti=',mti,' nk=',nk,' qi=',qi,
     &    ' lxs=',lxs,' loc=',lxs-kdlw+1
          if ((lxs+11).gt.nxss) then
            write(iou,*)' ERROR: Increase size XSS',nxss
            stop
          endif
          xss(kldlw+k-1)=lxs-kdlw+1
          xss(lxs)=0
          xss(lxs+1)=3
          xss(lxs+2)=lxs-kdlw+10
          xss(lxs+3)=0
          xss(lxs+4)=2
          l=nint(xss(klsig+k-1))
          j1=nint(xss(ksig+l-1))
          j2=nint(xss(ksig+l)+j1-1)
          xx=xss(kesz+j1-1)
          xss(lxs+5)=xx
          xss(lxs+6)=xss(kesz+j2-1)
          xss(lxs+7)=1.0d0
          xss(lxs+8)=1.0d0
          yy=(awr0+awi)/awr0
          zz=dabs(qi)*yy
          if (zz.gt.xx) then
            xss(lxs+9)=0.9999999d0*xx
          else
            xss(lxs+9)=zz
          endif
          xss(lxs+10)=1.0d0/(yy*yy)
          lxs=lxs+11
        endif
      enddo
      write(iou,*)' LDLW & DLW blocks prepared'
c
c      No photon transport
c
      do jx=12,20
        jxs(jx)=0
      enddo
      nxs(6)=0
      nxs(7)=0
      nxs(8)=0
      do jx=12,15
        nxs(jx)=0
      enddo
      nxs(16)=-1
c
c      Probability tables in the URR
c
      klunr=0
      if (mt153.gt.0) then
        write(iou,*)' Processing probability tables in the URR'
        klunr=lxs
        jxs(23)=klunr
        call findmt(nin,mat0,2,153,icod)
        call readcont(nin,c1,c2,l1,l2,iintt,nbin,mat,mf,mt,ns)
        call readlist(nin,tempz,c2,lssf,icomp,nw,nunr,b)
        mtabso=icomp/1000
        mtinel=icomp-mtabso*1000
        write(*,*)' PTABLE nbin=',nbin,' intunr=',iintt,' lssf=',lssf
        write(*,*)' mtinel=',mtinel,' mtabso=',mtabso,' temp=',tempz
        write(iou,*)' PTABLE nbin=',nbin,' intunr=',iintt,' lssf=',lssf
        write(iou,*)' mtinel=',mtinel,' mtabso=',mtabso,' temp=',tempz
        if (c2.gt.0.0d0.and.lssf.eq.0) then
          lssf0=1
          write(*,*)' PTABLE changed to self-shielding factors(LSSF=1)'
          write(iou,*)' PTABLE changed to self-shielding factors',
     &                '(LSSF=1)'
        else
          lssf0=0
        endif
        if (mtinel.ne.0) then
          ilf=mtinel
        else
          ilf=-1
        endif
        if (mtabso.ne.0) then
          ioa=mtabso
        else
          ioa=-1
        endif
        if (iintt.ne.2) then
          iintt=2
          write(*,*)' Interpolation law (intunr) changed to 2'
          write(iou,*)' Interpolation law (intunr) changed to 2'
        endif
        nbin2=nbin+nbin
        nbin3=nbin2+nbin
        nbin4=nbin3+nbin
        nbin5=nbin4+nbin
        nbin6=nbin5+nbin
        nrec=1+nbin6
        xx1=b(1)
        xx2=b(1+nrec)
        xx=edelta(xx1,1.0d0)
        if (xx.ge.xx2) xx=0.5d0*(xx1+xx2)
        b(1)=xx
        xss(lxs)=nunr
        xss(lxs+1)=nbin
        xss(lxs+2)=iintt
        xss(lxs+3)=ilf
        xss(lxs+4)=ioa
        if (lssf0.eq.1) then
          xss(lxs+5)=lssf0
        else
          xss(lxs+5)=lssf
        endif
        lxs=lxs+6
        do ie=1,nunr
          ie1=ie-1
          je=1+ie1*nrec
          xss(lxs+ie1)=b(je)*ev2mev
          lx=lxs+nunr-1+ie1*nbin6
          if (lssf0.eq.1) then
            do k=1,5
              x(k)=0.0d0
            enddo
            do ibin=1,nbin
              prob=b(je+ibin)
              x(1)=b(je+nbin+ibin)*prob+x(1)
              x(2)=b(je+nbin2+ibin)*prob+x(2)
              x(3)=b(je+nbin3+ibin)*prob+x(3)
              x(4)=b(je+nbin4+ibin)*prob+x(4)
              x(5)=b(je+nbin5+ibin)*prob+x(5)
            enddo
          endif
          sum=0.0d0
          do ibin=1,nbin
            sum=sum+b(je+ibin)
            xss(lx+ibin)=sum
            if (lssf0.eq.1) then
              if (x(1).ne.0.0d0) then
                xss(lx+nbin+ibin)=b(je+nbin+ibin)/x(1)
              else
                xss(lx+nbin+ibin)=1.0d0
              endif
              if (x(2).ne.0.0d0) then
                xss(lx+nbin2+ibin)=b(je+nbin2+ibin)/x(2)
              else
                xss(lx+nbin2+ibin)=1.0d0
              endif
              if (x(3).ne.0.0d0) then
                xss(lx+nbin3+ibin)=b(je+nbin3+ibin)/x(3)
              else
                xss(lx+nbin3+ibin)=1.0d0
              endif
              if (x(4).ne.0.0d0) then
                xss(lx+nbin4+ibin)=b(je+nbin4+ibin)/x(4)
              else
                xss(lx+nbin4+ibin)=1.0d0
              endif
              if (x(5).ne.0.0d0) then
                xss(lx+nbin5+ibin)=b(je+nbin5+ibin)/x(5)
              else
                xss(lx+nbin5+ibin)=1.0d0
              endif
            else
              xss(lx+nbin+ibin)=b(je+nbin+ibin)
              xss(lx+nbin2+ibin)=b(je+nbin2+ibin)
              xss(lx+nbin3+ibin)=b(je+nbin3+ibin)
              xss(lx+nbin4+ibin)=b(je+nbin4+ibin)
              if (lssf.eq.1) then
                xss(lx+nbin5+ibin)=b(je+nbin5+ibin)
              else
                xss(lx+nbin5+ibin)=b(je+nbin5+ibin)*ev2mev
              endif
            endif
          enddo
          xss(lx+nbin)=1.0d0
        enddo
        lxs=lxs+nunr*nrec
      endif
c
c      Delayed neutron blocks
c
      knud=0
      kdndat=0
      kldnd=0
      kdnd=0
      m5=iposm(nmf5,mf5,455)
      if (mt455.gt.0.and.m5.gt.0) then
c
c       delayed nu
c
        write(iou,*)' Processing delayed neutron blocks'
        knud=lxs
        jxs(24)=knud
        call findmt(nin,mat0,1,455,icod)
        call readcont(nin,c1,c2,ldg,lnu,n1,n2,mat,mf,mt,ns)
        if (ldg.ne.0.or.lnu.ne.2) then
          write(iou,*)' ERROR: LNU=',lnu,' LDG=',ldg,' for MT455'
          write(iou,*)'        LNU=2 and LDG=0 are expected'
          stop
        endif
        call readlist(nin,c1,c2,l1,l2,nnf,n2,y1)
        call readtab1(nin,c1,c2,l1,l2,nr,ne,nbt,ibt,x,y)
        call checklaw(nr,ibt,icod)
        if (icod.ne.0) then
          write(iou,*)' ERROR: MT455 linear interpolation expected'
          stop
        endif
        write(*,*)' Delayed nu-bar lnu=',lnu,' ldg=',ldg,' nnf=',nnf
        write(iou,*)' Delayed nu-bar lnu=',lnu,' ldg=',ldg,' nnf=',nnf
        xss(lxs)=lnu
        lxs=lxs+1
        xss(lxs)=0
        lxs=lxs+1
        xss(lxs)=ne
        do i=1,ne
          xss(lxs+i)=x(i)*ev2mev
          xss(lxs+ne+i)=y(i)
        enddo
        lxs=lxs+2*ne+1
c
c       delayed basic data (dndat)
c
        call findmt(nin,mat0,5,455,icod)
        call readcont(nin,c1,c2,l1,l2,nk,n2,mat,mf,mt,ns)
        if (nk.ne.nnf) then
          write(iou,*)' ERROR: Number of delayed groups not match'
          write(iou,*)' MF1/NNF=',nnf,' MF5/NK=',nk
          stop
        endif
        nxs(8)=nk
        kdndat=lxs
        jxs(25)=kdndat
        do k=1,nk
          call readtab1(nin,c1,c2,l1,law,nr,ne,nbt,ibt,x,y)
          if (law.ne.1) then
            write(iou,*)' LF=1 expected in MF5 ',' LF=',law
            stop
          endif
          call checklaw(nr,ibt,icod)
          if (icod.ne.0) then
            write(iou,*)' ERROR: MF5/MT455 linear interpolation expect'
            stop
          endif
          write(*,*)' Delayed basic data for group ',k,' decay=',y1(k)
          write(iou,*)' Delayed basic data for group ',k,' decay=',y1(k)
          xss(lxs)=y1(k)*s2shak
          lxs=lxs+1
          xss(lxs)=0
          lxs=lxs+1
          xss(lxs)=ne
          do i=1,ne
            xss(lxs+i)=x(i)*ev2mev
            xss(lxs+ne+i)=y(i)
          enddo
          lxs=lxs+2*ne+1
          call readtab2(nin,c1,z,l1,l2,nr,nz,nbt,ibt)
          do iz=1,nz
            call readtab1(nin,c1,c2,l1,l2,n1,n2,nbt,ibt,x,y)
          enddo
        enddo

c
c       kldnd and kdnd blocks
c
        kldnd=lxs
        jxs(26)=kldnd
        lxs=kldnd+nk
        kdnd=lxs
        jxs(27)=kdnd
        call findmt(nin,mat0,5,455,icod)
        call readcont(nin,c1,c2,l1,l2,nk,n2,mat,mf,mt,ns)
        do k=1,nk
          call readtab1(nin,c1,c2,l1,law,nr,npp,nbt,ibt,x,y)
          write(*,*)' Spectrum for group ',k,' lf=',law,' npp=',npp
          write(iou,*)' Spectrum for group ',k,' lf=',law,' npp=',npp
          xmin=x(1)*ev2mev
          xmax=x(npp)*ev2mev
          xss(kldnd+k-1)=lxs-kdnd+1
          xss(lxs)=0
          xss(lxs+1)=4
          xss(lxs+2)=lxs-kdnd+10
          xss(lxs+3)=0
          xss(lxs+4)=2
          xss(lxs+5)=xmin
          xss(lxs+6)=xmax
          xss(lxs+7)=1.0d0
          xss(lxs+8)=1.0d0
          call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
          xss(lxs+9)=0
          xss(lxs+10)=ne
          lxs=lxs+11
          ldat=lxs
          lxs=lxs+2*ne
          do ie=1,ne
            call readtab1(nin,c1,ep,l1,l2,nr,nep,nbt,ibt,x,y)
            xss(ldat+ie-1)=ep*ev2mev
            xss(ldat+ne+ie-1)=lxs-kdnd+1
            lep=ibt(1)
            if (lep.gt.2) lep=2
            call cdfcal(nep,x,y,lep,y2)
            xss(lxs)=lep
            lxs=lxs+1
            xss(lxs)=nep
            nep2=nep+nep
            do i=1,nep
              xss(lxs+i)=x(i)*ev2mev
              xss(lxs+nep+i)=y(i)/ev2mev
              xss(lxs+nep2+i)=y2(i)
            enddo
            lxs=lxs+1+3*nep
          enddo
        enddo
        write(iou,*)' Delayed neutrons blocks prepared'
      endif
      kend=lxs-1
      jxs(22)=kend
      nxs(1)=kend
      write(*,*)' end=',kend
      write(iou,*)' end=',kend
c
c      write the ace file
c
      close(nin)
      open (nou,file=fout)
      call change(nou,mcnpx)
      close(nou)
c
c      write *.xsd file for xsdir
c
      i=index(fout,'.',.true.)
      if (i.le.0) then
        i=len(trim(fout))
      else
        i=i-1
      endif
      write(fin1,'(a,a)')trim(fout(1:i)),'.xsd'
      fin1=trim(fin1)
      open(nou,file=fin1)
      if (klunr.gt.0) then
        if (mcnpx.eq.1) then
          write(nou,'(a13,f12.6,a,a,a7,i9,a4,1pe11.4,a,/,a12)')hz(1:13),
     &    awr0,' ',trim(fout),' 0 1 1 ',kend,' 0 0',tz,'  +','ptable'
        else
          write(nou,'(a10,f12.6,a,a,a7,i9,a4,1pe11.4,a7)')hz(1:10),
     &    awr0,' ',trim(fout),' 0 1 1 ',kend,' 0 0',tz,' ptable'
        endif
      else
        if (mcnpx.eq.1) then
          write(nou,'(a13,f12.6,a,a,a7,i9,a4,1pe11.4)')hz(1:13),
     &    awr0,' ',trim(fout),' 0 1 1 ',kend,' 0 0',tz
        else
          write(nou,'(a10,f12.6,a,a,a7,i9,a4,1pe11.4)')hz(1:10),
     &    awr0,' ',trim(fout),' 0 1 1 ',kend,' 0 0',tz
        endif
      endif
      close(nou)
      write(*,*)' DOACE ENDED'
      write(iou,*)' DOACE ENDED'
      close(iou)
      stop
      end
C======================================================================
C     Processing routines
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
      subroutine checkmon(nerr,n,x,icod)
c
c      Check if x is a non decreasing function
c
      implicit real*8 (a-h,o-z)
      dimension x(*)
      icod=0
      i=2
      do while (i.le.n)
        if (x(i-1).ge.x(i)) then
          write(*,*)' Non-decreasing value = ',x(i)
          if (nerr.gt.0) then
            write(nerr,*)' Non-decreasing value = ',x(i)
          endif
          icod=icod+1
        endif
        i=i+1
      enddo
      return
      end
C======================================================================
      subroutine checkeq(n1,x1,n2,x2,icod)
c
c     check if real arrays x1=x2
c
      implicit real*8 (a-h,o-z)
      dimension x1(*),x2(*)
      data eps/1.0d-10/
      icod=0
      if (n1.ne.n2) then
        icod=1
        return
      endif
      i=1
      do while (icod.eq.0.and.i.le.n1)
        xi=x1(i)
        if (dabs(xi-x2(i)).gt.dabs(eps*xi)) icod=1
        i=i+1
      enddo
      return
      end
C======================================================================
      subroutine checkdis(nerr,n,x,y,icod)
c
c     remove duplicate points & discontinuities from the
c     common energy grid
c
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      character*11 x0,x1
      icod=0
      call ff2chx(x(1),x0)
      read(x0,'(e11.0)')x(1)
      y0=y(1)
      j=1
      i=2
      do while (i.le.n)
        call ff2chx(x(i),x1)
        y1=y(i)
        if (x0.eq.x1.and.((y0.eq.y1).or.(i.eq.n.and.y1.eq.0.0d0))) then
          write(*,*)' Warning: duplicate point at x=',x1,' y=',y1,
     &      ' i=',i,' removed'
          if (nerr.gt.0) then
            write(nerr,*)' Warning: duplicate point at x=',x1,' y=',y1,
     &        ' i=',i,' removed'
          endif
          icod=icod+1
        else
          j=j+1
          read(x1,'(e11.0)')x(j)
          y(j)=y1
        endif
        x0=x1
        y0=y1
        i=i+1
      enddo
      n=j
      call ff2chx(x(1),x0)
      i=2
      j=1
      do while (i.le.n)
        call ff2chx(x(i),x1)
        if (x0.eq.x1) then
          if (i.lt.n) then
            do while(x0.eq.x1.and.i.lt.n)
              icod=icod+1
              i=i+1
              call ff2chx(x(i),x1)
            enddo
          endif
          if (i.lt.n) then
            i=i-1
            xx=edelta(x(i),1.0d0)
            if (xx.ge.x(i+1)) xx=0.5d0*(x(i)+x(i+1))
            call ff2chx(xx,x1)
            j=j+1
            x(j)=xx
            y(j)=y(i)
            write(*,*)' Warning: discontinuity at x=',x0,' i=',i,
     &        ' x(i) set to ',x1
            if (nerr.gt.0) then
              write(nerr,*)' Warning: discontinuity at x=',x0,' i=',i,
     &          ' x(i) set to ',x1
            endif
          else
            icod=icod+1
            write(*,*)' Warning: discontinuity at x=',x0,' i=',i,
     &        ' x(i) removed'
            if (nerr.gt.0) then
              write(nerr,*)' Warning: discontinuity at x=',x0,' i=',i,
     &          ' x(i) removed'
            endif
          endif
        else
          j=j+1
          x(j)=x(i)
          y(j)=y(i)
        endif
        x0=x1
        i=i+1
      enddo
      n=j
      return
      end
C======================================================================
      subroutine ff2chx(xx,strxx)
      implicit real*8 (a-h,o-z)
c
c     pack value into 11-character string (9 digit precision)
c
      character*11 strxx
      character*12 str12
      character*13 str13
      character*14 str14
      character*16 str16
      write(str16,'(1pe16.9)')xx
      read(str16,'(e16.0)')ff
      aff=dabs(ff)
      if (aff.lt.1.00000d-99) then
        strxx=' 0.0       '
      elseif (aff.lt.9.99999d+99) then
        read(str16,'(13x,i3)')iex
        select case (iex)
          case (-9,-8,-7,-6,-5,-4)
            if (ff.gt.0.0d0) then
              write(str14,'(1pe14.7)')ff
              strxx(1:9)=str14(2:10)
              strxx(10:10)=str14(12:12)
              strxx(11:11)=str14(14:14)
            else
              write(str13,'(1pe13.6)')ff
              strxx(1:9)=str13(1:9)
              strxx(10:10)=str13(11:11)
              strxx(11:11)=str13(13:13)
            endif
          case (-3,-2,-1)
            if (ff.gt.0.0d0) then
              write(str13,'(f13.10)')ff
              strxx(1:11)=str13(3:13)
            else
              write(str12,'(f12.9)')ff
              strxx(1:1)=str12(1:1)
              strxx(2:11)=str12(3:12)
            endif
          case (0)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.9)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.8)')ff
            endif
          case (1)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.8)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.7)')ff
            endif
          case (2)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.7)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.6)')ff
            endif
          case (3)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.6)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.5)')ff
            endif
          case (4)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.5)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.4)')ff
            endif
          case (5)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.4)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.3)')ff
            endif
          case (6)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.3)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.2)')ff
            endif
          case (7)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.2)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.1)')ff
            endif
          case (8)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.1)')ff
              strxx(1:11)=str12(2:12)
            else
              write(strxx,'(f11.0)')ff
            endif
          case (9)
            if (ff.gt.0.0d0) then
              write(str12,'(f12.0)')ff
              strxx(1:11)=str12(2:12)
            else
              write(str13,'(1pe13.6)')ff
              strxx(1:9)=str13(1:9)
              strxx(10:10)=str13(11:11)
              strxx(11:11)=str13(13:13)
            endif
          case default
            if (ff.gt.0.0d0) then
              write(str13,'(1pe13.6)')ff
              strxx(1:8)=str13(2:9)
              strxx(9:11)=str13(11:13)
            else
              write(str12,'(1pe12.5)')ff
              strxx(1:8)=str12(1:8)
              strxx(9:11)=str12(10:12)
            endif
        end select
      else
        if (ff.gt.0.0d0) then
          strxx='9.999999+99'
        else
          strxx='-9.99999+99'
        endif
      endif
      return
      end
C======================================================================
      subroutine chendf(ffin,str11)
      implicit real*8 (a-h,o-z)
c
c     pack value into 11-character string (8 digit precision)
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
      real*8 function edelta(ffin,fdig)
      implicit real*8 (a-h,o-z)
      character*16 str16
      afdig=dabs(fdig)
      if (afdig.eq.0.0d0) then
        edelta=ffin
      else
        if (afdig.gt.9.0d0) fdig=fdig/afdig*9.0d0
        write(str16,'(1pe16.9)')ffin
        read(str16,'(e16.0)')ff
        write(str16,'(1pe16.9)')ff
        read(str16,'(13x,i3)')iex
        select case (iex)
          case (-9,-8,-7,-6,-5,-4)
            if (ff.gt.0.0d0) then
              edelta=ff+(10.0d0**iex)*1.0d-7*fdig
            else
              edelta=ff+(10.0d0**iex)*1.0d-6*fdig
            endif
          case (-3,-2,-1)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-10*fdig
            else
              edelta=ff+1.0d-9*fdig
            endif
          case (0)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-9*fdig
            else
              edelta=ff+1.0d-8*fdig
            endif
          case (1)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-8*fdig
            else
              edelta=ff+1.0d-7*fdig
            endif
          case (2)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-7*fdig
            else
              edelta=ff+1.0d-6*fdig
            endif
          case (3)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-6*fdig
            else
              edelta=ff+1.0d-5*fdig
            endif
          case (4)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-5*fdig
            else
              edelta=ff+1.0d-4*fdig
            endif
          case (5)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-4*fdig
            else
              edelta=ff+1.0d-3*fdig
            endif
          case (6)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-3*fdig
            else
              edelta=ff+1.0d-2*fdig
            endif
          case (7)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-2*fdig
            else
              edelta=ff+1.0d-1*fdig
            endif
          case (8)
            if (ff.gt.0.0d0) then
              edelta=ff+1.0d-1*fdig
            else
              edelta=ff+fdig
            endif
          case (9)
            if (ff.gt.0.0d0) then
              edelta=ff+fdig
            else
              edelta=ff+(10.0d0**iex)*1.0d-6*fdig
            endif
          case default
            if (ff.gt.0.0d0) then
              edelta=ff+(10.0d0**iex)*1.0d-6*fdig
            else
              edelta=ff+(10.0d0**iex)*1.0d-5*fdig
            endif
        end select
      endif
      return
      end
C======================================================================
      function iposm(n,m,k)
c
c     find position of element k in the array m
c
      dimension m(*)
      ipos=0
      i=1
      do while (i.le.n.and.ipos.eq.0)
        if (m(i).eq.k) ipos=i
        i=i+1
      enddo
      iposm=ipos
      return
      end
C======================================================================
      function iposx(n,x,x0)
c
c     find the first element of array x greather x0
c
      implicit real*8 (a-h,o-z)
      dimension x(*)
      data eps/1.0d-10/
      epsx=eps*x0
      ipos=0
      i=n
      do while (i.ge.1.and.ipos.eq.0)
        xi=x(i)
        if (xi.gt.x0) then
          i=i-1
        elseif (dabs(xi-x0).lt.epsx) then
          ipos=i
        else
          ipos=i+1
        endif
      enddo
      if (ipos.eq.0) then
        iposx=1
      else
        iposx=ipos
      endif
      return
      end
C======================================================================
      function mtchkd(mtsum,mtlow,mtup,nmf3,mf3,nmf4,mf4,nmf6,mf6)
      dimension mf3(*),mf4(*),mf6(*)
      mts=0
      m3s=iposm(nmf3,mf3,mtsum)
      m4=iposm(nmf4,mf4,mtsum)
      m6=iposm(nmf6,mf6,mtsum)
      m46=m4+m6
      if (m3s.gt.0.and.m46.gt.0) mts=1
      mtl=0
      m3l=0
      do mt=mtlow,mtup
        m3=iposm(nmf3,mf3,mt)
        if (m3.gt.0) m3l=m3l+1
        m4=iposm(nmf4,mf4,mt)
        m6=iposm(nmf6,mf6,mt)
        m46=m4+m6
        if (m3.gt.0.and.m46.gt.0) mtl=mtl+1
      enddo
      if (mts.gt.0) then
        mtchkd=1
      elseif (mtl.gt.0) then
        mtchkd=0
      elseif (m3s.gt.0) then
         mtchkd=1
      else
         mtchkd=0
      endif
      return
      end
C======================================================================
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
      real*8 function fvalin(np,x,y,x0)
c
c      Return f(x0), value of function f evaluated at x0.
c      Function f is given by np linearly interpolable tabulated points.
c
        implicit real*8 (a-h, o-z)
        dimension x(*),y(*)
        character*11 cx0,cx1
        if (x0.lt.x(1).or.x0.gt.x(np)) then
          fvalin=0.0d0
        else
          call chendf(x0,cx0)
          i=2
          do while (i.le.np.and.x(i).lt.x0)
            i=i+1
          enddo
          i1=i-1
          call chendf(x(i1),cx1)
          if (cx0.eq.cx1) then
            fvalin=y(i1)
          elseif (x0.eq.x(i)) then
            fvalin=y(i)
          else
            call terp1m(x(i1),y(i1),x(i),y(i),x0,y0,2)
            fvalin=y0
          endif
        endif
        return
      end
C======================================================================
      function isconst(n,y)
c
c     check if array y is constant,integer or variable
c      if isconst=0, y array is not constant
c      if isconst=1, y array is constant
c      if isconst=2, y array is constant and integer
c
      implicit real*8 (a-h,o-z)
      dimension y(*)
      data epc/1.0d-10/,epi/1.0d-4/
      iconst=1
      iinteg=1
      yi=y(1)
      if (dabs(yi-nint(yi)).gt.epi) iinteg=0
      i=2
      do while (i.le.n.and.iconst.eq.1)
        yi=y(i)
        if (dabs(yi-y(i-1)).gt.epc) iconst=0
        if (dabs(yi-nint(yi)).gt.epi) iinteg=0
        i=i+1
      enddo
      if (iconst.eq.0) then
        isconst=0
      else
        if (iinteg.eq.1) then
          isconst=2
        else
          isconst=1
        endif
      endif
      return
      end
C======================================================================
      subroutine thresh(iou,mt,qi,awr,awi,nr,ne,nbt,ibt,x,y)
      implicit real*8 (a-h,o-z)
      dimension nbt(*),ibt(*),x(*),y(*)
      data eps/1.0d-10/
      if (qi.lt.0.0d0.and.mt.ne.5) then
        ethr=-qi*(awr+awi)/awr
        if (ethr.gt.x(1)) then
          write(iou,*)'  Threshold inconsistency found for MT=',mt
          write(iou,*)'  Threshold=',ethr,' E(1)=',x(1)
          epsx=eps*ethr
          ipos=0
          i=1
          do while (ipos.eq.0.and.i.le.ne)
            if ((ethr-x(i)).gt.epsx) then
              i=i+1
            else
              ipos=i
            endif
          enddo
          if (ipos.eq.0) then
            write(iou,*)'  MT=',mt,' Threshold=',ethr,' > E(ne)=',x(ne)
            write(iou,*)'  Cross section set to zero at all energies'
            nr=1
            nbt(1)=2
            ibt(1)=2
            xx=edelta(ethr,1.0d0)
            ne=2
            x(1)=xx
            x(2)=1.00001*xx
            y(1)=0.0d0
            y(2)=0.0d0
          else
            if((x(ipos)-ethr).gt.epsx.and.ipos.gt.1) ipos=ipos-1
            x(ipos)=edelta(ethr,1.0d0)
            y(ipos)=0.0d0
            ne=ne-ipos+1
            do i=1,ne
              x(i)=x(ipos+i-1)
              y(i)=y(ipos+i-1)
            enddo
            j=1
            do while (nbt(j).lt.ipos)
              j=j+1
            enddo
            nr=nr-j+1
            do i=1,nr
              nbt(i)=nbt(j+i-1)-ipos+1
              ibt(i)=ibt(j+i-1)
            enddo
          endif
          write(iou,*)'  E(1)corrected to ',x(1),' with ',y(1),' XS'
        endif
      endif
      return
      end
C======================================================================
      subroutine typr(nin,izai,mat0,mt,nmf4,mf4,nmf6,mf6,tyr,ntyr,
     &                nbt,ibt,x,y,b)
      implicit real*8 (a-h,o-z)
      dimension mf4(*),mf6(*),nbt(*),ibt(*)
      dimension x(*),y(*),b(*)
      data ityr/100/,epi/1.0d-4/,epc/1.0d-10/
      save ityr
      m4=0
      if (izai.eq.1) then
        m4=iposm(nmf4,mf4,mt)
        if (m4.gt.0) then
          call findmt(nin,mat0,4,mt,icod)
          call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
          call readcont(nin,c1,c2,l1,lct4,n1,n2,mat,mf,mt,ns)
        endif
      endif
      ntyr6=0
      m6=iposm(nmf6,mf6,mt)
      if (m6.gt.0) then
        call findmt(nin,mat0,6,mt,icod)
        call readcont(nin,c1,c2,l1,lct6,nk,n2,mat,mf,mt,ns)
        iconst=1
        iinteg=1
        yy=0
        do k=1,nk
          call readtab1(nin,zap,awp,l1,law,nr,ne,nbt,ibt,x,y)
          izap=nint(zap)
          if (izap.eq.izai) then
            ntyr6=ntyr6+1
            ii=isconst(ne,y)
            if (ii.eq.0) then
              iconst=0
              iinteg=0
            elseif (ii.eq.1) then
              iinteg=0
            endif
            if (iconst.eq.1) yy=yy+0.5d0*(y(1)+y(ne))
          endif
          call nextsub6(nin,law,nbt,ibt,x,b)
        enddo
      endif
      if (m4.gt.0.and.ntyr6.eq.0) then
        n1=numpart(izai,mt)
        tyr=(3-2*lct4)*n1
        ntyr=1
      elseif (ntyr6.gt.0) then
        if (m4.gt.0) then
          write(*,*)' MF4 and MF6 found for MT=',mt,' and IZAI=',izai
          write(*,*)' MF4 data ignored'
          mf4(m4)=0
        endif
        n1=numpart(izai,mt)
        y1=n1
        if (n1.eq.19) then
          yy=y1
        elseif (iconst.eq.1) then
          if (iinteg.eq.1.and.n1.lt.101) then
            if (dabs(yy-y1).gt.epi) then
              write(*,*)' Yields on MF6 not equal to ',n1,' for MT=',mt
            endif
          endif
        else
          n1=ityr+1
          yy=n1
        endif
        if (izai.le.2004.and.lct6.eq.3) lct6=2
        tyr=(3-2*lct6)*yy
        ntyr=ntyr6
      elseif (m4.eq.0.and.ntyr6.eq.0) then
        write(*,*)' ERROR: Angular distribution not found for MT=',mt
        close(nin)
        stop
      endif
      return
      end
C======================================================================
      function mtvalid(izai,mt)
      if (izai.eq.1) then
        select case (mt)
          case (1:4,5,11,16:18,19:21,38,22:25,28:30,32:37,41:42,44:45,
     &          51:91,102:109,111:117,152:200,201:207,301,444,452,455,
     &          456,600:649,650:699,700:749,750:799,800:849,875:891)
            mtvalid=1
          case default
            mtvalid=0
        end select
      endif
      return
      end
C======================================================================
      function mtrval(izai,mt)
      if (izai.eq.1) then
        select case (mt)
          case (5,11,16:18,19:21,38,22:25,28:30,32:37,41:42,44:45,51:91,
     &          102:109,111:117,152:200,201:207,444,600:649,650:699,
     &          700:749,750:799,800:849,875:891)
            mtrval=1
          case default
            mtrval=0
        end select
      endif
      return
      end
C======================================================================
      function mtprodi(izai,mt)
      if (izai.eq.1) then
        select case (mt)
          case (5,11,16:18,19:21,38,22:25,28:30,32:37,41:42,44:45,51:91,
     &        152:154,156:181,183:190,194:196,198:200,875:891)
            mtprodi=1
          case default
            mtprodi=0
        end select
      endif
      return
      end
C======================================================================
      function mtdisa(izai,mt)
      if (izai.eq.1) then
        select case (mt)
          case (102:109,111:117,155,182,191:193,197,600:649,650:699,
     &          700:749,750:799,800:849)
            mtdisa=1
          case default
            mtdisa=0
        end select
      endif
      return
      end
C======================================================================
      function numpart(izai,mt)
      if (izai.eq.1) then
        select case (mt)
          case (5)
            numpart=101
          case (18,19:21,38)
            numpart=19
          case (2,22,23,28,29,32:36,44,45,51:91,158,183:189,198)
            numpart=1
          case (11,16,24,30,41,154,159,176,190,875:891)
            numpart=2
          case (17,25,42,157,172,177,179:181,199)
            numpart=3
          case (37,156,165,169,173,178,194:196)
            numpart=4
          case (152,162,166,170,174,200)
            numpart=5
          case (153,163,167,171,175)
            numpart=6
          case (160,164,168)
            numpart=7
          case (161)
            numpart=8
          case default
            numpart=0
        end select
      endif
      return
      end
C======================================================================
      subroutine cdfcal(np,x,y,ilaw,cdf)
c
c     calculate cdf normalized to 1.0
c
c       lep=ilaw-(ilaw/10)*10
c              y        x
c       lep=1 constant
c       lep=2 lin      lin
c       lep=3 lin      log
c       lep=4 log      lin
c       lep=5 log      log
c
      implicit real*8 (a-h, o-z)
      dimension x(*),y(*),cdf(*)
      parameter (eps=1.0d-6, zero=0.0d0)
      ll=ilaw/10
      lep=ilaw-ll*10
      sum=0.0d0
      cdf(1)=sum
      if (lep.eq.1) then
        do j=2,np
          j1=j-1
          sum=sum+y(j1)*(x(j)-x(j1))
          cdf(j)=sum
        enddo
      elseif (lep.eq.2) then
        do j=2,np
          j1=j-1
          sum=sum+0.5d0*(y(j)+y(j1))*(x(j)-x(j1))
          cdf(j)=sum
        enddo
      elseif (lep.eq.4) then
        do j=2,np
          j1=j-1
          xr=x(j)-x(j1)
          if (xr.gt.zero) then
            y1=y(j1)
            y2=y(j)
            yr=y2-y1
            ay1=abs(y1)
            ayr=abs(yr)
            if (y1*y2.gt.zero.and.ayr.gt.eps*ay1) then
              sum=sum+xr*yr/log(y2/y1)
            else
              sum=sum+0.5d0*(y1+y2)*xr
            endif
          endif
          cdf(j)=sum
        enddo
      elseif (lep.eq.3) then
        do j=2,np
          j1=j-1
          x1=x(j1)
          x2=x(j)
          xr=x2-x1
          if (xr.gt.zero) then
            y1=y(j1)
            y2=y(j)
            ax1=abs(x1)
            if (x1*x2.gt.zero.and.xr.gt.eps*ax1) then
              sum=sum+x2*y2-x1*y1-xr*(y2-y1)/log(x2/x1)
            else
              sum=sum+0.5d0*(y1+y2)*xr
            endif
          endif
          cdf(j)=sum
        enddo
      elseif (lep.eq.5) then
        do j=2,np
          j1=j-1
          x1=x(j1)
          x2=y(j)
          xr=x2-x1
          if (xr.gt.zero) then
            y1=y(j1)
            y2=y(j)
            yr=y2-y1
            ax1=abs(x1)
            ay1=abs(y1)
            ayr=abs(yr)
            x12=x1*x2
            y12=y1*y2
            if (x12.gt.zero.and.xr.gt.eps*ax1.and.
     &          y12.gt.zero.and.ayr.gt.eps*ay1) then
              d=log(y2/y1)/log(x2/x1)+1.0d0
              sum=sum+y1*x1*(((x2/x1)**d)-1.0d0)/d
            elseif (y12.gt.zero.and.ayr.gt.eps*ay1) then
              sum=sum+xr*yr/log(y2/y1)
            elseif (x12.gt.zero.and.xr.gt.eps*ax1) then
              sum=sum+x2*y2-x1*y1-xr*yr/log(x2/x1)
            else
              sum=sum+0.5d0*(y1+y2)*xr
            endif
          endif
          cdf(j)=sum
        enddo
      else
        write(*,*)' ERROR in cdfcal: ilaw=',ilaw,' not coded'
        stop
      endif
      if (sum.ne.0.0d0) then
        if (sum.ne.1.0d0) then
          fn=1.0d0/sum
          do j=1,np
            y(j)=y(j)*fn
            cdf(j)=cdf(j)*fn
          enddo
        endif
      endif
      return
      end
C======================================================================
      subroutine cdfcal1(np,nd,x,y,ilaw,cdf)
c
c     calculate cdf normalized to 1.0 with nd discrete contributions
c
c       lep=ilaw-(ilaw/10)*10
c              y        x
c       lep=1 constant
c       lep=2 lin      lin
c       lep=4 log      lin
c
      implicit real*8 (a-h, o-z)
      dimension x(*),y(*),cdf(*)
      parameter (eps=1.0d-6, zero=0.0d0)
      ll=ilaw/10
      lep=ilaw-ll*10
      sum=0.0d0
      if (nd.gt.0) then
        do j=1,nd
          sum=sum+y(j)
          cdf(j)=sum
        enddo
      endif
      if (np.gt.nd) then
        j0=nd+1
        cdf(j0)=sum
        j0=j0+1
        if (lep.eq.1) then
          do j=j0,np
            j1=j-1
            sum=sum+y(j1)*(x(j)-x(j1))
            cdf(j)=sum
          enddo
        elseif (lep.eq.2) then
          do j=j0,np
            j1=j-1
            sum=sum+0.5d0*(y(j)+y(j1))*(x(j)-x(j1))
            cdf(j)=sum
          enddo
        elseif (lep.eq.4) then
          do j=j0,np
            j1=j-1
            xr=x(j)-x(j1)
            if (xr.gt.zero) then
              y1=y(j1)
              y2=y(j)
              yr=y2-y1
              ay1=abs(y1)
              ayr=abs(yr)
              if (y1*y2.gt.zero.and.ayr.gt.eps*ay1) then
                sum=sum+xr*yr/log(y2/y1)
              else
                sum=sum+0.5d0*(y1+y2)*xr
              endif
            endif
            cdf(j)=sum
          enddo
        elseif (lep.eq.3) then
          do j=j0,np
            j1=j-1
            x1=x(j1)
            x2=x(j)
            xr=x2-x1
            if (xr.gt.zero) then
              y1=y(j1)
              y2=y(j)
              ax1=abs(x1)
              if (x1*x2.gt.zero.and.xr.gt.eps*ax1) then
                sum=sum+x2*y2-x1*y1-xr*(y2-y1)/log(x2/x1)
              else
                sum=sum+0.5d0*(y1+y2)*xr
              endif
            endif
            cdf(j)=sum
          enddo
        elseif (lep.eq.5) then
          do j=j0,np
            j1=j-1
            x1=x(j1)
            x2=y(j)
            xr=x2-x1
            if (xr.gt.zero) then
              y1=y(j1)
              y2=y(j)
              yr=y2-y1
              ax1=abs(x1)
              ay1=abs(y1)
              ayr=abs(yr)
              x12=x1*x2
              y12=y1*y2
              if (x12.gt.zero.and.xr.gt.eps*ax1.and.
     &            y12.gt.zero.and.ayr.gt.eps*ay1) then
                d=log(y2/y1)/log(x2/x1)+1.0d0
                sum=sum+y1*x1*(((x2/x1)**d)-1.0d0)/d
              elseif (y12.gt.zero.and.ayr.gt.eps*ay1) then
                sum=sum+xr*yr/log(y2/y1)
              elseif (x12.gt.zero.and.xr.gt.eps*ax1) then
                sum=sum+x2*y2-x1*y1-xr*yr/log(x2/x1)
              else
                sum=sum+0.5d0*(y1+y2)*xr
              endif
            endif
            cdf(j)=sum
          enddo
        else
          write(*,*)' ERROR in cdfcal1: ilaw=',ilaw,' not coded'
          stop
        endif
      endif
      if (sum.ne.0.0d0) then
        if (sum.ne.1.0d0) then
          fn=1.0d0/sum
          do j=1,np
            y(j)=y(j)*fn
            cdf(j)=cdf(j)*fn
          enddo
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
      subroutine renorm(np,x,y,ilaw,c,fn)
c
c     renormalize tabulated data with interpolation law = ilaw to a
c     constant c. if c=0.0d0, then return the value of integral in c
c
c       lep=ilaw-(ilaw/10)*10
c              y        x
c       lep=1 constant
c       lep=2 lin      lin
c       lep=3 lin      log
c       lep=4 log      lin
c       lep=5 log      log
c
      implicit real*8 (a-h, o-z)
      parameter (eps=1.0d-6, zero=0.0d0)
      dimension x(*),y(*)
      ll=ilaw/10
      lep=ilaw-ll*10
      sum=0.0d0
      if (lep.eq.1) then
        do j=2,np
          j1=j-1
          sum=sum+y(j1)*(x(j)-x(j1))
        enddo
      elseif (lep.eq.2) then
        do j=2,np
          j1=j-1
          sum=sum+0.5d0*(y(j)+y(j1))*(x(j)-x(j1))
        enddo
      elseif (lep.eq.4) then
        do j=2,np
          j1=j-1
          xr=x(j)-x(j1)
          if (xr.gt.zero) then
            y1=y(j1)
            y2=y(j)
            yr=y2-y1
            ay1=abs(y1)
            ayr=abs(yr)
            if (y1*y2.gt.zero.and.ayr.gt.eps*ay1) then
              sum=sum+xr*yr/log(y2/y1)
            else
              sum=sum+0.5d0*(y1+y2)*xr
            endif
          endif
        enddo
      elseif (lep.eq.3) then
        do j=2,np
          j1=j-1
          x1=x(j1)
          x2=x(j)
          xr=x2-x1
          if (xr.gt.zero) then
            y1=y(j1)
            y2=y(j)
            ax1=abs(x1)
            if (x1*x2.gt.zero.and.xr.gt.eps*ax1) then
              sum=sum+x2*y2-x1*y1-xr*(y2-y1)/log(x2/x1)
            else
              sum=sum+0.5d0*(y1+y2)*xr
            endif
          endif
        enddo
      elseif (lep.eq.5) then
        do j=2,np
          j1=j-1
          x1=x(j1)
          x2=y(j)
          xr=x2-x1
          if (xr.gt.zero) then
            y1=y(j1)
            y2=y(j)
            yr=y2-y1
            ax1=abs(x1)
            ay1=abs(y1)
            ayr=abs(yr)
            x12=x1*x2
            y12=y1*y2
            if (x12.gt.zero.and.xr.gt.eps*ax1.and.
     &          y12.gt.zero.and.ayr.gt.eps*ay1) then
              d=log(y2/y1)/log(x2/x1)+1.0d0
              sum=sum+y1*x1*(((x2/x1)**d)-1.0d0)/d
            elseif (y12.gt.zero.and.ayr.gt.eps*ay1) then
              sum=sum+xr*yr/log(y2/y1)
            elseif (x12.gt.zero.and.xr.gt.eps*ax1) then
              sum=sum+x2*y2-x1*y1-xr*yr/log(x2/x1)
            else
              sum=sum+0.5d0*(y1+y2)*xr
            endif
          endif
        enddo
      else
        write(*,*)' ERROR in renorm: ilaw=',ilaw,' not coded'
        stop
      endif
      fn=1.0d0
      if (c.eq.zero) then
        c=sum
      else
        if (sum.ne.0.0d0) then
          fn=c/sum
          do j=1,np
            y(j)=y(j)*fn
          enddo
        else
          write(*,*)' Warning: Integral equal zero'
        endif
      endif
      return
      end
C======================================================================
C      General routines for ENDF-6 formatted files
C======================================================================
      subroutine readtext(nin,line,mat,mf,mt,ns)
      character*66 line
      read(nin,10)line,mat,mf,mt,ns
      return
   10 format(a66,i4,i2,i3,i5)
      end
C======================================================================
      subroutine readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
      implicit real*8 (a-h, o-z)
      read(nin,10)c1,c2,l1,l2,n1,n2,mat,mf,mt,ns
      return
   10 format(2e11.0,4i11,i4,i2,i3,i5)
      end
C======================================================================
      subroutine readlist(nin,c1,c2,l1,l2,npl,n2,b)
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
      implicit real*8 (a-h, o-z)
      dimension nbt(*),intp(*)
      read(nin,10)c1,z,l1,l2,nr,nz,mat,mf,mt,ns
      read(nin,20)(nbt(n),intp(n),n=1,nr)
      return
   10 format(2e11.0,4i11,i4,i2,i3,i5)
   20 format(6i11)
      end
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
      subroutine findnextmt(nin,mf,mt)
c
c     find next mt on mf file for current material
c
      character*66 line
      mt0=-1
      do while (mt0.ne.0)
        read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf0,mt0,ns
      enddo
      read(nin,'(a66,i4,i2,i3,i5)',err=10,end=10)line,mat0,mf0,mt,ns
      if (mt.ne.0.and.mf0.eq.mf) then
        backspace nin
      else
        mt=-1
      endif
      return
   10 mt=-2
      return
      end
C======================================================================
      subroutine nextsub6(nin,law,nbt,ibt,x,b)
c
c     find next subsection on MF6 section
c
      implicit real*8 (a-h,o-z)
      dimension nbt(*),ibt(*),x(*),b(*)
      if (law.eq.1.or.law.eq.2.or.law.eq.5) then
        call readtab2(nin,c1,c2,l1,l2,n1,ne,nbt,ibt)
        do i=1,ne
          call readlist(nin,c1,c2,l1,l2,n1,n2,b)
        enddo
      elseif (law.eq.6) then
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns)
      elseif (law.eq.7) then
        call readtab2(nin,c1,c2,l1,l2,n1,ne,nbt,ibt)
        do i=1,ne
          call readtab2(nin,c1,c2,l1,l2,n1,nmu,nbt,ibt)
          do j=1,nmu
            call readtab1(nin,c1,c2,l1,l2,n1,n2,nbt,ibt,x,b)
          enddo
        enddo
      elseif (law.lt.0.or.law.gt.7) then
       write(*,*)' ERROR: unknown LAW=',law,' on MF6'
       stop
      endif
      return
      end
C======================================================================
C     Get date and time
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
C     Legendre polynomials
C======================================================================
      subroutine leg2lin(na,b,np,xt,yt,epy,ymin,npmax)
      implicit real*8 (a-h, o-z)
      parameter (nmu0max=3, ns=25, h0=1.0d-8, yzero=1.0d-30)
      dimension b(*),xt(*),yt(*)
      dimension xs(ns),ys(ns),xmu0(nmu0max)
c
c     Starting convertion from Legendre representation to tabular data
c
      if (na.gt.0) then
        if (na.eq.1) then
          xmu0(1)=-1.0d0
          xmu0(2)=1.0d0
          nmu0=2
        else
          xmu0(1)=-1.0d0
          xmu0(2)=0.0d0
          xmu0(3)=1.0d0
          nmu0=3
        endif
        x0=xmu0(1)
        y0=yleg(x0,b,na)
        np=1
        xt(1)=x0
        yt(1)=y0
        do i=2,nmu0
          x1=xmu0(i)
          y1=yleg(x1,b,na)
          dymax=0.05d0*max(abs(y0),abs(y1))
          k=0
          nostop=1
          do while (nostop.eq.1)
            xm=0.5d0*(x0+x1)
            ymlin=0.5d0*(y0+y1)
            ymlaw=yleg(xm,b,na)
            dy=dabs(ymlin-ymlaw)
            ymabs=dabs(ymlaw)
            ylabs=dabs(ymlin)
            dx=dabs(xm-x0)
            if ((dy.le.epy*ymabs.and.dabs(y1-y0).le.(ylabs+dymax)).or.
     &          ymabs.le.ymin.or.dx.le.h0.or.k.eq.ns) then
              call inc(np,npmax)
              xt(np)=x1
              yt(np)=y1
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
              y1=ymlaw
            endif
          enddo
          x0=x1
          y0=y1
        enddo
        do i=1,np
          if (yt(i).lt.0.0d0) yt(i)=yzero
        enddo
        c=b(1)
        if (c.gt.0.0d0) then
          ilaw=2
          call renorm(np,xt,yt,ilaw,c,fn)
        endif
      else
        np=2
        xt(1)=-1.0d0
        xt(2)=1.0d0
        yt(1)=0.5d0
        yt(2)=0.5d0
      endif
      return
      end
C======================================================================
      real*8 function yleg(x,a,na)
c
c     calculate y(x) using a legendre expansion
c
      implicit real*8 (a-h,o-z)
      dimension a(*),p(65)
      call legndr(x,p,na)
      n=na+1
      yleg=0.0d0
      do l=1,n
        yleg=yleg+(dble(l)-0.5d0)*a(l)*p(l)
      enddo
      return
      end
C======================================================================
C     Linearization
C======================================================================
      subroutine linear(nr,nbt,ibt,np,x,y,epy,ymin,npmax)
c
c      Linearize x,y data
c
      implicit real*8(a-h, o-z)
      parameter (ns=32)
      dimension nbt(*),ibt(*),x(*),y(*)
      dimension xt(npmax),yt(npmax),xs(ns),ys(ns)
c
c     Check interpolation law
c
      call checklaw(nr,ibt,icod)
      if (icod.eq.0) then
        if (nr.ne.1) then
          nr=1
          nbt(nr)=np
          ibt(nr)=2
        endif
        return
      endif
c
c     Check for minimun space
c
      if (np.gt.npmax) then
        write(*,*)' === Increase dimension of arrays ',npmax
        write(*,*)' === Stop in linear'
        stop
      endif
c
c      Starting linearization
c
      ii=0
      do j=1,nr
        ilaw=ibt(j)
        ilaw=ilaw-(ilaw/10)*10
        iup=nbt(j)
        if (j.gt.1) then
          ilow=nbt(j-1)+1
        else
          ilow=1
        endif
        if(ilaw.eq.2) then
c
c         linear-linear interpolation law, just copy data
c
          do i=ilow,iup
            call inc(ii,npmax)
            xt(ii)=x(i)
            yt(ii)=y(i)
          enddo
        elseif (ilaw.eq.1) then
c
c         histogram interpolation
c         constant interpolation law, convert to linear by
c         given linear interpolation with zero slope
c
          call inc(ii,npmax)
          xt(ii)=x(ilow)
          yt(ii)=y(ilow)
          do i=ilow+1,iup
            call inc(ii,npmax)
            xt(ii)=x(i)
            yt(ii)=y(i-1)
            call inc(ii,npmax)
            xt(ii)=x(i)
            yt(ii)=y(i)
          enddo
        else
c
c         lin-log, log-lin, log-log
c
          x0=x(ilow)
          y0=y(ilow)
          call inc(ii,npmax)
          xt(ii)=x0
          yt(ii)=y0
          do i=ilow+1,iup
            x1=x(i)
            y1=y(i)
            dymax=0.1d0*max(dabs(y0),dabs(y1))
            k=0
            nostop=1
            do while (nostop.eq.1)
              xm=0.5d0*(x0+x1)
              ymlin=0.5d0*(y0+y1)
              call terp1m(x0,y0,x1,y1,xm,ymlaw,ilaw)
              dy=dabs(ymlin-ymlaw)
              ymabs=dabs(ymlaw)
              ylabs=dabs(ymlin)
              dx=dabs(xm-x0)
              h0=dabs(edelta(x0,5.0d0)-x0)
              if ((dy.le.epy*ymabs.and.dabs(y1-y0).le.(ylabs+dymax)).or.
     &            ymabs.le.ymin.or.dx.le.h0.or.k.eq.ns) then
                call inc(ii,npmax)
                xt(ii)=x1
                yt(ii)=y1
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
                y1=ymlaw
              endif
            enddo
            x0=x1
            y0=y1
          enddo
        endif
      enddo
c
c     back to original array
c
      np=ii
      do i=1,np
        x(i)=xt(i)
        y(i)=yt(i)
      enddo
      nr=1
      nbt(nr)=np
      ibt(nr)=2
      return
      end
C======================================================================
C      The following subroutines were taken from NJOY2016 and
C      adapted/modified by D. Lopez Aldama for ACEMAKER:
C       1. subroutine legndr
C       2. subroutine terp1 (renamed as terp1m)
C       3. real*8 function bachaa
C       4. subroutine change
C       5. subroutine typen
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
      subroutine legndr(x,p,np)
c     *****************************************************************
c       generate legendre polynomials at x by recursion.
c       place pl in p(l+1).
c     *****************************************************************
      implicit real*8 (a-h,o-z)
      dimension p(*)
      p(1)=1.0d0
      p(2)=x
      if (np.lt.2) return
      m1=np-1
      do i=1,m1
         g=x*p(i+1)
         h=g-p(i)
         p(i+2)=h+g-h/(i+1)
      enddo
      return
      end
C======================================================================
      subroutine terp1m(x1,y1,x2,y2,x,y,i)
c
c      interpolate one point using ENDF-6 interpolation laws
c      (x1,y1) and (x2,y2) are the end points
c      (x,y) is the interpolated point
c      i is the interpolation law (1-5)
c
c      Adapted by D. L. Aldama for ACEMAKER
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
      real*8 function bachaa(iza1i,iza2,izat,e,ep)
c     Adapted by D. L. Aldama for ACEMAKER
c     ******************************************************************
c     compute the kalbach a parameter
c     ******************************************************************
      implicit real*8 (a-h,o-z)
      real*8 nc,nb
      data amassn/1.00866491578d0/
      data amu/1.6605402d-24/
      data ev/1.60217733d-12/
      data clight/2.99792458d10/
      data third,twoth,fourth/.333333333d0,.666666667d0,1.33333333d0/
      data c1,c2,c3,c4,c5,c6/15.68d0,-28.07d0,-18.56d0,33.22d0,
     &  -0.717d0,1.211d0/
      data s2,s3,s4,s5/2.22d0,8.48d0,7.72d0,28.3d0/
      data b1,b2,b3/0.04d0,1.8d-6,6.7d-7/
      data d1/9.3d0/
      data ea1,ea2/41.d0,130.d0/
      data emev/1.d6/
      emc2=amassn*amu*clight*clight/ev/emev
c
      iza1=iza1i
      if (iza1i.eq.0) iza1=1
      iza=izat
      if (iza.eq.6000) iza=6012
      if (iza.eq.12000) iza=12024
      if (iza.eq.14000) iza=14028
      if (iza.eq.16000) iza=16032
      if (iza.eq.17000) iza=17035
      if (iza.eq.19000) iza=19039
      if (iza.eq.20000) iza=20040
      if (iza.eq.22000) iza=22048
      if (iza.eq.23000) iza=23051
      if (iza.eq.24000) iza=24052
      if (iza.eq.26000) iza=26056
      if (iza.eq.28000) iza=28058
      if (iza.eq.29000) iza=29063
      if (iza.eq.31000) iza=31069
      if (iza.eq.40000) iza=40090
      if (iza.eq.42000) iza=42096
      if (iza.eq.48000) iza=48112
      if (iza.eq.49000) iza=49115
      if (iza.eq.50000) iza=50120
      if (iza.eq.63000) iza=63151
      if (iza.eq.72000) iza=72178
      if (iza.eq.74000) iza=74184
      if (iza.eq.82000) iza=82208
      aa=mod(iza,1000)
      if (aa.eq.0.) then
         write(*,*)' Dominant isotope not known for ',iza,' in bachaa'
         stop
      endif
      za=int(iza/1000)
      ac=aa+mod(iza1,1000)
      zc=za+int(iza1/1000)
      ab=ac-mod(iza2,1000)
      zb=zc-int(iza2/1000)
      na=nint(aa-za)
      nb=nint(ab-zb)
      nc=nint(ac-zc)
      sa=c1*(ac-aa)
     &  +c2*((nc-zc)**2/ac-(na-za)**2/aa)
     &  +c3*(ac**twoth-aa**twoth)
     &  +c4*((nc-zc)**2/ac**fourth-(na-za)**2/aa**fourth)
     &  +c5*(zc**2/ac**third-za**2/aa**third)
     &  +c6*(zc**2/ac-za**2/aa)
      if (iza1.eq.1002) sa=sa-s2
      if (iza1.eq.1003) sa=sa-s3
      if (iza1.eq.2003) sa=sa-s4
      if (iza1.eq.2004) sa=sa-s5
      sb=c1*(ac-ab)
     &  +c2*((nc-zc)**2/ac-(nb-zb)**2/ab)
     &  +c3*(ac**twoth-ab**twoth)
     &  +c4*((nc-zc)**2/ac**fourth-(nb-zb)**2/ab**fourth)
     &  +c5*(zc**2/ac**third-zb**2/ab**third)
     &  +c6*(zc**2/ac-zb**2/ab)
      if (iza2.eq.1002) sb=sb-s2
      if (iza2.eq.1003) sb=sb-s3
      if (iza2.eq.2003) sb=sb-s4
      if (iza2.eq.2004) sb=sb-s5
      ecm=aa*e/ac
      ea=ecm+sa
      eb=ep*ac/ab+sb
      x1=eb
      if (ea.gt.ea2) x1=ea2*eb/ea
      x3=eb
      if (ea.gt.ea1) x3=ea1*eb/ea
      fa=1
      if (iza1.eq.2004) fa=0
      fb=1
      if (iza2.eq.1) fb=fb/2
      if (iza2.eq.2004) fb=2
      bb=b1*x1+b2*x1**3+b3*fa*fb*x3**4
      if (iza1i.eq.0) then
         fact=d1
         if (ep.ne.0.) fact=fact/sqrt(ep)
         test=1
         if (fact.lt.test) fact=test
         test=4
         if (fact.gt.test) fact=test
         bb=bb*sqrt(e/(2*emc2))*fact
      endif
      bachaa=bb
      return
      end
C======================================================================
      subroutine change(nout,mcnpx)
c     Adapted by D. Lopez Aldama for ACEMAKER
c     ******************************************************************
c     change ace data fields from integer to real or vice versa.
c     if nout.gt.1, the results are written in type 1 format
c        (all fields are assumed to contain real numbers).
c     if nout.eq.0, real fields are changed to integers in memory
c        (all fields are assumed to contain real numbers).
c     if nout.eq.1, integer fields are changed to real in memory
c        (fields are assumed to contain mixed reals and integers).
c     ******************************************************************
      implicit real*8 (a-h,o-z)
      integer esz,sig,and,tyr,dlw,gpd,fis,sigp,andp,dlwp,yp,end
      integer dndat,dnd,ptype,ploct
      character*10 hd,hm
      character*13 hz
      character*70 hk
      common/acetxt/hz,hd,hm,hk
      common/acecte/awr0,tz,awn(16),izn(16)
      common/acepnt/len2,mzaid,nes,ntr,nrx,ntrp,ntype,ndnf,
     &              liso,izt,iat,nxsd(5),
     &              esz,nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,gpd,
     &              mtrp,lsigp,sigp,landp,andp,ldlwp,dlwp,yp,fis,end,
     &              iurpt,nud,dndat,ldnd,dnd,jxsd(2),ptype,ntro,ploct
      common/acedat/xss(200000000),nxss
      dimension nxs(16),jxs(32)
      dimension iss(1)
      equivalence (nxs(1),len2)
      equivalence (jxs(1),esz)
      equivalence (xss(1),iss(1))
      external typen
      if (nout.gt.1) then
       if (mcnpx.eq.1) then
         write(nout,'(a13,f12.6,1pe12.4,1x,a10)')hz(1:13),awr0,tz,hd
       else
         write(nout,'(a10,f12.6,1pe12.4,1x,a10)')hz(1:10),awr0,tz,hd
       endif
       write(nout,'(a70,a10)')hk,hm
       write(nout,'(4(i7,f11.0))')(izn(i),awn(i),i=1,16)
       write(nout,'(8i9)')(nxs(i),i=1,16)
       write(nout,'(8i9)')(jxs(i),i=1,32)
      endif
c
c     ***write or convert esz block
      n=5*nes
      do i=1,n
         call typen(i,nout,2)
      enddo
c
c     ***nu block
      if (nu.ne.0) then
         l=nu
         if (nout.ne.1) lnu=nint(xss(l))
         if (nout.eq.1) lnu=iss(l)
         if (lnu.gt.0) then
           m=1
         else
           m=2
           call typen(l,nout,1)
           l=l+1
         endif
         do i=1,m
            if (nout.ne.1) lnu=nint(xss(l))
            if (nout.eq.1) lnu=iss(l)
            call typen(l,nout,1)
            l=l+1
            if (lnu.ne.2) then
               if (nout.ne.1) nc=nint(xss(l))
               if (nout.eq.1) nc=iss(l)
               call typen(l,nout,1)
               l=l+1
               do j=1,nc
                  call typen(l,nout,2)
                  l=l+1
               enddo
            else
               if (nout.ne.1) nr=nint(xss(l))
               if (nout.eq.1) nr=iss(l)
               call typen(l,nout,1)
               l=l+1
               if (nr.ne.0) then
                  n=2*nr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               if (nout.ne.1) ne=nint(xss(l))
               if (nout.eq.1) ne=iss(l)
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
            endif
         enddo
      endif
c
c     ***cross section data
      if (ntr.ne.0) then
c
c        ***mtr block
         l=mtr
         do i=1,ntr
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***lqr block
         l=lqr
         do i=1,ntr
            call typen(l,nout,2)
            l=l+1
         enddo
c
c        ***tyr block
         l=tyr
         do i=1,ntr
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***lsig block
         l=lsig
         do i=1,ntr
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***sig block
         l=sig
         do i=1,ntr
            call typen(l,nout,1)
            l=l+1
            if (nout.ne.1) ne=nint(xss(l))
            if (nout.eq.1) ne=iss(l)
            call typen(l,nout,1)
            l=l+1
            do j=1,ne
               call typen(l,nout,2)
               l=l+1
            enddo
         enddo
      endif
c
c     ***land block
      n=nrx+1
      l=land
      li=l-1
      do i=1,n
         call typen(l,nout,1)
         l=l+1
      enddo
c
c     ***and block
      l=and
      do i=1,n
         if (nout.ne.0) nn=nint(xss(li+i))
         if (nout.eq.0) nn=iss(li+i)
         if (nn.gt.0) then
            if (nout.ne.1) ne=nint(xss(l))
            if (nout.eq.1) ne=iss(l)
            call typen(l,nout,1)
            l=l+1
            do j=1,ne
               call typen(l,nout,2)
               l=l+1
            enddo
            ll=l-1
            do j=1,ne
               call typen(l,nout,1)
               l=l+1
            enddo
            do j=1,ne
               if (nout.ne.0) nn=nint(xss(ll+j))
               if (nout.eq.0) nn=iss(ll+j)
               if (nn.ne.0) then
                  if (nn.ge.0) then
                     do k=1,33
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  else
                     if (nout.ne.0) np=nint(xss(iabs(nn)+and))
                     if (nout.eq.0) np=iss(iabs(nn)+and)
                     call typen(l,nout,1)
                     l=l+1
                     call typen(l,nout,1)
                     l=l+1
                     nw=3*np
                     do k=1,nw
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  endif
               endif
            enddo
         endif
      enddo
c
c     ***distributions
      if (nrx.ne.0) then
c
c        ***ldlw block
         l=ldlw
         do i=1,nrx
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***dlw block
         l=dlw
         do i=1,nrx
            if (nout.ne.0) ly=nint(xss(tyr+i-1))
            if (nout.eq.0) ly=iss(tyr+i-1)
            ly=iabs(ly)
            if (ly.gt.100) then
               l=ly-100+dlw-1
               if (nout.ne.1) nr=nint(xss(l))
               if (nout.eq.1) nr=iss(l)
               call typen(l,nout,1)
               l=l+1
               if (nr.gt.0) then
                  n=2*nr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               if (nout.ne.1) ne=nint(xss(l))
               if (nout.eq.1) ne=iss(l)
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
            endif
c
c           ***loop over laws
            lnw=1
            do while (lnw.gt.0)
               if (nout.ne.1) lnw=nint(xss(l))
               if (nout.eq.1) lnw=iss(l)
               call typen(l,nout,1)
               l=l+1
               if (nout.ne.1) law=nint(xss(l))
               if (nout.eq.1) law=iss(l)
               call typen(l,nout,1)
               l=l+1
               call typen(l,nout,1)
               l=l+1
               if (nout.ne.1) nr=nint(xss(l))
               if (nout.eq.1) nr=iss(l)
               call typen(l,nout,1)
               l=l+1
               if (nr.gt.0) then
                  n=2*nr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               if (nout.ne.1) ne=nint(xss(l))
               if (nout.eq.1) ne=iss(l)
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
c
c              ***law 1
               if (law.eq.1) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  if (nout.ne.1) net=nint(xss(l))
                  if (nout.eq.1) net=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     do k=1,net
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  enddo
c
c              ***law 3
               else if (law.eq.3) then
                  call typen(l,nout,2)
                  l=l+1
                  call typen(l,nout,2)
                  l=l+1
c
c              ***law 4
               else if (law.eq.4) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                     if (nout.ne.1) np=nint(xss(l))
                     if (nout.eq.1) np=iss(l)
                     call typen(l,nout,1)
                     l=l+1
                     n=3*np
                     do k=1,n
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  enddo
c
c              ***law 5
               else if (law.eq.5) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  n=2*ne
                  do j=1,n
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  if (nout.ne.1) net=nint(xss(l))
                  if (nout.eq.1) net=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,net
                        call typen(l,nout,2)
                     l=l+1
                  enddo
c
c              ***law 7 or law 9
               else if (law.eq.7.or.law.eq.9) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  n=2*ne
                  do j=1,n
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  call typen(l,nout,2)
                  l=l+1
c
c              ***law 11
               else if (law.eq.11) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  n=2*ne
                  do j=1,n
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  n=2*ne
                  do j=1,n
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  call typen(l,nout,2)
                  l=l+1
c
c              ***law 44
               else if (law.eq.44) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                     if (nout.ne.1) np=nint(xss(l))
                     if (nout.eq.1) np=iss(l)
                     call typen(l,nout,1)
                     l=l+1
                     n=5*np
                     do k=1,n
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  enddo
c
c              ***law 61
               else if (law.eq.61) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.ne.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                     do j=1,ne
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                     if (nout.ne.1) np=nint(xss(l))
                     if (nout.eq.1) np=iss(l)
                     call typen(l,nout,1)
                     l=l+1
                     n=3*np
                     do k=1,n
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                     do k=1,np
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                     do k=1,np
                        call typen(l,nout,1)
                        l=l+1
                        if (nout.ne.1) nmu=nint(xss(l))
                        if (nout.eq.1) nmu=iss(l)
                        call typen(l,nout,1)
                        l=l+1
                        nw=3*nmu
                        do kk=1,nw
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                     enddo
                  enddo
c
c              ***law 66
               else if (law.eq.66) then
                  call typen(l,nout,1)
                  l=l+1
                  call typen(l,nout,2)
                  l=l+1
c                  call typen(l,nout,1)
c                  l=l+1
c                  if (nout.ne.1) nn=nint(xss(l))
c                  if (nout.eq.1) nn=iss(l)
c                  n=3*nn
c                  call typen(l,nout,1)
c                  l=l+1
c                  do k=1,n
c                     call typen(l,nout,2)
c                     l=l+1
c                  enddo
c
c              ***law 67
               else if (law.eq.67) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                     if (nout.ne.1) nmu=nint(xss(l))
                     if (nout.eq.1) nmu=iss(l)
                     call typen(l,nout,1)
                     l=l+1
                     do k=1,nmu
                        call typen(l,nout,2)
                        l=l+1
                        enddo
                     do k=1,nmu
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                     do k=1,nmu
                        call typen(l,nout,1)
                        l=l+1
                        if (nout.ne.1) nep=nint(xss(l))
                        if (nout.eq.1) nep=iss(l)
                        call typen(l,nout,1)
                        l=l+1
                        nn=3*nep
                        do n=1,nn
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif
c
c     ***unresolved-range probability-table block
      if (iurpt.gt.0) then
         l=iurpt
         if (nout.ne.1) nure=nint(xss(l))
         if (nout.eq.1) nure=iss(l)
         call typen(l,nout,1)
         l=l+1
         if (nout.ne.1) nurb=nint(xss(l))
         if (nout.eq.1) nurb=iss(l)
         call typen(l,nout,1)
         l=l+1
         call typen(l,nout,1)
         l=l+1
         call typen(l,nout,1)
         l=l+1
         call typen(l,nout,1)
         l=l+1
         call typen(l,nout,1)
         l=l+1
         n=nure*(1+6*nurb)
         do i=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
      endif
c
c
c     ***delayed neutron block
      if (ndnf.ne.0) then
c        delayed nubar
         l=nud
         if (nout.ne.1) lnu=nint(xss(l))
         if (nout.eq.1) lnu=iss(l)
         call typen(l,nout,1)
         l=l+1
         if (nout.ne.1) nr=nint(xss(l))
         if (nout.eq.1) nr=iss(l)
         call typen(l,nout,1)
         l=l+1
         if (nr.ne.0) then
            n=2*nr
            do j=1,n
               call typen(l,nout,1)
               l=l+1
            enddo
         endif
         if (nout.ne.1) ne=nint(xss(l))
         if (nout.eq.1) ne=iss(l)
         call typen(l,nout,1)
         l=l+1
         n=2*ne
         do j=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
c        precursor data
         l=dndat
         do i=1,ndnf
            call typen(l,nout,2)
            l=l+1
            if (nout.ne.1) nr=nint(xss(l))
            if (nout.eq.1) nr=iss(l)
            call typen(l,nout,1)
            l=l+1
            if (nr.ne.0) then
               n=2*nr
               do j=1,n
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
            if (nout.ne.1) ne=nint(xss(l))
            if (nout.eq.1) ne=iss(l)
            call typen(l,nout,1)
            l=l+1
            n=2*ne
            do j=1,n
               call typen(l,nout,2)
               l=l+1
            enddo
         enddo
c        precursor energy distribution locators
         do i=1,ndnf
            call typen(l,nout,1)
            l=l+1
         enddo
c        precursor energy distributions
         do i=1,ndnf
            call typen(l,nout,1)
            l=l+1
            call typen(l,nout,1)
            l=l+1
            call typen(l,nout,1)
            l=l+1
            if (nout.ne.1) nr=nint(xss(l))
            if (nout.eq.1) nr=iss(l)
            call typen(l,nout,1)
            l=l+1
            if (nr.ne.0) then
               n=2*nr
               do j=1,n
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
            if (nout.ne.1) ne=nint(xss(l))
            if (nout.eq.1) ne=iss(l)
            call typen(l,nout,1)
            l=l+1
            n=2*ne
            do j=1,n
               call typen(l,nout,2)
               l=l+1
            enddo
c           law=4 data
            if (nout.ne.1) nr=nint(xss(l))
            if (nout.eq.1) nr=iss(l)
            call typen(l,nout,1)
            l=l+1
            if (nr.gt.0) then
               n=2*nr
               do j=1,n
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
            if (nout.ne.1) ne=nint(xss(l))
            if (nout.eq.1) ne=iss(l)
            call typen(l,nout,1)
            l=l+1
            do j=1,ne
               call typen(l,nout,2)
               l=l+1
            enddo
            do j=1,ne
               call typen(l,nout,1)
               l=l+1
            enddo
            do j=1,ne
               call typen(l,nout,1)
               l=l+1
               if (nout.ne.1) np=nint(xss(l))
               if (nout.eq.1) np=iss(l)
               call typen(l,nout,1)
               l=l+1
               n=3*np
               do k=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
            enddo
         enddo
      endif
c     ***gpd block
      if (gpd.ne.0) then
         l=gpd
         do i=1,nes
            call typen(l,nout,2)
            l=l+1
         enddo
         if (mtrp.gt.gpd+nes) then
            n=20*30
            do i=1,n
               call typen(l,nout,2)
               l=l+1
            enddo
         endif
      endif
c
c     ***detailed photon production
      if (ntrp.ne.0) then
c
c        ***mtrp block
         l=mtrp
         do i=1,ntrp
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***lsigp block
         l=lsigp
         do i=1,ntrp
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***sigp block
         l=sigp
         do i=1,ntrp
            if (nout.ne.1) mftype=nint(xss(l))
            if (nout.eq.1) mftype=iss(l)
            call typen(l,nout,1)
            l=l+1
            if (mftype.ne.12.and.mftype.ne.16) then
               call typen(l,nout,1)
               l=l+1
               if (nout.ne.1) ne=nint(xss(l))
               if (nout.eq.1) ne=iss(l)
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
            else
               call typen(l,nout,1)
               l=l+1
               if (nout.ne.1) nr=nint(xss(l))
               if (nout.eq.1) nr=iss(l)
               call typen(l,nout,1)
               l=l+1
               if (nr.gt.0) then
                  n=2*nr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               if (nout.ne.1) ne=nint(xss(l))
               if (nout.eq.1) ne=iss(l)
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
            endif
         enddo
c
c        ***landp block
         l=landp
         li=l-1
         do i=1,ntrp
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***andp block
         l=andp
         do i=1,ntrp
            if (nout.ne.0) nn=nint(xss(li+i))
            if (nout.eq.0) nn=iss(li+i)
            if (nn.gt.0) then
               if (nout.ne.1) ne=nint(xss(l))
               if (nout.eq.1) ne=iss(l)
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
               ll=l-1
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
               enddo
               do j=1,ne
                  if (nout.ne.0) nn=nint(xss(ll+j))
                  if (nout.eq.0) nn=iss(ll+j)
                  if (nn.gt.0) then
                     do k=1,33
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  endif
               enddo
            endif
         enddo
c
c        ***ldlwp block
         l=ldlwp
         do i=1,ntrp
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***dlwp block
         l=dlwp
         do i=1,ntrp
            lnw=1
            do while (lnw.ne.0)
               if (nout.ne.1) lnw=nint(xss(l))
               if (nout.eq.1) lnw=iss(l)
               call typen(l,nout,1)
               l=l+1
               if (nout.ne.1) law=nint(xss(l))
               if (nout.eq.1) law=iss(l)
               call typen(l,nout,1)
               l=l+1
               call typen(l,nout,1)
               l=l+1
               if (nout.ne.1) nr=nint(xss(l))
               if (nout.eq.1) nr=iss(l)
               call typen(l,nout,1)
               l=l+1
               if (nr.gt.0) then
                  n=2*nr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               if (nout.ne.1) ne=nint(xss(l))
               if (nout.eq.1) ne=iss(l)
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
c
c              ***law 1
               if (law.eq.1) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  if (nout.ne.1) net=nint(xss(l))
                  if (nout.eq.1) net=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     do k=1,net
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  enddo
c
c              ***law 2
               else if (law.eq.2) then
                  call typen(l,nout,1)
                  l=l+1
                  call typen(l,nout,2)
                  l=l+1
c
c              ***law 4 and law 44
               else if (law.eq.4.or.law.eq.44) then
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do j=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  if (nout.ne.1) ne=nint(xss(l))
                  if (nout.eq.1) ne=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                  enddo
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                     if (nout.ne.1) np=nint(xss(l))
                     if (nout.eq.1) np=iss(l)
                     call typen(l,nout,1)
                     l=l+1
                     n=3*np
                     do k=1,n
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  enddo
               endif
            enddo
         enddo
c
c        ***yp block
         l=yp
         if (nout.ne.1) nyp=nint(xss(l))
         if (nout.eq.1) nyp=iss(l)
         call typen(l,nout,1)
         l=l+1
         do i=1,nyp
            call typen(l,nout,1)
            l=l+1
         enddo
      endif
c
c     ***particle production blocks
      if (ntype.gt.0) then
c
c        ***ptype, ntro, and ixs arrays
         do i=1,ntype
            call typen(l,nout,1)
            l=l+1
         enddo
         ntro=l
         do i=1,ntype
            call typen(l,nout,1)
            l=l+1
         enddo
         nw=10*ntype
         do i=1,nw
            call typen(l,nout,1)
            l=l+1
         enddo
c
c        ***loop over particle types
         do i=1,ntype
c
c           ***hpd block
            call typen(l,nout,1)
            l=l+1
            if (nout.ne.0) ne=nint(xss(l))
            if (nout.eq.0) ne=iss(l)
            call typen(l,nout,1)
            l=l+1
            if (ne.ne.0) then
               nw=2*ne
               do j=1,nw
                  call typen(l,nout,2)
                  l=l+1
               enddo
               ntrh=nint(xss(ntro-1+i))
c
c              ***mtrh block
               do k=1,ntrh
                  call typen(l,nout,1)
                  l=l+1
               enddo
c
c                 ***tyrh block
               do k=1,ntrh
                  call typen(l,nout,1)
                  l=l+1
               enddo
c
c              ***lsigh block
               do k=1,ntrh
                  call typen(l,nout,1)
                  l=l+1
               enddo
c
c              ***sigh block
               do j=1,ntrh
                  call typen(l,nout,1)
                  l=l+1
                  call typen(l,nout,1)
                  l=l+1
                  if (nout.ne.1) nr=nint(xss(l))
                  if (nout.eq.1) nr=iss(l)
                  call typen(l,nout,1)
                  l=l+1
                  if (nr.gt.0) then
                     n=2*nr
                     do jj=1,n
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  ne=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  nw=2*ne
                  do k=1,nw
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               enddo
c
c              ***landh block
               li=l-1
               do k=1,ntrh
                  call typen(l,nout,1)
                  l=l+1
               enddo
c
c              ***andh block
               do ir=1,ntrh
                  if (nout.ne.0) nn=nint(xss(li+ir))
                  if (nout.eq.0) nn=iss(li+ir)
                  if (nn.gt.0) then
                     if (nout.ne.1) ne=nint(xss(l))
                     if (nout.eq.1) ne=iss(l)
                     call typen(l,nout,1)
                     l=l+1
                     do j=1,ne
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                     ll=l-1
                     do j=1,ne
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                     do j=1,ne
                        if (nout.ne.0) nn=nint(xss(ll+j))
                        if (nout.eq.0) nn=iss(ll+j)
                        if (nn.gt.0) then
                           do k=1,33
                              call typen(l,nout,2)
                              l=l+1
                           enddo
                        else if (nn.lt.0) then
                           call typen(l,nout,1)
                           l=l+1
                           if (nout.ne.0) np=nint(xss(l))
                           if (nout.eq.0) np=iss(l)
                           call typen(l,nout,1)
                           l=l+1
                           nw=3*np
                           do k=1,nw
                              call typen(l,nout,2)
                              l=l+1
                           enddo
                        endif
                     enddo
                  endif
               enddo
c
c              ***ldlwh block
               li=l-1
               do k=1,ntrh
                  call typen(l,nout,1)
                  l=l+1
               enddo
c
c              ***dlwh block
               do ii=1,ntrh
                  if (nout.ne.0) nn=nint(xss(li+ii))
                  if (nout.eq.0) nn=iss(li+ii)
                  if (nn.gt.0) then
                     call typen(l,nout,1)
                     l=l+1
                     law=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     call typen(l,nout,1)
                     l=l+1
                     nr=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     if (nr.ne.0) then
                        nw=2*nr
                        do k=1,nw
                           call typen(l,nout,1)
                           l=l+1
                        enddo
                     endif
                     ne=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     nw=2*ne
                     do k=1,nw
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                     if (law.eq.4) then
                        nr=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        ne=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        do k=1,ne
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                        do k=1,ne
                           call typen(l,nout,1)
                           l=l+1
                        enddo
                        do k=1,ne
                           call typen(l,nout,1)
                           l=l+1
                           np=nint(xss(l))
                           call typen(l,nout,1)
                           l=l+1
                           nw=3*np
                           do kk=1,nw
                              call typen(l,nout,2)
                              l=l+1
                           enddo
                        enddo
                     else if (law.eq.44) then
                        nr=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        ne=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        do k=1,ne
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                        do k=1,ne
                           call typen(l,nout,1)
                           l=l+1
                        enddo
                        do j=1,ne
                           call typen(l,nout,1)
                           l=l+1
                           np=nint(xss(l))
                           call typen(l,nout,1)
                           l=l+1
                           nw=5*np
                           do k=1,nw
                              call typen(l,nout,2)
                              l=l+1
                           enddo
                        enddo
                     else if (law.eq.33) then
                        call typen(l,nout,2)
                        l=l+1
                        call typen(l,nout,2)
                        l=l+1
                     else if (law.eq.66) then
                        call typen(l,nout,1)
                        l=l+1
                        call typen(l,nout,2)
                        l=l+1
                        call typen(l,nout,1)
                        l=l+1
                        if (nout.ne.1) nn=nint(xss(l))
                        if (nout.eq.1) nn=iss(l)
                        n=3*nn
                        call typen(l,nout,1)
                        l=l+1
                        do k=1,n
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                     else if (law.eq.61) then
                        if (nout.ne.1) nr=nint(xss(l))
                        if (nout.eq.1) nr=iss(l)
                        call typen(l,nout,1)
                        l=l+1
                        if (nr.ne.0) then
                           n=2*nr
                           do j=1,n
                              call typen(l,nout,1)
                              l=l+1
                           enddo
                        endif
                        if (nout.ne.1) ne=nint(xss(l))
                        if (nout.eq.1) ne=iss(l)
                        call typen(l,nout,1)
                        l=l+1
                        do j=1,ne
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                        do j=1,ne
                           call typen(l,nout,1)
                           l=l+1
                        enddo
                        do j=1,ne
                           call typen(l,nout,1)
                           l=l+1
                           if (nout.ne.1) np=nint(xss(l))
                           if (nout.eq.1) np=iss(l)
                           call typen(l,nout,1)
                           l=l+1
                           n=3*np
                           do k=1,n
                              call typen(l,nout,2)
                              l=l+1
                           enddo
                           do k=1,np
                              call typen(l,nout,1)
                              l=l+1
                           enddo
                           do k=1,np
                              call typen(l,nout,1)
                              l=l+1
                              if (nout.ne.1) nmu=nint(xss(l))
                              if (nout.eq.1) nmu=iss(l)
                              call typen(l,nout,1)
                              l=l+1
                              nw=3*nmu
                              do kk=1,nw
                                 call typen(l,nout,2)
                                 l=l+1
                              enddo
                           enddo
                        enddo
                     else if (law.eq.67) then
                        do j=1,ne
                           call typen(l,nout,1)
                           l=l+1
                           if (nout.ne.1) nmu=nint(xss(l))
                           if (nout.eq.1) nmu=iss(l)
                           call typen(l,nout,1)
                           l=l+1
                           do k=1,nmu
                              call typen(l,nout,2)
                              l=l+1
                           enddo
                           do k=1,nmu
                              call typen(l,nout,1)
                              l=l+1
                           enddo
                           do k=1,nmu
                              call typen(l,nout,1)
                              l=l+1
                              if (nout.ne.1) nep=nint(xss(l))
                              if (nout.eq.1) nep=iss(l)
                              call typen(l,nout,1)
                              l=l+1
                              nn=3*nep
                              do n=1,nn
                                 call typen(l,nout,2)
                                 l=l+1
                              enddo
                           enddo
                        enddo
                     endif
                  endif
               enddo
c
c              ***yh block
               nyh=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do k=1,nyh
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
         enddo
      endif
      call typen(0,nout,3)
c
      return
      end
C======================================================================
      subroutine typen(l,nout,iflag)
c     ******************************************************************
c     write an integer or a real number to a type-1 ace file
c     or (if nout=0) convert real to integer for type-3 output
c     or (if nout=1) convert integer to real for type-3 input
c     use iflag.eq.1 to write an integer (i20)
c     use iflag.eq.2 to write a real number (1pe20.11)
c     use iflag.eq.3 to write partial line at end of file
c
c     Adapted by D. Lopez Aldama for ACEMAKER
c
c     ******************************************************************
      implicit real*8 (a-h,o-z)
      common/acedat/xss(200000000),nxss
      dimension iss(1)
      equivalence (xss(1),iss(1))
      character*20 hl(4)
      save hl,i
c
      if (iflag.eq.3.and.nout.gt.1.and.i.lt.4) then
         write(nout,'(4a20)') (hl(j),j=1,i)
      else if (nout.eq.0) then
         if (iflag.eq.1) iss(l)=nint(xss(l))
      else if (nout.eq.1) then
         if (iflag.eq.1) xss(l)=iss(l)
      else
         i=mod(l-1,4)+1
         if (iflag.eq.1) write(hl(i),'(i20)') nint(xss(l))
         if (iflag.eq.2) write(hl(i),'(1p,e20.11)') xss(l)
         if (i.eq.4) write(nout,'(4a20)') (hl(j),j=1,i)
      endif
      return
      end
C======================================================================