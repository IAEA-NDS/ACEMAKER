      program dophn
c     version 1.0
c
c     prepare a photonuclear continuous energy ACE-formatted file from
c     ENDF-6 formatted data for MCNP and other Monte Carlo codes.
c
c     To run DOPHN all endf data should be pre-processed by ACEMAKER.
c     The ACEMAKER driver uses LINEAR,LEGEND,SPECTRA,MERGE and DICTIN
c     from PREPRO. The driver also invokes SIXLIN and GAMLIN to get
c     prepared a pendf tape ready for use the module DOPHN.
c
c     DOPHN calculate the heating taking into account the relativistic
c     convertion of the average outgoing energy from the CM to the LAB
c     system if required
c
c     INPUT data for DOPHN (for ACEMAKER see ACEMAKER.FOR):
c
c     Input data options should be entered on the DOPHN.INP text file.
c
c     line 1:       isel       imon      mcnpx      idneu   idos  (4i11)
c     line 2:       tol        ymin                             (2e11.4)
c     line 3:       input PENDF filename                           (a72)
c     line 4:       output ACE-formatted filename                  (a72)
c     line 5:       imat                                           (i11)
c     line 6:       suff                                            (a4)
c
c     where
c       isel: selection criterium      (0/1) = ZA/MAT   Default=0
c       imon: Monitor printing trigger (0/1) = min/max  Default=0
c      mcnpx: MCNP trigger (0/1) = MCNP/MCNPX  Default=0
c      idneu: Delayed neutrons data trigger (-1/0/1) = no/yes/yes+corr.
c             (Default=1)
c       idos: Production=(dosimetry/activation) trigger (0/1) = no/yes
c             (default=0)
c        tol: Linearization/reconstruction tolerance
c       ymin: Minimum value for linearization/reconstruction
c       imat: Selected material number
c       suff: ZAID suffix for ACE-formatted file
c             (Examples: .00, .32, .80, .067, the dot '.' is required)
c
c     Example of DOPHN.INP:
c
c               0          0          0
c           0.001    1.0E-30
c     \PENDF\U235.PENDF
c     \MC\U235.ACEF
c            9228
c     .00
c
c     Retrieve material 9228 (U-235) from PENDF tape \PENDF\U235.PENDF
c     and generate ACE-formatted file U235.ACEF on \MC\ sub-directory.
c     ZAID for U235 should be 92235.00c and minimum printout will be
c     produced during DOPHN processing. Use a tolerance of 0.1% and a
c     minimum value of 1.0E-30 for reconstruction and linearization.
c
      implicit real*8 (a-h, o-z)
      parameter (nnxc=450,npmx=7,ndos=1000, nprodmax=504)
      parameter (npmax=2000000,nbmax=2000000,nxsmax=200000000,nnux=5000)
      parameter (emnc2=939.56542052539d+6, bk=8.6173303d-11)
      parameter (ev2mev=1.0d-6, ktvn=1)
      character*1 line1(80)
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
      dimension nbt(20),ibt(20),nbt1(20),ibt1(20),nbt2(20),ibt2(20)
      dimension mf3(nnxc),mf4(nnxc),mf5(nnxc),mf6(nnxc),mf86(ndos)
      dimension mf10(nnxc),mf12(nnxc),mf14(nnxc),ethr0(nnxc),n6ry(nnxc)
      dimension kfs86(ndos),lzap86(ndos),lfs86(ndos)
      dimension lfnd(npmx),nmtp(npmx),iptype(npmx),jzap(npmx),thrp(npmx)
      dimension xnu(nnux),ynu(nnux)
      dimension x(npmax),x1(npmax),x2(npmax),xe(npmax)
      dimension y(npmax),y1(npmax),y2(npmax),ye(npmax)
      dimension b(nbmax)
      allocatable eek(:),avepk(:),x6(:),y6(:)
      data line1/80*'-'/,suff/'    '/
      data in1/2/,nin/3/,iou/11/,nou/12/,nscr/20/,nfi40/40/,nfi60/60/
      data iptype/1,2,9,31,32,33,34/
      data jzap/1,0,1001,1002,1003,2003,2004/
c
c      open input file DOPHN.INP and list file DOPHN.LST
c
      open (in1, file='DOPHN.INP')
      open (iou, file='DOPHN.LST')
      write(*,*)' ============================================='
      write(*,*)' PROGRAM DOPHN: Prepare photonuclear ACE-files'
      write(*,*)' ============================================='
      write(*,*)
      write(iou,*)' ============================================='
      write(iou,*)' PROGRAM DOPHN: Prepare photonuclear ACE-files'
      write(iou,*)' ============================================='
      write(iou,*)
c
c      read input data from DOPHN.INP
c
      isel=0
      imon=0
      mcnpx=0
      read(in1,'(5i11)')isel,imon,mcnpx,idneu,idos
      if (isel.ne.0) isel=1
      if (imon.ne.0) imon=1
      if (mcnpx.ne.0)mcnpx=1
      if (idneu.lt.0) then
        idneu=-1
      elseif (idneu.gt.0) then
        idneu=1
      endif
      if (idos.le.0) then
        idos=0
      elseif (idos.gt.2) then
        idos=1
      endif
      read(in1,'(2e11.0)')tol,ymin
      read(in1,'(a72)')fin1
      read(in1,'(a72)')fout
      read(in1,'(i11)')nsel
      if (nsel.lt.1) nsel=0
      read(in1,'(a)')suff
      close(in1)
c
c      Printing input data
c
      open (nin, file=fin1)
      write(iou,*)' Input parameters'
      write(iou,'(80a1)')line1
      write(iou,*)' MAT/ZA selection        =',isel
      write(iou,*)' Printing option         =',imon
      write(iou,*)' MCNP trigger            =',mcnpx
      write(iou,*)' Delayed neutrons option =',idneu
      write(iou,*)' Dosimetry data option   =',idos
      write(iou,*)' Tolerance               =',tol
      write(iou,*)' Minimun value           =',ymin
      write(iou,*)' Input PENDF-file name   = ',fin1
      write(iou,*)' Output ACE-file name    = ',fout
      if (nsel.eq.0) then
        write(iou,*)' Selected material       = first material'
      else
        write(iou,*)' Selected material       =',nsel
      endif
      write(iou,'(a,a)')' Suffix                  = ',suff
      write(iou,'(80a1)')line1
c
c      Initialize pointers arrays for ACE-formatted file
c
      do i=1,16
        nxs(i)=0
        izn(i)=0
        awn(i)=0.0d0
      enddo
      nxs(16)=ktvn
      do i=1,32
        jxs(i)=0
      enddo
      do i=1,70
        hk(i:i)=' '
      enddo
      nxss=nxsmax
      do i=1,nxss
        xss(i)=0.0d0
      enddo
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
      if (nsub.eq.0.and.awi.eq.0.0d0) then
        zai=0.0d0
        izai=0
      else
        write(iou,*)'  === Error: non a photo-nuclear library'
        write(iou,*)'  === NSUB=',nsub,' IPART=',(nsub/10),' AWI=',awi
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
        write(hz,'(i6,a4,a3)')izaid,suff(1:4),'nu '
      else
        write(hz,'(i6,a3,a4)')izaid,suff(1:3),'u   '
      endif
      write(hm,'(a6,i4)')'   mat',mat0
      tz=bk*temp
      write(*,*)' Material=',mat0,' found'
      write(iou,*)
      write(iou,'(a,i5,a,i7,a,1p,e15.8,a,e15.8)')' Material=',mat0,
     &  ' ZA=',matza,' AWR=',awr0,' ELIS [eV]=',elis
      write(iou,'(a,a,a,i4,a,i7,a,1pe15.8)')' SYM=',zsymam,' LFI=',lfi,
     &  ' ZAI=',izai,' AWI=',awi
      write(iou,'(a,1p,e15.8,a,e13.6,a,e13.6,a)')' EMAX [eV]=', emax,
     &  ' Temperature=',temp,' K = ',tz,' MeV'
      nxs(2)=matza
      nxs(9)=liso
      nxs(10)=matza/1000
      nxs(11)=mod(matza,1000)
c
c      Open ACE summary information file for the final user
c
      i=index(fout,'.',.true.)
      if (i.le.0) then
        i=len(trim(fout))
      else
        i=i-1
      endif
      write(fin1,'(a,a)')trim(fout(1:i)),'.txt'
      fin1=trim(fin1)
      open(in1,file=fin1)
      write(in1,*)
      write(in1,'(80a1)')line1
      write(in1,'(a)')trim(hk)
      write(in1,'(3a,i4,a,i6)')' ZAID=',trim(hz),'  MAT=',mat0,
     &  '  ZA=',matza
      write(in1,'(a,f11.6,a,i3,a,f8.3,a,2(a,i3))')
     &  ' AWR=',awr0,'  STA=',nint(sta),'  ELIS=',elis*ev2mev,' MeV',
     &  '  LIS=',lis,'  LISO=',liso
      write(in1,'(80a1)')line1
c
c      Explore ENDF content (MF3,MF4,MF5,MF8,MF6,MF12,MF14)
c      Prepare a union grid for incident energies from MF3
c
      nmf3=0
      call findmf(nin,mat0,3,icod)
      if (icod.eq.0) then
c
c       Get first cross section data
c
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat1,mf,mt,ns1)
        call readtab1(nin,c1,c2,l1,lr,nr,ne1,nbt,ibt,xe,ye)
        nmf3=nmf3+1
        mf3(nmf3)=mt
        call getthr(ne1,xe,ye)
        call checkdis(iou,ne1,xe,ye,icod)
        if (icod.gt.0) then
          write(iou,*)' Warning: MF3/MT=',mt,', ',icod,' points removed'
        endif
        ethr0(nmf3)=xe(1)
        np2=ne1
        do ie=1,ne1
          x2(ie)=xe(ie)
        enddo
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat1,mf,mt,ns1)
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat1,mf,mt,ns1)
c
c       Explore MF3
c
        do while (mf.eq.3)
          call readtab1(nin,c1,c2,l1,lr,nr,np,nbt,ibt,x,y)
          ival=mtvalid(izai,mt)
          if (ival.gt.0) then
            nmf3=nmf3+1
            mf3(nmf3)=mt
            call getthr(np,x,y)
            call checkdis(iou,np,x,y,icod)
            if (icod.gt.0) then
              write(iou,*)' Warning: MF3/MT=',mt,', ',icod,
     &          ' points removed'
            endif
            ethr0(nmf3)=x(1)
            call union(np,x,np2,x2,ne1,xe,npmax)
            np2=ne1
            do ie=1,ne1
              x2(ie)=xe(ie)
            enddo
          endif
          call readcont(nin,c1,c2,l1,l2,n1,n2,mat1,mf,mt,ns1)
          call readcont(nin,c1,c2,l1,l2,n1,n2,mat1,mf,mt,ns1)
        enddo
        write(iou,*)' MF3 sections:',nmf3
      endif
c
c      Explore MF4
c
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
        write(iou,*)' MF4 sections:',nmf4
      endif
c
c      Explore MF5
c
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
        write(iou,*)' MF5 sections:',nmf5
      endif
c
c      Search for MF8 data (products with yields on MF6)
c
      nmf86=0
      if (idos.gt.0) then
        call findmf(nin,mat0,8,icod8)
        if (icod8.eq.0) then
          backspace(nin)
          mt=1000
          do while (mt.gt.0)
            call findnextmt(nin,8,mt)
            ival=mtvalid(izai,mt)
            if (ival.gt.0) then
              call readcont(nin,c1,c2,l1,l2,nfs,n0,mat1,mf1,mt1,ns1)
              kzap=-99999
              do k=1,nfs
                if (n0.eq.0) then
                  call readlist(nin,zap,ee,lmf,lfs,nw,n2,b)
                else
                  call readcont(nin,zap,ee,lmf,lfs,nw,n2,
     &              mat1,mf1,mt1,ns1)
                endif
                lzap=nint(zap)
                if (lmf.eq.6.and.lzap.gt.2004) then
                  nmf86=nmf86+1
                  mf86(nmf86)=mt
                  lzap86(nmf86)=lzap
                  lfs86(nmf86)=lfs
                  if (kzap.ne.lzap) then
                    kk=0
                    kzap=lzap
                  else
                    kk=kk+1
                  endif
                  kfs86(nmf86)=kk
                endif
              enddo
            endif
          enddo
          write(iou,*)' Number of production cross-sections',
     &      ' taken from MF6/MF8:',nmf86
        endif
      endif
c
c      Explore MF6 yields
c
      idos6=0
      if (idos.gt.0.and.icod8.ne.0) idos6=1
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
            n6ry(nmf6)=0
            m3=iposm(nmf3,mf3,mt)
            if (m3.gt.0) then
              thresh=ethr0(m3)
              call readcont(nin,c1,awr,jp6,lct,nk,n2,mat1,mf1,mt1,ns1)
              kzap=-99999
              do k=1,nk
                call readtab1(nin,zap,awp,lip,law,nr,np,nbt,ibt,x,y)
                izap=nint(zap)
                if (izap.le.2004.or.idos6.gt.0) then
                  call getthr(np,x,y)
                  call checkdis(iou,np,x,y,icod)
                  if (icod.gt.0) then
                    write(iou,*)' Warning: MF6/MT=',mt,' ZAP=',
     &                nint(zap),' LIP=',lip,', ',icod,' points removed'
                  endif
                  i0=np
                  do i=1,np-1
                    if (x(i).ge.thresh) then
                      i0=i
                      exit
                    endif
                  enddo
                  np=np-i0+1
                  if (x(i0).lt.thresh) x(i0)=thresh
                  call union(np,x(i0),np2,x2,ne1,xe,npmax)
                  np2=ne1
                  do ie=1,ne1
                    x2(ie)=xe(ie)
                  enddo
                  if (idos6.gt.0.and.izap.gt.2004) then
                    nmf86=nmf86+1
                    mf86(nmf86)=mt
                    lzap86(nmf86)=izap
                    lfs86(nmf86)=lip
                    if (kzap.ne.izap) then
                      kk=0
                      kzap=izap
                    else
                      kk=kk+1
                    endif
                    kfs86(nmf86)=kk
                  endif
                endif
                call nextsub6(nin,law,nbt,ibt,x,b)
              enddo
            endif
          endif
        enddo
        write(iou,*)' MF6 sections:',nmf6
        if (idos6.gt.0) then
          write(iou,*)' Warning: Number of production cross-sections',
     &      ' taken from MF6(No MF8 data):',nmf86
        endif
      endif
c
c       Search for MF10 production reactions
c
      nmf10=0
      ndmf10=0
      if (idos.gt.0) then
        call findmf(nin,mat0,10,icod)
        if (icod.eq.0) then
          backspace(nin)
          mt=1000
          do while (mt.gt.0)
            call findnextmt(nin,10,mt)
            ival=mtvalid(izai,mt)
            if (ival.gt.0) then
              nmf10=nmf10+1
              mf10(nmf10)=mt
              m3=iposm(nmf3,mf3,mt)
              if (m3.gt.0) then
                thresh=ethr0(m3)
              else
                thresh=1.0d-5
              endif
              call readcont(nin,c1,c2,l1,l2,nfs,n2,mat1,mf1,mt1,ns1)
              ndmf10=ndmf10+nfs
              do k=1,nfs
                call readtab1(nin,c1,c2,izap,lfs,nr,np,nbt,ibt,x,y)
                call getthr(np,x,y)
                call checkdis(iou,np,x,y,icod)
                if (icod.gt.0) then
                  write(iou,*)' Warning: MF10/MT=',mt,' ZAP=',izap,
     &              ' LFS=',lfs,', ',icod,' points removed'
                endif
                i0=np
                do i=1,np-1
                  if (x(i).ge.thresh) then
                    i0=i
                    exit
                  endif
                enddo
                np=np-i0+1
                if (x(i0).lt.thresh) x(i0)=thresh
                call union(np,x(i0),np2,x2,ne1,xe,npmax)
                np2=ne1
                do ie=1,ne1
                  x2(ie)=xe(ie)
                enddo
              enddo
            endif
          enddo
          write(iou,*)' MF10 sections:',nmf10,' Products:',ndmf10
        endif
      endif
c
c      Explore MF12
c
      nmf12=0
      call findmf(nin,mat0,12,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,12,mt)
          ival=mtvalid(izai,mt)
          if (ival.gt.0) then
            nmf12=nmf12+1
            mf12(nmf12)=mt
          endif
        enddo
        write(iou,*)' MF12 sections:',nmf12
      endif
c
c      Explore MF14
c
      nmf14=0
      call findmf(nin,mat0,14,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,14,mt)
          ival=mtvalid(izai,mt)
          if (ival.gt.0) then
            nmf14=nmf14+1
            mf14(nmf14)=mt
          endif
        enddo
        write(iou,*)' MF14 sections:',nmf14
      endif
      if (nmf3.eq.0.or.(nmf4.eq.0.and.nmf6.eq.0.and.nmf14.eq.0)) then
        write(iou,*)' ERROR: No angular distribution found'
        close(nin)
        close(iou)
        stop
      endif
c
c       Get photo-fission nubar and merge prompt and delayed
c       angle-energy distributions on MF6 if delayed data are available
c
      m18=iposm(nmf3,mf3,18)
      if (m18.gt.0) then
        call photofis(nin,iou,nfi,awi,mat0,nmf4,mf4,nmf5,mf5,
     &    nmf6,mf6,nnu,xnu,ynu,idneu)
        if (nnu.gt.0) then
          thresh=ethr0(m18)
          i0=nnu
          do i=1,nnu-1
            if (xnu(i).ge.thresh) then
              i0=i
              exit
            endif
          enddo
          np=nnu-i0+1
          call union(np,xnu(i0),np2,x2,ne1,xe,npmax)
             write(iou,*)
             write(iou,*)' NU MT=',mf3(m18),' NE1=',ne1
             write(iou,'(1p4e20.10)')(xe(ie),ie=1,ne1)
        endif
      else
        nnu=-1
        nfi=-1
      endif
c
c       Set MF6 distributions for outgoing photons of discrete
c       reactions given on MF12/MF14
c
      if (nmf12.gt.0) then
        call disgam(nin,iou,ngi,mat0,nmf3,mf3,nmf12,mf12,nmf14,mf14)
      else
        ngi=-1
      endif
c
c      check if common grid is monotonic
c
      call checkmon(iou,ne1,xe,icod)
      if (icod.ne.0) then
        write(iou,*)' ERROR: ',icod,' non monotonic energies found'
        write(iou,'(1p5e16.9)')(xe(k),k=1,ne1)
        stop
      endif
c
c      Check use of summa versus partials
c
      mx=iposm(nmf3,mf3,16)
      if (mx.gt.0)call mtchkd(16,nmf3,mf3,nmf4,mf4,nmf6,mf6,nmf14,mf14,
     &    mt16,mt16g)
      mx=iposm(nmf3,mf3,103)
      if (mx.gt.0)call mtchkd(103,nmf3,mf3,nmf4,mf4,nmf6,mf6,nmf14,mf14,
     &    mt103,mt103g)
      mx=iposm(nmf3,mf3,104)
      if (mx.gt.0)call mtchkd(104,nmf3,mf3,nmf4,mf4,nmf6,mf6,nmf14,mf14,
     &    mt104,mt104g)
      mx=iposm(nmf3,mf3,105)
      if (mx.gt.0)call mtchkd(105,nmf3,mf3,nmf4,mf4,nmf6,mf6,nmf14,mf14,
     &    mt105,mt105g)
      mx=iposm(nmf3,mf3,106)
      if (mx.gt.0)call mtchkd(106,nmf3,mf3,nmf4,mf4,nmf6,mf6,nmf14,mf14,
     &    mt106,mt106g)
      mx=iposm(nmf3,mf3,107)
      if (mx.gt.0)call mtchkd(107,nmf3,mf3,nmf4,mf4,nmf6,mf6,nmf14,mf14,
     &    mt107,mt107g)
c
c      Printing mt control parameters
c
      write(iou,'(80a1)')line1
      write(iou,*)' lfi    =',lfi
      write(iou,*)' nfi    =',nfi
      write(iou,*)' nnu    =',nnu
      write(iou,*)' ngi    =',ngi
      write(iou,*)' mt16   =',mt16
      write(iou,*)' mt16g  =',mt16g
      write(iou,*)' mt103  =',mt103
      write(iou,*)' mt103g =',mt103g
      write(iou,*)' mt104  =',mt104
      write(iou,*)' mt104g =',mt104g
      write(iou,*)' mt105  =',mt105
      write(iou,*)' mt105g =',mt105
      write(iou,*)' mt106  =',mt106
      write(iou,*)' mt106g =',mt106
      write(iou,*)' mt107  =',mt107
      write(iou,*)' mt107g =',mt107
      write(iou,*)' nmf3   =',nmf3
      write(iou,*)' nmf4   =',nmf4
      write(iou,*)' nmf5   =',nmf5
      write(iou,*)' nmf6   =',nmf6
      write(iou,*)' nmf8(6)=',nmf86
      write(iou,*)' nmf10  =',nmf10
      write(iou,*)' nmf14  =',nmf14
      write(iou,'(80a1)')line1
c
c     set number of partial cross sections
c     remove redundant, if required
c
      ntr=nmf3
      mx=iposm(nmf3,mf3,1)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,3)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,4)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      nred3=nmf3-ntr
      mx=iposm(nmf3,mf3,16)
      if (mx.gt.0.and.mt16.eq.1.and.mt16g.eq.1) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,103)
      if (mx.gt.0.and.mt103.eq.1.and.mt103g.eq.1) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,104)
      if (mx.gt.0.and.mt104.eq.1.and.mt104g.eq.1) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,105)
      if (mx.gt.0.and.mt105.eq.1.and.mt105g.eq.1) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,106)
      if (mx.gt.0.and.mt106.eq.1.and.mt106g.eq.1) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,107)
      if (mx.gt.0.and.mt107.eq.1.and.mt107g.eq.1) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,201)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,202)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,203)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,204)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,205)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,206)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      mx=iposm(nmf3,mf3,207)
      if (mx.gt.0) then
        ntr=ntr-1
        mf3(mx)=-mf3(mx)
      endif
      nprod=nmf86+ndmf10
      if (idos.gt.0) then
        if (nprod.gt.nprodmax) then
          write(iou,*)
          write(iou,*)' Error: too many production reactions',
     &      ' from MF=10&MF6*MF3. nprod=',nprod,' nprodmax=',nprodmax
          stop
        endif
        ntdr=nmf3-nred3+nprod
      else
        ntdr=ntr
      endif
      nxs(4)=ntdr
c
c      Set final common energy grid
c
      thrh=xe(ne1)
      do i=1,nmf3
        mti=mf3(i)
        if (mti.gt.0.or.(ntdr.gt.ntr.and.mti.ne.-1.and.mti.ne.-3.and.
     &    mti.ne.-4)) then
          if (ethr0(i).lt.thrh) thrh=ethr0(i)
        endif
      enddo
      ie0=1
      do ie=1,ne1
        if(abs(xe(ie)-thrh).le.1.0d-8*thrh) then
          ie0=ie
          exit
        endif
      enddo
c
c     set ESZ block and triggers
c
      write(iou,*)
      write(iou,*)' Common incident energy grid'
      lesz=1
      do ie=ie0,ne1
        xss(lesz+ie-ie0)=xe(ie)
        write(iou,*)' i=',ie-ie0+1,' e=',xe(ie)
      enddo
      ne1=ne1-ie0+1
      emax=xss(lesz+ne1-1)
      nxs(3)=ne1
      jxs(1)=lesz
      write(iou,*)' ESZ block done (ESZ, NES, emin, emax)'
      write(iou,'(i5,i9,1p2e17.9)')lesz,ne1,xss(lesz),emax
      write(iou,'(80a1)')line1
c
c     set block triggers
c
      ltot=lesz+ne1
      iela=iposm(nmf3,mf3,2)
      if (iela.gt.0) then
        lnon=ltot+ne1
        lels=lnon+ne1
        lthn=lels+ne1
      else
        lnon=ltot
        lels=0
        lthn=ltot+ne1
      endif
      lxs=lthn+ne1
      jxs(2)=ltot
      jxs(3)=lnon
      jxs(4)=lels
      jxs(5)=lthn
c
c     set locators for partials
c
      lmtr=lxs
      lqr=lmtr+ntdr
      lsig=lqr+ntdr
      ksig=lsig+ntdr
      lxs=ksig
      jxs(6)=lmtr
      jxs(7)=lqr
      jxs(8)=lsig
      jxs(9)=ksig
c
c      Set TOT, NON, ELS, MTR, LQR, LSIG and SIG blocks from MF3 data
c
      write(iou,*)' Total number of reactions     : ',ntdr
      write(iou,*)' Number of transport reactions : ',ntr
      write(iou,*)' Number of production reactions: ',ntdr-ntr
      write(iou,*)' List of reactions (transport+production)'
      write(iou,*)' ----------------------------------------'
      write(iou,'(9a)')'   i','     MT/MTD','     Q [MeV]',
     &  '     T [MeV]',' LMF','  MT','    ZAP',' LFS','  Warning'
      write(in1,'(a,i4,a,i4,a,i4,a,i8,a,f8.3,a)')' NTOTR=',ntdr,
     &  ' NTRAN=',ntr,' NPROD=',ntdr-ntr,' NE=',ne1,
     &  ' EMAX=',emax*ev2mev,' MeV'
      write(in1,'(8a)')'   i','     MT/MTD','     Q [MeV]',
     &  '     T [MeV]',' LMF','  MT','    ZAP',' LFS'
      k=0
      do i=1,nmf3
        mti=mf3(i)
        if (mti.gt.0) then
          call findmt(nin,mat0,3,mti,icod)
          call readxs(nin,c1,c2,qm,qi,lr,nr,nbt,ibt,np,x,y)
c
c         set energy index for reaction mti
c
          ethri=ethr0(i)
          ie0=1
          do ie=1,ne1
            if(abs(xss(lesz+ie-1)-ethri).le.1.0d-8*ethri) then
              ie0=ie
              exit
            endif
          enddo
          ethri=xss(lesz+ie0-1)
c
c         save reaction parameters in MTR,LQR,LSIG blocks
c
          k=k+1
          xss(lmtr+k-1)=mti
          xss(lqr+k-1)=qi*ev2mev
          xss(lsig+k-1)=lxs-ksig+1
          xss(lxs)=ie0
          lxs=lxs+1
          xss(lxs)=ne1-ie0+1
          lxs=lxs+1
          if (qi.gt.0.0d0.and.ethri.gt.1.0d-5) then
            write(iou,'(i4,i11,2f12.6,2i4,11x,a)')
     &        k,mti,qi*ev2mev,ethri*ev2mev,3,mti,'  Q>0,T>1.0E-11'
          elseif (qi.eq.0.0d0.and.ethri.gt.1.0d-5) then
            write(iou,'(i4,i11,2f12.6,2i4,11x,a)')
     &        k,mti,qi*ev2mev,ethri*ev2mev,3,mti,'  Q=0,T>1.0E-11'
          elseif (qi.lt.0.0d0.and.ethri*1.00000001d0.lt.(-qi)) then
            write(iou,'(i4,i11,2f12.6,2i4,11x,a)')
     &        k,mti,qi*ev2mev,ethri*ev2mev,3,mti,'  Q<0,T<(-Q)'
          else
            write(iou,'(i4,i11,2f12.6,2i4)')
     &        k,mti,qi*ev2mev,ethri*ev2mev,3,mti
          endif
          write(in1,'(i4,i11,2f12.6,2i4)')
     &      k,mti,qi*ev2mev,ethri*ev2mev,3,mti
c
c         save cross section value (xsv)
c         calculate total, save elastic and
c         calculate non elastic, if required
c
          do ie=ie0,ne1
             e=xss(lesz+ie-1)
             xsv=fvalue(nr,nbt,ibt,np,x,y,e)
c
c            SIG block
c
             xss(lxs)=xsv
             lxs=lxs+1
c
c              TOT, NON and ELS blocks
c
             if((mti.eq.16.and.(mt16.eq.1.or.mt16g.eq.1)).or.
     &          (mti.eq.103.and.(mt103.eq.1.or.mt103g.eq.1)).or.
     &          (mti.eq.104.and.(mt104.eq.1.or.mt104g.eq.1)).or.
     &          (mti.eq.105.and.(mt105.eq.1.or.mt105g.eq.1)).or.
     &          (mti.eq.106.and.(mt106.eq.1.or.mt106g.eq.1)).or.
     &          (mti.eq.107.and.(mt107.eq.1.or.mt107g.eq.1))) then
c
c              to protect again possible double-counting if outgoing
c              distributions for neutrons and gammas are not given
c              using the the same set of cross sections
c
               cycle
             else
c
c              store total, elastic and non elastic if required
c
               xss(ltot+ie-1)=xss(ltot+ie-1)+xsv
               if (mti.eq.2) then
                 xss(lels+ie-1)=xsv
               elseif (lnon.ne.ltot) then
                 xss(lnon+ie-1)=xss(lnon+ie-1)+xsv
               endif
             endif
          enddo
        endif
      enddo
c
c      Include production cross sections (idos>0 ==> ntdr>ntr)
c
      if (ntdr.gt.ntr) then
c
c       Production reactions from MF3
c
        if ((ntdr-ntr-nprod).gt.0) then
          do i=1,nmf3
            mti=mf3(i)
            if (mti.lt.-5) then
              mti=abs(mti)
              call findmt(nin,mat0,3,mti,icod)
              call readxs(nin,c1,c2,qm,qi,lr,nr,nbt,ibt,np,x,y)
c
c             set energy index for reaction mti
c
              ethri=ethr0(i)
              ie0=1
              do ie=1,ne1
                if(abs(xss(lesz+ie-1)-ethri).le.1.0d-8*ethri) then
                  ie0=ie
                  exit
                endif
              enddo
              ethri=xss(lesz+ie0-1)
c
c             save reaction parameters in MTR,LQR,LSIG blocks
c
              k=k+1
              xss(lmtr+k-1)=mti
              xss(lqr+k-1)=qi*ev2mev
              xss(lsig+k-1)=lxs-ksig+1
              xss(lxs)=ie0
              lxs=lxs+1
              xss(lxs)=ne1-ie0+1
              lxs=lxs+1
              if (qi.gt.0.0d0.and.ethri.gt.1.0d-5) then
                write(iou,'(i4,i11,2f12.6,2i4,11x,a)')
     &            k,mti,qi*ev2mev,ethri*ev2mev,3,mti,'  Q>0,T>1.0E-11'
              elseif (qi.eq.0.0d0.and.ethri.gt.1.0d-5) then
                write(iou,'(i4,i11,2f12.6,2i4,11x,a)')
     &            k,mti,qi*ev2mev,ethri*ev2mev,3,mti,'  Q=0,T>1.0E-11'
              elseif (qi.lt.0.0d0.and.ethri*1.00000001d0.lt.(-qi)) then
                write(iou,'(i4,i11,2f12.6,2i4,11x,a)')
     &            k,mti,qi*ev2mev,ethri*ev2mev,3,mti,'  Q<0,T<(-Q)'
              else
                write(iou,'(i4,i11,2f12.6,2i4)')
     &            k,mti,qi*ev2mev,ethri*ev2mev,3,mti
              endif
              write(in1,'(i4,i11,2f12.6,2i4)')
     &          k,mti,qi*ev2mev,ethri*ev2mev,3,mti
c
c             save cross section value (xsv)
c
              do ie=ie0,ne1
                e=xss(lesz+ie-1)
                xsv=fvalue(nr,nbt,ibt,np,x,y,e)
c
c               SIG block
c
                xss(lxs)=xsv
                lxs=lxs+1
              enddo
            endif
          enddo
        endif
c
c       Production reactions from MF10
c
        imtd=0
        if (ndmf10.gt.0) then
          do i=1,nmf10
            mti=mf10(i)
            call findmt(nin,mat0,10,mti,icod)
            call readcont(nin,c1,c2,l1,l2,nfs,n2,mat1,mf1,mt1,ns1)
            do kk=1,nfs
              call readtab1(nin,qm,qi,izap,lfs,nr,np,nbt,ibt,x,y)
              mt=mtdos(idos,mti,izap,lfs,imtd)
              ithr0=np
              do ie=1,np-1
                if (y(ie).gt.0.0d0) then
                  ithr0=ie
                  exit
                endif
              enddo
              if (ithr0.gt.1) ithr0=ithr0-1
              ethri=x(ithr0)
              ie0=1
              do ie=1,ne1
                if(abs(xss(lesz+ie-1)-ethri).le.1.0d-8*ethri) then
                  ie0=ie
                  exit
                endif
              enddo
              ethri=xss(lesz+ie0-1)
              k=k+1
              xss(lmtr+k-1)=mt
              xss(lqr+k-1)=qi*ev2mev
              xss(lsig+k-1)=lxs-ksig+1
              xss(lxs)=ie0
              lxs=lxs+1
              xss(lxs)=ne1-ie0+1
              lxs=lxs+1
              if (qi.gt.0.0d0.and.ethri.gt.1.0d-5) then
                write(iou,'(i4,i11,2f12.6,2i4,i7,i4,a)')k,mt,qi*ev2mev,
     &            ethri*ev2mev,10,mti,izap,lfs,'  Q>0,T>1.0E-11'
              elseif (qi.eq.0.0d0.and.ethri.gt.1.0d-5) then
                write(iou,'(i4,i11,2f12.6,2i4,i7,i4,a)')k,mt,qi*ev2mev,
     &            ethri*ev2mev,10,mti,izap,lfs,'  Q=0,T>1.0E-11'
              elseif (qi.lt.0.0d0.and.ethri*1.00000001d0.lt.(-qi)) then
                write(iou,'(i4,i11,2f12.6,2i4,i7,i4,a)')k,mt,qi*ev2mev,
     &            ethri*ev2mev,10,mti,izap,lfs,'  Q<0,T<(-Q)'
              else
                write(iou,'(i4,i11,2f12.6,2i4,i7,i4)')
     &            k,mt,qi*ev2mev,ethri*ev2mev,10,mti,izap,lfs
              endif
              write(in1,'(i4,i11,2f12.6,2i4,i7,i4)')
     &          k,mt,qi*ev2mev,ethri*ev2mev,10,mti,izap,lfs
c
c             save cross section value (xsv)
c
              do ie=ie0,ne1
                e=xss(lesz+ie-1)
                xsv=fvalue(nr,nbt,ibt,np,x,y,e)
c
c               SIG block
c
                xss(lxs)=xsv
                lxs=lxs+1
              enddo
            enddo
          enddo
        endif
c
c        Production reaction from MF6 yields and MF3 cross sections
c
        if (nmf86.gt.0) then
          i=1
          do while (i.le.nmf86)
            mti=mf86(i)
            izap86=lzap86(i)
            lip86=lfs86(i)
            kip86=kfs86(i)
            call findmt(nin,mat0,3,mti,icod)
            call readxs(nin,c1,c2,qm,qi,lr,nr,nbt,ibt,np,x,y)
            call findmt(nin,mat0,6,mti,icod)
            call readcont(nin,za,awr,jp,lct,nk6,n2,mat1,mf1,mt1,ns1)
            do j=1,nk6
              call readtab1(nin,zap,awp,lip,law,nr1,np1,nbt1,ibt1,x1,y1)
              izap=nint(zap)
              if (mt.eq.18.and.izap.eq.1.and.nnu.gt.0) then
                np1=nnu
                do ii=1,nnu
                  x1(ii)=xnu(i)
                  y1(ii)=ynu(i)
                enddo
                nr1=1
                nbt1(1)=nnu
                ibt1(1)=2
              endif
              if (izap.eq.izap86.and.
     &          (lip.eq.lip86.or.lip.eq.kip86))then
                do ie=1,ne1
                   x0=xss(lesz+ie-1)
                   xsv=fvalue(nr,nbt,ibt,np,x,y,x0)
                   yld=fvalue(nr1,nbt1,ibt1,np1,x1,y1,x0)
                   y2(ie)=xsv*yld
                enddo
                ie0=ne1
                do ie=1,ne1-1
                  if (y2(ie).gt.0.0d0) then
                    ie0=ie
                    exit
                  endif
                enddo
                if (ie0.gt.1)ie0=ie0-1
                ethri=xss(lesz+ie0-1)
                mt=mtdos(idos,mti,izap86,lip86,imtd)
                k=k+1
                xss(lmtr+k-1)=mt
                xss(lqr+k-1)=qi*ev2mev
                xss(lsig+k-1)=lxs-ksig+1
                xss(lxs)=ie0
                lxs=lxs+1
                xss(lxs)=ne1-ie0+1
                lxs=lxs+1
                if (qi.gt.0.0d0.and.ethri.gt.1.0d-5) then
                  write(iou,'(i4,i11,2f12.6,2i4,i7,i4,a)')k,mt,
     &              qi*ev2mev,ethri*ev2mev,6,mti,izap86,lip86,
     &              '  Q>0,T>1.0E-11'
                elseif (qi.eq.0.0d0.and.ethri.gt.1.0d-5) then
                  write(iou,'(i4,i11,2f12.6,2i4,i7,i4,a)')k,mt,
     &              qi*ev2mev,ethri*ev2mev,6,mti,izap86,lip86,
     &              '  Q=0,T>1.0E-11'
                elseif (qi.lt.0.0d0.and.ethri*1.00000001d0.lt.(-qi))then
                  write(iou,'(i4,i11,2f12.6,2i4,i7,i4,a)')k,mt,
     &              qi*ev2mev,ethri*ev2mev,6,mti,izap86,lip86,
     &              '  Q<0,T<(-Q)'
                else
                  write(iou,'(i4,i11,2f12.6,2i4,i7,i4)')
     &              k,mt,qi*ev2mev,ethri*ev2mev,6,mti,izap86,lip86
                endif
                write(in1,'(i4,i11,2f12.6,2i4,i7,i4)')
     &            k,mt,qi*ev2mev,ethri*ev2mev,6,mti,izap86,lip86
c
c               save cross section value (xsv)
c
                do ie=ie0,ne1
                  e=xss(lesz+ie-1)
c
c                 SIG block
c
                  xss(lxs)=y2(ie)
                  lxs=lxs+1
                enddo
                i=i+1
                if (mf86(i).eq.mti) then
                  izap86=lzap86(i)
                  lip86=lfs86(i)
                  kip86=kfs86(i)
                else
                  i=i-1
                  exit
                endif
              endif
              call nextsub6(nin,law,nbt1,ibt1,x1,b)
            enddo
            i=i+1
          enddo
        endif
      endif
      write(iou,*)' TOT, NON, ELS, MTR, LQR, LSIG and SIG blocks done'
      write(iou,'(80a1)')line1
c
c      explore distributions on MF4, MF6 & MF14 to count the different
c      particles produced, production thresholds, and to accumulate the
c      heating from recoils (MF6)
c
c       1  iptype=1   jzap=1      neutron
c       2  iptype=2   jzap=0      photon
c       3  iptype=9   jzap=1001   proton
c       4  iptype=31  jzap=1002   deuteron
c       5  iptype=32  jzap=1003   triton
c       6  iptype=33  jzap=2003   He-3
c       7  iptype=34  jzap=2004   alpha(He-4)
c
      do i=1,npmx
        nmtp(i)=0
        thrp(i)=emax
      enddo
c
c      search for neutron production on MF4
c      MF4 can be used by exception in some old photonuclear libraries
c      for describing the neutron angular distribution
c      The recommendation for photo-nuclear data is to use MF6.
c
      if (nmf4.gt.0) then
        do i=1,nmf4
          mti=mf4(i)
          nn=numpart(1,mti)
          if (nn.gt.0) then
            nmtp(1)=nmtp(1)+1
            mtt=0
            k=0
            do while (mtt.ne.mti.and.k.le.ntr)
              k=k+1
              mtt=nint(xss(lmtr+k-1))
            enddo
            if (mtt.eq.mti) then
              ik=nint(xss(lsig+k-1))+ksig-1
              ie=nint(xss(ik))
              thresh=xss(lesz+ie-1)
              if (thresh.lt.thrp(1)) thrp(1)=thresh
            else
              write(iou,*)' ERROR: MF4 inconsistency found for MT=',mti
              stop
            endif
          else
            write(iou,*)' ERROR: MT=',mti,' not allowed on MF4 (NLIB=0)'
            stop
          endif
        enddo
      endif
c
c      search for discrete gamma production on MF12/MF14 ==> MF6
c      It is also an exceptional feature used by IAEA/PD-2019 library
c      MF6 representation is the ENDF-6 format recommendation.
c
      if (nmf14.gt.0) then
        do i=1,nmf14
          mti=mf14(i)
          nmtp(2)=nmtp(2)+1
          mtt=0
          k=0
          do while (mtt.ne.mti.and.k.le.ntr)
            k=k+1
            mtt=nint(xss(lmtr+k-1))
          enddo
          if (mtt.eq.mti) then
            ik=nint(xss(lsig+k-1))+ksig-1
            ie=nint(xss(ik))
            thresh=xss(lesz+ie-1)
            if (thresh.lt.thrp(2)) thrp(2)=thresh
          else
            write(iou,*)' ERROR: MF14 inconsistency found for MT=',mti
            stop
          endif
        enddo
      endif
c
c      search for particles on MF6
c
      if (nmf6.gt.0) then
        do i=1,nmf6
          mti=mf6(i)
c
c         Find reaction MT parameters in the xss array
c
          mtt=0
          ii=0
          do while (mtt.ne.mti.and.ii.le.ntr)
            ii=ii+1
            mtt=nint(xss(lmtr+ii-1))
          enddo
          if (mtt.eq.mti) then
            q=xss(lqr+ii-1)/ev2mev
            iik=nint(xss(lsig+ii-1))+ksig-1
            ie0=nint(xss(iik))
            thresh=xss(lesz+ie0-1)
          else
            write(iou,*)' ERROR: MF6 inconsistency found for MT=',mti
            stop
          endif
c
c         Explore MF6 for MT and recoil nucleus given by LAW=4
c
          do j=1,npmx
            lfnd(j)=0
          enddo
          call findmt(nin,mat0,6,mti,icod)
          call readcont(nin,c1,awr,jp6,lct,nk,n2,mat1,mf1,mt1,ns1)
          if (jp6.ne.0.and.mti.eq.18) then
            write(iou,*)' JP=',jp6,' not coded for photo-fission (MT18)'
            stop
          endif
          j2bd=id2body(mti)
          ilaw2=-1
          do k=1,nk
            call readtab1(nin,zap,awp,lip,law,nr,np,nbt,ibt,x,y)
            izap=nint(zap)
            if (izap.le.2004) then
c
c             particles: izap = [1,0,1001,1002,1003,2003,2004]
c
              j=iposm(npmx,jzap,izap)
              if (j.gt.0) then
                lfnd(j)=lfnd(j)+1
                if (lfnd(j).eq.1) nmtp(j)=nmtp(j)+1
                if (izap.eq.1.and.mti.eq.18.and.nnu.gt.0) then
                  np=nnu
                  do iy=1,nnu
                    x(iy)=xnu(iy)
                    y(iy)=ynu(iy)
                  enddo
                  nr=1
                  nbt(1)=nnu
                  ibt(1)=2
                endif
                ee=x(np-1)
                do iy=2,np
                  if (y(iy).gt.0.0d0) then
                    ee=x(iy-1)
                    exit
                  endif
                enddo
                if (ee.lt.thresh) ee=thresh
                if (ee.lt.thrp(j)) thrp(j)=ee
              endif
              if (law.eq.2.and.j2bd.gt.0.and.k.lt.nk) then
c
c               Calculate in advance the recoil heating for two body
c               reactions to be used if LAW=4(two-body recoil) is given
c
c               The treatment should be relativistic for photo-nuclear
c               data, but it has not been yet implemented in MCNP.
c               Here a rather crude approximation is applied, which is
c               based on incident neutron kinematics. It could work
c               better for massive particles at low energies.
c
                ilaw2=1
                lep=1
                nd=1
                nep=1
                awp4=awr0-awp
                izap4=matza+izai-izap
                call readtab2(nin,c1,c2,l1,l2,nr2,ne,nbt2,ibt2)
                np2=ne
                do iz=1,np2
                  call readlist(nin,c1,ee,lang,l2,na,nmu,b)
                  call recoil(lang,nmu,b)
                  do jl=na,1,-1
                    b(jl+2)=b(jl)
                  enddo
                  b(1)=awp*(ee+q)/awr0
                  b(2)=1.0d0
                  x2(iz)=ee
                  y2(iz)=aveep(ee,izai,awi,matza,awr0,izap4,awp4,
     &                   lct,lang,lep,nd,na,nep,b)
                enddo
              elseif (law.eq.3.and.j2bd.gt.0.and.k.lt.nk) then
                ilaw2=-1
              elseif (law.eq.4) then
                n6ry(i)=n6ry(i)+1
              else
                call nextsub6(nin,law,nbt,ibt,x,b)
              endif
            else
c
c             check for heavy recoil nucleus (zap>2004)
c             heating contribution for recoil nucleus is calculated
c             for allowed LAW=4(two-body) or LAW=1(continuous)
c
              if (law.eq.4.or.law.eq.1) then
                if (law.eq.4.and.ilaw2.le.0) then
                  write(iou,*)' Warning: MF6/MT=',mti,' ZAP=',awp,
     &              ' LAW=4 discrete two body recoil found without',
     &              ' corresponding LAW=2 for the outgoing particle'
                  write(iou,*)' Isotropic distribution used (LAW=3)',
     &              ' for outgoing particle and recoil nucleus'
                  lang=0
                  lep=1
                  nd=1
                  na=0
                  nep=1
                  np2=ne1-ie0+1
                  do ie=ie0,ne1
                    ee=xss(lesz+ie-1)
                    iz=ie-ie0+1
                    b(1)=(awr0-awp)*(ee+q)/awr0
                    b(2)=1.0d0
                    x2(iz)=ee
                    y2(iz)=aveep(ee,izai,awi,matza,awr0,izap,awp,
     &                     lct,lang,lep,nd,na,nep,b)
                  enddo
                else
                  call readtab2(nin,c1,c2,lang,lep,nr2,np2,nbt2,ibt2)
                  do iz=1,np2
                    call readlist(nin,c1,ee,nd,na,npl,nep,b)
                    x2(iz)=ee
                    y2(iz)=aveep(ee,izai,awi,matza,awr0,izap,awp,
     &                      lct,lang,lep,nd,na,nep,b)
                  enddo
                endif
                do ie=ie0,ne1
                  e=xss(lesz+ie-1)
                  yld=fvalue(nr,nbt,ibt,np,x,y,e)
                  hh=fvalin(np2,x2,y2,e)
                  h=yld*hh*ev2mev*xss(2+iik+ie-ie0)
                  if (xss(ltot+ie-1).ne.0.0d0) then
                    h=h/xss(ltot+ie-1)
                  else
                    h=0.0d0
                  endif
                  xss(lthn+ie-1)=xss(lthn+ie-1)+h
                enddo
                n6ry(i)=n6ry(i)+1
              elseif (law.eq.0) then
                write(iou,*)' LAW=0: No heating info for recoil zap=',
     &           izap,' in MF6/MT=',mti
              else
                write(iou,*)' LAW=',law,' not expected for recoil zap=',
     &            izap,' in MF6/MT=',mti
                call nextsub6(nin,law,nbt,ibt,x,b)
              endif
            endif
          enddo
        enddo
      endif
c
c      add in e+(q=0) heating for capture (MT=102) when recoil is not
c      given in MF6.
c
      m3=iposm(nmf3,mf3,102)
      mrec=0
      if (nmf6.gt.0) then
        m6=iposm(nmf6,mf6,102)
        if (m6.gt.0) mrec=n6ry(m6)
      endif
      if (m3.gt.0.and.mrec.eq.0) then
        do ii=1,ntr
          mtt=nint(xss(lmtr+ii-1))
          if (mtt.eq.102) then
            ik=nint(xss(lsig+ii-1))+ksig-1
            j=nint(xss(ik))
            do ie=j,ne1
              e=xss(lesz+ie-1)*ev2mev
              if (xss(ltot+ie-1).ne.0.0d0) then
                xss(lthn+ie-1)=xss(lthn+ie-1)+
     &            e*xss(2+ik+ie-j)/xss(ltot+ie-1)
              endif
            enddo
            exit
          endif
        enddo
      endif
c
c       add in e+q for MTs in MF6 with no recoil given. Later,
c       the particle heatings are subtracted, just to get the recoil
c       heating and finally the charged particles heating is added to
c       get the kerma values
c
      if (nmf6.gt.0) then
        do ii=1,ntr
          mtt=nint(xss(lmtr+ii-1))
          m6=iposm(nmf6,mf6,mtt)
          if (m6.gt.0) then
            if (m6.gt.0.and.n6ry(m6).eq.0) then
              q=xss(lqr+ii-1)
              ik=nint(xss(lsig+ii-1))+ksig-1
              j=nint(xss(ik))
              do ie=j,ne1
                e=xss(lesz+ie-1)*ev2mev
                if (xss(ltot+ie-1).ne.0.0d0) then
                  xss(lthn+ie-1)=xss(lthn+ie-1)+
     &              (e+q)*xss(2+ik+ie-j)/xss(ltot+ie-1)
                endif
              enddo
            endif
          endif
        enddo
      endif
c
c       prepare ixs array block
c
      write(iou,'(5a)')'   p','   ZAP',' IP','   NR','    Tp [MeV]'
      write(in1,*)
      write(in1,'(5a)')'   p','   ZAP',' IP','   NR','    Tp [MeV]'
      npixs=2
      neixs=12
      ixsa=lxs
      ntype=0
      do i=1,npmx
        if(nmtp(i).gt.0) then
          iptr=ixsa+neixs*ntype
          xss(iptr)=iptype(i)
          xss(iptr+1)=nmtp(i)
          ntype=ntype+1
          write(iou,'(i4,i6,i3,i5,f12.6)')ntype,jzap(i),iptype(i),
     &      nmtp(i),thrp(i)*ev2meV
          write(in1,'(i4,i6,i3,i5,f12.6)')ntype,jzap(i),iptype(i),
     &      nmtp(i),thrp(i)*ev2meV
        endif
      enddo
      lxs=ixsa+neixs*ntype
      ixs=lxs
      nxs(5)=ntype
      nxs(6)=npixs
      nxs(7)=neixs
      jxs(10)=ixsa
      jxs(11)=ixs
      write(iou,*)' IXS block (ntype, npixs, neixs, ixsa, ixs) done'
      write(iou,'(5i10)')ntype,npixs,neixs,ixsa,ixs
      write(iou,'(80a1)')line1
      write(in1,'(80a1)')line1
      write(in1,*)
      close(in1)
c
c      production data block for each particle
c
      write(iou,*)' Production data block by particle type'
      do itype=1,ntype
       ityp=ixsa+neixs*(itype-1)
       ipt=nint(xss(ityp))
       ntrp=nint(xss(ityp+1))
       jpt=iposm(npmx,iptype,ipt)
       ip=jzap(jpt)
       thresh=thrp(jpt)
       it=1
       do while(xss(lesz+it-1).le.thresh.and.it.le.ne1)
         it=it+1
       enddo
       if (it.gt.1) it=it-1
       lpxs=lxs
       xss(ityp+2)=lpxs
       xss(lpxs)=it
       xss(lpxs+1)=ne1-it+1
       lphn=lpxs+2+ne1-it+1
       xss(ityp+3)=lphn
       xss(lphn)=it
       xss(lphn+1)=ne1-it+1
       lmtrp=lphn+2+ne1-it+1
       xss(ityp+4)=lmtrp
       ltyrp=lmtrp+ntrp
       xss(ityp+5)=ltyrp
       lsigp=ltyrp+ntrp
       xss(ityp+6)=lsigp
       ksigp=lsigp+ntrp
       xss(ityp+7)=ksigp
       lxs=ksigp
       jp=0
       write(iou,*)' Block No.=',itype,' ACE-particle=',ipt,' ZAP=',ip
       write(iou,'(80a1)')line1
       write(iou,*)' number of production reactions=',ntrp
       write(iou,*)' ithr=',it,'  production threshold=',xss(lesz+it-1)
       write(iou,*)' lpxs=',lpxs,' lphn=',lphn,' lmtrp=',lmtrp
       write(iou,*)' ltyrp=',ltyrp,' lsigp=',lsigp,' ksigp=',ksigp
       write(iou,*)' Processed production yields by MT'
c
c      Find production cross sections
c
c      neutrons originally given on mf4/mf5
c
       if (ip.eq.1.and.nmf4.gt.0) then
c
c        only neutrons are allowed for mf4/mf5 representation
c
         do i=1,nmf4
           mt=mf4(i)
           mtt=0
           ii=0
           do while (mtt.ne.mt)
             ii=ii+1
             mtt=nint(xss(lmtr+ii-1))
             q=xss(lqr+ii-1)/ev2meV
             ik=nint(xss(lsig+ii-1))+ksig-1
             ie0=nint(xss(ik))
           enddo
           if (mtt.eq.mt) then
             call findmt(nin,mat0,4,mt,icod)
             read(nin,*)line
             read(nin,'(33x,i11)')lct
             jp=jp+1
             xss(lmtrp+jp-1)=mt
             if (lct.eq.1) then
               xss(ltyrp+jp-1)=1
             else
               xss(ltyrp+jp-1)=-1
             endif
             xss(lsigp+jp-1)=lxs-ksigp+1
             if (nmf6.gt.0) then
               m6=iposm(nmf6,mf6,mt)
             else
               m6=0
             endif
             write(iou,*)' MT=',mt,' jp=',jp,' ie0=',ie0,' MF4=yes',
     &         ' MF6=',m6
             if (mt.eq.18) then
c
c              Photo-fission (neutron yields are energy dependent)
c
               ite0=max(ie0,it)
               do ie=ite0,ne1
                 e=xss(lesz+ie-1)
                 yld=fvalin(nnu,xnu,ynu,e)
                 xx=xss(2+ik+ie-ie0)
                 xss(lpxs+2+ie-it)=xss(lpxs+2+ie-it)+yld*xx
                 if (xss(ltot+ie-1).ne.0.0d0.and.m6.le.0) then
                   xss(lthn+ie-1)=xss(lthn+ie-1)+
     &                            xx*ev2mev*(e+q)/xss(ltot+ie-1)
                 endif
               enddo
               xss(lxs)=12
               xss(lxs+1)=mt
               xss(lxs+2)=0
               xss(lxs+3)=nnu
               lxsx=lxs+3
               lxsy=lxsx+nnu
               do ie=1,nnu
                 xss(lxsx+ie)=xnu(ie)*ev2mev
                 xss(lxsy+ie)=ynu(ie)
               enddo
               lxs=lxs+4+2*nnu
             else
c
c              Other production reactions (constant neutron yield)
c
               nn=numpart(1,mt)
               yld=float(nn)
               ite0=max(ie0,it)
               do ie=ite0,ne1
                 e=xss(lesz+ie-1)
                 xx=xss(2+ik+ie-ie0)
                 xss(lpxs+2+ie-it)=xss(lpxs+2+ie-it)+yld*xx
                 if (xss(ltot+ie-1).ne.0.0d0.and.m6.le.0) then
                   xss(lthn+ie-1)=xss(lthn+ie-1)+
     &                            xx*ev2mev*(e+q)/xss(ltot+ie-1)
                 endif
               enddo
               xss(lxs)=12
               xss(lxs+1)=mt
               xss(lxs+2)=0
               xss(lxs+3)=2
               xss(lxs+4)=xss(lesz+ie0-1)*ev2mev
               xss(lxs+5)=xss(lesz+ne1-1)*ev2mev
               xss(lxs+6)=yld
               xss(lxs+7)=yld
               lxs=lxs+8
             endif
           endif
         enddo
       endif
c
c       gammas in mf12/mf14 converted to mf6
c
       if (ip.eq.0.and.nmf14.gt.0) then
c
c        only discrete two body reactions are allowed in mf12/mf14
c
         do i=1,nmf14
           mt=mf14(i)
           mtt=0
           ii=0
           do while (mtt.ne.mt)
             ii=ii+1
             mtt=nint(xss(lmtr+ii-1))
             q=xss(lqr+ii-1)/ev2mev
             ik=nint(xss(lsig+ii-1))+ksig-1
             ie0=nint(xss(ik))
           enddo
           if (mtt.eq.mt) then
             call findmt(ngi,mat0,6,mt,icod)
             call readcont(ngi,c1,c2,l1,l2,n1,n2,mat1,mf1,mt1,ns1)
             jp=jp+1
             xss(lmtrp+jp-1)=mt
             xss(ltyrp+jp-1)=1
             xss(lsigp+jp-1)=lxs-ksigp+1
             if (nmf6.gt.0) then
               m6=iposm(nmf6,mf6,mt)
             else
               m6=0
             endif
             if (nmf4.gt.0) then
               m4=iposm(nmf4,mf4,mt)
             else
               m4=0
             endif
             write(iou,*)' MT=',mt,' jp=',jp,' ie0=',ie0,' MF14=yes',
     &         ' MF4=',m4,' MF6=',m6
             call readtab1(ngi,c1,c2,l1,l2,nr,np,nbt,ibt,x,y)
             ite0=max(ie0,it)
             do ie=ite0,ne1
               e=xss(lesz+ie-1)
               yld=fvalue(nr,nbt,ibt,np,x,y,e)
               xx=xss(2+ik+ie-ie0)
               xss(lpxs+2+ie-it)=xss(lpxs+2+ie-it)+yld*xx
               if (xss(ltot+ie-1).ne.0.0d0.and.mt.ne.102.and.
     &             m6.le.0.and.m4.le.0) then
                 xss(lthn+ie-1)=xss(lthn+ie-1)+
     &                          xx*ev2mev*(e+q)/xss(ltot+ie-1)
               endif
             enddo
             xss(lxs)=6
             xss(lxs+1)=mt
             xss(lxs+2)=0
             xss(lxs+3)=np
             lxsx=lxs+3
             lxsy=lxsx+np
             do ie=1,np
               xss(lxsx+ie)=x(ie)*ev2mev
               xss(lxsy+ie)=y(ie)
             enddo
             lxs=lxs+4+2*np
           endif
         enddo
       endif
c
c       any particle in mf6
c
       if (nmf6.gt.0) then
         do i=1,nmf6
           mt=mf6(i)
           mf12(i)=0
           mtt=0
           ii=0
           do while (mtt.ne.mt)
             ii=ii+1
             mtt=nint(xss(lmtr+ii-1))
             q=xss(lqr+ii-1)/ev2mev
             iik=nint(xss(lsig+ii-1))+ksig-1
             ie0=nint(xss(iik))
             thresh=xss(lesz+ie0-1)
           enddo
           call findmt(nin,mat0,6,mt,icod)
           call readcont(nin,c1,c2,l1,lct,nk,n2,mat1,mf1,mt1,ns1)
           call yld6(nin,iou,b,nnu,xnu,ynu,ip,mt,nk,thresh,
     &               kp,np,x,y,nr,nbt,ibt)
           if (kp.gt.0) then
             jp=jp+1
             xss(lmtrp+jp-1)=mt
             if (lct.eq.1) then
               xss(ltyrp+jp-1)=1
             else
               xss(ltyrp+jp-1)=-1
             endif
             xss(lsigp+jp-1)=lxs-ksigp+1
             write(iou,*)' MT=',mt,' jp=',jp,' ie0=',ie0,' MF6=yes'
             ite0=max(ie0,it)
             do ie=ite0,ne1
               e=xss(lesz+ie-1)
               yld=fvalue(nr,nbt,ibt,np,x,y,e)
               xx=xss(2+iik+ie-ie0)
               xss(lpxs+2+ie-it)=xss(lpxs+2+ie-it)+yld*xx
             enddo
             xss(lxs)=6
             xss(lxs+1)=mt
             xss(lxs+2)=0
             xss(lxs+3)=np
             lxsx=lxs+3
             mf12(i)=lxsx
             lxsy=lxsx+np
             ny0=0
             do ie=1,np
               xss(lxsx+ie)=x(ie)*ev2mev
               xss(lxsy+ie)=y(ie)
               if (y(ie).gt.0.0d0) ny0=ny0+1
             enddo
             lxs=lxs+4+2*np
             if (ny0.eq.0) then
               write(iou,*)'  Warning: Yields are zero at all incident',
     &           ' energies for particle ',ip, '  MT=',mt
             endif
           endif
         enddo
       endif
       write(iou,*)' Particle production cross-section data done'
       write(iou,'(80a1)')line1
c
c      (LANDP and ANDP blocks) angular distributions
c
       kland=lxs
       xss(ityp+8)=kland
       kand=kland+ntrp
       xss(ityp+9)=kand
       lxs=kand
       write(iou,*)' Angular distributions for production MT'
       do i=1,ntrp
         mt=nint(xss(lmtrp+i-1))
         mtt=0
         ii=0
         do while (mtt.ne.mt)
           ii=ii+1
           mtt=nint(xss(lmtr+ii-1))
           q=xss(lqr+ii-1)/ev2mev
           iik=nint(xss(lsig+ii-1))+ksig-1
           ie0=nint(xss(iik))
         enddo
         if (nmf4.le.0.or.(mt.eq.18.and.nfi.gt.0)) then
           m4=0
         else
           m4=iposm(nmf4,mf4,mt)
         endif
         if (nmf14.gt.0) then
           m14=iposm(nmf14,mf14,mt)
         else
           m14=0
         endif
         if (nmf6.gt.0) then
           m6=iposm(nmf6,mf6,mt)
         else
           m6=0
         endif
         if (m6.gt.0) then
           mrec=n6ry(m6)
         else
           mrec=0
         endif
         j2bd=id2body(mt)
         jp=0
         if (ip.eq.1.and.m4.gt.0) then
c
c         neutron angular distribution given in mf4
c
          jp=i
          call findmt(nin,mat0,4,mt,icod)
          call readcont(nin,c1,c2,l1,ltt,n1,n2,mat1,mf1,mt1,ns1)
          call readcont(nin,c1,c2,li,lct,n1,n2,mat1,mf1,mt1,ns1)
          lep=1
          nd=1
          nep=1
          awp=1.0d0
          izap=1
          if (ltt.eq.0.and.li.eq.1) then
c
c           Isotropic distribution on MF4
c
            xss(kland+jp-1)=0
            write(iou,*)' MT=',mt,' jp=',jp,' MF4 isotropic',
     &        ' LCT=',lct,' LTT=',ltt,' MREC=',mrec,' J2BD=',j2bd
            if (j2bd.gt.0) then
c
c             Heating calculation for isotropic discrete two body
c             reactions because they do not need MF5 because the
c             outgoing energies are defined by kinematics.
c
              lang=0
              na=0
              ite0=max(ie0,it)
              do ie=ite0,ne1
                e=xss(lesz+ie-1)
                b(1)=(awr0-awp)*(e+q)/awr0
                b(2)=1.0d0
                h=aveep(e,izai,awi,matza,awr0,izap,awp,
     &             lct,lang,lep,nd,na,nep,b)
                hh=h*ev2mev*xss(2+iik+ie-ie0)
                xss(lphn+2+ie-it)=xss(lphn+2+ie-it)+hh
                if (xss(ltot+ie-1).ne.0.0d0.and.mrec.eq.0) then
                  xss(lthn+ie-1)=xss(lthn+ie-1)-hh/xss(ltot+ie-1)
                endif
              enddo
            endif
          elseif (ltt.eq.2) then
c
c           Angular distribution is given in MF4
c
            write(iou,*)' MT=',mt,' jp=',jp,' MF4 pdf given',
     &        ' LCT=',lct,' LTT=',ltt,' MREC=',mrec,' J2BD=',j2bd
            call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
            call checklaw(nr,ibt,icod)
            if (icod.ne.0.and.icod.ne.1) then
              write(iou,*)' Incident energy grid is not linearly',
     &        ' interpolable. Linear interpolation assumed'
            endif
            xss(kland+jp-1)=lxs-kand+1
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
                write(iou,*)' Run LEGEND/ACEMAKER'
                stop
              endif
              xss(le+j)=e*ev2mev
              if (imon.gt.0) then
                write(iou,*)' Incident energy j=',j,' E=',xss(le+j),
     &            ' MeV'
              endif
              xss(lc+j)=-(lxs-kand+1)
              call cdfcal(np,x,y,jj,y2)
              xss(lxs)=jj
              imu=lxs+1
              xss(imu)=np
              ipdf=imu+np
              icdf=ipdf+np
              if (imon.gt.0) then
                write(iou,*)' Angular distribution intt=',jj,' np=',np
                write(iou,*)'  i  ',' xmu                ',
     &          ' pdf                ',' cdf                '
              endif
              do ii=1,np
                xss(imu+ii)=x(ii)
                xss(ipdf+ii)=y(ii)
                xss(icdf+ii)=y2(ii)
                if (imon.gt.0) then
                  write(iou,'(i4,1p3e20.11)')ii,xss(imu+ii),
     &            xss(ipdf+ii),xss(icdf+ii)
                endif
              enddo
              lxs=lxs+2+3*np
              if (j2bd.gt.0) then
                lang=10+jj
                na=2*np
                b(1)=(awr0-awp)*(e+q)/awr0
                b(2)=1.0d0
                jb=2
                do ii=1,np
                  jb=jb+1
                  b(jb)=x(ii)
                  jb=jb+1
                  b(jb)=y(ii)
                enddo
                xe(j)=e
                ye(j)=aveep(e,izai,awi,matza,awr0,izap,awp,
     &               lct,lang,lep,nd,na,nep,b)
              endif
            enddo
            if (j2bd.gt.0) then
c
c             Heating calculation for non isotropic discrete two body
c             reactions because they do not need MF5 because the
c             outgoing energies are defined by kinematics.
c
              ite0=max(ie0,it)
              do ie=ite0,ne1
                e=xss(lesz+ie-1)
                h=fvalin(ne,xe,ye,e)
                hh=h*ev2mev*xss(2+iik+ie-ie0)
                xss(lphn+2+ie-it)=xss(lphn+2+ie-it)+hh
                if (xss(ltot+ie-1).ne.0.0d0.and.mrec.eq.0) then
                  xss(lthn+ie-1)=xss(lthn+ie-1)-hh/xss(ltot+ie-1)
                endif
              enddo
            endif
          else
            write(iou,*)' ERROR: ltt=',ltt,' only tabular data allowed'
            write(iou,*)' Run LEGEND/ACEMAKER'
            stop
          endif

         elseif (ip.eq.1.and.mt.eq.18.and.nfi.eq.nfi40) then
c
c          photo-fission merged neutron distribution stored on nfi=40
c          using MF6/LAW=1 prepared from MF4/MT18+MF5/MT18 prompt- and
c          MF5/MT455 delayed- distributions
c
           jp=i
           write(iou,*)' MT=',mt,' jp=',jp,' photo-fission neutron',
     &       ' total distribution prepared from MF4/MT18+MF5/MT18',
     &       ' and MF5/MT455 distributions'
         elseif (ip.eq.1.and.mt.eq.18.and.nfi.eq.nfi60) then
c
c          photo-fission merged neutron distribution stored on nfi=60
c          using MF6/LAW=1 prepared from MF6/MT18 prompt- and
c          MF5/MT455 delayed- distributions
c
           jp=i
           write(iou,*)' MT=',mt,' jp=',jp,' photo-fission neutron',
     &       ' total distribution prepared from MF6/MT18 and MF5/MT455',
     &       ' distributions'
         elseif (ip.eq.0.and.m14.gt.0) then
c
c          gamma angular distribution stored on ngi using MF6/LAW=1
c          (ND=NEP, NA=0) from data on MF12/MF14(LI=1)
c
           jp=i
           write(iou,*)' MT=',mt,' jp=',jp,' gamma angular energy',
     &       ' distribution converted from MF12/MF14 to MF6/LAW=1'
         elseif (m6.gt.0) then
c
c          search for angular distributions in mf6 using (LAW=2/3 & 4)
c
           call findmt(nin,mat0,6,mt,icod)
           call readcont(nin,c1,c2,l1,lct,nk,n2,mat1,mf1,mt1,ns1)
           open(nscr,file='law2.tmp')
           ilaw2=-1
           do ik=1,nk
             call readtab1(nin,zap,awp,lip,law,nr,np,nbt,ibt,x,y)
             izap=nint(zap)
             if (izap.le.2004.and.izap.ne.ip.and.j2bd.gt.0.and.
     &          (law.eq.2.or.law.eq.3).and.ik.lt.nk) then
c
c              save data for possible use if two body recoil is given
c
               rewind(nscr)
               ilaw2=1
               ns1=0
               if (law.eq.2) then
                call readtab2(nin,c1,c2,l1,l2,nr2,ne,nbt2,ibt2)
                call wrtab2(nscr,mat1,mf1,mt1,ns1,zap,awp,lip,law,
     &           nr2,ne,nbt2,ibt2)
                 do j=1,ne
                   call readlist(nin,c1,e,lang,l2,nw,np,b)
                   call wrtlist(nscr,mat1,mf1,mt1,ns1,c1,e,lang,l2,
     &               nw,np,b)
                 enddo
               else
                 nr=1
                 call wrtab2(nscr,mat1,mf1,mt1,ns1,zap,awp,lip,law,
     &             nr,np,nbt,ibt)
               endif
             elseif (izap.eq.ip.and.(law.eq.2.or.law.eq.3.or.
     &              (law.eq.4.and.ilaw2.gt.0.and.j2bd.gt.0))) then
c
c              process discrete two body reactions
c
               jp=i
               rewind(nscr)
               if (law.eq.2.or.law.eq.3) then
                 ilaw2=1
                 ns1=0
                 if (law.eq.2) then
                   if (lct.eq.1) then
                     write(iou,*)' Fatal error: Law=2 lct=1 (LAB) for',
     &                 ' MF6/MT=',mt
                     stop
                   endif
                   call readtab2(nin,c1,c2,l1,l2,nr2,ne,nbt2,ibt2)
                   call wrtab2(nscr,mat1,mf1,mt1,ns1,zap,awp,lip,law,
     &             nr2,ne,nbt2,ibt2)
                 else
                   nr=1
                   call wrtab2(nscr,mat1,mf1,mt1,ns1,zap,awp,lip,law,
     &               nr,np,nbt,ibt)
                 endif
               else
                 call readtab2(nscr,za2,aw2,li2,law42,nr2,ne,nbt2,ibt2)
               endif
               lep=1
               nd=1
               nep=1
               if (law.eq.2.or.(law.eq.4.and.law42.eq.2)) then
                 if (j2bd.gt.0) then
                   write(iou,*)' MT=',mt,' jp=',jp,' MF6 pdf given',
     &              ' LCT=',lct,' LAW=',law,' MREC=',mrec,' J2BD=',j2bd
                 else
                   write(iou,*)' MT=',mt,' jp=',jp,' MF6 pdf given',
     &              ' LCT=',lct,' LAW=',law,' MREC=',mrec,' J2BD=',j2bd,
     &              ' non two-body MT!'
                 endif
                 xss(kland+jp-1)=lxs-kand+1
                 le=lxs
                 xss(le)=ne
                 lc=le+ne
                 lxs=lc+ne+1
                 do ie=1,ne
                   if (law.eq.2) then
                     call readlist(nin,c1,e,lang,l2,nw,nl,b)
                     call wrtlist(nscr,mat1,mf1,mt1,ns1,c1,e,lang,l2,
     &                 nw,nl,b)
                   else
                     call readlist(nscr,c1,e,lang,l2,nw,nl,b)
                   endif
                   if (law.eq.4) call recoil(lang,nl,b)
                   if (lang.eq.0) then
                     do jl=nl,1,-1
                       b(jl+1)=b(jl)
                     enddo
                     b(1)=1.0d0
                     call leg2lin(nl,b,nmu,x1,y1,tol,ymin,npmax)
                     nl=nmu
                   else
                     j=0
                     do jl=1,nl
                       j=j+1
                       x1(jl)=b(j)
                       j=j+1
                       y1(jl)=b(j)
                     enddo
                     if (lang.ne.12) then
                      nr1=1
                      nbt1(1)=nl
                      ibt1(1)=lang-10
                      call linear(nr1,nbt1,ibt1,nl,x1,y1,tol,ymin,npmax)
                      do jl=1,nl
                        if (y1(jl).lt.0.0d0) y1(jl)=1.0d-30
                      enddo
                      ilaw=2
                      c=1.0d0
                      call renorm(nl,x1,y1,ilaw,c,fn)
                     endif
                   endif
                   intlaw=2
                   call cdfcal(nl,x1,y1,intlaw,y2)
                   xss(le+ie)=e*ev2mev
                   xss(lc+ie)=-(lxs-kand+1)
                   xss(lxs)=intlaw
                   imu=lxs+1
                   xss(imu)=nl
                   ipdf=imu+nl
                   icdf=ipdf+nl
                   if (imon.gt.0) then
                     write(iou,*)' Angular distribution intt=',intlaw,
     &                 ' np=',nl, ' e=',e
                     write(iou,*)'  i  ',' xmu                ',
     &                 ' pdf                ',' cdf                '
                   endif
                   do ii=1,nl
                     xss(imu+ii)=x1(ii)
                     xss(ipdf+ii)=y1(ii)
                     xss(icdf+ii)=y2(ii)
                     if (imon.gt.0) then
                       write(iou,'(i4,1p3e20.11)')ii,xss(imu+ii),
     &                 xss(ipdf+ii),xss(icdf+ii)
                     endif
                   enddo
                   lxs=lxs+2+3*nl
                   na=2*nl
                   lang=12
                   b(1)=(awr0-awp)*(e+q)/awr0
                   b(2)=1.0d0
                   j=2
                   do ii=1,nl
                     j=j+1
                     b(j)=x1(ii)
                     j=j+1
                     b(j)=y1(ii)
                   enddo
                   xe(ie)=e
                   ye(ie)=aveep(e,izai,awi,matza,awr0,izap,awp,
     &                     lct,lang,lep,nd,na,nep,b)
                 enddo
               else
                 xss(kland+jp-1)=0
                 write(iou,*)' MT=',mt,' jp=',jp,' MF6 isotropic',
     &             ' LCT=',lct,' LAW=',law,' MREC=',mrec,' J2BD=',j2bd
                 na=0
                 lang=0
                 ne=ne1-ie0+1
                 do ie=ie0,ne1
                   e=xss(lesz+ie-1)
                   b(1)=(awr0-awp)*(e+q)/awr0
                   b(2)=1.0d0
                   iz=ie-ie0+1
                   xe(iz)=e
                   ye(iz)=aveep(e,izai,awi,matza,awr0,izap,awp,
     &                     lct,lang,lep,nd,na,nep,b)
                 enddo
               endif
               if (x(1).lt.0.0d0) then
                 x(1)=xss(lesz+ie0-1)
                 if (x(2).le.x(1)) then
                   y(1)=1.0d0
                   x(2)=xss(lesz+ne1-1)
                   y(2)=1.0d0
                   np=2
                   nr=1
                   nbt(1)=2
                   ibt(1)=2
                 endif
               endif
               ite0=max(ie0,it)
               do ie=ite0,ne1
                 e=xss(lesz+ie-1)
                 h=fvalin(ne,xe,ye,e)
                 yld=fvalue(nr,nbt,ibt,np,x,y,e)
                 ss=yld*h*ev2mev*xss(2+iik+ie-ie0)
                 xss(lphn+2+ie-it)=xss(lphn+2+ie-it)+ss
                 if (xss(ltot+ie-1).ne.0.0d0.and.mrec.eq.0) then
                   xss(lthn+ie-1)=xss(lthn+ie-1)-ss/xss(ltot+ie-1)
                 endif
               enddo
             else
               if (izap.eq.ip) then
                 write(iou,*)' MT=',mt,' jp=',i,' MF6 (correlated)',
     &             ' LCT=',lct,' LAW=',law,' MREC=',mrec,' J2BD=',j2bd
               endif
               call nextsub6(nin,law,nbt,ibt,x,b)
             endif
           enddo
           close(nscr,status='DELETE')
         endif
       enddo
       write(iou,*)' LANDP and ANDP blocks done'
       write(iou,'(80a1)')line1
c
c      LDLWP and DLWP blocks (energy distributions)
c
       kldlw=lxs
       xss(ityp+10)=kldlw
       kdlw=kldlw+ntrp
       xss(ityp+11)=kdlw
       lxs=kdlw
       write(iou,*)' Energy or correlated angle-energy distributions',
     &   ' by production MT'
       do i=1,ntrp
         mt=nint(xss(lmtrp+i-1))
         mtt=0
         ii=0
         do while (mtt.ne.mt)
           ii=ii+1
           mtt=nint(xss(lmtr+ii-1))
           q=xss(lqr+ii-1)
           iik=nint(xss(lsig+ii-1))+ksig-1
           ie0=nint(xss(iik))
         enddo
         if (nmf4.le.0.or.(mt.eq.18.and.nfi.gt.0)) then
           m4=0
         else
           m4=iposm(nmf4,mf4,mt)
         endif
         if (nmf14.gt.0) then
           m14=iposm(nmf14,mf14,mt)
         else
           m14=0
         endif
         if (nmf6.gt.0) then
           m6=iposm(nmf6,mf6,mt)
         else
           m6=0
         endif
         if (m6.gt.0) then
           mrec=n6ry(m6)
         else
           mrec=0
         endif
         j2bd=id2body(mt)
         jp=0
         if (ip.eq.1.and.m4.gt.0) then
c
c          outgoing neutron distributions given by mf4/mf5
c
           if (mt.ge.50.and.mt.le.90) then
c
c            outgoing energy distribution for discrete level two body
c            reactions (ACE-LAW=33)
c            heating already calculated
c
             jp=i
             write(iou,*)' MT=',mt,' jp=',jp,' MF4 two-body reaction',
     &         ' LAW=33',' MREC=',mrec,' J2BD=',j2bd
             awp=1.0d0
             xss(kldlw+jp-1)=lxs-kdlw+1
             lxs0=lxs
             xss(lxs)=0
             xss(lxs+1)=33
             lxs=lxs+3
             xss(lxs)=0
             xss(lxs+1)=2
             xss(lxs+2)=xss(lesz+ie0-1)*ev2mev
             xss(lxs+3)=xss(lesz+ne1-1)*ev2mev
             xss(lxs+4)=1.0d0
             xss(lxs+5)=1.0d0
             lxs=lxs+6
             xss(lxs0+2)=lxs-kdlw+1
             xss(lxs)=-q
             xss(lxs+1)=(awr0-awp)/awr0
             lxs=lxs+2
           else
c
c            Process outgoing neutron energy distribution in MF5
c            (ACE-LAW=4)
c
             if (nmf5.gt.0) then
               m5=iposm(nmf5,mf5,mt)
             else
               m5=0
             endif
             if (m5.gt.0) then
               jp=i
               call findmt(nin,mat0,5,mt,icod)
               call readcont(nin,c1,c2,l1,l2,nk,n2,mat1,mf1,mt1,ns1)
               xss(kldlw+jp-1)=lxs-kdlw+1
               izap=1
               awp=1.0d0
               lct=1
               lang=12
               nd=0
               na=0
               do ie=1,ne1
                 ye(ie)=0.0d0
               enddo
               write(iou,*)' MT=',mt,' jp=',jp,' MF5 distribution',
     &         ' nk=',nk,' MREC=',mrec,' J2BD=',j2bd
               do ik=1,nk
                 call readtab1(nin,c1,c2,l1,lf,nr1,np1,nbt1,ibt1,x1,y1)
                 if (lf.ne.1) then
                   write(iou,*)' ERROR: LF=',lf,' not allowed'
                   write(iou,*)' Run SPECTRA'
                   stop
                 endif
                 if ((lxs+2*(nr1+np1)+3).gt.nxss) then
                   write(iou,*)' ERROR: Increase size xss ',nxss
                    stop
                 endif
                 if (ik.gt.1) xss(lxs0)=lxs-kdlw+1
                 lxs0=lxs
                 xss(lxs)=0
                 xss(lxs+1)=4
                 write(iou,*)'   LAW=4 for k=',ik
                 call checklaw(nr1,ibt1,icod)
                 if (icod.ne.0) then
                   j1=lxs+3
                   xss(j1)=nr1
                   j2=j1+nr1
                   if (imon.gt.0) then
                     write(iou,*)'  Probability interpolation nr=',nr1
                     write(iou,*)'  == nbt ==  == ibt == '
                   endif
                   do j=1,nr1
                     xss(j1+j)=nbt1(j)
                     xss(j2+j)=ibt1(j)
                     if (imon.gt.0) then
                       write(iou,'(1x,i10,1x,i10)')nbt1(j),ibt1(j)
                     endif
                   enddo
                   lxs=j2+nr1+1
                 else
                   if (imon.gt.0) then
                     write(iou,*)'  LIN-LIN probability interpolation',
     &                 ' nr=0'
                   endif
                   xss(lxs+3)=0
                   lxs=lxs+4
                 endif
                 xss(lxs)=np1
                 if (imon.gt.0) then
                   write(iou,*)'  Law probability npp=',np1
                   write(iou,*)'  Energy [MeV]              Probability'
                 endif
                 do j=1,np1
                   xss(lxs+j)=x1(j)*ev2mev
                   xss(lxs+np1+j)=y1(j)
                   if (imon.gt.0) then
                     write(iou,'(1p2e21.11)')xss(lxs+j),xss(lxs+np1+j)
                   endif
                 enddo
                 lxs=lxs+2*np1+1
                 xss(lxs0+2)=lxs-kdlw+1
                 call readtab2(nin,c1,c2,l1,l2,nr2,ne,nbt2,ibt2)
                 if((lxs+2*(nr2+ne)+3).gt.nxss)  then
                   write(iou,*)' ERROR: Increase size xss ',nxss
                   stop
                 endif
                 call checklaw(nr2,ibt2,icod)
                 if (icod.eq.0) then
                   if (imon.gt.0) then
                     write(iou,*)'  LIN-LIN energy interpolation'
                   endif
                   xss(lxs)=0
                   lxs=lxs+1
                 else
                   xss(lxs)=nr2
                   if (imon.gt.0) then
                     write(iou,*)'  Energy interpolation nr=',nr2
                     write(iou,*)'  == nbt ==  == ibt == '
                   endif
                   do j=1,nr2
                     ibt2(j)=min(mod(ibt2(j),10),2)
                     xss(lxs+j)=nbt2(j)
                     xss(lxs+nr2+j)=ibt2(j)
                     if (imon.gt.0) then
                       write(iou,'(1x,i10,1x,i10)')nbt2(j),ibt2(j)
                     endif
                   enddo
                   lxs=lxs+2*nr2+1
                 endif
                 xss(lxs)=ne
                 lxsn=lxs+ne
                 lxsd=lxsn+ne+1
                 allocate(eek(ne),avepk(ne))
                 do j=1,ne
                   call readtab1(nin,c1,e,l1,l2,nr,nep,nbt,ibt,x,y)
                   if (nep.gt.npmax) then
                    write(iou,*)' ERROR: Increase TAB1 arrays size ',
     &                npmax
                    stop
                   endif
                   if ((lxs+2+3*nep).gt.nxss) then
                     write(iou,*)' ERROR: Increase size xss ',nxss
                     stop
                   endif
                   xss(lxs+j)=e*ev2mev
                   xss(lxsn+j)=lxsd-kdlw+1
                   call checklaw(nr,ibt,icod)
                   if (icod.lt.0) then
                     write(iou,*)' ERROR: E'' data are not lin-lin or',
     &                 ' constant interpolable for MF5/MT=',mt
                     write(iou,*)' Use SPECTRA'
                     stop
                   endif
                   if (icod.eq.0)icod=2
                   call cdfcal(nep,x,y,icod,y2)
                   xss(lxsd)=icod
                   lxsd=lxsd+1
                   xss(lxsd)=nep
                   if (imon.gt.0) then
                     write(iou,*)'  Incident energy E=',xss(lxs+j),
     &                 ' MeV'
                     write(iou,*)'  E'' intt=',icod,' nep=',nep
                     write(iou,*)'   i  ','   Outgoing energy  ',
     &                 ' pdf                ',' cdf                '
                   endif
                   nep2=nep+nep
                   kk=0
                   do iep=1,nep
                     xss(lxsd+iep)=x(iep)*ev2mev
                     xss(lxsd+nep+iep)=y(iep)/ev2mev
                     xss(lxsd+nep2+iep)=y2(iep)
                     kk=kk+1
                     b(kk)=x(iep)
                     kk=kk+1
                     b(kk)=y(iep)
                     if (imon.gt.0) then
                       write(iou,'(i5,1p3e20.11)')iep,xss(lxsd+iep),
     &                   xss(lxsd+nep+iep),xss(lxsd+nep2+iep)
                     endif
                   enddo
                   lxsd=lxsd+3*nep+1
                   eek(j)=e
                   avepk(j)=aveep(e,izai,awi,matza,awr0,izap,awp,lct,
     &                            lang,icod,nd,na,nep,b)
                 enddo
                 lxs=lxsd
                 do ie=ie0,ne1
                   ee=xss(lesz+ie-1)
                   pe=fvalue(nr1,nbt1,ibt1,np1,x1,y1,ee)
                   ep=fvalin(ne,eek,avepk,ee)
                   ye(ie)=ye(ie)+pe*ep
                 enddo
                 deallocate(eek,avepk)
               enddo
               ite0=max(ie0,it)
               do ie=ite0,ne1
                 ee=xss(lesz+ie-1)
                 if (mt.eq.18) then
                   yld=fvalin(nnu,xnu,ynu,ee)
                 else
                   nyld=numpart(1,mt)
                   yld=float(nyld)
                 endif
                 hh=yld*ye(ie)*ev2mev*xss(2+iik+ie-ie0)
                 xss(lphn+2+ie-it)=xss(lphn+2+ie-it)+hh
                 if (xss(ltot+ie-1).ne.0.0d0.and.mrec.eq.0) then
                   xss(lthn+ie-1)=xss(lthn+ie-1)-hh/xss(ltot+ie-1)
                 endif
               enddo
             endif
           endif
         elseif (ip.eq.1.and.mt.eq.18.and.
     &           (nfi.eq.nfi40.or.nfi.eq.nfi60)) then
c
c          angle-energy distribution for photo-fission neutrons
c          including delayed neutron data
c            prompt neutrons: ACE-LAW=44, if ENDF-6/LAW=1, LANG=2
c                             ACE-LAW=61, otherwise
c           delayed neutrons: ACE-LAW=61, from original MF5/MT455
c                                         converted to MF6/LAW=1,LANG=12
c
           jp=i
           call findmt(nfi,mat0,6,mt,icod)
           call readcont(nfi,c1,c2,l1,lct,nk,n2,mat1,mf1,mt1,ns1)
           if (nfi.eq.nfi40) then
             write(iou,*)' MT=',mt,' jp=',jp,' Merged MF6 distribution',
     &        ' from (MF4/MT18+MF5/MT18) + MF5/MT455',' lct=',lct,
     &        ' nk=',nk,' MREC=',mrec
           else
             write(iou,*)' MT=',mt,' jp=',jp,' Merged MF6 distribution',
     &        ' from MF6/MT18 + MF5/MT455',' lct=',lct,
     &        ' nk=',nk,' MREC=',mrec
           endif
           xss(kland+jp-1)=-1
           xss(kldlw+jp-1)=lxs-kdlw+1
           do ie=1,ne1
             ye(ie)=0.0d0
           enddo
           do ik=1,nk
             call readtab1(nfi,zap,awp,lip,law,nr1,np1,nbt1,ibt1,x1,y1)
             izap=nint(zap)
             if (izap.ne.1) then
               write(iou,*)' ERROR: IZAP=',izap,' not allowed for MT=18'
               write(iou,*)' See ACEMAKER documentation'
               stop
             endif
             if (law.ne.1) then
               write(iou,*)' ERROR: LAW=',law,' not allowed for MT=18'
               write(iou,*)' Run ACEMAKER'
               stop
             endif
             if ((lxs+2*(nr1+np1)+3).gt.nxss) then
               write(iou,*)' ERROR: Increase size xss ',nxss
               stop
             endif
             if (ik.gt.1) xss(lxs0)=lxs-kdlw+1
             lxs0=lxs
             xss(lxs)=0
             lawxs=lxs+1
             ldatxs=lxs+2
             xss(lxs+3)=0
             lxs=lxs+4
             xss(lxs)=np1
             do iy=1,np1
               xss(lxs+iy)=x1(iy)*ev2mev
               xss(lxs+np1+iy)=y1(iy)
             enddo
             lxs=lxs+2*np1+1
             xss(ldatxs)=lxs-kdlw+1
             call readtab2(nfi,c1,c2,lang,lep,nr2,ne,nbt2,ibt2)
             if (lang.eq.2)then
               if (lct.eq.1) then
                 write(iou,*)' Fatal error: Law=1 lang=2 lct=1 (LAB)',
     &             ' for MF6/MT=',mt
                 stop
               endif
               xss(lawxs)=44
             else
               xss(lawxs)=61
               if (lang.eq.11) then
                 intmu=1
               else
                 intmu=2
               endif
             endif
             write (iou,*)'   LAW=',nint(xss(lawxs)), ' LEP=',lep,
     &         ' for k=',ik
             call checklaw(nr2,ibt2,icod)
             if (icod.eq.0) then
               xss(lxs)=0
               lxs=lxs+1
             else
               xss(lxs)=nr2
               do j=1,nr2
                 xss(lxs+j)=nbt2(j)
                 xss(lxs+nr2+j)=min(mod(ibt2(j),10),2)
               enddo
               lxs=lxs+2*nr2+1
             endif
             xss(lxs)=ne
             lxsn=lxs+ne
             lxsd=lxsn+ne+1
             allocate(eek(ne),avepk(ne))
             do j=1,ne
               call readlist(nfi,c1,e,nd,na,nw,nep,b)
               if (nw.gt.nbmax) then
                 write(iou,*)' ERROR: Increase the size of list arrays'
                 write(iou,*)' Set nbmax greater than ',nbmax
                 stop
               endif
               eemev=e*ev2mev
               xss(lxs+j)=eemev
               xss(lxsn+j)=lxsd-kdlw+1
               na2=na+2
               do iep=1,nep
                 i0=na2*(iep-1)
                 x(iep)=b(i0+1)
                 y(iep)=b(i0+2)
               enddo
               eek(j)=e
               avepk(j)=aveep(e,izai,awi,matza,awr0,izap,awp,lct,
     &                  lang,lep,nd,na,nep,b)
               call cdfcal1(nep,nd,x,y,lep,y2)
               if ((lxs+5*nep).gt.nxss) then
                 write(iou,*)' Increase size of XSS array ',nxss
                 stop
               endif
               xss(lxsd)=lep+10*nd
               lxsd=lxsd+1
               xss(lxsd)=nep
               nep2=nep+nep
               nep3=nep2+nep
               nep4=nep3+nep
               do iep=1,nep
                 xss(lxsd+iep)=x(iep)*ev2mev
                 if (iep.gt.nd) then
                   xss(lxsd+nep+iep)=y(iep)/ev2mev
                 else
                   xss(lxsd+nep+iep)=y(iep)
                 endif
                 xss(lxsd+nep2+iep)=y2(iep)
               enddo
               if (lang.eq.2) then
c
c                ACE-LAW=44
c
                 do iep=1,nep
                   if (na.eq.0) then
                     xss(lxsd+nep3+iep)=0.0d0
                     xss(lxsd+nep4+iep)=1.0d-20
                   else
                     i0=na2*(iep-1)
                     xss(lxsd+nep3+iep)=b(i0+3)
                     if (na.eq.2) then
                       xss(lxsd+nep4+iep)=b(i0+4)
                     else
                       epmev=x(iep)*ev2mev
                       aa=bachaa(izai,izap,matza,eemev,epmev)
                       xss(lxsd+nep4+iep)=aa
                     endif
                   endif
                 enddo
                 lxsd=lxsd+5*nep+1
               else
c
c                ACE-LAW=61
c
                 lxscd=lxsd+nep4+1
                 do iep=1,nep
                   xss(lxsd+nep3+iep)=lxscd-kdlw+1
                   i0=na2*(iep-1)+2
                   if (na.gt.0.and.b(i0).gt.0.0d0) then
                     i0=na2*(iep-1)+2
                     if (lang.eq.1) then
                       call leg2lin(na,b(i0),nmu,x,y,tol,ymin,npmax)
                     else
                       im0=i0
                       nmu=na/2
                       do im=1,nmu
                         im0=im0+1
                         x(im)=b(im0)
                         im0=im0+1
                         y(im)=2.0d0*b(im0)*b(i0)
                       enddo
                       if (lang.gt.12) then
                         nr=1
                         nbt(1)=nl
                         ibt(1)=lang-10
                         call linear(nr,nbt,ibt,nmu,x,y,tol,ymin,npmax)
                         do im=1,nmu
                           if (y(im).lt.0.0d0) y(im)=1.0d-30
                         enddo
                         ilaw=2
                         c=b(i0)
                         call renorm(nmu,x,y,ilaw,c,fn)
                       endif
                     endif
                     call cdfcal(nmu,x,y,intmu,y2)
                     xss(lxscd)=intmu
                     lxscd=lxscd+1
                     xss(lxscd)=nmu
                     nmu2=nmu+nmu
                     do im=1,nmu
                       xss(lxscd+im)=x(im)
                       xss(lxscd+nmu+im)=y(im)
                       xss(lxscd+nmu2+im)=y2(im)
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
               endif
             enddo
             lxs=lxsd
             do ie=ie0,ne1
               ee=xss(lesz+ie-1)
               pe=fvalue(nr1,nbt1,ibt1,np1,x1,y1,ee)
               ep=fvalin(ne,eek,avepk,ee)
               ye(ie)=ye(ie)+pe*ep
             enddo
             deallocate(eek,avepk)
           enddo
           ite0=max(ie0,it)
           do ie=ite0,ne1
             ee=xss(lesz+ie-1)
             yld=fvalin(nnu,xnu,ynu,ee)
             hh=yld*ye(ie)*ev2mev*xss(2+iik+ie-ie0)
             xss(lphn+2+ie-it)=xss(lphn+2+ie-it)+hh
             if (xss(ltot+ie-1).ne.0.0d0.and.mrec.eq.0) then
               xss(lthn+ie-1)=xss(lthn+ie-1)-hh/xss(ltot+ie-1)
             endif
           enddo
         elseif (ip.eq.0.and.m14.gt.0) then
c
c          outgoing photon distribution for discrete two body reactions
c          ACE-LAW=61, prepared from original MF12/MF14, converted to
c          MF6/LAW=1,LANG=12
c
           jp=i
           call findmt(ngi,mat0,6,mt,icod)
           call readcont(ngi,c1,c2,l1,lct,nk,n2,mat1,mf1,mt1,ns1)
           write(iou,*)' MT=',mt,' jp=',jp,' Merged MF6 distribution',
     &       ' from MF12/MF14',' lct=',lct,' nk=',nk,' MREC=',mrec,
     &       ' J2BD=',j2bd
           if (nk.ne.1) then
             write(iou,*)' ERROR: nk=',nk,' only nk=1 is allowed'
             stop
           endif
           xss(kland+jp-1)=-1
           call readtab1(ngi,zap,awp,lip,law,nr1,np1,nbt1,ibt1,x1,y1)
           izap=nint(zap)
           if (izap.ne.0) then
             write(iou,*)' ERROR: IZAP=',izap,' not allowed for MT=',mt
             write(iou,*)' Run ACEMAKER'
             stop
           endif
           if (law.ne.1) then
             write(iou,*)' ERROR: LAW=',law,' not allowed for MT=',mt
             write(iou,*)' Run ACEMAKER'
             stop
           endif
           if ((lxs+10).gt.nxss) then
             write(iou,*)' ERROR: Increase size xss ',nxss
             stop
           endif
           xss(kldlw+jp-1)=lxs-kdlw+1
           xss(lxs)=0
           xss(lxs+1)=61
           write(iou,*)'   LAW=61 for k=',nk
           ldatxs=lxs+2
           xss(lxs+3)=0
           xss(lxs+4)=2
           xss(lxs+5)=x1(1)*ev2mev
           xss(lxs+6)=x1(np1)*ev2mev
           xss(lxs+7)=1.0d0
           xss(lxs+8)=1.0d0
           lxs=lxs+9
           xss(ldatxs)=lxs-kdlw+1
           call readtab2(ngi,c1,c2,lang,lep,nr2,ne,nbt2,ibt2)
           if (lang.ne.11.and.lang.ne.12) then
             write(iou,*)' ERROR: LANG=',lang,' not allowed for MT=',mt
             write(iou,*)' Run GAMLIN/ACEMAKER'
             stop
           endif
           call checklaw(nr2,ibt2,icod)
           if (icod.eq.0) then
             xss(lxs)=0
             lxs=lxs+1
           else
             xss(lxs)=nr2
             do j=1,nr2
               xss(lxs+j)=nbt2(j)
               xss(lxs+nr2+j)=min(mod(ibt2(j),10),2)
             enddo
             lxs=lxs+2*nr2+1
           endif
           xss(lxs)=ne
           lxsn=lxs+ne
           lxsd=lxsn+ne+1
           allocate(eek(ne),avepk(ne))
           do j=1,ne
             call readlist(ngi,c1,e,nd,na,nw,nep,b)
             if (nw.gt.nbmax) then
               write(iou,*)' ERROR: Increase b array(list) size ',nbmax
               stop
             endif
             xss(lxs+j)=e*ev2mev
             xss(lxsn+j)=lxsd-kdlw+1
             na2=na+2
             do iep=1,nep
               i0=na2*(iep-1)
               x(iep)=b(i0+1)
               y(iep)=b(i0+2)
             enddo
             eek(j)=e
             avepk(j)=aveep(e,izai,awi,matza,awr0,izap,awp,lct,
     &                lang,lep,nd,na,nep,b)
             call cdfcal1(nep,nd,x,y,lep,y2)
             if ((lxs+5*nep).gt.nxss) then
               write(iou,*)' Increase size of XSS array ',nxss
               stop
             endif
             xss(lxsd)=lep+10*nd
             lxsd=lxsd+1
             xss(lxsd)=nep
             nep2=nep+nep
             nep3=nep2+nep
             nep4=nep3+nep
             lxscd=lxsd+nep4+1
             do iep=1,nep
               xss(lxsd+iep)=x(iep)*ev2mev
               if (iep.gt.nd) then
                 xss(lxsd+nep+iep)=y(iep)/ev2mev
               else
                 xss(lxsd+nep+iep)=y(iep)
               endif
               xss(lxsd+nep2+iep)=y2(iep)
             enddo
             do iep=1,nep
               xss(lxsd+nep3+iep)=lxscd-kdlw+1
               if (na.eq.0) then
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
               else
                 write(iou,*)' Fatal error: Non-isotropic angular',
     &             ' distribution for outgoing photons given in MF14'
                 write(iou,*)' Option not expected in DOPHN for ',
     &             ' discrete two-body reaction MT=',mt
                 stop
               endif
               lxscd=lxscd+3*nmu+1
             enddo
             lxsd=lxscd
           enddo
           lxs=lxsd
           ite0=max(ie0,it)
           do ie=ite0,ne1
             ee=xss(lesz+ie-1)
             yld=fvalue(nr1,nbt1,ibt1,np1,x1,y1,ee)
             ep=fvalin(ne,eek,avepk,ee)
             hh=yld*ep*ev2mev*xss(2+iik+ie-ie0)
             xss(lphn+2+ie-it)=xss(lphn+2+ie-it)+hh
             if (xss(ltot+ie-1).ne.0.0d0.and.mrec.eq.0) then
               xss(lthn+ie-1)=xss(lthn+ie-1)-hh/xss(ltot+ie-1)
             endif
           enddo
           deallocate(eek,avepk)
         elseif (m6.gt.0) then
c
c          any particle given on MF6
c          ACE-LAW=44 is prepared from ENDF-6/LAW=1/LANG=2,
c          ACE-LAW=61 is prepared from ENDF-6/LAW=1/LANG=0,11-15
c          ACE-LAW=33 is prepared for two-body reaction from
c                     ENDF-6/LAW=2,3,4
c
           jp=i
           call findmt(nin,mat0,6,mt,icod)
           call readcont(nin,c1,c2,l1,lct,nk,n2,mat1,mf1,mt1,ns1)
           xss(kldlw+jp-1)=lxs-kdlw+1
           do ie=1,ne1
             ye(ie)=0.0d0
           enddo
           lyt=mf12(m6)
           nty=nint(xss(lyt))
           allocate (x6(nty),y6(nty))
           do iy=1,nty
             x6(iy)=xss(lyt+iy)/ev2mev
             y6(iy)=xss(lyt+nty+iy)
           enddo
           kk=0
           do ik=1,nk
             call readtab1(nin,zap,awp,lip,law,nr1,np1,nbt1,ibt1,x1,y1)
             izap=nint(zap)
             if (ip.eq.izap) then
               if (law.lt.1.or.law.gt.4) then
                 write(iou,*)' ERROR: IZAP=',izap,' with LAW=',law,
     &             ' not recommended for photonuclear data (NLIB=0)'
                 stop
               endif
               if ((lxs+2*(nr1+nty)+3).gt.nxss) then
                 write(iou,*)' ERROR: Increase size xss ',nxss
                 stop
               endif
               if (x1(1).ge.0.0d0) then
                 write(iou,*)' MT=',mt,' jp=',jp,' MF6 distribution',
     &             ' law=',law,' lct=',lct,' k=',ik,' MREC=',mrec,
     &             ' J2BD=',j2bd
               else
                 if (x1(2).gt.xss(lesz+ie0-1)) then
                   write(iou,*)' MT=',mt,' jp=',jp,' MF6 distribution',
     &               ' law=',law,' lct=',lct,' k=',ik,' MREC=',mrec,
     &               ' J2BD=',j2bd
                   write(iou,*)' Negative yield threshold T=',x1(1),
     &               ' found. Set to ',xss(lesz+ie0-1)
                   x1(1)=xss(lesz+ie0-1)
                 else
                   write(iou,*)' MT=',mt,' jp=',jp,' MF6 distribution',
     &               ' law=',law,' lct=',lct,' k=',ik,' MREC=',mrec,
     &               ' J2BD=',j2bd
                   write(iou,*)' Negative yield threshold T=',x1(1),
     &               ' found. Fatal error.'
                   stop
                 endif
               endif
               kk=kk+1
               if (kk.gt.1) xss(lxs0)=lxs-kdlw+1
               lxs0=lxs
               xss(lxs)=0
               lawxs=lxs+1
               ldatxs=lxs+2
               xss(lxs+3)=0
               lxs=lxs+4
               if (ip.eq.1.and.mt.eq.18.and.nnu.gt.0) then
                 ii=isconst(np1,y1)
                 if (ii.gt.0) then
                   np1=nnu
                   do iy=1,nnu
                     x1(iy)=xnu(iy)
                     y1(iy)=ynu(iy)
                   enddo
                   nr1=1
                   nbt1(1)=nnu
                   ibt1(1)=2
                 endif
               endif
               n1=0
               do iy=1,nty
                 xx=x6(iy)
                 x2(iy)=xx
                 yyk=fvalue(nr1,nbt1,ibt1,np1,x1,y1,xx)
                 yy=y6(iy)
                 if (yy.eq.yyk) then
                   n1=n1+1
                   y2(iy)=1.0d0
                 elseif (yy.gt.0.0d0) then
                   y2(iy)=min(yyk/yy,1.0d0)
                 else
                   y2(iy)=0.0d0
                 endif
               enddo
               if (n1.eq.nty) then
                 np1=2
                 x1(1)=x2(1)
                 x1(2)=x2(nty)
                 y1(1)=1.0d0
                 y1(2)=1.0d0
               else
                 np1=nty
                 do iy=1,nty
                   x1(iy)=x2(iy)
                   y1(iy)=y2(iy)
                 enddo
               endif
               nr1=1
               nbt1(1)=np1
               ibt1(1)=2
               xss(lxs)=np1
               do iy=1,np1
                 xss(lxs+iy)=x1(iy)*ev2mev
                 xss(lxs+np1+iy)=y1(iy)
               enddo
               lxs=lxs+2*np1+1
               xss(ldatxs)=lxs-kdlw+1
               if (law.eq.1) then
c
c                ENDF-6 LAW=1 (LANG=2)       ==> ACE-LAW=44
c                ENDF-6 LAW=1 (LANG=1,11-15) ==> ACE-LAW=61
c
                 xss(kland+jp-1)=-1
                 call readtab2(nin,c1,c2,lang,lep,nr2,ne,nbt2,ibt2)
                 if (lang.eq.2)then
                   if (lct.eq.1) then
                     write(iou,*)' Fatal error: Law=1 lang=2 lct=1',
     &                 ' (LAB) for MF6/MT=',mt
                     stop
                   endif
                   xss(lawxs)=44
                 else
                   xss(lawxs)=61
                   if (lang.eq.11) then
                     intmu=1
                   else
                     intmu=2
                   endif
                 endif
                 write (iou,*)'   LAW=',nint(xss(lawxs)), ' LEP=',lep
                 call checklaw(nr2,ibt2,icod)
                 if (icod.eq.0) then
                   xss(lxs)=0
                   lxs=lxs+1
                 else
                   xss(lxs)=nr2
                   do j=1,nr2
                     xss(lxs+j)=nbt2(j)
                     xss(lxs+nr2+j)=min(mod(ibt2(j),10),2)
                   enddo
                   lxs=lxs+2*nr2+1
                 endif
                 xss(lxs)=ne
                 lxsn=lxs+ne
                 lxsd=lxsn+ne+1
                 allocate(eek(ne),avepk(ne))
                 do j=1,ne
                   call readlist(nin,c1,e,nd,na,nw,nep,b)
                   if (nw.gt.nbmax) then
                     write(iou,*)' ERROR: Increase size of list arrays'
                     write(iou,*)' Set nbmax greater than ',nbmax
                     stop
                   endif
                   eemev=e*ev2mev
                   xss(lxs+j)=eemev
                   xss(lxsn+j)=lxsd-kdlw+1
                   na2=na+2
                   do iep=1,nep
                     i0=na2*(iep-1)
                     x(iep)=b(i0+1)
                     y(iep)=b(i0+2)
                   enddo
                   eek(j)=e
                   avepk(j)=aveep(e,izai,awi,matza,awr0,izap,awp,lct,
     &                      lang,lep,nd,na,nep,b)
                   call cdfcal1(nep,nd,x,y,lep,y2)
                   if ((lxs+5*nep).gt.nxss) then
                     write(iou,*)' Increase size of XSS array ',nxss
                     stop
                   endif
                   xss(lxsd)=lep+10*nd
                   lxsd=lxsd+1
                   xss(lxsd)=nep
                   nep2=nep+nep
                   nep3=nep2+nep
                   nep4=nep3+nep
                   do iep=1,nep
                     xss(lxsd+iep)=x(iep)*ev2mev
                     if (iep.gt.nd) then
                       xss(lxsd+nep+iep)=y(iep)/ev2mev
                     else
                       xss(lxsd+nep+iep)=y(iep)
                     endif
                     xss(lxsd+nep2+iep)=y2(iep)
                   enddo
                   if (lang.eq.2) then
c
c                    ACE-LAW=44
c
                     do iep=1,nep
                       if (na.eq.0) then
                         xss(lxsd+nep3+iep)=0.0d0
                         xss(lxsd+nep4+iep)=1.0d-20
                       else
                         i0=na2*(iep-1)
                         xss(lxsd+nep3+iep)=b(i0+3)
                         if (na.eq.2) then
                           xss(lxsd+nep4+iep)=b(i0+4)
                         else
                           epmev=x(iep)*ev2mev
                           aa=bachaa(izai,izap,matza,eemev,epmev)
                           xss(lxsd+nep4+iep)=aa
                         endif
                       endif
                     enddo
                     lxsd=lxsd+5*nep+1
                   else
c
c                    ACE-LAW=61
c
                     lxscd=lxsd+nep4+1
                     do iep=1,nep
                       xss(lxsd+nep3+iep)=lxscd-kdlw+1
                       i0=na2*(iep-1)+2
                       if (na.gt.0.and.b(i0).gt.0.0d0) then
                         if (lang.eq.1) then
                           call leg2lin(na,b(i0),nmu,x,y,tol,ymin,npmax)
                         else
                           im0=i0
                           nmu=na/2
                           do im=1,nmu
                             im0=im0+1
                             x(im)=b(im0)
                             im0=im0+1
                             y(im)=2.0d0*b(im0)*b(i0)
                           enddo
                           if (lang.gt.12) then
                             nr=1
                             nbt(1)=nl
                             ibt(1)=lang-10
                             call linear(nr,nbt,ibt,nmu,x,y,tol,ymin,
     &                         npmax)
                             do im=1,nmu
                               if (y(im).lt.0.0d0) y(im)=1.0d-30
                             enddo
                             ilaw=2
                             c=b(i0)
                             call renorm(nmu,x,y,ilaw,c,fn)
                           endif
                         endif
                         call cdfcal(nmu,x,y,intmu,y2)
                         xss(lxscd)=intmu
                         lxscd=lxscd+1
                         xss(lxscd)=nmu
                         nmu2=nmu+nmu
                         do im=1,nmu
                           xss(lxscd+im)=x(im)
                           xss(lxscd+nmu+im)=y(im)
                           xss(lxscd+nmu2+im)=y2(im)
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
                   endif
                 enddo
                 lxs=lxsd
                 ite0=max(ie0,it)
                 do ie=ite0,ne1
                   ee=xss(lesz+ie-1)
                   pe=fvalue(nr1,nbt1,ibt1,np1,x1,y1,ee)
                   ep=fvalin(ne,eek,avepk,ee)
                   yld=fvalin(nty,x6,y6,ee)
                   hh=yld*pe*ep*ev2mev*xss(2+iik+ie-ie0)
                   xss(lphn+2+ie-it)=xss(lphn+2+ie-it)+hh
                   if (xss(ltot+ie-1).ne.0.0d0.and.mrec.eq.0) then
                     xss(lthn+ie-1)=xss(lthn+ie-1)-hh/xss(ltot+ie-1)
                   endif
                 enddo
                 deallocate(eek,avepk)
               else
c
c                ENDF-6 LAW=2,3,4 => ACE-LAW=33 for photonuclear data
c                Heating already calculated
c
                 if (j2bd.gt.0) then
                   write (iou,*)'   LAW= 33 (two body law)'
                 else
                   write (iou,*)'   LAW= 33 (two body law) used.',
     &               ' Warning: It is not a two body reaction'
                 endif
                 xss(lawxs)=33
                 xss(lxs)=-q
                 xss(lxs+1)=(awr0-awp)/awr0
                 lxs=lxs+2
                 if (law.eq.2) call nextsub6(nin,law,nbt,ibt,x,b)
               endif
             else
               call nextsub6(nin,law,nbt,ibt,x,b)
             endif
           enddo
           deallocate(x6,y6)
         endif
       enddo
c
c      divide the heating contribution by the total cross section
c      to get the heating numbers in Mev/collisions.
c      Add it to total heating for charged particles (ip>1), because
c      it was subtracted or not included during MF6 processing
c
       do ie=it,ne1
         if (xss(ltot+ie-1).ne.0.0d0) then
            xss(lphn+2+ie-it)=xss(lphn+2+ie-it)/xss(ltot+ie-1)
         endif
         if (ip.gt.1) then
           xss(lthn+ie-1)=xss(lthn+ie-1)+xss(lphn+2+ie-it)
         endif
       enddo
       write(iou,*)' LDLWP and DLWP block done'
       write(iou,'(80a1)')line1
c
c      ending cycle for production production data block of particle ip
c
      enddo
      kend=lxs-1
      nxs(1)=kend
      write(*,*)' end=',kend
      write(iou,*)' end=',kend
      write(iou,*)
c
c      convert main incident energy grid to MeV
c      checking negative heating values
c
      nneg=0
      do ie=1,ne1
        xss(lesz+ie-1)=xss(lesz+ie-1)*ev2mev
        if (xss(lthn+ie-1).lt.0.0d0) then
          write(iou,*)' Warning: negative heating number found at E=',
     &    xss(lesz+ie-1),' heat=',xss(lthn+ie-1)
          xss(lthn+ie-1)=0.0d0
          nneg=nneg+1
        endif
      enddo
      if (nneg.gt.0) then
        write(iou,*)' Warning: ',nneg,' negative heating numbers found',
     &    ' and set to zero'
      endif
c
c      write the ace file
c
      close(nin)
      open (nou,file=fout)
      call change(iou,nou,mcnpx)
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
      if (mcnpx.eq.1) then
        write(nou,'(a13,f12.6,a,a,a7,i9,a4,1pe11.4)')hz(1:13),
     &  awr0,' ',trim(fout),' 0 1 1 ',kend,' 0 0',tz
      else
        write(nou,'(a10,f12.6,a,a,a7,i9,a4,1pe11.4)')hz(1:10),
     &  awr0,' ',trim(fout),' 0 1 1 ',kend,' 0 0',tz
      endif
      close(nou)
      if (nfi.gt.0) close(nfi, status='DELETE')
      if (ngi.gt.0) close(ngi, status='DELETE')
      write(*,*)' DOPHN ENDED'
      write(iou,*)
      write(iou,*)' DOPHN ENDED'
      close(iou)
      stop
      end
C======================================================================
      subroutine getthr(np,x,y)
      implicit real*8 (a-h, o-z)
      dimension x(*),y(*)
      i0=np
      do i=1,np-1
        if (y(i).gt.0.0d0) then
          i0=i
          exit
        endif
      enddo
      if (i0.gt.1) i0=i0-1
      if (i0.gt.1) then
        do i=i0,np
          ii=i-i0+1
          x(ii)=x(i)
          y(ii)=y(i)
        enddo
        np=np-i0+1
      endif
      return
      end
C======================================================================
      subroutine photofis(nin,iou,nfi,awi,mat0,nmf4,mf4,nmf5,mf5,
     &  nmf6,mf6,nnu,xnu,ynu,idneu)
c
c       Get nubar if fissionable
c       Prepare unified MF6/MT18 file if delayed data (MF1/MT455 and
c       MF5/MT455) are available.
c
      implicit real*8(a-h,o-z)
      parameter (emnc2=939.56542052539d+6, nbmax=2000000, nmax=50000)
      dimension mf4(*),mf5(*),mf6(*)
      dimension xnu(*),ynu(*)
      dimension ibt(20),nbt(20),ibt1(20),nbt1(20),ibt2(20),nbt2(20)
      dimension xnup(nmax),ynup(nmax),xnud(nmax),ynud(nmax)
      dimension xnut(nmax),ynut(nmax)
      dimension xe4(nmax),xe5(nmax),xe(nmax),ye(nmax)
      dimension x(nmax),y(nmax),u(nmax),v(nmax)
      allocatable b(:)
      allocate(b(nbmax))
      data nf4/46/,nf5/56/
      nnu=0
      nut=0
      nup=0
      nud=0
      call findmt(nin,mat0,1,452,icod)
      if (icod.eq.0) then
        call readcont(nin,c1,c2,l1,lnu,n1,n2,mat,mf,mt,ns1)
        if (lnu.ne.2) then
          write(iou,*)' ERROR: LNU not equal 2 for MT452'
          close(nin)
          close(iou)
          stop
        endif
        call readtab1(nin,c1,c2,l1,l2,nr,nut,nbt,ibt,xnut,ynut)
      endif
      call findmt(nin,mat0,1,456,icod)
      if (icod.eq.0) then
        call readcont(nin,c1,c2,l1,lnu,n1,n2,mat,mf,mt,ns1)
        if (lnu.ne.2) then
          write(iou,*)' ERROR: LNU not equal 2 for MT456'
          close(nin)
          close(iou)
          stop
        endif
        call readtab1(nin,c1,c2,l1,l2,nr1,nup,nbt1,ibt1,xnup,ynup)
      endif
      call findmt(nin,mat0,1,455,icod)
      if (icod.eq.0) then
        call readcont(nin,c1,c2,ldg,lnu,n1,n2,mat,mf,mt,ns1)
        if (ldg.ne.0.or.lnu.ne.2) then
          write(iou,*)' ERROR: LNU=',lnu,' LDG=',ldg,' for MT455'
          write(iou,*)'        LNU=2 and LDG=0 are expected'
          stop
        endif
        call readlist(nin,c1,c2,l1,l2,ndf,n2,b)
        call readtab1(nin,c1,c2,l1,l2,nr2,nud,nbt2,ibt2,xnud,ynud)
      endif
      m4=0
      m5=0
      m6=0
      m455=0
      if (nmf4.gt.0) m4=iposm(nmf4,mf4,18)
      if (nmf6.gt.0) m6=iposm(nmf6,mf6,18)
      if (nmf5.gt.0) then
         m5=iposm(nmf5,mf5,18)
         m455=iposm(nmf5,mf5,455)
      endif
      if ((nup.gt.0.or.nut.gt.0).and.(m4*m5.gt.0.or.m6.gt.0).and.
     &    (nud.gt.0.and.m455.gt.0).and.idneu.ge.0) then
c
c       Prepare common fision mf6 file on tape nfi with extensions from
c        (mf4/mt18 + mf5/mt18 + mf5/mt455)  or  (mf6/mt18 + mf5/mt455)
c
        zap=1.0d0
        awp=1.0d0
        ewi=awi*emnc2
        law=1
        lang=12
        lep=2
        nd=0
        lip=0
        ns=0
        if (nut.gt.0.and.nup.gt.0) then
c
c         total, prompt and delayed nubar given
c
          call union(nud,xnud,nut,xnut,nx,x,nmax)
          call union(nup,xnup,nx,x,nnu,xnu,nmax)
          do i=1,nnu
            ynu(i)=fvalue(nr,nbt,ibt,nut,xnut,ynut,xnu(i))
          enddo
          call checkdis(iou,nnu,xnu,ynu,icod)
          if (icod.gt.0) then
            write(iou,*)' Warning: nubar ',icod,' points removed'
          endif
        elseif (nut.gt.0) then
c
c         total and delayed nubar given. Prompt nubar is calculated
c
          call union(nud,xnud,nut,xnut,nnu,xnu,nmax)
          do i=1,nnu
            ynu(i)=fvalue(nr,nbt,ibt,nut,xnut,ynut,xnu(i))
          enddo
          call checkdis(iou,nnu,xnu,ynu,icod)
          if (icod.gt.0) then
            write(iou,*)' Warning: nubar ',icod,' points removed'
          endif
          nup=nnu
          do i=1,nup
            xx=xnu(i)
            yy=ynu(i)-fvalue(nr2,nbt2,ibt2,nud,xnud,ynud,xx)
            xnup(i)=xx
            if (yy.gt.0) then
              ynup(i)=yy
            else
              ynup(i)=0.0d0
            endif
          enddo
          nr1=1
          nbt1(1)=nup
          ibt1(1)=2
        elseif (nup.gt.0) then
c
c         prompt and delayed nubar given. Total nubar is calculated
c
          call union(nud,xnud,nup,xnup,nnu,xnu,nmax)
          do i=1,nnu
            xx=xnu(i)
            ynu(i)=fvalue(nr1,nbt1,ibt1,nup,xnup,ynup,xx)+
     &             fvalue(nr2,nbt2,ibt2,nud,xnud,ynud,xx)
          enddo
          call checkdis(iou,nnu,xnu,ynu,icod)
          if (icod.gt.0) then
            write(iou,*)' Warning: nubar ',icod,' points removed'
          endif
        endif
c
c       angle-energy distributions
c
        call findmt(nin,mat0,5,455,icod)
        call readcont(nin,c1,c2,l1,l2,nkd,n2,mat1,mf1,mt1,ns1)
        if (m4.gt.0) then
          nfi=40
          open(nfi,file='MF6MT18T.dat')
          write(nfi,'(a10,56x,i4,i2,i3)')' MF6MT18T ',7000,0,0
c
c         mf4/mf5 representation of prompt neutrons
c
          call findmt(nin,mat0,4,18,icod)
          call readcont(nin,c1,c2,l1,ltt,n1,n2,mat1,mf1,mt1,ns1)
          call readcont(nin,c1,c2,li,lct,n1,n2,mat1,mf1,mt1,ns1)
          if (ltt.eq.0.and.li.eq.1) then
            iso=1
          elseif (ltt.eq.2) then
            iso=0
            call readtab2(nin,c1,c2,l1,l2,nr,ne4,nbt,ibt)
            call checklaw(nr,ibt,icod)
            if (icod.ne.0) then
              write(iou,*)' Incident energy grid is not linearly',
     &        ' interpolable. Run LEGEND'
              stop
            endif
            open(nf4,file='MF4MT18P.dat',form='UNFORMATTED')
            do j=1,ne4
              call readtab1(nin,c1,xe4(j),l1,l2,nr,np,nbt,ibt,x,y)
              call checklaw(nr,ibt,icod)
              if (icod.ne.0) then
                write(iou,*)' Angular grid is not linearly',
     &            ' interpolable. Run LEGEND/ACEMAKER'
                stop
              endif
              x0=x(1)
              dx=x(np)-x0
              write(nf4)xe4(j),np,x0,dx
              write(nf4)(x(i),y(i),i=1,np)
            enddo
          else
            write(iou,*)' ERROR: ltt=',ltt,' only tabular data allowed'
            write(iou,*)' Run LEGEND/ACEMAKER'
            stop
          endif
          call findmt(nin,mat0,5,18,icod)
          call readcont(nin,c1,awr,l1,l2,nkp,n2,mat1,mf1,mt1,ns1)
          ewr=awr*emnc2
          nk=nkp+nkd
          call wrtcont(nfi,mat0,6,18,ns,c1,awr,0,lct,nk,0)
          do ik=1,nkp
            call readtab1(nin,c1,c2,l1,lf,nr,np,nbt,ibt,x,y)
            if (lf.ne.1) then
              write(iou,*)' ERROR: LF=',lf,' not allowed'
              write(iou,*)' Run SPECTRA/ACEMAKER'
              stop
            endif
            do iy=1,nnu
              xx=xnu(iy)
              yx=ynu(iy)
              if (yx.gt.0.0d0) then
                yy=fvalue(nr1,nbt1,ibt1,nup,xnup,ynup,xx)
                yy=yy*fvalue(nr,nbt,ibt,np,x,y,xx)
                ye(iy)=yy/yx
              else
                ye(iy)=0.0d0
              endif
            enddo
            nr=1
            nbt(1)=nnu
            ibt(1)=2
            call wrtab1(nfi,mat0,6,18,ns,zap,awp,lip,law,
     &        nr,nbt,ibt,nnu,xnu,ye)
            lip=lip+1
            call readtab2(nin,c1,c2,l1,l2,nr,ne5,nbt,ibt)
            call checklaw(nr,ibt,icod)
            if (icod.ne.0) then
              write(iou,*)' Incident energy grid is not linearly',
     &          ' interpolable. Not allowed for LF=1 in MF5'
              stop
            endif
            open(nf5,file='MF5MT18P.dat',form='UNFORMATTED')
            do ie=1,ne5
              call readtab1(nin,xe5(ie),c2,l1,l2,nr,np,nbt,ibt,x,y)
              call checklaw(nr,ibt,icod)
              if (icod.ne.0) then
                write(iou,*)' Secondary energy grid is not linearly',
     &            ' interpolable. Use SPECTRA/ACEMAKER'
                stop
              endif
              x0=x(1)
              dx=x(np)-x0
              write(nf5)xe5(ie),np,x0,dx
              write(nf5)(x(i),y(i),i=1,np)
            enddo
            if (iso.eq.0) then
              call union(ne4,xe4,ne5,xe5,ne,xe,nmax)
            else
              ne=ne5
              do ie=1,ne
                xe(ie)=xe5(ie)
              enddo
            endif
            nr=1
            nbt(1)=ne
            ibt(1)=2
            call wrtab2(nfi,mat0,6,18,ns,
     &        0.0d0,0.0d0,lang,lep,nr,ne,nbt,ibt)
            e=-1.0d0
            call fep(e,nep,x,y,nf5)
            if (iso.eq.0) then
              e=-1.0d0
              call fmu(e,na0,u,v,nf4)
            endif
            do ie=1,ne
              e=xe(ie)
              call fep(e,nep,x,y,nf5)
              if (iso.eq.0) then
                call fmu(e,na0,u,v,nf4)
              else
                na0=0
              endif
              if (lct.ne.1.and.idneu.gt.0) then
c
c               Energy distribution of prompt neutrons given in the LAB
c               system and uncorrelated angular distribution given in
c               the CM system: Full conversion from LAB to CM is not
c               applied. A first order correction is used considering
c               that for incident photons with energies below 200 MeV
c               the CM and the LAB systems are almost the same for the
c               important heavy photo fissionable materials with
c               delayed data available.
c               For the correction, it is assumed that the angular
c               distribution in the LAB system is not too anisotropic,
c               so the average cosine is close to zero.
c               Warning: It can be a bad aproximation for very
c                        anisotropic photofission.
c
                if (ewi.eq.0.0d0) then
                  beta=e/(e+ewr)
                else
                  beta=sqrt(e*(e+2.0d0*ewi))/(e+ewi+ewr)
                endif
                gam=1.0d0/sqrt(1.0d0-beta*beta)
                e0=(gam-1.0d0)*emnc2
              else
c
c               angular and energy distribution given in the LAB system
c
                gam=1.0d0
                e0=0.0d0
              endif
              na=2*na0
              nw=(na+2)*nep
              allocate(b(nw))
              j=0
              do iep=1,nep
                j=j+1
                b(j)=gam*x(iep)+e0
                f0=y(iep)/gam
                j=j+1
                if (f0.gt.0.0d0) then
                  b(j)=f0
                else
                  b(j)=0.0d0
                endif
                if (na0.gt.0) then
                  do kk=1,na0
                    j=j+1
                    b(j)=u(kk)
                    j=j+1
                    b(j)=0.5d0*v(kk)
                  enddo
                endif
                call wrtlist(nfi,mat0,6,18,ns,0.0d0,e,nd,na,nw,nep,b)
              enddo
              deallocate(b)
            enddo
            close(nf5,status='DELETE')
          enddo
          if (iso.eq.0) close(nf4,status='DELETE')
c
c         Now prompt neutron distribution for MT=18 is on MF6/LAW=1
c
          i0=m5+1
          do i=i0,nmf5
            mf5(i-1)=mf5(i)
          enddo
          nmf5=nmf5-1
        else
c
c         mf6 representation of prompt neutrons
c
          nfi=60
          open(nfi,file='MF6MT18T.dat')
          write(nfi,'(a10,56x,i4,i2,i3)')' MF6MT18T ',7000,0,0
          call findmt(nin,mat0,6,18,icod)
          call readcont(nin,za,awr,jp,lct,nk6,n2,mat1,mf1,mt1,ns1)
          if (jp.ne.0) then
            write(iou,*)' JP=',jp,' not coded for photo-fission (MT18)'
            stop
          endif
          do ik=1,nk6
            call readtab1(nin,zap6,awp6,l1,law6,nr,np,nbt,ibt,x,y)
            izap6=nint(zap6)
            if (izap6.eq.1) then
              if (law6.ne.1) then
                write(iou,*)' LAW=',law6,' not coded/allowed for MT=18',
     &            ' Run SIXLIN'
                stop
              endif
              nk=nkd+1
              call wrtcont(nfi,mat0,6,18,ns,za,awr,jp,lct,nk,0)
              do iy=1,nnu
                xx=xnu(iy)
                yx=ynu(iy)
                if (yx.gt.0.0d0) then
                  yy=fvalue(nr1,nbt1,ibt1,nup,xnup,ynup,xx)
                  ye(iy)=yy/yx
                else
                  ye(iy)=0.0d0
                endif
              enddo
              nr=1
              nbt(1)=nnu
              ibt(1)=2
              call wrtab1(nfi,mat0,6,18,ns,zap,awp,lip,law6,
     &          nr,nbt,ibt,nnu,xnu,ye)
              lip=lip+1
              call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
              call wrtab2(nfi,mat0,6,18,ns,c1,c2,l1,l2,nr,ne,nbt,ibt)
              do ie=1,ne
                call readlist(nin,c1,c2,l1,l2,nw,n2,b)
                call wrtlist(nfi,mat0,6,18,ns,c1,c2,l1,l2,nw,n2,b)
              enddo
              exit
            else
              call nextsub6(nin,law6,nbt,ibt,x,b)
            endif
          enddo
        endif
c
c       add isotropic delayed neutron distributions
c
        call findmt(nin,mat0,5,455,icod)
        call readcont(nin,za,awr,l1,l2,nkd,n2,mat1,mf1,mt1,ns1)
        ewr=awr*emnc2
        do ik=1,nkd
          call readtab1(nin,c1,c2,l1,lf,nr,np,nbt,ibt,x,y)
          if (lf.ne.1) then
            write(iou,*)' ERROR: LF=',lf,' not allowed'
            write(iou,*)' Run SPECTRA'
            stop
          endif
          do iy=1,nnu
            xx=xnu(iy)
            yx=ynu(iy)
            if (yx.gt.0.0d0) then
              yy=fvalue(nr2,nbt2,ibt2,nud,xnud,ynud,xx)
              yy=yy*fvalue(nr,nbt,ibt,np,x,y,xx)
              ye(iy)=yy/yx
            else
              ye(iy)=0.0d0
            endif
          enddo
          nr=1
          nbt(1)=nnu
          ibt(1)=2
          call wrtab1(nfi,mat0,6,18,ns,zap,awp,lip,law,
     &      nr,nbt,ibt,nnu,xnu,ye)
          lip=lip+1
          call readtab2(nin,c1,c2,l1,l2,nr,ne5,nbt,ibt)
          call checklaw(nr,ibt,icod)
          if (icod.ne.0) then
            write(iou,*)' Incident energy grid is not linearly',
     &        ' interpolable. Not allowed for LF=1 in MF5'
            stop
          endif
          c1=0.0d0
          c2=0.0d0
          call wrtab2(nfi,mat0,6,18,ns,c1,c2,lang,lep,nr,ne5,nbt,ibt)
          do ie=1,ne5
            call readtab1(nin,c1,e,l1,l2,nr,np,nbt,ibt,x,y)
            call checklaw(nr,ibt,icod)
            if (icod.ne.0) then
              write(iou,*)' Secondary energy grid is not linearly',
     &          ' interpolable. Use SPECTRA'
              stop
            endif
            if (lct.ne.1.and.idneu.gt.0) then
c
c             Full conversion from LAB to CM is not applied. A first
c             order correction is used considering that the fraction
c             of delayed neutrons is small (<1%) and the angular
c             distribution of the delayed neutrons is considered
c             isotropic in the LAB system (average cosine=0).
c             Furthermore, for incident photons with energies below
c             200 MeV, the CM and the LAB systems are almost the same
c             for the important heavy photo fissionable isotopes with
c             delayed data available.
c             The applied correction conserves the average energy
c             of the delayed neutrons in the LAB system.
c
              if (ewi.eq.0.0d0) then
                beta=e/(e+ewr)
              else
                beta=sqrt(e*(e+2.0d0*ewi))/(e+ewi+ewr)
              endif
              gam=1.0d0/sqrt(1.0d0-beta*beta)
              e0=(gam-1.0d0)*emnc2
              do iep=1,np
                x(iep)=gam*x(iep)+e0
                y(iep)=y(iep)/gam
              enddo
            endif
            na=0
            nw=2*np
            j=0
            do iep=1,np
              j=j+1
              b(j)=x(iep)
              j=j+1
              b(j)=y(iep)
            enddo
            call wrtlist(nfi,mat0,6,18,ns,0.0d0,e,nd,na,nw,np,b)
          enddo
        enddo
        call wrtsend(nfi,mat0,6,ns)
        call wrtfend(nfi,mat0,ns)
        call wrtmend(nfi,ns)
        call wrtend(nfi,ns)
        rewind(nfi)
      elseif (nup.gt.0) then
c
c       prompt nubar given.
c
        nnu=nup
        do i=1,nup
          xnu(i)=xnup(i)
          ynu(i)=ynup(i)
        enddo
        nfi=-1
      elseif (nut.gt.0) then
c
c       total nubar given
c
        nnu=nut
        do i=1,nut
          xnu(i)=xnut(i)
          ynu(i)=ynut(i)
        enddo
        nfi=-1
      else
        nnu=-1
        nfi=-1
      endif
      deallocate(b)
      return
      end
C======================================================================
      subroutine fep(z,n,x,y,inp)
c
c       E' interpolation fep=f(E,E')
c       using unit-base linear interpolation
c
      implicit real*8 (a-h, o-z)
      parameter (nmax=50000)
      dimension x(*),y(*)
      dimension x1(nmax),y1(nmax),x2(nmax),y2(nmax)
      save z1,n1,x10,dx1,x1,y1,z2,n2,x20,dx2,x2,y2,dlow,dupp,dz
      if (z.lt.0.0d0) then
c
c       initio
c
        rewind(inp)
        read(inp)z1,n1,x10,dx1
        read(inp)(x1(i),y1(i),i=1,n1)
        read(inp)z2,n2,x20,dx2
        read(inp)(x2(i),y2(i),i=1,n2)
        dlow=x20-x10
        dupp=dlow+dx2-dx1
        dz=z2-z1
        na0=0
        do i=1,n1
          x(i)=0.0d0
          y(i)=0.0d0
        enddo
        return
      endif
      if (z.gt.z2) then
c
c        load new panel
c
         do while (z.gt.z2)
           z1=z2
           n1=n2
           x10=x20
           dx1=dx2
           do i=1,n2
             x1(i)=x2(i)
             y1(i)=y2(i)
           enddo
           read(inp,end=10)z2,n2,x20,dx2
           read(inp,end=10)(x2(i),y2(i),i=1,n2)
         enddo
         dlow=x20-x10
         dupp=dlow+dx2-dx1
         dz=z2-z1
      elseif (z.lt.z1) then
c
c       out of range (below lower limit)
c
        n=n1
        do i=1,n
          x(i)=x1(i)
          y(i)=0.0d0
        enddo
        return
      endif
c
c       Interpolation between z1 and z2
c
      if (z.eq.z1) then
        n=n1
        do i=1,n
          x(i)=x1(i)
          y(i)=y1(i)
        enddo
      elseif (z.eq.z2) then
        n=n2
        do i=1,n
          x(i)=x2(i)
          y(i)=y2(i)
        enddo
      else
        dz1=(z-z1)/dz
        dz2=(z2-z)/dz
        xmin=x10+dlow*dz1
        xmax=x10+dx1+dupp*dz1
        dx=xmax-xmin
        dx1j=dx1/dx
        dx2j=dx2/dx
        if (dz1.lt.dz2) then
          n=n1
          do i=1,n
            xx1=x1(i)
            xxx=xmin+(xx1-x10)/dx1j
            xx2=x20+(xxx-xmin)*dx2j
            yy1=fvalin(n1,x1,y1,xx1)*dx1j
            yy2=fvalin(n2,x2,y2,xx2)*dx2j
            x(i)=xxx
            y(i)=(yy2-yy1)*dz1+yy1
          enddo
        else
          n=n2
          do i=1,n
            xx2=x2(i)
            xxx=xmin+(xx2-x20)/dx2j
            xx1=x10+(xxx-xmin)*dx1j
            yy1=fvalin(n1,x1,y1,xx1)*dx1j
            yy2=fvalin(n2,x2,y2,xx2)*dx2j
            x(i)=xxx
            y(i)=(yy2-yy1)*dz1+yy1
          enddo
        endif
      endif
      return
c
c       out of range (above upper limit)
c
   10 n=n2
      do i=1,n
        x(i)=x2(i)
        y(i)=0.0d0
      enddo
      return
      end
C======================================================================
      subroutine fmu(z,n,x,y,inp)
c
c       Cosine interpolation fmu=f(E,u)
c
      implicit real*8 (a-h, o-z)
      parameter (nmax=50000)
      dimension x(*),y(*)
      dimension x1(nmax),y1(nmax),x2(nmax),y2(nmax)
      save z1,n1,x10,dx1,x1,y1,z2,n2,x20,dx2,x2,y2,dlow,dupp,dz
      if (z.lt.0.0d0) then
c
c       initio
c
        rewind(inp)
        read(inp)z1,n1,x10,dx1
        read(inp)(x1(i),y1(i),i=1,n1)
        read(inp)z2,n2,x20,dx2
        read(inp)(x2(i),y2(i),i=1,n2)
        dlow=x20-x10
        dupp=dlow+dx2-dx1
        dz=z2-z1
        na0=0
        do i=1,n1
          x(i)=0.0d0
          y(i)=0.0d0
        enddo
        return
      endif
      if (z.gt.z2) then
c
c        load new panel
c
         do while (z.gt.z2)
           z1=z2
           n1=n2
           x10=x20
           dx1=dx2
           do i=1,n2
             x1(i)=x2(i)
             y1(i)=y2(i)
           enddo
           read(inp,end=10)z2,n2,x20,dx2
           read(inp,end=10)(x2(i),y2(i),i=1,n2)
         enddo
         dlow=x20-x10
         dupp=dlow+dx2-dx1
         dz=z2-z1
      elseif (z.lt.z1) then
c
c       out of range (below lower limit)
c
        n=n1
        do i=1,n
          x(i)=x1(i)
          y(i)=0.0d0
        enddo
        return
      endif
c
c       Interpolation between z1 and z2
c
      if (z.eq.z1) then
        n=n1
        do i=1,n
          x(i)=x1(i)
          y(i)=y1(i)
        enddo
      elseif (z.eq.z2) then
        n=n2
        do i=1,n
          x(i)=x2(i)
          y(i)=y2(i)
        enddo
      else
        dz1=(z-z1)/dz
        dz2=(z2-z)/dz
        xmin=x10+dlow*dz1
        xmax=x10+dx1+dupp*dz1
        dx=xmax-xmin
        dx1j=dx1/dx
        dx2j=dx2/dx
        if (dz1.lt.dz2) then
          n=n1
          do i=1,n
            xx1=x1(i)
            xxx=xmin+(xx1-x10)/dx1j
            xx2=x20+(xxx-xmin)*dx2j
            yy1=fvalin(n1,x1,y1,xx1)*dx1j
            yy2=fvalin(n2,x2,y2,xx2)*dx2j
            x(i)=xxx
            y(i)=(yy2-yy1)*dz1+yy1
          enddo
        else
          n=n2
          do i=1,n
            xx2=x2(i)
            xxx=xmin+(xx2-x20)/dx2j
            xx1=x10+(xxx-xmin)*dx1j
            yy1=fvalin(n1,x1,y1,xx1)*dx1j
            yy2=fvalin(n2,x2,y2,xx2)*dx2j
            x(i)=xxx
            y(i)=(yy2-yy1)*dz1+yy1
          enddo
        endif
      endif
      return
c
c       out of range (above upper limit)
c
   10 n=n2
      do i=1,n
        x(i)=x2(i)
        y(i)=0.0d0
      enddo
      return
      end
C======================================================================
      subroutine disgam(nin,iou,ngi,mat0,nmf3,mf3,nmf12,mf12,nmf14,mf14)
c
c      Prepare MF6/LAW=1 data for photon production of two-body
c      reactions given in MF12/MF14
c
      implicit real*8 (a-h,o-z)
      dimension mf3(*),mf12(*),mf14(*)
      dimension e(5000),yt(5000)
      dimension nbt(20),ibt(20)
      allocatable x(:),y(:),ep(:),yp(:,:),b(:)
      ngi=14
      open (ngi,file='MF6GAM.DAT')
      write(ngi,*)' FILE MF6 from MF12/MF14 outgoing gammas'
      ns=0
      do i=1,nmf12
        mt=mf12(i)
        ngam=numpart(0,mt)
        m3=iposm(nmf3,mf3,mt)
        if (ngam.gt.0.and.m3.gt.0) then
          nk14=-1
          m14=0
          if (nmf14.gt.0) m14=iposm(nmf14,mf14,mt)
          if (m14.gt.0) then
            call findmt(nin,mat0,14,mt,icod)
            call readcont(nin,c1,c2,li,l2,nk14,n2,mat1,mf1,mt1,ns1)
            if (li.ne.1) then
              write(iou,*)' Warning: LI=',li,' on MF14 for MT=',mt
              write(iou,*)' Option not coded for photonuclear data'
              write(iou,*)' Isotropic angular distribution will be',
     &          ' assumed for all outgoing photons'
            endif
          else
            write(iou,*)' Warning: MF14 not found for MT=',mt
            write(iou,*)' Isotropic angular distribution will be',
     &        ' assumed for all outgoing photons specified on MF12'
          endif
          call findmt(nin,mat0,12,mt,icod)
          call readcont(nin,za,awr,lo,l2,nk,n2,mat1,mf1,mt1,ns1)
          if (lo.eq.1) then
            if (nk14.gt.0.and.nk.ne.nk14) then
              write(iou,*)' Warning: Number of outgoing photons(nk) on',
     &         ' MF12 and MF14 are not the same',nk,nk14
              write(iou,*)' Assumed nk=',nk
            endif
            call readtab1(nin,eg,es,lp,lf,nr,ne,nbt,ibt,e,yt)
            nn=2*ne
            allocate (x(nn),y(nn),ep(nk),yp(nk,ne))
            if (nk.gt.1) then
              do ie=1,ne
                yt(ie)=0.0d0
              enddo
              jk=0
              do ik=1,nk
                call readtab1(nin,eg,es,lp,lf,nr,np,nbt,ibt,x,y)
                if (lf.eq.1.or.eg.eq.0.0d0) then
                  write(iou,*)' Warning: Continuous energy',
     &              ' distribution not recommended for discrete',
     &              ' photoreaction',' MT=',mt,' ik=',ik,' data ignored'
                else
                  jk=jk+1
                  if (lp.eq.2) eg=-eg
                  ep(jk)=eg
                  do ie=1,ne
                    ee=e(ie)
                    yp(jk,ie)=fvalue(nr,nbt,ibt,np,x,y,ee)
                    yt(ie)=yt(ie)+yp(jk,ie)
                  enddo
                endif
              enddo
              nk=jk
            else
              if (lf.eq.1.or.eg.eq.0.0d0) then
                write(iou,*)' Warning: Continuous energy distribution',
     &            ' is not recommended for discrete photoreaction',
     &            ' MT=',mt,'ik=',1,' eg=',eg,' Data ignored'
                nk=0
              else
                if (lp.eq.2) eg=-eg
                ep(1)=eg
                do ie=1,ne
                  yp(1,ie)=yt(ie)
                enddo
              endif
            endif
            if (nk.gt.0) then
              jp=0
              lct=1
              call wrtcont(ngi,mat0,6,mt,ns,za,awr,jp,lct,1,0)
              zap=0.0d0
              awp=0.0d0
              lip=0
              law=1
              nr=1
              nbt(1)=ne
              ibt(1)=2
              call wrtab1(ngi,mat0,6,mt,ns,zap,awp,lip,law,
     &          nr,nbt,ibt,ne,e,yt)
              c1=0.0d0
              c2=0.0d0
              lang=12
              lep=1
              call wrtab2(ngi,mat0,6,mt,ns,c1,c2,lang,lep,nr,ne,nbt,ibt)
              na=0
              nw=2*nk
              allocate(b(nw))
              do ie=1,ne
                ee=e(ie)
                j=0
                do ik=1,nk
                  epp=ep(ik)
                  if (epp.lt.0.0d0) epp=abs(epp)+ee
                  j=j+1
                  b(j)=epp
                  j=j+1
                  if (yt(ie).gt.0.0d0) then
                    b(j)=yp(ik,ie)/yt(ie)
                  else
                    b(j)=0.0d0
                  endif
                enddo
                call wrtlist(ngi,mat0,6,mt,ns,c1,ee,nk,na,nw,nk,b)
              enddo
              call wrtsend(ngi,mat0,6,ns)
              mf12(i)=-mt
              deallocate (b,x,y,ep,yp)
            else
c
c             Non discrete photon
c
              write(iou,*)' Warning: Non discrete outgoing photon',
     &          ' found for MT=',mt,' Section will be ignored'
            endif
          else
c
c           lo not equal 1
c
            write(iou,*)' Warning: LO=',lo,' is not equal 1',
     &        '(multiplicities) on MF12 for MT=',mt,' Run GAMLIN.',
     &        ' Section will be ignored'
          endif
        elseif (ngam.gt.0.and.m3.le.0) then
          write(iou,*)' Warning: MT=',mt,' missing on MF3, but given',
     &      ' on MF12. Reaction will be ignored'
        else
          write(iou,*)' Warning: MT=',mt,' not allowed on MF14 for',
     &      ' photonuclear data (NLIB=0). Section will be ignored'
        endif
      enddo
      nmf14=0
      do i=1,nmf12
        mt=mf12(i)
        if (mt.lt.0) then
          nmf14=nmf14+1
          mf14(nmf14)=-mt
          mf12(i)=-mt
        endif
      enddo
      if (nmf14.gt.0) then
        call wrtfend(ngi,mat0,ns)
        call wrtmend(ngi,ns)
        call wrtend(ngi,ns)
        rewind (ngi)
      else
        close(ngi,status='DELETE')
        ngi=-1
      endif
      return
      end
C======================================================================
      subroutine yld6(nin,iou,scr,nnu,xnu,ynu,ip,mt,nk,thresh,
     &                kk,ny,xt,yt,nr,nbt,ibt)
c
c      Get the total yield of the particle ip for reaction MT
c      given on MF6
c
      parameter (nymax=50000)
      implicit real*8 (a-h, o-z)
      dimension nbt(*),ibt(*)
      dimension scr(*),xnu(*),ynu(*),xt(*),yt(*)
      dimension nbt1(20),ibt1(20)
      dimension xu(nymax),xx(nymax),yy(nymax)
c
c     calculate total yield for particle izap=ip from MF6/MT=mt
c       xt: incident energy grid
c       yt: yield
c       kk: number of ip particles found
c       linear interpolation law is required
c
      kk=0
      if (mt.eq.18.and.ip.eq.1.and.nnu.gt.0) then
        ny=nnu
        do i=1,nnu
          xt(i)=xnu(i)
          yt(i)=ynu(i)
        enddo
        nr=1
        nbt(1)=nnu
        ibt(1)=2
        kk=1
      else
        do ik=1,nk
          call readtab1(nin,zap,awp,lip,law,nr,npy,nbt,ibt,xx,yy)
          izap=nint(zap)
          if (izap.eq.ip) then
            if (xx(1).lt.0.0d0) then
              if (xx(2).gt.thresh) then
                xx(1)=thresh
                write(iou,*)' Negative yield threshold found for',
     &            ' MF6/MT=',mt,'. Set to ',thresh
              else
                write(iou,*)' Negative yield threshold found for',
     &            ' MF6/MT=',mt,'. Fatal error. Stop at subroutine yld6'
                stop
              endif
            endif
            if (kk.eq.0) then
              call checklaw(nr,ibt,icod1)
              ny=npy
              do i=1,npy
                xt(i)=xx(i)
                yt(i)=yy(i)
              enddo
              nr1=nr
              do i=1,nr
                nbt1(i)=nbt(i)
                ibt1(i)=ibt(i)
              enddo
              kk=1
            else
              call checklaw(nr,ibt,icod)
              if (icod.ne.0.or.icod1.ne.0) then
                write(iou,*)' Warning: MF6/MT=',mt,' yields are not',
     &          ' linearly interpolable'
                write(iou,*)' Use ACEMAKER/SIXLIN'
                write(iou,*)' Linear interpolation will be assumed'
              endif
              call union(npy,xx,ny,xt,npu,xu,nymax)
              do i=1,npu
                xxx=xu(i)
                yyy=fvalin(npy,xx,yy,xxx)
                yyt=fvalin(ny,xt,yt,xxx)
                scr(i)=yyy+yyt
              enddo
              ny=npu
              do i=1,npu
                xt(i)=xu(i)
                yt(i)=scr(i)
              enddo
              kk=kk+1
            endif
          endif
          call nextsub6(nin,law,nbt,ibt,xx,scr)
        enddo
        if (kk.eq.1) then
          nr=nr1
          do i=1,nr1
            nbt(i)=nbt1(i)
            ibt(i)=ibt1(i)
          enddo
        else
          nr=1
          nbt(1)=ny
          ibt(1)=2
        endif
      endif
      return
      end
C======================================================================
      subroutine recoil(lang,nmu,b)
      implicit real*8 (a-h, o-z)
      dimension b(*)
      allocatable x(:), y(:)
      allocate (x(nmu),y(nmu))
      if (lang.eq.0) then
        do jl=1,nmu,2
          b(jl)=-b(jl)
        enddo
      else
        il=0
        do jl=1,nmu
          il=il+1
          x(jl)=b(il)
          il=il+1
          y(jl)=b(il)
        enddo
        il=0
        do jl=nmu,1,-1
          il=il+1
          b(il)=-x(jl)
          il=il+1
          b(il)=y(jl)
        enddo
      endif
      deallocate (x,y)
      return
      end
C======================================================================
      real*8 function aveep(ee,izai,awi,matza,awr,izap,awp,
     &  lct,lang,lep,nd,na,nep,b)
c
c      Calculate average energy taking into account the reference
c      system (LAB or CM). Relativistic transformation is applied.
c
      implicit real*8 (a-h,o-z)
      parameter (emnc2=939.56542052539d+6, ev2mev=1.0d-6)
      dimension b(*)
      allocatable ep(:),f0(:),uave(:),fu(:),xnu(:)
      allocate(ep(nep),f0(nep),uave(nep))
      eemev=ee*ev2mev
      ewi=awi*emnc2
      ewr=awr*emnc2
      ewp=awp*emnc2
      if (izai.eq.0) then
        beta=ee/(ee+ewr)
      else
        beta=sqrt(ee*(ee+2.0d0*ewi))/(ee+ewi+ewr)
      endif
      if (lct.eq.2.or.(lct.eq.3.and.awp.lt.4.0d0).or.lct.eq.4) then
c
c       CM system
c
c       Warning:
c        lct=4 will be treated as lct=2. It is not correct for break-up
c        products. This feature has not be used for photonuclear data.
c
        lcm=1
      else
c
c       LAB system
c
        lcm=0
      endif
      ncyc=na+2
      do iep=1,nep
        j=(iep-1)*ncyc+1
        ep(iep)=abs(b(j))
        f0(iep)=b(j+1)
        uave(iep)=0.0d0
        if (lcm.gt.0.and.na.gt.0) then
          if (lang.eq.0.or.lang.eq.1) then
c
c           Legendre coefficient representation
c
            if(f0(iep).ne.0.0d0) then
              uave(iep)=b(j+2)/f0(iep)
            endif
          elseif (lang.eq.2) then
c
c           Kalbach-Mann representation
c
            r=b(j+2)
            if (na.eq.1) then
              epmev=ep(iep)*ev2mev
              aa=bachaa(izai,izap,matza,eemev,epmev)
            else
              aa=b(j+3)
            endif
            if (abs(aa).gt.1.0d-6) then
              xplus=exp(aa)
              xmin=exp(-aa)
              uave(iep)=r*((xplus+xmin)/(xplus-xmin)-1.0d0/aa)
            else
              uave(iep)=r*aa/3.0d0
            endif
          else
c
c           Tabulated data
c
            nu=na/2
            allocate(fu(nu),xnu(nu))
            j=j+2
            do iu=1,nu
              xnu(iu)=b(j)
              fu(iu)=b(j+1)
              j=j+2
            enddo
            uave(iep)=avecos(nu,xnu,fu,lang)
            deallocate(fu,xnu)
          endif
          if (uave(iep).lt.-1.0d0) then
            uave(iep)=-1.0d0
          elseif (uave(iep).gt.1.0d0) then
            uave(iep)=1.0d0
          endif
        endif
      enddo
      sumep=0.0d0
      sumf0=0.0d0
      if (nd.gt.0) then
        do iep=1,nd
          epp=ep(iep)
          eplab=elab(epp,uave(iep),ewp,lcm,beta)
          sumep=sumep+eplab*f0(iep)
          sumf0=sumf0+f0(iep)
        enddo
      endif
      if (nep.gt.nd) then
        iepc=nd+1
        nepc=nep-1
        do iep=iepc,nepc
          iep1=iep+1
          call ef0int(epfint,fint,ep(iep),ep(iep1),uave(iep),uave(iep1),
     &      f0(iep),f0(iep1),ewp,beta,lcm,lep)
          sumep=sumep+epfint
          sumf0=sumf0+fint
        enddo
      endif
      deallocate(ep,f0,uave)
      if (sumf0.ne.0.0d0) then
        aveep=sumep/sumf0
      else
        aveep=0.0d0
      endif
      return
      end
C======================================================================
      real*8 function avecos(nu,x,y,intlaw)
c
c      calculate the average cosine for tabulated data
c
      implicit real*8 (a-h,o-z)
      parameter (c03=0.33333333333333d0, eps=1.0d-8, nmax=25)
      dimension xy0(nmax),xy1(nmax), y0(nmax),y1(nmax)
      dimension x(*),y(*)
      ll=intlaw/10
      ilaw=intlaw-10*ll
      sumxy=0.0d0
      sumy=0.0d0
      u1=x(1)
      f1=y(1)
      if (ilaw.eq.1) then
        do i=2,nu
          u2=x(i)
          f2=y(i)
          du=u2-u1
          if (du.gt.0) then
            su=u2+u1
            fdu=f1*du
            sumxy=sumxy+0.5d0*fdu*su
            sumy=sumy+fdu
          endif
          u1=u2
          f1=f2
        enddo
      elseif (ilaw.eq.2) then
        do i=2,nu
          u2=x(i)
          f2=y(i)
          du=u2-u1
          if (du.ne.0.0d0) then
            su=u2+u1
            slope=(f2-f1)/du
            b=f1-slope*u1
            sumxy=sumxy+du*(c03*slope*(u2*u2+u2*u1+u1*u1)+0.5d0*b*su)
            sumy=sumy+du*(0.5d0*slope*su+b)
          endif
          u1=u2
          f1=f2
        enddo
      else
        do l=2,nu
          u2=x(l)
          f2=y(l)
          h=u2-u1
          if (h.gt.0.0d0) then
            h2=0.5d0*h
            xy0(1)=h2*(u1*f1+u2*f2)
            y0(1)=h2*(f1+f2)
            h=h2
            jmax=1
            do i=1,nmax-1
              sumxy1=0.5d0*xy0(1)
              sumy1=0.5d0*y0(1)
              do j=1,jmax
                uj=u1+(2*j-1)*h
                call terp1m(u1,f1,u2,f2,uj,fj,ilaw)
                sumxy1=sumxy1+h*uj*fj
                sumy1=sumy1+h*fj
              enddo
              xy1(1)=sumxy1
              y1(1)=sumy1
              jmax=2*jmax
              pow=4.0d0
              do k=1,i
                den=pow-1.0d0
                xy1(k+1)=xy1(k)+(xy1(k)-xy0(k))/den
                y1(k+1)=y1(k)+(y1(k)-y0(k))/den
                pow=4.0d0*pow
              enddo
              i1=i+1
              if ((abs(xy1(i1)-xy1(i)).le.eps*abs(xy1(i1)).and.
     &             abs(y1(i1)-y1(i)).le.eps*abs(y1(i1))).or.
     &             i1.eq.nmax) then
                 exit
              else
                h=0.5d0*h
                do k=1,i1
                  xy0(k)=xy1(k)
                  y0(k)=y1(k)
                enddo
              endif
            enddo
            sumxy=sumxy+xy1(i1)
            sumy=sumy+y1(i1)
          endif
          u1=u2
          f1=f2
        enddo
      endif
      if (sumy.ne.0.0d0) then
        avecos=sumxy/sumy
      else
        avecos=0.0d0
      endif
      return
      end
C======================================================================
      real*8 function elab(tp,up,ewp,lcm,beta)
c
c       Convert secondary energy from CM to LAB for averaging
c
      implicit real*8 (a-h,o-z)
      if (lcm.gt.0) then
        gam=1.0d0/(sqrt(1.0d0-beta*beta))
        elab=gam*(tp+beta*up*sqrt(tp*(tp+2.0d0*ewp)))+(gam-1.0d0)*ewp
      else
        elab=tp
      endif
      return
      end
C======================================================================
      subroutine ef0int(epfint,fint,e1,e2,u1,u2,f1,f2,ewp,beta,lcm,lep)
c
c       integration of e'*f(e,e') and f(e,e') between e1 and e2
c
      implicit real*8 (a-h,o-z)
      parameter (eps=1.0d-8, nmax=25)
      dimension y0(nmax),y1(nmax)
      h=e2-e1
      if (h.le.0.0d0) then
        epfint=0.0d0
        fint=0.0d0
      else
        el1=elab(e1,u1,ewp,lcm,beta)
        el2=elab(e2,u2,ewp,lcm,beta)
        h2=0.5d0*h
        if (lep.eq.1) then
          fint=h*f1
          y0(1)=h2*(el1+el2)*f1
        else
          fint=h2*(f1+f2)
          y0(1)=h2*(el1*f1+el2*f2)
        endif
        h=h2
        jmax=1
        do i=1,nmax-1
          sumy=0.5d0*y0(1)
          do j=1,jmax
            ej=e1+(2*j-1)*h
            call terp1m(e1,u1,e2,u2,ej,uj,lep)
            elj=elab(ej,uj,ewp,lcm,beta)
            call terp1m(e1,f1,e2,f2,ej,fj,lep)
            sumy=sumy+h*elj*fj
          enddo
          y1(1)=sumy
          jmax=2*jmax
          pow=4.0d0
          do k=1,i
            den=pow-1.0d0
            y1(k+1)=y1(k)+(y1(k)-y0(k))/den
            pow=4.0d0*pow
          enddo
          i1=i+1
          if (abs(y1(i1)-y1(i)).le.eps*abs(y1(i1)).or.i1.eq.nmax) then
            exit
          else
            h=0.5d0*h
            do k=1,i1
              y0(k)=y1(k)
            enddo
          endif
        enddo
        epfint=y1(i1)
      endif
      return
      end
C======================================================================
C     General routines for processing
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
        if (x0.eq.x1.and.y0.eq.y1) then
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
c     pack value into 11-character string
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
      real*8 function edelta(ffin,fdig)
c
c      Generate the difference fdig in the least significant digit
c      assuming 11 character representation (ENDF-6)
c
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
      subroutine mtchkd(mt,nmf3,mf3,nmf4,mf4,nmf6,mf6,nmf14,mf14,
     &  mtneu,mtgam)
c
c      Check posible use of summation vs partials
c
      dimension mf3(*),mf4(*),mf6(*),mf14(*)
      if (mt.eq.16) then
        mtlow=875
        mthigh=891
      elseif (mt.eq.103) then
        mtlow=600
        mthigh=649
      elseif (mt.eq.104) then
        mtlow=650
        mthigh=699
      elseif (mt.eq.105) then
        mtlow=700
        mthigh=749
      elseif (mt.eq.106) then
        mtlow=750
        mthigh=799
      elseif (mt.eq.107) then
        mtlow=800
        mthigh=849
      else
        write(*,*)' ERROR: MT=',mt,' not allowed as summa reaction'
        stop
      endif
      mtln=0
      mtlg=0
      do mti=mtlow,mthigh
        m3=iposm(nmf3,mf3,mti)
        m4=iposm(nmf4,mf4,mti)
        m6=iposm(nmf6,mf6,mti)
        m14=iposm(nmf14,mf14,mti)
        m46=m4+m6
        m146=m14+m6
        if (m3.gt.0) then
          if (m46.gt.0) mtln=mtln+1
          if (m146.gt.0) mtlg=mtlg+1
        endif
      enddo
      if (mtln.gt.0) then
        mtneu=1
      else
        mtneu=0
      endif
      if (mtlg.gt.0) then
        mtgam=1
      else
        mtgam=0
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
      function mtvalid(izai,mt)
c
c      Valid MT numbers for explicit photonuclear data processing
c
      if (izai.eq.0) then
        select case (mt)
          case (1:4,5,11,16:18,22:25,28:30,32:37,41:42,44:45,
     &          50:91,102:109,111:119,152:207,452,455,456,
     &          600:649,650:699,700:749,750:799,800:849,875:891)
            mtvalid=1
          case default
            mtvalid=0
        end select
      endif
      return
      end
C======================================================================
      function numpart(izp,mt)
      if (izp.eq.1) then
c
c       Neutron yield (izp=1)
c
        select case (mt)
          case (18)
            numpart=19
          case (2,5,22,23,28,29,32:36,44,45,50:91,158,183:189,198)
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
      elseif (izp.eq.0) then
        select case (mt)
          case (50:90,600:648,650:698,700:748,750:798,800:848,875:890)
            numpart=19
          case default
            numpart=0
        end select
      else
        numpart=0
      endif
      return
      end
C======================================================================
      function id2body(mt)
c
c      Allowed two-body reactions for photonuclear data processing
c
        select case (mt)
          case (2,50:90,102:107,600:648,650:698,700:748,750:798,800:848)
            id2body=1
          case default
            id2body=0
        end select
      return
      end
C======================================================================
      function mtdos(idos,mt,izap,lfs,i)
      parameter (nprodmax=504)
c
c     Assign production reactions MT-numbers
c
c       nprodmax=504 Maximum production reactions from MF10&MF6*MF3
c
c       MTD=219-451 (N1=232) i=  1..232 =>  MTD=218+i
c       MTD=461-599 (N2=139) i=233..371 =>  MTD=228+i
c       MTD=850-874 (N3=25)  i=372..396 =>  MTD=478+i
c       MTD=892-999 (N4=108) i=397..504 =>  MTD=495+i
c
      if (idos.eq.2) then
        i=i+1
        if (i.le.232) then
          mtdos=218+i
        elseif (i.le.371) then
          mtdos=228+i
        elseif (i.le.396) then
          mtdos=478+i
        elseif (i.le.nprodmax) then
          mtdos=495+i
        else
          write(*,*)' Error: too many production reactions ',i,
     &      ' > ',nprodmax
          stop
        endif
      else
        if (mt.eq.5) then
          mtdos=1000000*(50+lfs)+izap
        elseif (mt.eq.18) then
          mtdos=1000000*(80+lfs)+izap
        else
          mtdos=1000*(10+lfs)+mt
        endif
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
      subroutine wrtcont(lib,mat,mf,mt,ns,c1,c2,l1,l2,n1,n2)
c
c     write CONT record
c
      implicit real*8 (a-h, o-z)
      character*11 sc1,sc2
      call chendf(c1,sc1)
      call chendf(c2,sc2)
      ns=ns+1
      if (ns.gt.99999) ns=0
      write(lib,10)sc1,sc2,l1,l2,n1,n2,mat,mf,mt,ns
      return
   10 format(2a11,4i11,i4,i2,i3,i5)
      end
C======================================================================
      subroutine wrtsend(lib,mat,mf,ns)
c
c     write SECTION END record
c
      ns=0
      write(lib,'(66x,i4,i2,i3,i5)')mat,mf,0,99999
      return
      end
C======================================================================
      subroutine wrtfend(lib,mat,ns)
c
c     Write FILE END record
c
      ns=0
      write(lib,'(66x,i4,i2,i3,i5)')mat,0,0,ns
      return
      end
C======================================================================
      subroutine wrtmend(lib,ns)
c
c     write MATERIAL END record
c
      ns=0
      write(lib,'(66x,i4,i2,i3,i5)')0,0,0,ns
      return
      end
C======================================================================
      subroutine wrtend(lib,ns)
c
c     write TAPE END record
c
      ns=0
      write(lib,'(66x,i4,i2,i3,i5)')-1,0,0,ns
      return
      end
C======================================================================
      subroutine wrtlist(lib,mat,mf,mt,ns,c1,c2,l1,l2,npl,n2,b)
c
c     write LIST record
c
      implicit real*8 (a-h, o-z)
      character*11 rec(6)
      dimension b(*)
      call wrtcont(lib,mat,mf,mt,ns,c1,c2,l1,l2,npl,n2)
      nlm=npl/6
      nlr=npl-nlm*6
      if (nlm.gt.0) then
        do i=1,nlm
          do j=1,6
            k=6*(i-1)+j
            call chendf(b(k),rec(j))
          enddo
          ns=ns+1
          if (ns.gt.99999) ns=0
          write(lib,10)(rec(j),j=1,6),mat,mf,mt,ns
        enddo
      endif
      if (nlr.ne.0) then
        do j=1,nlr
          k=6*nlm+j
          call chendf(b(k),rec(j))
        enddo
        i=nlr+1
        do j=i,6
          rec(j)='           '
        enddo
        ns=ns+1
        if (ns.gt.99999) ns=0
        write(lib,10)(rec(j),j=1,6),mat,mf,mt,ns
      endif
      return
   10 format(6a11,i4,i2,i3,i5)
      end
C======================================================================
      subroutine wrtab1(lib,mat,mf,mt,ns,c1,c2,l1,l2,nr,nbt,inr,np,x,y)
c
c     write TAB1 record
c
      implicit real*8 (a-h, o-z)
      character*11 strx(3),stry(3)
      dimension nbt(*),inr(*),x(*),y(*)
      dimension nn(3),ii(3)
      call wrtcont(lib,mat,mf,mt,ns,c1,c2,l1,l2,nr,np)
      nlm=nr/3
      nlr=nr-3*nlm
      if (nlm.gt.0) then
        do i=1,nlm
          do j=1,3
            k=3*(i-1)+j
            nn(j)=nbt(k)
            ii(j)=inr(k)
          enddo
          ns=ns+1
          if (ns.gt.99999) ns=0
          write(lib,10)(nn(j),ii(j),j=1,3),mat,mf,mt,ns
        enddo
      endif
      if (nlr.ne.0) then
        do j=1,nlr
          k=3*nlm+j
          nn(j)=nbt(k)
          ii(j)=inr(k)
        enddo
        ns=ns+1
        if (ns.gt.99999) ns=0
        if (nlr.eq.1) then
          write(lib,11)(nn(j),ii(j),j=1,nlr),mat,mf,mt,ns
        else
          write(lib,12)(nn(j),ii(j),j=1,nlr),mat,mf,mt,ns
        endif
      endif
      nlm=np/3
      nlr=np-3*nlm
      if (nlm.gt.0) then
        do i=1,nlm
          do j=1,3
            k=3*(i-1)+j
            call chendf(x(k),strx(j))
            call chendf(y(k),stry(j))
          enddo
          ns=ns+1
          if (ns.gt.99999) ns=0
          write(lib,20)(strx(j),stry(j),j=1,3),mat,mf,mt,ns
        enddo
      endif
      if (nlr.ne.0) then
        do j=1,nlr
          k=3*nlm+j
          call chendf(x(k),strx(j))
          call chendf(y(k),stry(j))
        enddo
        i=nlr+1
        do j=i,3
          strx(j)='           '
          stry(j)='           '
        enddo
        ns=ns+1
        if (ns.gt.99999) ns=0
        write(lib,20)(strx(j),stry(j),j=1,3),mat,mf,mt,ns
      endif
      return
   10 format(6i11,i4,i2,i3,i5)
   11 format(2i11,44x,i4,i2,i3,i5)
   12 format(4i11,22x,i4,i2,i3,i5)
   20 format(6a11,i4,i2,i3,i5)
      end
C======================================================================
      subroutine wrtab2(lib,mat,mf,mt,ns,c1,z,l1,l2,nr,nz,nbt,inr)
c
c     write TAB2 record
c
      implicit real*8 (a-h, o-z)
      dimension nbt(*),inr(*)
      dimension nn(3),ii(3)
      call wrtcont(lib,mat,mf,mt,ns,c1,z,l1,l2,nr,nz)
      nlm=nr/3
      nlr=nr-3*nlm
      if (nlm.gt.0) then
        do i=1,nlm
          do j=1,3
            k=3*(i-1)+j
            nn(j)=nbt(k)
            ii(j)=inr(k)
          enddo
          ns=ns+1
          if (ns.gt.99999) ns=0
          write(lib,10)(nn(j),ii(j),j=1,3),mat,mf,mt,ns
        enddo
      endif
      if (nlr.ne.0) then
        do j=1,nlr
          k=3*nlm+j
          nn(j)=nbt(k)
          ii(j)=inr(k)
        enddo
        ns=ns+1
        if (ns.gt.99999) ns=0
        if (nlr.eq.1) then
          write(lib,11)(nn(j),ii(j),j=1,nlr),mat,mf,mt,ns
        else
          write(lib,12)(nn(j),ii(j),j=1,nlr),mat,mf,mt,ns
        endif
      endif
   10 format(6i11,i4,i2,i3,i5)
   11 format(2i11,44x,i4,i2,i3,i5)
   12 format(4i11,22x,i4,i2,i3,i5)
      return
      end
C=============================================================================
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
C=============================================================================
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
      elseif (law.lt.-15.or.law.gt.7) then
       write(*,*)' ERROR: unknown LAW=',law,' on MF6'
       stop
      endif
      return
      end
C======================================================================
      subroutine readxs(lib,za,awr,qm,qi,lr,nr,nbt,ibt,np,x,y)
c
c      read a cross section(MT section) on MF3
c
      implicit real*8 (a-h, o-z)
      dimension nbt(*),ibt(*),x(*),y(*)
      call readcont(lib,za,awr,l1,l2,n1,n2,mat0,mf0,mt0,ns0)
      call readtab1(lib,qm,qi,l1,lr,nr,np,nbt,ibt,x,y)
      return
      end
C======================================================================
C       Get date and time
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
c       2. subroutine terp1 (renamed as terp1m)
C       3. real*8 function bachaa
C       4. subroutine phnout (renamed as change)
C       5. subroutine advance_to_locator
C       6. subroutine write_integer
C       7. subroutine write_real
C       8. subroutine write_integer_list
C       9. subroutine write_real_list
C      10. subroutine typen
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
      real*8 function bachaa(iza1i,iza2,izat,e,ep)
c     ******************************************************************
c     compute the kalbach a parameter (adapted for ACEMAKER)
c     ******************************************************************
c     (adapted by D. Lopez Aldama for ACEMAKER)
c
      implicit real*8 (a-h,o-z)
      real*8 nc,nb
      parameter (emnc2=939.56542052539d+6)
      data third,twoth,fourth/.333333333d0,.666666667d0,1.33333333d0/
      data c1,c2,c3,c4,c5,c6/15.68d0,-28.07d0,-18.56d0,33.22d0,
     &  -0.717d0,1.211d0/
      data s2,s3,s4,s5/2.22d0,8.48d0,7.72d0,28.3d0/
      data b1,b2,b3/0.04d0,1.8d-6,6.7d-7/
      data d1/9.3d0/
      data ea1,ea2/41.d0,130.d0/
      data emev/1.d6/
      emc2=emnc2/emev
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
      subroutine change(iou,nout,mcnpx)
c     (adapted by D. Lopez Aldama for ACEMAKER)
      implicit real*8 (a-h,o-z)
      external typen
      character*10 hd,hm
      character*13 hz
      character*70 hk
      integer esz, tot, els, thn, sig
      integer pxs, phn, tyrp, sigp, andp, dlwp
      integer rlocator, plocator, ielocator, oelocator
      common/acetxt/hz,hd,hm,hk
      common/acecte/awr0,tz,awn(16),izn(16)
      common/acepnt/nxs(16),jxs(32)
      common/acedat/xss(200000000),nxss
c
c     nout>1: write type 1 ACE-formatted file on unit nout
c
c     header block
c
      if (mcnpx.eq.1) then
        write(nout,'(a13,f12.6,1pe12.4,1x,a10)')hz(1:13),awr0,tz,hd
      else
        write(nout,'(a10,f12.6,1pe12.4,1x,a10)')hz(1:10),awr0,tz,hd
      endif
      write(nout,'(a70,a10)')hk,hm
      write(nout,'(4(i7,f11.0))')(izn(i),awn(i),i=1,16)
      write(nout,'(8i9)')(nxs(i),i=1,16)
      write(nout,'(8i9)')(jxs(i),i=1,32)
c
c     assign flags and triggers
c
      nes=nxs(3)
      ntr=nxs(4)
      ntype=nxs(5)
      npixs=nxs(6)
      neixs=nxs(7)
      esz=jxs(1)
      tot=jxs(2)
      non=jxs(3)
      els=jxs(4)
      thn=jxs(5)
      mtr=jxs(6)
      lqr=jxs(7)
      lsig=jxs(8)
      sig=jxs(9)
      ixsa=jxs(10)
c
c       -- esz block
c
      l=1
      call advance_to_locator(iou,nout,l,esz)
      call write_real_list(nout,l,2*nes)
c
c       -- non, els blocks
c
      if (els.ne.0) then
c       -- non block
       call advance_to_locator(iou,nout,l,non)
       call write_real_list(nout,l,nes)
c       -- els block
       call advance_to_locator(iou,nout,l,els)
       call write_real_list(nout,l,nes)
      endif
c
c       -- thn block
c
      if (thn.ne.0) then
        call advance_to_locator(iou,nout,l,thn)
        call write_real_list(nout,l,nes)
      endif
c
c       -- mtr block
c
      call advance_to_locator(iou,nout,l,mtr)
      call write_integer_list(nout,l,ntr)
c
c       -- lqr block
c
      call advance_to_locator(iou,nout,l,lqr)
      call write_real_list(nout,l,ntr)
c
c       -- lsig block
c
      call advance_to_locator(iou,nout,l,lsig)
      rlocator=l
      call write_integer_list(nout,l,ntr)
c
c      --sig block
c
      call advance_to_locator(iou,nout,l,sig)
      do i=1,ntr
        call advance_to_locator(iou,nout,l,sig+nint(xss(rlocator))-1)
        call write_integer(nout,l)
        ne=nint(xss(l))
        call write_integer(nout,l)
        call write_real_list(nout,l,ne)
        rlocator=rlocator+1
      enddo
c
c      -- particle production blocks
c
      if (ntype.gt.0) then
c
c       -- ixs arrays
c
        call advance_to_locator(iou,nout,l,ixsa)
        plocator=l
        call write_integer_list(nout,l,neixs*ntype)
c
c       -- loop over the ntype productions
c
        do ip=1,ntype
c
c         -- IXS array entries
c
          ipt=nint(xss(plocator))
          ntrp=nint(xss(plocator+1))
          pxs=nint(xss(plocator+2))
          phn=nint(xss(plocator+3))
          mtrp=nint(xss(plocator+4))
          tyrp=nint(xss(plocator+5))
          lsigp=nint(xss(plocator+6))
          sigp=nint(xss(plocator+7))
          landp=nint(xss(plocator+8))
          andp=nint(xss(plocator+9))
          ldlwp=nint(xss(plocator+10))
          dlwp=nint(xss(plocator+11))
c
c         -- pxs block
c
          call advance_to_locator(iou,nout,l,pxs)
          call write_integer(nout,l)
          ne=nint(xss(l))
          call write_integer(nout,l)
          call write_real_list(nout,l,ne)
c
c         -- phn block
c
          call advance_to_locator(iou,nout,l,phn)
          call write_integer(nout,l)
          ne=nint(xss(l))
          call write_integer(nout,l)
          call write_real_list(nout,l,ne)
c
c         -- mtrp block
c
          call advance_to_locator(iou,nout,l,mtrp)
          call write_integer_list(nout,l,ntrp)
c
c         --tyrp block
c
          call advance_to_locator(iou,nout,l,tyrp)
          call write_integer_list(nout,l,ntrp)
c
c         --lsigp block
c
          call advance_to_locator(iou,nout,l,lsigp)
          rlocator=l
          call write_integer_list(nout,l,ntrp)
c
c         --sigp block
c
          call advance_to_locator(iou,nout,l,sigp)
          do i=1,ntrp
            irsigp=sigp+nint(xss(rlocator))-1
            call advance_to_locator(iou,nout,l,irsigp)
            mftype=nint(xss(l))
            call write_integer(nout,l)
            if (mftype.eq.13) then
c
c             MFTYPE=13
c
              call write_integer(nout,l)
              ne=nint(xss(l))
              call write_integer(nout,l)
              call write_real_list(nout,l,ne)
            else
c
c             MFTYPE=6,12
c
              call write_integer(nout,l)
              nr=nint(xss(l))
              call write_integer(nout,l)
              if (nr.gt.0) then
                 call write_integer_list(nout,l,2*nr)
              endif
              ne=nint(xss(l))
              call write_integer(nout,l)
              call write_real_list(nout,l,2*ne)
            endif
            rlocator=rlocator+1
          enddo
c
c         -- landp block
c
          call advance_to_locator(iou,nout,l,landp)
          rlocator=l
          call write_integer_list(nout,l,ntrp)
c
c         -- andp block
c
          call advance_to_locator(iou,nout,l,andp)
          do i=1,ntrp
            nn=nint(xss(rlocator))
            if (nn.gt.0) then
              call advance_to_locator(iou,nout,l,andp+nn-1)
              ne=nint(xss(l))
              call write_integer(nout,l)
              call write_real_list(nout,l,ne)
              ielocator=l
              call write_integer_list(nout,l,ne)
              do j=1,ne
                nn=nint(xss(ielocator))
                if (nn.ne.0) then
                  call advance_to_locator(iou,nout,l,andp+abs(nn)-1)
                  if (nn.gt.0) then
                    call write_real_list(nout,l,33)
                  else if (nn.lt.0) then
                    call write_integer(nout,l)
                    np=nint(xss(l))
                    call write_integer(nout,l)
                    call write_real_list(nout,l,3*np)
                  endif
                endif
                ielocator=ielocator+1
              enddo
            endif
            rlocator=rlocator+1
          enddo
c
c         -- ldlwp block
c
          call advance_to_locator(iou,nout,l,ldlwp)
          rlocator=l
          call write_integer_list(nout,l,ntrp)
c
c         -- dlwp block
c
          call advance_to_locator(iou,nout,l,dlwp)
          do i=1,ntrp
             nn=nint(xss(rlocator))
             if (nn.gt.0) then
               call advance_to_locator(iou,nout,l,dlwp+nn-1)
               lnw=1
               do while (lnw.ne.0)
c
c                law header
c
                 lnw=nint(xss(l))
                 call write_integer(nout,l)
                 law=nint(xss(l))
                 call write_integer(nout,l)
                 call write_integer(nout,l)
                 nr=nint(xss(l))
                 call write_integer(nout,l)
                 if (nr.gt.0) then
                   call write_integer_list(nout,l,2*nr)
                 endif
                 ne=nint(xss(l))
                 call write_integer(nout,l)
                 call write_real_list(nout,l,2*ne)
c
c                -- law 4
c
                 if (law.eq.4) then
                   nr=nint(xss(l))
                   call write_integer(nout,l)
                   if (nr.gt.0) then
                     call write_integer_list(nout,l,2*nr)
                   endif
                   ne=nint(xss(l))
                   call write_integer(nout,l)
                   call write_real_list(nout,l,ne)
                   ielocator=l
                   call write_integer_list(nout,l,ne)
                   do j=1,ne
                     iedlwp=dlwp+nint(xss(ielocator))-1
                     call advance_to_locator(iou,nout,l,iedlwp)
                     call write_integer(nout,l)
                     np=nint(xss(l))
                     call write_integer(nout,l)
                     call write_real_list(nout,l,3*np)
                     ielocator=ielocator+1
                   enddo
c
c                -- law 44
c
                 else if (law.eq.44) then
                   nr=nint(xss(l))
                   call write_integer(nout,l)
                   if (nr.gt.0) then
                     call write_integer_list(nout,l,2*nr)
                   endif
                   ne=nint(xss(l))
                   call write_integer(nout,l)
                   call write_real_list(nout,l,ne)
                   ielocator=l
                   call write_integer_list(nout,l,ne)
                   do j=1,ne
                     iedlwp=dlwp+nint(xss(ielocator))-1
                     call advance_to_locator(iou,nout,l,iedlwp)
                     call write_integer(nout,l)
                     np=nint(xss(l))
                     call write_integer(nout,l)
                     call write_real_list(nout,l,5*np)
                     ielocator=ielocator+1
                   enddo
c
c                -- law 61
c
                 else if (law.eq.61) then
                   nrr=nint(xss(l))
                   call write_integer(nout,l)
                   if (nrr.gt.0) then
                     call write_integer_list(nout,l,2*nrr)
                   endif
                   ne=nint(xss(l))
                   call write_integer(nout,l)
                   call write_real_list(nout,l,ne)
                   ielocator=l
                   call write_integer_list(nout,l,ne)
                   do j=1,ne
                     iedlwp=dlwp+nint(xss(ielocator))-1
                     call advance_to_locator(iou,nout,l,iedlwp)
                     call write_integer(nout,l)
                     np=nint(xss(l))
                     call write_integer(nout,l)
                     call write_real_list(nout,l,3*np)
                     oelocator=l
                     call write_integer_list(nout,l,np)
                     do k=1,np
                       jedlwp=dlwp+nint(xss(oelocator))-1
                       call advance_to_locator(iou,nout,l,jedlwp)
                       call write_integer(nout,l)
                       nmu=nint(xss(l))
                       call write_integer(nout,l)
                       call write_real_list(nout,l,3*nmu)
                       oelocator=oelocator+1
                     enddo
                     ielocator=ielocator+1
                   enddo
c
c                -- law 7 or 9
c
                 else if (law.eq.7.or.law.eq.9) then
                   nr=nint(xss(l))
                   call write_integer(nout,l)
                   if (nr.gt.0) then
                      call write_integer_list(nout,l,2*nr)
                   endif
                   ne=nint(xss(l))
                   call write_integer(nout,l)
                   call write_real_list(nout,l,2*ne)
                   call write_real(nout,l)
c
c                -- law 33
c
                 else if (law.eq.33) then
                   call write_real(nout,l)
                   call write_real(nout,l)
c
c                -- unknown law
c
                 else
                   write(iou,*)' Undefined law=',law,' for NLIB=0',
     &               ' (photonuclear data)'
                   stop
                 endif
               enddo
             endif
             rlocator=rlocator+1
          enddo
          plocator=plocator+neixs
c
c         -- continue loop over productions
c
        enddo
      else
        write(iou,*)' ERROR: None outgoing particle ntype=0. Check data'
        stop
      endif
c
c     Clear typen buffer
c
      call typen(0,nout,3)
      return
      end
C======================================================================
      subroutine advance_to_locator(iou,nout,l,locator)
c     (adapted by D. Lopez Aldama for ACEMAKER)
c     -----------------------------------------------------------------
c     Advance to the next locator position from the current position l.
c     If the current position is not equal to the locator position, the
c     function will advance l until it is equal to the locator position.
c     It will write the values in the xss array while advancing to the
c     new position.
c     -----------------------------------------------------------------
      if (l.lt.locator) then
        write(iou,'(a,i6,a,i6)')' Warning: expected xss index ',
     &    locator,' greater than current index ',l
        write(iou,'(a)')' xss array was padded accordingly'
        do while (l.lt.locator)
          call typen(l,nout,1)
          l=l+1
        enddo
      else if (l.gt.locator) then
        write(iou,'(a,i6,a,i6)')' Warning: expected xss index ',
     &    locator,' less than current index ',l
        write(iou,'(a)')' This may be a serious problem'
        stop
      endif
      return
      end
C======================================================================
      subroutine write_integer(nout,l)
c     -----------------------------------------------------------------
c     Write an integer value at the position l, and advance l to the
c     next position
c     -----------------------------------------------------------------
      call typen(l,nout,1)
      l=l+1
      return
      end
C======================================================================
      subroutine write_real(nout,l)
c     -----------------------------------------------------------------
c     Write a real value at the position l, and advance l to the
c     next position
c     -----------------------------------------------------------------
      call typen(l,nout,2)
      l=l+1
      return
      end
C======================================================================
      subroutine write_integer_list(nout,l,n)
c     -----------------------------------------------------------------
c     Write n integer values from position l, and advance l to the
c     next position
c     -----------------------------------------------------------------
      do i=1,n
        call typen(l,nout,1)
        l=l+1
      enddo
      return
      end
C======================================================================
      subroutine write_real_list(nout,l,n)
c     -----------------------------------------------------------------
c     Write n real values from position l, and advance l to the
c     next position
c     -----------------------------------------------------------------
      do i=1,n
        call typen(l,nout,2)
        l=l+1
      enddo
      return
      end
C======================================================================
      subroutine typen(l,nout,iflag)
c     (adapted by D. Lopez Aldama for ACEMAKER)
c     *****************************************************************
c     write an integer or a real number to a type-1 ace file
c     or (if nout=0) convert real to integer for type-3 output
c     or (if nout=1) convert integer to real for type-3 input
c     use iflag.eq.1 to write an integer (i20)
c     use iflag.eq.2 to write a real number (1pe20.11)
c     use iflag.eq.3 to write partial line at end of file
c     *****************************************************************
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
