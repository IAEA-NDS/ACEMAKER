      program dodos
c     version 2.0
c
c     prepare dosimetry ace-formatted files for MCNP
c     the code belong to the ACEMAKER code system
c
c     Input data:
c
c     line 1:  input endf-6 formatted filename                     (A72)
c     line 2:  output dosimetry ace-formatted filename             (A72)
c     line 3:        mat       imon       idos6                   (3I11)
c     line 4:        tol       ymin                             (2E11.0)
c     line 5:       suff      mcnpx                              (7X,A4)
c
c     where,
c        mat: Requested material
c       imon: Monitor printing trigger (0/1/2) = min/max/max+plt
c             (Default: imon=0)
c      idos6: Force to produce dosimetry data from MF6 yields even if
c             MF8/LMF=6 data are not available
c             (default: idos6=0) No process if MF8/LMF=6 not available
c        tol: linearization tolerance (Default: tol=0.001)
c       ymin: minimum cross section (Default: ymin=1.0e-20)
c       suff: ZAID suffix for ACE-formatted file
c             Examples: .00, .32, .80, .067, the dot '.' is required
c             (Default: suff=.00)
c      mcnpx: MCNP trigger (0/1) = MCNP/MCNPX (Default mcnpx=0)
c
c     Output files:
c       1. Dosimetry ACE-formatted file
c       2. DODOS.LST listing file (fix name)
c      if imon=2
c       3. DODOS.PLT PLOTTAB input option file
c       4. DODOS.CUR PLOTTAB curve file
c
c     MTD: MT numbers for the MTR block in the dosimetry ACE-file
c       If the reaction is given in MF3, then the MT number remains the
c       same (MTD=MT). If the dosimetry reaction is given in MF10 or
c       multiplicities are supplied in MF9 or MF6, then the MTD number
c       is computed as:
c         if MT==5 then
c           MTD=1000000*(50+lfs)+zap
c         elseif MT==18
c           if (izap==-1) then
c             MTD=18
c           else
c             MTD=1000000*(80+lfs)+zap
c           endif
c         else
c           MTD=1000*(10+lfs)+MT
c         endif
c
c       where,
c         MTD: reaction identifier in the dosimetry ACE-formatted file
c         zap: ZA number of the product nuclide
c         lfs: level number of the product nuclide
c         MT:  reaction identifier according to the ENDF-6 format
c
c       For MT=5 and (MT=18 with zap!=-1) the zap is codified in MTD,
c       but for the rest of reactions the recoil is implicitly given
c       according to the reaction MT number
c
c
c     Example of input
c
c     \IRDFF-II\IRDFF-II.endf
c     \DOS\ZA013027.acef
c             1325          0          1
c            0.001 1.0000E-20
c              .00          0
c
c     In this example the material 1325 (Al-27) is retrived from the
c     \IRDFF-II\IRDFF-II.endf tape and the dosimetry ace-formatted file
c     \DOS\ZA013027.acef is generated. Minimun printout is reqquested.
c     The linearization tolerance is 0.1% and the minimum cross section
c     allowed is 1.0E-20 barn. They are applied in cases where the files
c     MF3/MF9 or MF3/MF6 have been used for describing a dosimetry
c     reaction. The ZAID suffix will be .00y taking into account that
c     the input values of suff and mcnpx are .00 and 0 respectively.
c     Dosimetry data from MF6 yields are generated even if MF8/LMF=6
c     data are not available.
c
      implicit real*8 (a-h, o-z)
      parameter (nnxc=450,nnx6=2000,npmax=2000000,nxssmx=50000000)
      parameter (bk=8.6173303d-11)
      character*1 ch
      character*4 suff
      character*10 hd,hm,cdate
      character*11 zsymam,str11,ctime
      character*13 hz
      character*66 line
      character*70 hk
      character*72 fin1,fout
      common/acetxt/hz,hd,hm,hk
      common/acecte/awr0,tz,awn(16),izn(16)
      common/acepnt/nxs(16),jxs(32)
      common/acedat/xss(50000000),nxss
      dimension mf3(nnxc),mf9(nnxc),mf10(nnxc)
      dimension lzap6(nnx6),lip6(nnx6),mt6(nnx6),mf6(nnx6),kip6(nnx6)
      dimension nzap6(nnx6),nzap(nnx6)
      dimension nbt(20),ibt(20),nbty(20),ibty(20)
      dimension x(npmax),xy(npmax)
      dimension y(npmax),yy(npmax)
      data ev2mev/1.0d-6/
      data suff/'    '/
      data in1/2/,nin/3/,ntp/10/,nou/20/,lou/21/,ncur/30/,nplt/31/
c
c      Initialize pointers and triggers
c
      nxss=nxssmx
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
c
c     write heading and open files DODOS.INP and DODOS.LST
c
      open(in1, file='DODOS.INP')
      open(lou, file='DODOS.LST')
      call getdtime(cdate,ctime)
      write(lou,'(a)')' PROGRAM DODOS v2.0: Prepare dosimetry ACE-files'
      write(lou,'(a)')' ==============================================='
      write(lou,*)
      write(lou,'(a,a11,a,a10)')' Started at ',ctime,' on ',cdate
      write(lou,*)
      write(*,'(a)')' PROGRAM DODOS v2.0: Prepare dosimetry ACE-files'
      write(*,'(a)')' ==============================================='
      write(*,'(a,a11,a,a10)')' Started at ',ctime,' on ',cdate
      write(*,'(a)')
c
c     read input data from DODOS.INP
c
      read(in1,'(a)')fin1
      read(in1,'(a)')fout
      read(in1,'(3i11)')mat,imon,idos6
      if (imon.lt.0) then
        imon=0
      elseif (imon.gt.2) then
        imon=2
      endif
      if (idos6.gt.0) then
        idos6=1
      else
        idos6=0
      endif
      read(in1,'(2e11.0)')tol,ymin
      if (tol.le.0.0d0) tol=1.0d-3
      if (ymin.le.0.0d0) ymin=1.0d-30
      read(in1,'(e11.0,i11)')xsuff,mcnpx
      if (xsuff.lt.0.0d0.or.xsuff.gt.1.0d0) xsuff=0.0d0
      if (mcnpx.ne.1) mcnpx=0
      write(suff,'(a1,i3)')'.',nint(1000*xsuff)
      if (suff(2:2).eq.' ') suff(2:2)='0'
      if (suff(3:3).eq.' ') suff(3:3)='0'
      close(in1)
c
c     print out input data
c
      write(lou,'(a)')' Input parameters'
      write(lou,'(a)')' ================'
      write(lou,'(a,a)')' Input file name   = ',fin1
      write(lou,'(a,a)')' Output file name  = ',fout
      write(lou,'(a,i5)')' ENDF material(MAT)=',mat
      write(lou,'(a,i2)')' Printing option   =',imon
      write(lou,'(a,i2)')' Force use MF6 data=',idos6
      write(lou,'(a,1pe13.6)')' Tolerance [%]     =',tol*100.0d0
      write(lou,'(a,1pe13.6)')' Minimum XS [barn] =',ymin
      if (mcnpx.eq.1) then
        write(lou,'(a,a4)')' ZAID suffix(suff) = ',suff(1:4)
      else
        write(lou,'(a,a3)')' ZAID suffix(suff) = ',suff(1:3)
      endif
      write(lou,'(a,i2)')' MCNP trigger      =',mcnpx
c
c      Read tape header and find selected material
c      mat=0 implies first material on tape
c
      open (nin, file=fin1)
      read(nin,'(a66)')line
      write(lou,*)
      write(lou,'(a,a66)')' TAPE ID. = ',line
      if (mat.ne.0) then
        call findmat(nin,mat,icod)
        if (icod.ne.0) then
          write(*,*)' Material =',mat,' not found on tape ',fin1
          write(lou,*)' Material =',mat,' not found on tape ',fin1
          close(nin)
          close(lou)
          stop
        endif
      endif
c
c      Material to be processed. Reading MF1/MT451
c
      call readcont(nin,za0,awr0,lrp,lfi,nlib,n2,mat0,mf0,mt0,ns0)
      if (mat.eq.0) mat=mat0
      call readcont(nin,elis,sta,lis,liso,n1,nfor,mat0,mf0,mt0,ns0)
      call readcont(nin,awi,emax,lrel,l2,nsub,nver,mat0,mf0,mt0,ns0)
      call readcont(nin,temp,c2,ldrv,l2,nwd,nxc,mati,mf,mt,nsi)
c
c      Check for incident particle type (library type)
c
      if (nsub.eq.10.or.nsub.eq.19) then
        zai=1.0d0
        izai=1
        ch='n'
      elseif (nsub.eq.10010) then
        zai=1001.0d0
        izai=1001
        ch='h'
      elseif (nsub.eq.10020) then
        zai=1002.0d0
        izai=1002
        ch='o'
      elseif (nsub.eq.10030) then
        zai=1003.0d0
        izai=1003
        ch='r'
      elseif (nsub.eq.20030) then
        zai=2003.0d0
        izai=2003
        ch='s'
      elseif (nsub.eq.20040) then
        zai=2004.0d0
        izai=2004
        ch='a'
      elseif (nsub.eq.0) then
        zai=0.0d0
        izai=0
        ch='u'
      elseif (nsub.eq.3) then
        zai=0.0d0
        izai=0
        ch='p'
      else
        write(lou,*)'  === Error: incident particle is not coded'
        write(lou,*)'  === NSUB=',nsub
        close(nin)
        close(lou)
        stop
      endif
c
c     Prepare ACE-file heading information
c
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
      call readtext(nin,line,mat0,mf0,mt0,ns0)
      zsymam=line(1:11)
      call readtext(nin,line,mat0,mf0,mt0,ns0)
      call getdtime(hd,str11)
      str11=' '
      call readtext(nin,line,mat0,mf0,mt0,ns0)
      hk(1:11)=zsymam
      hk(12:14)=' T='
      if (temp.lt.1.0d6) then
        write(str11,'(f11.2)')temp
      else
        write(str11,'(1pe11.4)')temp
      endif
      i=index(str11,' ',.true.)
      if (i.eq.0) i=1
      k=11-i
      hk(15:15+k)=str11(i:11)
      hk(16+k:23+k)=' K from '
      hk(24+k:41+k)=line(5:22)
      hk(42+k:52+k)=' (ACEMAKER)'
      hk(63:69)=' ipart='
      hk(70:70)=ch
      str11=' '
      if (mcnpx.eq.1) then
        write(hz,'(i6,a4,a1,a2)')izaid,suff(1:4),ch,'y '
      else
        write(hz,'(i6,a3,a4)')izaid,suff(1:3),'y   '
      endif
      write(hm,'(a6,i4)')'   mat',mat
      tz=bk*temp
      write(*,*)' Material=',mat
      write(lou,*)
      write(lou,'(a,i7,a,1pe15.8)')' ZAI=',izai,' AWI=',awi
      write(lou,'(a,i5,a,i7,a,i3,a,i7,a,1pe15.8)')' Material=',mat,
     &  ' ZA=',matza,' LISO=',liso,' ZAID=',izaid,' AWR=',awr0
      write(lou,'(a,a,a,i3,a,1pe15.8)')' SYM=',zsymam,
     &  ' LFI=',lfi,' ELIS=',elis
      write(lou,'(a,1p,e15.8,a,e13.6,a,e13.6,a)')' EMAX=', emax,
     &  ' Temperature=',temp,' K = ',tz,' MeV'
      write(lou,*)
c
c     screening MF3
c
      nmf3=0
      call findmf(nin,mat,3,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,3,mt)
          if (mt.gt.2) then
            nmf3=nmf3+1
            mf3(nmf3)=mt
          endif
        enddo
        if (nmf3.gt.nnxc) then
          write(lou,'(2(a,i6))')' Number of reactions greater than',
     &      nnxc,' in MF3. Increase the value of nnxc in ',nmf3-nnxc
          write(*,'(2(a,i6))')' Number of reactions greater than',
     &      nnxc,' in MF3. Increase the value of nnxc in ',nmf3-nnxc
        else
          write(lou,'(1x,i6,a)')nmf3,' reactions found on MF3'
          write(*,'(1x,i6,a)')nmf3,' reactions found on MF3'
        endif
      endif
c
c     screening MF10
c
      nmf10=0
      call findmf(nin,mat,10,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,10,mt)
          if (mt.gt.2) then
            call readcont(nin,c1,c2,l1,l2,nfs,n2,mat0,mf0,mt0,ns0)
            do i=1,nfs
              call readtab1(nin,c1,c2,izap,lfs,nr,ne,nbt,ibt,x,y)
              nmf10=nmf10+1
              mf10(nmf10)=mtdos(mt,izap,lfs)
            enddo
          endif
        enddo
        if (nmf10.gt.nnxc) then
          write(lou,'(2(a,i6))')' Number of reactions greater than',
     &      nnxc,' in MF10. Increase the value of nnxc in ',nmf10-nnxc
          write(*,'(2(a,i6))')' Number of reactions greater than',
     &      nnxc,' in MF10. Increase the value of nnxc in ',nmf10-nnxc
        else
          write(lou,'(1x,i6,a)')nmf10,' reactions found on MF10'
          write(*,'(1x,i6,a)')nmf10,' reactions found on MF10'
        endif
      endif
c
c     screening MF9 and saving yields
c
      open(ntp,file='yields.tmp')
      nmf9=0
      call findmf(nin,mat,9,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        do while (mt.gt.0)
          call findnextmt(nin,9,mt)
          if (mt.gt.2) then
            if (nmf3.gt.0) then
              i3=iposm(nmf3,mf3,mt)
            else
              i3=0
            endif
            if (i3.gt.0) then
              call readcont(nin,c1,c2,l1,l2,nfs,n2,mat0,mf0,mt0,ns0)
              do i=1,nfs
                call readtab1(nin,c1,c2,izap,lfs,nr,ne,nbt,ibt,x,y)
                mt9=mtdos(mt,izap,lfs)
                if (nmf10.gt.0) then
                  i10=iposm(nmf10,mf10,mt9)
                else
                  i10=0
                endif
                if (i10.eq.0) then
                  nmf9=nmf9+1
                  mf9(nmf9)=mt9
                  nst=0
                  c1=dble(mt)
                  c2=dble(mt9)
                  call wrtab1(ntp,mat,9,mt,nst,c1,c2,izap,lfs,
     &              nr,nbt,ibt,ne,x,y)
                else
                  write(lou,*)' MF9 yield for mt=',mt,' izap=',izap,
     &              ' lfs=',lfs,' was found on MF9 and MF10.',
     &              ' Data on MF9 ignored.'
                endif
              enddo
            else
              write(lou,*)' MF9 yield for mt=',mt,' was found on MF9,',
     &          ' but not on MF3. Reaction ignored.'
            endif
          endif
        enddo
        if (nmf9.gt.0) then
          if (nmf9.gt.nnxc) then
            write(lou,'(2(a,i6))')' Number of reactions greater than',
     &        nnxc,' in MF9. Increase the value of nnxc in ',nmf9-nnxc
            write(*,'(2(a,i6))')' Number of reactions greater than',
     &        nnxc,' in MF9. Increase the value of nnxc in ',nmf9-nnxc
          else
            write(lou,'(1x,i6,a)')nmf9,' set of yields found on MF9'
            write(*,'(1x,i6,a)')nmf9,' set of yields found on MF9'
          endif
        endif
      endif
c
c     screening MF8 for LMF=6 (MF6 yields)
c
      nmf68=0
      call findmf(nin,mat,8,icod)
      if (icod.eq.0) then
        backspace(nin)
        mt=1000
        ii=0
        do while (mt.gt.0)
          call findnextmt(nin,8,mt)
          if (mt.gt.2) then
            call readcont(nin,c1,c2,l1,l2,nfs,n2,mat0,mf0,mt0,ns0)
            kzap=-99999
            do i=1,nfs
              if (n2.eq.1) then
                call readcont(nin,zap,c2,lmf,lfs,n3,n4,mat0,mf0,mt0,ns0)
              else
                call readlist(nin,zap,c2,lmf,lfs,n3,n4,y)
              endif
              if (lmf.eq.6) then
                if (nmf3.gt.0) then
                  i3=iposm(nmf3,mf3,mt)
                else
                  i3=0
                endif
                if (i3.gt.0) then
                  izap=nint(zap+1.0d-6)
                  mt8=mtdos(mt,izap,lfs)
                  if (nmf10.gt.0) then
                    i10=iposm(nmf10,mf10,mt8)
                  else
                    i10=0
                  endif
                  if (nmf9.gt.0) then
                    i9=iposm(nmf9,mf9,mt8)
                  else
                    i9=0
                  endif
                  if (i9.eq.0.and.i10.eq.0) then
                    nmf68=nmf68+1
                    mt6(nmf68)=mt
                    lip6(nmf68)=lfs
                    lzap=izap
                    lzap6(nmf68)=lzap
                    if (kzap.ne.lzap) then
                      kk=0
                      kzap=lzap
                      ii=ii+1
                    else
                      kk=kk+1
                    endif
                    kip6(nmf68)=kk
                    nzap6(ii)=kk+1
                  endif
                else
                  write(lou,*)' MF6 yield for mt=',mt,' was specified',
     &              ' on MF8, but not on MF3. Data ignored.'
                endif
              endif
            enddo
          endif
        enddo
        if (nmf68.gt.0) then
          if (nmf68.gt.nnx6) then
            write(lou,'(2(a,i6))')' Number of reactions greater than',
     &        nnx6,' in MF6. Increase the value of nnx6 in ',nmf68-nnx6
            write(*,'(2(a,i6))')' Number of reactions greater than',
     &        nnx6,' in MF6. Increase the value of nnx6 in ',nmf68-nnx6
          else
            write(lou,'(1x,i6,a)')nmf68,
     &        ' set of yields found on MF8 with LMF=6'
            write(*,'(1x,i6,a)')nmf68,
     &        ' set of yields found on MF8 with LMF=6'
          endif
        endif
      endif
c
c     Force dosimetry data from MF6 yields if idos6>0 and none MF8/LMF=6
c
      if (idos6.gt.0.and.nmf68.eq.0) then
        call findmf(nin,mat,6,icod)
        if (icod.eq.0) then
          backspace(nin)
          mt=1000
          ii=0
          do while (mt.gt.0)
            call findnextmt(nin,6,mt)
            if (mt.gt.2) then
              if (nmf3.gt.0) then
                i3=iposm(nmf3,mf3,mt)
              else
                i3=0
              endif
              if (i3.gt.0) then
                call readcont(nin,c1,awr,jp6,lct,nk,n2,mat1,mf1,mt1,ns1)
                kzap=-99999
                do k=1,nk
                  call readtab1(nin,zap,awp,lfs,law,nr,np,nbt,ibt,x,y)
                  izap=nint(zap+1.0d-6)
                  if (izap.gt.2004) then
                    mt8=mtdos(mt,izap,lfs)
                    if (nmf10.gt.0) then
                      i10=iposm(nmf10,mf10,mt8)
                    else
                      i10=0
                    endif
                    if (nmf9.gt.0) then
                      i9=iposm(nmf9,mf9,mt8)
                    else
                      i9=0
                    endif
                    if (i9.eq.0.and.i10.eq.0) then
                      nmf68=nmf68+1
                      mt6(nmf68)=mt
                      lip6(nmf68)=lfs
                      lzap=izap
                      lzap6(nmf68)=lzap
                      if (kzap.ne.lzap) then
                        kk=0
                        kzap=lzap
                        ii=ii+1
                      else
                       kk=kk+1
                      endif
                      kip6(nmf68)=kk
                      nzap6(ii)=kk+1
                    endif
                  endif
                  call nextsub6(nin,law,nbt,ibt,y,yy)
                enddo
              else
                write(lou,*)' MF6 yield for mt=',mt,' was specified',
     &            ' but no MF3 data availble. Data ignored.'
              endif
            endif
          enddo
          if (nmf68.gt.0) then
            if (nmf68.gt.nnx6) then
              write(lou,'(2(a,i6))')' Number of reactions greater than',
     &         nnx6,' in MF6. Increase the value of nnx6 in ',nmf68-nnx6
              write(*,'(2(a,i6))')' Number of reactions greater than',
     &         nnx6,' in MF6. Increase the value of nnx6 in ',nmf68-nnx6
            else
              write(lou,'(1x,i6,a)')nmf68,
     &         ' set of yields found on MF6 (no data available on MF8)'
              write(*,'(1x,i6,a)')nmf68,
     &         ' set of yields found on MF6 (no data available on MF8)'
            endif
          endif
        endif
      endif
c
c     searching and saving MF6 yields, if required
c
      nmf6=0
      if (nmf68.gt.0) then
        kzap=-99999
        mtk=-999
        ii=0
        do i=1,nmf68
          mti=mt6(i)
          lzapi=lzap6(i)
          if (mti.ne.mtk.or.lzapi.ne.kzap) then
            mtk=mti
            kzap=lzapi
            call findmt(nin,mat,6,mti,icod)
            if (icod.eq.0) then
              call readcont(nin,c1,c2,jp,lct,nk,n3,mat0,mf0,mt0,ns0)
              kk=0
              do k=1,nk
                call readtab1(nin,zap,c2,lip,law,nr,ne,nbt,ibt,x,y)
                izap=nint(zap+1.0d-6)
                if (izap.eq.kzap) kk=kk+1
                call nextsub6(nin,law,nbt,ibt,x,y)
              enddo
              ii=ii+1
              nzap(ii)=kk
            else
              write(lou,*)' MF6 yield for mt=',mti,' zap=',lzapi,
     &          ' lip=',lip6(i),' was specified on MF8, but was not',
     &          ' found on MF6. Data ignored.'
            endif
          endif
        enddo
c
c       save yields in temporary file
c
        kzap=-99999
        mtk=-999
        ii=0
        do i=1,nmf68
          mti=mt6(i)
          lzap6i=lzap6(i)
          lip6i=lip6(i)
          if (mti.ne.mtk.or.lzap6i.ne.kzap) then
            ii=ii+1
            if (nzap(ii).eq.nzap6(ii)) then
              ktest=1
            else
              ktest=0
            endif
            write(lou,*)' nzap=',nzap(ii),' nzap6=',nzap6(ii),' ktest=',
     &       ktest
            kk=0
            mtk=mti
            kzap=lzap6i
          endif
          write(lou,*)' mt=',mti,' zap=',lzap6i,' lfs=',lip6i,' kfs=',
     &     kip6(i)
          call findmt(nin,mat,6,mti,icod)
          if (icod.eq.0) then
            call readcont(nin,c1,c2,jp,lct,nk,n3,mat0,mf0,mt0,ns0)
            ifound=0
            do k=1,nk
              call readtab1(nin,zap,c2,lip,law,nr,ne,nbt,ibt,x,y)
              izap=nint(zap+1.0d-6)
              if (izap.eq.lzap6i.and.
     &          (lip.eq.lip6i.or.(ktest.eq.1.and.kk.eq.kip6(i)))) then
                mti6=mtdos(mti,izap,lip6i)
                nmf6=nmf6+1
                mf6(nmf6)=mti6
                nst=0
                c1=dble(mti)
                c2=dble(mti6)
                call wrtab1(ntp,mat,6,mti,nst,c1,c2,izap,lip6i,
     &            nr,nbt,ibt,ne,x,y)
                ifound=1
                kk=kk+1
                exit
              endif
              call nextsub6(nin,law,nbt,ibt,x,y)
            enddo
            if (ifound.eq.0) then
              write(lou,*)' MF6 yield for mt=',mti,' zap=',lzap6i,
     &          ' lip=',lip6i,' was specified on MF8, but was not',
     &          ' found on MF6. Data ignored.'
            endif
          endif
        enddo
      endif
c
c     Setting triggers and main pointers for ACE-formatted file
c
      nmtr=nmf3+nmf10+nmf9+nmf6
      write(lou,*)
      write(lou,'(1x,i6,a)')nmtr,' dosimetry reactions found'
      write(lou,*)
      write(*,*)
      write(*,'(1x,i6,a)')nmtr,' dosimetry reactions found'
      write(*,*)
      nxs(2)=matza
      nxs(4)=nmtr
      nxs(9)=liso
      nxs(10)=matza/1000
      nxs(11)=mod(matza,1000)
      jxs(1)=1
      lmt=1
      jxs(3)=lmt
      lsig=lmt+nmtr
      jxs(6)=lsig
      lxsd=lsig+nmtr
      jxs(7)=lxsd
c
c     open files DODOS.CUR and DODOS.PLT for PLOTTAB, if imon=2
c
      if (imon.eq.2) then
        open(ncur,file='DODOS.CUR')
        open(nplt,file='DODOS.PLT')
      endif
c
c      load dosimetry data from MF3
c      XSD=XS(MF3)
c
      imt=0
      iloc=1
      if (nmf3.gt.0) then
        call findmf(nin,mat,3,icod)
        backspace(nin)
        i=1
        do while (i.le.nmf3)
          call findnextmt(nin,3,mti)
          if (mti.eq.mf3(i)) then
            call readcont(nin,za,awr,l1,l2,n1,n2,mat0,mf0,mt0,ns0)
            call readtab1(nin,qm,qi,l1,lr,nr,ne,nbt,ibt,x,y)
            imt=imt+1
            xss(lmt+imt-1)=mti
            xss(lsig+imt-1)=iloc
            call checklaw(nr,ibt,icod)
            if (icod.eq.0) then
              call remdup(lou,ne,x,y,irem)
              if (irem.gt.0) then
                write(lou,*)' Warning: ',irem,' points',
     &            ' removed from dosimetry reaction ',mti
              endif
              nr=0
            endif
            ll=lxsd+iloc-1
            xss(ll)=nr
            if (nr.gt.0) then
              do j=1,nr
                xss(ll+j)=nbt(j)
                xss(ll+nr+j)=ibt(j)
              enddo
            endif
            ll=ll+2*nr+1
            xss(ll)=ne
            do j=1,ne
              xss(ll+j)=x(j)*ev2mev
              xss(ll+ne+j)=y(j)
            enddo
            iloc=iloc+2*(1+nr+ne)
            izap=0
            lfs=0
            if (imon.gt.0) then
              call prtxsd(lou,imt,mti,3,nr,nbt,ibt,
     &          ne,xss(ll+1),xss(ll+ne+1),mti,izap,lfs)
            else
              write(lou,'(i4,a,i9,a)')imt,'. Dosimetry reaction mtd=',
     &          mti,' from MF3'
            endif
            if (imon.eq.2) then
              call plotxsd(ncur,nplt,hk,nmtr,imt,mti,3,
     &          ne,xss(ll+1),xss(ll+ne+1),mti,izap,lfs)
            endif
            write(*,'(i4,a,i9,a)')imt,'. Dosimetry reaction mtd=',
     &        mti,' from MF3'
            i=i+1
          endif
        enddo
      endif
c
c      load dosimetry data from MF10
c      XSD=XS(MF10)
c
      if (nmf10.gt.0) then
        call findmf(nin,mat,10,icod)
        backspace(nin)
        i=1
        do while (i.le.nmf10)
          call findnextmt(nin,10,mt)
          call readcont(nin,c1,c2,l1,l2,nfs,n2,mat0,mf0,mt0,ns0)
          do k=1,nfs
            call readtab1(nin,c1,c2,izap,lfs,nr,ne,nbt,ibt,x,y)
            mti=mtdos(mt,izap,lfs)
            if (mti.eq.mf10(i)) then
              imt=imt+1
              xss(lmt+imt-1)=mti
              xss(lsig+imt-1)=iloc
              call checklaw(nr,ibt,icod)
              if (icod.eq.0) then
                call remdup(lou,ne,x,y,irem)
                if (irem.gt.0) then
                write(lou,*)' Warning: ',irem,' points',
     &            ' removed from dosimetry reaction ',mti
                endif
                nr=0
              endif
              ll=lxsd+iloc-1
              xss(ll)=nr
              if (nr.gt.0) then
                do j=1,nr
                  xss(ll+j)=nbt(j)
                  xss(ll+nr+j)=ibt(j)
                enddo
              endif
              ll=ll+2*nr+1
              xss(ll)=ne
              do j=1,ne
                xss(ll+j)=x(j)*ev2mev
                xss(ll+ne+j)=y(j)
              enddo
              iloc=iloc+2*(1+nr+ne)
              if (imon.gt.0) then
                call prtxsd(lou,imt,mti,10,nr,nbt,ibt,
     &            ne,xss(ll+1),xss(ll+ne+1),mt,izap,lfs)
              else
                write(lou,'(i4,a,i9,a)')imt,'. Dosimetry reaction mtd=',
     &            mti,' from MF10'
              endif
              if (imon.eq.2) then
                call plotxsd(ncur,nplt,hk,nmtr,imt,mti,10,
     &            ne,xss(ll+1),xss(ll+ne+1),mt,izap,lfs)
              endif
              write(*,'(i4,a,i9,a)')imt,'. Dosimetry reaction mtd=',
     &          mti,' from MF10'
              i=i+1
            endif
          enddo
        enddo
      endif
c
c      load dosimetry data from MF9
c      XSD=YLD(MF9)*XS(MF3)
c
      rewind(ntp)
      if (nmf9.gt.0) then
        do i=1,nmf9
          call readtab1(ntp,c1,c2,izap,lfs,nry,ney,nbty,ibty,xy,yy)
          mt3=nint(c1+1.0d-6)
          mti=nint(c2+1.0d-6)
          if (mti.eq.mf9(i)) then
            call findmt(nin,mat,3,mt3,icod)
            call readcont(nin,c1,c2,l1,l2,n1,n2,mat0,mf0,mt0,ns0)
            call readtab1(nin,c1,c2,l1,l2,nr,ne,nbt,ibt,x,y)
            call setxsd(nr,nbt,ibt,ne,x,y,nry,nbty,ibty,ney,xy,yy,
     &        tol,ymin,npmax)
            imt=imt+1
            xss(lmt+imt-1)=mti
            xss(lsig+imt-1)=iloc
            nr=0
            ll=lxsd+iloc-1
            xss(ll)=nr
            ll=ll+1
            xss(ll)=ne
            do j=1,ne
              xss(ll+j)=x(j)*ev2mev
              xss(ll+ne+j)=y(j)
            enddo
            iloc=iloc+2*(1+ne)
            if (imon.gt.0) then
              call prtxsd(lou,imt,mti,9,nr,nbt,ibt,
     &          ne,xss(ll+1),xss(ll+ne+1),mt3,izap,lfs)
            else
              write(lou,'(i4,a,i9,a)')imt,'. Dosimetry reaction mtd=',
     &          mti,' from MF9*MF3'
            endif
            if (imon.eq.2) then
              call plotxsd(ncur,nplt,hk,nmtr,imt,mti,9,
     &          ne,xss(ll+1),xss(ll+ne+1),mt3,izap,lfs)
            endif
            write(*,'(i4,a,i9,a)')imt,'. Dosimetry reaction mtd=',
     &        mti,' from MF9*MF3'
          else
            write(lou,*)' === Fatal error processing MF9 data for mtd=',
     &        mf9(i),' mt=',mt3
            write(*,*)' === Fatal error processing MF9 data for mtd=',
     &        mf9(i),' mt=',mt3
            stop
          endif
        enddo
      endif
c
c      load dosimetry data from MF6
c      XSD=YLD(MF6)*XS(MF3)
c
      if (nmf6.gt.0) then
        do i=1,nmf6
          call readtab1(ntp,c1,c2,izap,lfs,nry,ney,nbty,ibty,xy,yy)
          mt3=nint(c1+1.0d-6)
          mti=nint(c2+1.0d-6)
          if (mti.eq.mf6(i)) then
            call findmt(nin,mat,3,mt3,icod)
            call readcont(nin,c1,c2,l1,l2,n1,n2,mat0,mf0,mt0,ns0)
            call readtab1(nin,c1,c2,l1,l2,nr,ne,nbt,ibt,x,y)
            call setxsd(nr,nbt,ibt,ne,x,y,nry,nbty,ibty,ney,xy,yy,
     &        tol,ymin,npmax)
            imt=imt+1
            xss(lmt+imt-1)=mti
            xss(lsig+imt-1)=iloc
            nr=0
            ll=lxsd+iloc-1
            xss(ll)=nr
            ll=ll+1
            xss(ll)=ne
            do j=1,ne
              xss(ll+j)=x(j)*ev2mev
              xss(ll+ne+j)=y(j)
            enddo
            iloc=iloc+2*(1+ne)
            if (imon.gt.0) then
              call prtxsd(lou,imt,mti,6,nr,nbt,ibt,
     &          ne,xss(ll+1),xss(ll+ne+1),mt3,izap,lfs)
            else
              write(lou,'(i4,a,i9,a)')imt,'. Dosimetry reaction mtd=',
     &          mti,' from MF6*MF3'
            endif
            if (imon.eq.2) then
              call plotxsd(ncur,nplt,hk,nmtr,imt,mti,6,
     &          ne,xss(ll+1),xss(ll+ne+1),mt3,izap,lfs)
            endif
            write(*,'(i4,a,i9,a)')imt,'. Dosimetry reaction mtd=',
     &        mti,' from MF6*MF3'
          else
            write(lou,*)' === Fatal error processing MF6 data for mtd=',
     &        mf6(i),' mt=',mt3
            write(*,*)' === Fatal error processing MF6 data for mtd=',
     &        mf6(i),' mt=',mt3
            stop
          endif
        enddo
      endif
      nxs(1)=iloc+2*nmtr-1
      jxs(22)=nxs(1)
      close(ntp,status='DELETE')
      close(nin)
c
c     write out ACE and xsdir files
c
      call dosout(fout,mcnpx)
c
c     print out summary
c
      if (imon.eq.0) write(lou,*)
      write(lou,'(a)')' ACE-formatted file heading lines:'
      if (mcnpx.eq.1) then
        write(lou,'(a13,f12.6,1x,1pe11.4,1x,a10)')hz(1:13),awr0,tz,hd
      else
        write(lou,'(a10,f12.6,1x,1pe11.4,1x,a10)')hz(1:10),awr0,tz,hd
      endif
      write(lou,'(a70,a10)')hk,hm
      write(lou,*)
      write(lou,'(a,i11)')' len2 = Lenght of XSS array:       ',nxs(1)
      write(lou,'(a,i11)')' za   = ZA number:                 ',nxs(2)
      write(lou,'(a,i11)')' ntr  = no. of dosimetry reactions:',nxs(4)
      write(lou,'(a,i11)')' lone = Location of first word:    ',jxs(1)
      write(lou,'(a,i11)')' mtr  = Location of MTR  block:    ',jxs(3)
      write(lou,'(a,i11)')' lsig = Location of LSIG block:    ',jxs(6)
      write(lou,'(a,i11)')' sigd = Location of SIGD block:    ',jxs(7)
      write(lou,*)
      write(lou,'(1x,a,2x,a)')'    MTD    ','    LOC    '
      write(lou,'(1x,a,2x,a)')'===========','==========='
      do i=1,nmtr
        mtd=nint(xss(lmt+i-1)+1.0d-10)
        loc=nint(xss(lsig+i-1)+1.0d-10)
        write(lou,'(i11,2x,i11)')mtd,loc
      enddo
      write(lou,*)
      write(lou,'(a,a)')' Output ACE-formatted file: ',fout
      write(*,*)
      write(*,'(a,a)')' Output ACE-formatted file: ',fout
c
c     The end
c
      call getdtime(cdate,ctime)
      write(lou,*)
      write(lou,'(a,a11,a,a10)')' DODOS ended at ',ctime,' on ',cdate
      write(*,*)
      write(*,'(a,a11,a,a10)')' DODOS ended at ',ctime,' on ',cdate
      if (imon.eq.2) then
        close(ncur)
        close(nplt)
      endif
      close(lou)
      stop
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
      subroutine remdup(nerr,n,x,y,irem)
c
c     remove duplicate information
c
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      character*11 x0,x1
      irem=0
      call ff2chx(x(1),x0)
      read(x0,'(e11.0)')x(1)
      y0=y(1)
      j=1
      i=2
      n0=n
      do while (i.le.n)
        call ff2chx(x(i),x1)
        y1=y(i)
        if (x0.eq.x1.and.((y0.eq.y1).or.(i.eq.n.and.y1.eq.0.0d0))) then
          if (nerr.gt.0) then
            write(nerr,*)' Warning: duplicate point at x=',x1,' y=',y1,
     &        ' i=',i,' removed'
          endif
        elseif (x0.ne.x1.and.y0.eq.y1) then
          i0=i
          do while(x0.ne.x1.and.y0.eq.y1.and.i.lt.n)
            x0=x1
            y0=y1
            i=i+1
            call ff2chx(x(i),x1)
            y1=y(i)
          enddo
          jrem=i-i0
          if (i.lt.n.and.((x0.eq.x1).or.(y0.ne.y1))) then
            j=j+1
            read(x0,'(e11.0)')x(j)
            y(j)=y0
            jrem=jrem-1
          endif
          j=j+1
          read(x1,'(e11.0)')x(j)
          y(j)=y1
          if (jrem.gt.0.and.nerr.gt.0) then
            if (i.lt.n.and.((x0.eq.x1).or.(y0.ne.y1))) then
              write(nerr,*)' Warning: ',jrem,' points removed between',
     &          x(j-2),' and ',x(j-1)
            else
              write(nerr,*)' Warning: ',jrem,' points removed between',
     &          x(j-1),' and ',x(j)
            endif
          endif
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
      irem=n0-n
      return
      end
C======================================================================
      subroutine setxsd(nr,nbt,ibt,ne,x,y,nry,nbty,ibty,ney,xy,yy,
     &  tol,ymin,nnmax)
c
c      calculate a linearly interpolable dosimetry reaction as the
c      product of the yield by the cross section
c
c      nr,nbt,ibt,ne,x,y : TAB1 data for representing the cross section
c      nry,nbty,ibty,ney,xy,yy: TAB1 data for representing the yield
c      tol: linearization tolerance
c      ymin: minimum cros section allowed for linearization
c      nnmax: array dimension (maximum number of energy points)
c
      implicit real*8 (a-h,o-z)
      parameter (ns=20, tole=1.0d-6)
      dimension nbt(*),ibt(*),nbty(*),ibty(*)
      dimension x(*),y(*),xy(*),yy(*)
      dimension xu(nnmax),xl(nnmax),yl(nnmax)
      dimension xs(ns),ys(ns)
c
c     Prepare union grid
c
      call checkdis(0,ne,x,y,icod)
      call checkdis(0,ney,xy,yy,icod)
      call union(ne,x,ney,xy,neu,xu,nnmax)
c
c      linearize the product yy(E)*y(E)=yld(E)*xsd(E)
c
      x0=xu(1)
      yld=fvalue(nry,nbty,ibty,ney,xy,yy,x0)
      xsd=fvalue(nr,nbt,ibt,ne,x,y,x0)
      y0=yld*xsd
      j=1
      xl(1)=x0
      yl(1)=y0
      do i=2,neu
        x1=xu(i)
        yld=fvalue(nry,nbty,ibty,ney,xy,yy,x1)
        xsd=fvalue(nr,nbt,ibt,ne,x,y,x1)
        y1=yld*xsd
        k=0
        nostop=1
        do while (nostop.eq.1)
          xm=0.5d0*(x0+x1)
          ym=0.5d0*(y0+y1)
          yld=fvalue(nry,nbty,ibty,ney,xy,yy,xm)
          xsd=fvalue(nr,nbt,ibt,ne,x,y,xm)
          yf=yld*xsd
          if (abs(yf-ym).le.abs(tol*yf).or.(xm-x0).le.tole*xm.or.
     &        yf.lt.ymin.or.k.eq.ns) then
            call inc(j,nnmax)
            xl(j)=x1
            yl(j)=y1
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
            y1=yf
          endif
        enddo
        x0=x1
        y0=y1
      enddo
      ne=j
      call remdup(0,ne,xl,yl,irem)
      if (irem.gt.0) then
        write(*,*)' Warning: ',irem,' points',
     &    ' removed from dosimetry reaction'
      endif
      do i=1,ne
        x(i)=xl(i)
        y(i)=yl(i)
      enddo
      nr=1
      nbt(1)=ne
      ibt(1)=2
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
c     find next subsection on MF6 section according to the value of law
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
      elseif (law.gt.7) then
        write(*,*)' ERROR: unknown LAW=',law,' on MF6'
        stop
      endif
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
C     Printing and plotting routines for dosimetry cross sections
C======================================================================
      subroutine prtxsd(lou,imtd,mtd,mfd,nr,nbt,ibt,np,x,y,mt,izap,l)
c
c     print the dosimetry cross section mtd from file mfd
c
      implicit real*8 (a-h, o-z)
      dimension nbt(*),ibt(*),x(*),y(*)
      character*15 ii,ni,li,ui
      character*18 ei,xsd,uu
      data ii/'    interval   '/,ni/'     nbt(i)    '/
      data li/'     int(i)    '/
      data ei/'      ENERGY      '/,xsd/' DOSIMETRY X-SEC. '/
      data ui/'==============='/,uu/'=================='/
      if (mfd.eq.10) then
        write(lou,'(i4,a,i9,a,i4,a,i8,a,i4,a)')imtd,
     &    '. Dosimetry reaction mtd=',mtd,' from MF10. ( MT=',
     &    mt,', ZAP=',izap,', LFS=',l,' )'
      elseif (mfd.eq.9) then
        write(lou,'(i4,a,i9,a,i4,a,i8,a,i4,a)')imtd,
     &    '. Dosimetry reaction mtd=',mtd,' from MF9*MF3. ( MT=',
     &    mt,', ZAP=',izap,', LFS=',l,' )'
      elseif (mfd.eq.6) then
        write(lou,'(i4,a,i9,a,i4,a,i8,a,i4,a)')imtd,
     &    '. Dosimetry reaction mtd=',mtd,' from MF6*MF3. ( MT=',
     &    mt,', ZAP=',izap,', LFS=',l,' )'
      else
        write(lou,'(i4,a,i9,a)')imtd,
     &    '. Dosimetry reaction mtd=',mtd,' from MF3.'
      endif
      write(lou,*)
      write(lou,'(1x,a,i9)')' number of energy points: ',np
      if (nr.gt.0) then
        write(lou,'(1x,a,3(1x,a15))')' interpolation law:',ii,ni,li
        write(lou,'(1x,a,3(1x,a15))')'===================',ui,ui,ui
        do i=1,nr
          write(lou,'(21x,i15,1x,i15,1x,i15)')i,nbt(i),ibt(i)
        enddo
      else
        write(lou,'(1x,a)')' interpolation law: LIN-LIN'
        write(lou,'(1x,a)')'============================'
      endif
      write(lou,*)
      write(lou,'(3(1x,a18,1x,a18))')(ei,xsd,i=1,3)
      write(lou,'(3(1x,a18,1x,a18))')(uu,uu,i=1,3)
      write(lou,'(3(1x,1pe18.11,1x,e18.11))')(x(i),y(i),i=1,np)
      write(lou,*)
      return
      end
C======================================================================
      subroutine plotxsd(ncur,nplt,hk,nmtr,imtd,mtd,mfd,np,x,y,
     &                   mt,izap,l)
c
c     Add dosimetry cross section mtd from mfd to PLOTTAB files
c
      implicit real*8 (a-h, o-z)
      parameter (maxcur=10)
      dimension x(*),y(*)
      character*11 chx,chy
      character*70 hk

      if (mfd.eq.10) then
        write(ncur,'(a9,i9,11x,i4,i8,i4)')'MF10/MTD=',mtd,mt,izap,l
      elseif (mfd.eq.9) then
        write(ncur,'(a9,i9,11x,i4,i8,i4)')'MF09/MTD=',mtd,mt,izap,l
      elseif (mfd.eq.6) then
        write(ncur,'(a9,i9,11x,i4,i8,i4)')'MF06/MTD=',mtd,mt,izap,l
      else
        write(ncur,'(a9,i9,11x,i4)')'MF03/MTD=',mtd,mt
      endif
      do i=1,np
        call chendf(x(i),chx)
        call chendf(y(i),chy)
        write(ncur,'(2a11)')chx,chy
      enddo
      write(ncur,*)
      jmt=imtd/maxcur
      mmt=imtd-jmt*maxcur
      if (mmt.eq.0.or.imtd.eq.nmtr) then
        if (imtd.gt.maxcur) write(nplt,*)
        write(nplt,'(1p,4e11.4,2i11,0p,f4.2)')0.0,13.5,0.0,10.0,1,1,1.5
        if (mmt.eq.0) then
          write(nplt,'(6i11,i4)')maxcur,0,1,0,0,0,0
        else
          write(nplt,'(6i11,i4)')mmt,0,1,0,0,0,0
        endif
        write(nplt,'(a15,25x,a3)')'Incident Energy','MeV'
        write(nplt,'(a13,27x,a4)')'Cross section','barn'
        write(nplt,'(a70)')hk
        write(nplt,'(a)')'Dosimetry cross sections'
        write(nplt,'(22x,4i11)')0,0,0,0
        write(nplt,'(22x,4i11)')0,0,0,0
      endif
      return
      end
C======================================================================
      function mtdos(mt,izap,lfs)
      if (mt.eq.5) then
        mtdos=1000000*(50+lfs)+izap
      elseif (mt.eq.18) then
        if (izap.eq.-1) then
          mtdos=18
        else
          mtdos=1000000*(80+lfs)+izap
        endif
      else
        mtdos=1000*(10+lfs)+mt
      endif
      return
      end
C======================================================================
      function mtdchk(nmf68,mt6,lzap6,lip6,mt8)
      dimension mt6(*),lzap6(*),lip6(*)
      mtdchk=0
      if (nmf68.gt.0) then
        do i=1,nmf68
          mt=mt6(i)
          izap=lzap6(i)
          lfs=lip6(i)
          mt0=mtdos(mt,izap,lfs)
          if (mt0.eq.mt8) then
            mtdchk=1
            exit
          endif
        enddo
      endif
      return
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
C      The following subroutines were taken from NJOY and
C      adapted/modified by D. Lopez Aldama for ACEMAKER:
C       1. subroutine terp1 (renamed as terp1m)
c       2. subroutine dosout
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
      subroutine dosout(fout,mcnpx)
c
c      Write out a type 1 dosimetry ACE-formatted file.
c
c      Adapted by D. Lopez Aldama for ACEMAKER
c
      implicit real*8 (a-h, o-z)
      character*10 hd,hm
      character*13 hz
      character*70 hk
      character*72 fout,fdir
      common/acetxt/hz,hd,hm,hk
      common/acecte/awr0,tz,awn(16),izn(16)
      common/acepnt/nxs(16),jxs(32)
      common/acedat/xss(50000000),nxss
      data nout/30/,ndir/31/
c
c       open output ACE-formatted file
c
      open(nout,file=fout)
c
c      Type 1
c
      if (mcnpx.eq.1) then
        write(nout,'(a13,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')
     &    hz(1:13),awr0,tz,hd,hk,hm
      else
        write(nout,'(a10,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')
     &    hz(1:10),awr0,tz,hd,hk,hm
      endif
      write(nout,'(4(i7,f11.0))') (izn(i),awn(i),i=1,16)
      len2=nxs(1)
      ntr=nxs(4)
      mtr=jxs(3)
      lsig=jxs(6)
      lxsd=jxs(7)
      write(nout,'(8i9)')(nxs(i),i=1,16)
      write(nout,'(8i9)')(jxs(i),i=1,32)
c
c     mtr block
c
      l=mtr
      do i=1,ntr
        call typen(l,nout,1)
        l=l+1
      enddo
c
c     lsig block
c
      l=lsig
      do i=1,ntr
        call typen(l,nout,1)
        l=l+1
      enddo
c
c     sigd block
c
      l=lxsd
      do i=1,ntr
        nr=nint(xss(l))
        call typen(l,nout,1)
        l=l+1
        if (nr.ne.0) then
          n=2*nr
          do j=1,n
            call typen(l,nout,1)
            l=l+1
          enddo
        endif
        ne=nint(xss(l))
        call typen(l,nout,1)
        l=l+1
        n=2*ne
        do j=1,n
          call typen(l,nout,2)
          l=l+1
        enddo
      enddo
      call typen(0,nout,3)
      close(nout)
      nern=0
      lrec=0
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
     &    hz(1:13),awr0,trim(fout),'0 1 1 ',len2,lrec,nern,tz
      else
        write(ndir,'(a10,f12.6,1x,a,1x,a,i9,2i3,1pe11.4)')
     &    hz(1:10),awr0,trim(fout),'0 1 1 ',len2,lrec,nern,tz
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
      common/acedat/xss(50000000),nxss
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