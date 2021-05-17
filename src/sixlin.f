      program sixlin
c     version 2.0
c
c     Process and linearize MF6 data
c
c     SIXLIN converts continuum energy-angle distribution given
c     by Legendre coefficients (LAW=1/LANG=1) and laboratory angle-
c     energy distribution (LAW=7) to continuum energy-angle distribution
c     given by linearly interpolable tabulated data (LAW=1/LANG=12).
c     Additionally, SIXLIN linearizes continuum energy-angle
c     distribution given by Kalbach-Mann systematics or by non-linearly
c     interpolable tabulated data on MF6.
c
c     The code belongs to ACEMAKER code package to produce continuous
c     energy ACE-formatted file for Monte Carlo calculations
c
c     INPUT data:
c
c     Input data options should be entered on the SIXLIN.INP text file.
c
c     line 1:       isel       imon       ymin     (2i11,e11.0)
c     line 2:       input  ENDF filename           (a72)
c     line 3:       output ENDF filename           (a72)
c     line 4:       ilow       iupp                (2i11)
c     line 5:        tol       dxmu                (2e11.0)
c
c     where
c       isel: selection criterium      (0/1) = ZA/MAT   (Default = 0)
c       imon: Monitor printing trigger (0/1) = min/max  (Default = 0)
c       ymin: Minimum value alowable for linearization/reconstruction
c             of angle-energy distributions (1.0d-30)
c       ilow: Lower material number requested
c       iupp: Upper material number requested
c            (Default ilow=iupp=0,first material on the input ENDF tape)
c        tol: Fraction tolerance for linearization/reconstruction
c             (Default = 0.01)
c       dxmu: Maximum cosine interval for linearization/reconstruction
c             (Default = 0.001)
c
c     Example of SIXLIN.INP:
c
c               0          0    1.0D-20
c     \PENDF\U235.PENDF
c     \LIN\U235.SIX
c            9228       9228
c         1.0D-3      5.0D-3
c
c     Retrieve material 9228 (U-235) from tape \PENDF\U235.PENDF
c     and process/linearize MF6 file. Write results on \LIN\U235.SIX
c     tape. The minimum allowable value for angle-energy reconstruction
c     will be 1.e-20. Use a linearizaation tolerance of 0.1% (0.001) and
c     a maximum cosine interval of 0.005.
c
c     Output files:
c      Output ENDF formatted file
c      SIXLIN.LST listing file
c
      implicit real*8(a-h, o-z)
      parameter (npmax=500000, nbmax=2000000, npmumax=5001)
      parameter (etol0=0.01d0,dmumax0=0.001d0,ymin0=1.0d-30)
      character*11 snum
      character*66 line
      character*72 fin1,fout
      dimension nepmu(npmumax)
      dimension nbt(20),ibt(20),nbtm(20),ibtm(20),nbte(20),ibte(20)
      dimension x(npmax),y(npmax)
      dimension x0(npmax),y0(npmax),x1(npmax),y1(npmax)
      dimension b(nbmax),b0(nbmax)
      data in1/2/,nin/3/,jou/10/,iou/11/,nou/12/
c
c      open input file SIXLIN.INP and list file SIXLIN.LST
c
      open (in1, file='SIXLIN.INP')
      open (iou, file='SIXLIN.LST')
      write(iou,*)' PROGRAM SIXLIN: Process and linearize File 6 data'
      write(iou,*)' ================================================='
      write(iou,*)
c
c      read input data from SIXLIN.INP
c
      isel=0
      imon=0
      llower=0
      lupper=0
      read(in1,'(2i11,e11.0)')isel,imon,ymin
      if (isel.ne.1) isel=0
      if (imon.ne.1) imon=0
      if (ymin.le.0.0d0) ymin=ymin0
      read(in1,'(a72)')fin1
      read(in1,'(a72)')fout
      read(in1,'(2i11)')llower,lupper
      if (llower.lt.0) llower=0
      if (lupper.lt.llower) lupper=llower
      read(in1,'(2e11.0)')tolmax,dmumax
      etol=tolmax
      if (etol.le.0.0d0) etol=etol0
      if (dmumax.le.0.0d0) dmumax=dmumax0
      close(in1)
c
c     printing input data and open input and output files
c
      open (nin, file=fin1)
      open (nou, file=fout)
      write(iou,*)' Input parameters'
      write(iou,*)' ================'
      write(iou,*)' MAT/ZA selection =',isel
      write(iou,*)' Printing option  =',imon
      write(iou,'(1x,a,1pe11.5)')' Minimum Y-value  =',ymin
      write(iou,*)' Input file name  =',fin1
      write(iou,*)' Output file name =',fout
      if (llower.eq.0.and.lupper.eq.0) then
        write(iou,*)' Selection range  = first material'
      else
        write(iou,*)' Selection range  =',llower, ' to ',lupper
      endif
      write(iou,'(1x,a,1pe11.5)')' Tolerance [%]    =',etol*100.0d0
      write(iou,'(1x,a,1pe11.5)')' Cosine interval  =',dmumax
      write(iou,*)
c
c     read tape header
c
      ns=-1
      call readtext(nin,line,mat,mf,mt,nsi)
      call wrtext(nou,mat,mf,mt,ns,line)
      call readcont(nin,za,awr,l1,l2,n1,n2,mat,mf,mt,nsi)
      if (llower.eq.0.and.lupper.eq.0) then
        llower=mat
        lupper=mat
      endif
      call checksel(isel,mat,za,llower,lupper,nsel)
      do while (nsel.ne.-1)
       if (nsel.eq.0) then
c
c        Skip this material and find the next one
c
         call findnextmat(nin)
       else
c
c        Material should be processed
c
c        processing  mf=1/mt=451
c
         call wrtcont(nou,mat,mf,mt,ns,za,awr,l1,l2,n1,n2)
         call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,nsi)
         call wrtcont(nou,mat,mf,mt,ns,c1,c2,l1,l2,n1,n2)
         call readcont(nin,awi,emax,lrel,l2,nsub,nver,mat,mf,mt,nsi)
         call wrtcont(nou,mat,mf,mt,ns,awi,emax,lrel,l2,nsub,nver)
         call readcont(nin,temp,c2,ldrv,l2,nwd,nxc,mat,mf,mt,nsi)
         nwd1=nwd+3
         call wrtcont(nou,mat,mf,mt,ns,temp,c2,ldrv,l2,nwd1,nxc)
         do i=1,nwd
           call readtext(nin,line,mat,mf,mt,nsi)
           call wrtext(nou,mat,mf,mt,ns,line)
         enddo
         line( 1:33)=' ***************** Program SIXLIN'
         line(34:66)=' (VERSION 2017-1) ***************'
         call wrtext(nou,mat,mf,mt,ns,line)
         write(line,'(a26,1pe11.4,a24)')' For All Data Greater than',
     &     ymin,' in Absolute Value'
         call wrtext(nou,mat,mf,mt,ns,line)
         write(line,'(a42,f10.7,a9)')
     &     ' Data Linearized to Within an Accuracy of ',
     &     etol*100.d0, ' per-cent'
         call wrtext(nou,mat,mf,mt,ns,line)
         do i=1,nxc
           call readtext(nin,line,mat,mf,mt,nsi)
           call wrtext(nou,mat,mf,mt,ns,line)
         enddo
         call readtext(nin,line,mat,mf,mt,nsi)
         call wrtsend(nou,mat,mf,ns)
c
c        Check for incident particle type (library type)
c
         if (nsub.eq.0) then
           zai=0.0d0
         elseif (nsub.eq.5.or.nsub.eq.10) then
           zai=1.0d0
         elseif (nsub.eq.10010) then
           zai=1001.0d0
         elseif (nsub.eq.10020) then
           zai=1002.0d0
         elseif (nsub.eq.10030) then
           zai=1003.0d0
         elseif (nsub.eq.20030) then
           zai=2003.0d0
         elseif (nsub.eq.20040) then
           zai=2004.0d0
         else
           write(iou,*)'  === Error: incident particle not allowed'
           write(iou,*)'  === NSUB=',nsub,' IPART=',(nsub/10)
           stop
         endif
         write(iou,*)' === MAT=',mat,' ZA=',int(za),' AWR=',awr,
     &     ' ZAI=',int(zai), '==='
c
c         Copy MF1 to MF5 to output file
c
         write(*,*)' MF1 to MF5 copied to ouput tape'
         write(iou,*)' MF1 to MF5 copied to ouput tape'
         call readtext(nin,line,mat,mf,mt,nsi)
         do while (mf.lt.6)
           if (mt.eq.0.and.mf.ne.0) then
             call wrtsend(nou,mat,mf,ns)
           elseif (mf.eq.0) then
             call wrtfend(nou,mat,ns)
           else
             call wrtext(nou,mat,mf,mt,ns,line)
           endif
           call readtext(nin,line,mat,mf,mt,nsi)
         enddo
         backspace(nin)
c
c          Process MF6 data, if any
c
         if (mf.eq.6) then
           write(*,*)' Processing MF6'
           call readcont(nin,za,awr,jp,lct,nk,n2,mat,mf,mt,nsi)
           do while(mf.eq.6)
            if ((mt.ne.18).or.(mt.eq.18.and.jp.eq.0)) then
             call wrtcont(nou,mat,mf,mt,ns,za,awr,l1,lct,nk,n2)
             write(iou,*)' Processing MF6/MT=',mt,' LCT=',lct,' NK=',nk
             write(iou,*)' ============================================'
             do kk=1,nk
               call readtab1(nin,zap,awp,lip,law,nr,np,nbt,ibt,x,y)
               write(iou,*)' Particle ',kk,' ZAP=',zap,' AWP=',awp
               write(iou,*)' LAW=',law
c
c              Check interpolation law for product yields
c              call linear if non linear
c
               call checklaw(nr,ibt,icod)
               if (icod.ne.0) then
                 call linear(nr,nbt,ibt,np,x,y,etol,ymin,npmax)
               endif
c
c              for law 1, 2 & 7 data is processed
c
               if (law.eq.1) then
c
c                law=1
c
                 call wrtab1(nou,mat,mf,mt,ns,
     &                       zap,awp,lip,law,nr,nbt,ibt,np,x,y)
                 call readtab2(nin,c1,c2,lang,lep,nr,ne,nbt,ibt)
                 if (lang.eq.2.or.lang.eq.11.or.lang.eq.12) then
                   lang2=lang
                   write(iou,*)' LANG=',lang
                 else
                   lang2=12
                   write(iou,*)' LANG=',lang,
     &             ' Angular parameters will be converted to LANG=12'
                 endif
                 if (lep.gt.2) then
                   lep2=2
                   write(iou,*)' LEP =',lep,
     &             ' Secundary energies will be linearized'
                 else
                   lep2=lep
                   write(iou,*)' LEP =',lep
                 endif
                 call wrtab2(nou,mat,mf,mt,ns,c1,c2,lang2,lep2,
     &            nr,ne,nbt,ibt)
                 nr=1
                 ibt(1)=lep
                 nrm=1
                 ibtm(1)=lang-10
                 do ie=1,ne
                   call readlist(nin,c1,e,nd,na,nw,nep,b)
                   write(iou,*)' Incident energy =',e
                   na2=na+2
                   nepc=nep-nd
                   i0=nd*na2
                   nbt(1)=nepc
                   if (lep.gt.2) then
c
c                    linearize list data and prepare linearly
c                    interpolable data with new nep,
c                    new nw=nep(na+2),
c                    nep-nd=nepc (continuous),
c                    where nd is the number of discrete lines.
c
c                    saving original list record
c
                     nepc0=nepc
                     do j=1,nw
                       b0(j)=b(j)
                     enddo
c
c                    saving original secundary energy grid
c
                     j=i0+1
                     do i=1,nepc0
                       ep=b0(j)
                       x0(i)=ep
                       x(i)=ep
                       j=j+na2
                     enddo
c
c                    linearize angular parameters on E' grid
c
                     do k=2,na2
                       j=i0+k
                       do i=1,nepc0
                         y(i)=b0(j)
                         j=j+na2
                       enddo
                       do i=1,nepc
                         b(i)=fvalue(nr,nbt,ibt,npec0,x0,y,x(i))
                       enddo
                       call linear(nr,nbt,ibt,nepc,x,b,etol,ymin,npmax)
                     enddo
c
c                    Prepare linearized list record
c                     1. Calculate new nep and nw
c                     2. Process discrete data
c                     3. Set new secondary energy grid for lin-lin
c                     4. Update lineaely interpolable angular parameters
c                     5. Check renormalization
c
                     nep=nepc+nd
                     nw=nep*na2
                     if (nw.gt.nbmax) then
                       write(iou,*)' ERROR: increase size of arrays,',
     &                 ' nbmax=',nbmax,'is too small ',nw,' required'
                       stop
                     endif
                     if (nd.gt.0) then
                       do i=1,i0
                         b(i)=b0(i)
                       enddo
                     endif
                     j=i0+1
                     do i=1,nepc
                       b(j)=x(i)
                       j=j+na2
                     enddo
                     do k=2,na2
                       j=i0+k
                       do i=1,nepc0
                         y(i)=b0(j)
                         j=j+na2
                       enddo
                       j=i0+k
                       do i=1,nepc
                         b(j)=fvalue(nr,nbt,ibt,npec0,x0,y,x(i))
                         j=j+na2
                       enddo
                     enddo
c
c                    Check normalization of f0 over E'
c
                     sum=0.0d0
                     if (nd.gt.0) then
                       j=2
                       do i=1,nd
                         sum=sum+b(j)
                         j=j+na2
                       enddo
                     endif
                     j=i0+1
                     do i=1,nepc
                       x(i)=b(j)
                       f0=b(j+1)
                       if (f0.ge.0.0d0) then
                         y(i)=f0
                       else
                         y(i)=0.0d0
                       endif
                       j=j+na2
                     enddo
                     c=0.0d0
                     call renorm(nepc,x,y,lep,c,fn)
                     call renorm(nepc,x,y,lep2,c,fn)
                     if (imon.gt.0) then
                      write(iou,*)' Integral of f0 over E'' at E = ',e
                      write(iou,*)' Discrete   contribution      = ',sum
                      write(iou,*)' Continuous contribution      = ',c
                      write(iou,*)' Renormalization  factor      = ',fn
                      write(iou,*)' ======================= '
                     endif
                     j=i0+2
                     do i=1,nepc
                       b(j)=y(i)
                       j=j+na2
                     enddo
                   endif
c
c                  1. for (na=0) or (lang=2) or (lep=1/2 and lang=11/12)
c                     data are ready to be written.
c                  2. for lang =13,14,15 angular data should be
c                     linearized and renormalized.
c                  3. for lang=11 or 12 with lep>2 data should be
c                     renormalized
c                  4. for lang=1 (Legendre) data should be converted to
c                     LAW=1 LANG=12
c
                   if ((na.eq.0).or.(lang.eq.2).or.
     &                ((lep.le.2).and.(lang.eq.11.or.lang.eq.12))) then
c
c                    write list record
c
                     call wrtlist(nou,mat,mf,mt,ns,c1,e,nd,na,nw,nep,b)
                   elseif (lang.ge.11.and.lang.le.15) then
                     nmu=na/2
                     if ((na-2*nmu).ne.0) then
                       write(iou,*)' ERROR: NA is not even'
                       stop
                     endif
                     if (lang.gt.12) then
c
c                      linearize angular distribution
c                        1. prepare a union grid
c                        2. linearize angular data for lang>12
c
c                        prepare union grid for cosines
c
                       nmu0=nmu
                       nbtm(1)=nmu0
                       do j=1,nw
                         b0(j)=b(j)
                       enddo
                       j=3
                       do i=1,nmu
                         x(i)=b0(j)
                         j=j+2
                       enddo
                       do k=2,nep
                         j=(k-1)*na2+3
                         do i=1,nmu0
                           x0(i)=b0(j)
                           j=j+2
                         enddo
                         call union(nmu0,x0,nmu,x,nmu1,y,npmax)
                         nmu=nmu1
                         do i=1,nmu
                           x(i)=y(i)
                         enddo
                       enddo
c
c                      linearize angular data for lang>12
c
                       do k=1,nepc
                         j=i0+(k-1)*na2+3
                         do i=1,nmu0
                           x0(i)=b0(j)
                           y(i)=b0(j+1)
                           j=j+2
                         enddo
                         do i=1,nmu
                           b(i)=fvalue(nrm,nbtm,ibtm,nmu0,x0,y,x(i))
                         enddo
                         call linear(nrm,nbtm,ibtm,nmu,x,b,
     &                    etol,ymin,npmax)
                       enddo
c
c                      reconstruction of the list record
c
                       na0=na
                       na20=na2
                       i00=i0
                       na=2*nmu
                       na2=na+2
                       nw=nep*na2
                       i0=nd*na2
                       if (nw.gt.nbmax) then
                         write(iou,*)' ERROR: increase arrays size, ',
     &                   ' nbmax=',nbmax,'is too small ',nw,' required'
                         stop
                       endif
c
c                      E',f0
c
                       j=1
                       i=1
                       do k=1,nep
                         b(j)=b0(i)
                         b(j+1)=b0(i+1)
                         i=i+na20
                         j=j+na2
                       enddo
c
c                      discrete distribution, if any
c
                       if (nd.gt.0) then
                         do k=1,nd
                           j=(k-1)*na20+3
                           do i=1,nmu0
                             x0(i)=b0(j)
                             y0(i)=b0(j+1)
                             j=j+2
                           enddo
                           j=(k-1)*na2+3
                           do i=1,nmu
                             xnu=x(i)
                             b(j)=xnu
                             b(j+1)=bvalue(xnu,nmu0,x0,y0)
                             j=j+2
                           enddo
                         enddo
                       endif
c
c                      continuous distribution
c
                       do k=1,nepc
                         j=i00+(k-1)*na20+3
                         do i=1,nmu0
                           x0(i)=b0(j)
                           y0(i)=b0(j+1)
                           j=j+2
                         enddo
                         j=i0+(k-1)*na2+3
                         do i=1,nmu
                           xnu=x(i)
                           b(j)=xnu
                           b(j+1)=fvalue(nrm,nbtm,ibtm,nmu0,x0,y0,xnu)
                           j=j+2
                         enddo
                       enddo
                     endif
c
c                    Normalize angular distribution
c
                     do k=1,nepc
                       j=i0+(k-1)*na2+2
                       c=b(j)
                       j=j+1
                       do i=1,nmu
                         x(i)=b(j)
                         f0=b(j+1)
                         if (f0.ge.0.0d0) then
                           y(i)=f0
                         else
                           y(i)=0.0d0
                         endif
                         j=j+2
                       enddo
                       call renorm(nmu,x,y,lang2,c,fn)
                       if (imon.gt.0) then
                       write(iou,*)' Renormalized f(u) at E'' number ',k
                       write(iou,*)' Integral value               = ',c
                       write(iou,*)' Renormalization factor       = ',fn
                       write(iou,*)' ======================'
                       endif
                       if (c.ne.0.0d0) then
                         j=i0+(k-1)*na2+4
                         do i=1,nmu
                           b(j)=0.5d0*y(i)/c
                           j=j+2
                         enddo
                       endif
                     enddo
                     call wrtlist(nou,mat,mf,mt,ns,c1,e,nd,na,nw,nep,b)
                   elseif (lang.eq.1) then
c
c                    Legendre coefficient representation should be
c                    converted to linearly interpolable tabulated
c                    data (law=1, lang=12)
c
c                    prepare a common grid (cosine)
c
                     do i=1,nw
                       b0(i)=b(i)
                     enddo
                     y(1)=-1.0d0
                     y(2)=1.0d0
                     nmu0=2
                     if (na.gt.1) then
                       n1=(npmax-1)/2 + 1
                       call setxmu(x1,nmumax,dmumax,n1)
                     else
                       nmumax=2
                     endif
                     imu0=nmumax
                     npmax1=nmumax-1
                     npmax0=2*npmax1+1
                     nl=na
                     nl1=nl+1
                     nl2=nl1+1
                     j=1
                     do k=1,nep
                       do i=1,nl1
                         b(i)=b0(j+i)
                       enddo
                       xmu=x1(1)
                       x0(imu0)=xmu
                       y0(imu0)=yleg(xmu,b,nl)
                       do i=1,npmax1
                         xmu=x1(i+1)
                         call yleg2(xmu,b,nl,yplus,yminus)
                         iminus=imu0-i
                         x0(iminus)=-xmu
                         y0(iminus)=yminus
                         iplus=imu0+i
                         x0(iplus)=xmu
                         y0(iplus)=yplus
                       enddo
                       call thinxs(x0,y0,npmax0,x,b,nmu,etol)
                       call union(nmu0,y,nmu,x,nmu1,x0,npmax)
                       nmu0=nmu1
                       do i=1,nmu0
                         y(i)=x0(i)
                       enddo
                       j=j+nl2
                     enddo
                     if(nmu0.gt.2) then
                       x(1)=-1.0d0
                       x(2)=0.0d0
                       x(3)=1.0d0
                       nmu=3
                       call union(nmu,x,nmu1,y,nmu0,x0,npmax)
                     endif
                     write(iou,*)' Common cosine grid for mt=',mt,
     &                           ' at E=',e,' nmu=',nmu0
c
c                    Convert data to lang=12 using the common grid
c
                     nmu=nmu0
                     j=1
                     i=1
                     do k=1,nep
                       b(i)=b0(j)
                       write(iou,*)' Secundary energy E'' =',b(i)
                       j=j+1
                       i=i+1
                       c=b0(j)
                       b(i)=c
                       y0(1)=c
                       do l=1,nl
                         j=j+1
                         y0(l+1)=b0(j)
                       enddo
                       do l=1,nmu
                         xnu=x0(l)
                         f0=yleg(xnu,y0,nl)
                         if (f0.ge.0.0d0) then
                           y(l)=f0
                         else
                           y(l)=0.0d0
                         endif
                       enddo
                       call renorm(nmu,x0,y,lang2,c,fn)
                       if (imon.gt.0) then
                       write(iou,*)' Renormalized f(u) at E'' number ',k
                       write(iou,*)' Integral value               = ',c
                       write(iou,*)' Renormalization factor       = ',fn
                       write(iou,*)' ======================'
                       endif
                       do l=1,nmu
                         i=i+1
                         b(i)=x0(l)
                         i=i+1
                         if(c.ne.0.0d0) then
                           b(i)=0.5d0*y(l)/c
                         else
                           b(i)=0.0d0
                         endif
                       enddo
                       j=j+1
                       i=i+1
                     enddo
                     na=2*nmu
                     na2=na+2
                     nw=nep*na2
                     call wrtlist(nou,mat,mf,mt,ns,c1,e,nd,na,nw,nep,b)
                   else
                     write(iou,*)' Error: LAW=1 LANG=',lang,' is wrong'
                     stop
                   endif
                 enddo
               elseif (law.eq.2) then
c
c                convert lang=0 and lang=14 to lang 12
c
                 call wrtab1(nou,mat,mf,mt,ns,
     &            zap,awp,lip,law,nr,nbt,ibt,np,x,y)
                 call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
                 call wrtab2(nou,mat,mf,mt,ns,c1,c2,l1,l2,nr,ne,nbt,ibt)
                 do ie=1,ne
                   call readlist(nin,c1,e,lang,l2,nw,nl,b)
                   if (lang.eq.0) then
c
c                    convert legendre coefficient (lang=0) to linearly
c                    interpolable tabulated data  (lang=12)
c
                     write(iou,*)' LANG=0 converted to LANG=12'
                     nl1=nl+1
                     b0(1)=0.5d0
                     j=2
                     do i=1,nl
                       b0(j)=b(i)
                       j=j+1
                     enddo
                     if (nl.gt.1) then
                       n1=(npmax-1)/2 + 1
                       dmumax=dmumax0
                       call setxmu(x1,nmumax,dmumax,n1)
                     else
                       nmumax=2
                     endif
                     imu0=nmumax
                     npmax1=nmumax-1
                     npmax0=2*npmax1+1
                     xmu=x1(1)
                     x0(imu0)=xmu
                     y0(imu0)=yleg(xmu,b0,nl)
                     do i=1,npmax1
                       xmu=x1(i+1)
                       call yleg2(xmu,b0,nl,yplus,yminus)
                       iminus=imu0-i
                       x0(iminus)=-xmu
                       y0(iminus)=yminus
                       iplus=imu0+i
                       x0(iplus)=xmu
                       y0(iplus)=yplus
                     enddo
                     call thinxs(x0,y0,npmax0,x,y,nmu,etol)
                     if(nmu.gt.2) then
                       y1(1)=-1.0d0
                       y1(2)=0.0d0
                       y1(3)=1.0d0
                       nmu1=3
                       call union(nmu1,y1,nmu,x,nmu0,x0,npmax)
                       nmu=nmu0
                     else
                       do i=1,nmu
                         x0(i)=x(i)
                       enddo
                     endif
                     write(iou,*)' Common cosine grid for mt=',
     &                            mt,' at E=',e,' nmu=',nmu
c
c                    prepare data to lang=12 using the common grid
c
                     do j=1,nmu
                       xnu=x0(j)
                       f0=yleg(xnu,b0,nl)
                       if (f0.ge.0.0d0) then
                         y0(j)=f0
                       else
                         y0(j)=0.0d0
                       endif
                     enddo
                     lang=12
                     c=1.0d0
                     call renorm(nmu,x0,y0,lang,c,fn)
                     if (imon.gt.0) then
                      write(iou,*)' Renormalized f(u,E) at E= ',e
                      write(iou,*)' Integral value          = ',c
                      write(iou,*)' Renormalization factor  = ',fn
                      write(iou,*)' ======================'
                     endif
                     nl=nmu
                     nw=2*nmu
                     if (nw.le.nbmax) then
                       j=0
                       do i=1,nmu
                         j=j+1
                         b(j)=x0(i)
                         j=j+1
                         b(j)=y0(i)
                       enddo
                     else
                       write(iou,*)' ERROR: increase size of arrays,',
     &                 ' nbmax=',nbmax,'is too small ',nw,' required'
                       stop
                     endif
                   elseif (lang.eq.14) then
c
c                    convert lang=14 to lang=12
c
                     write(iou,*)' LANG=14 converted to LANG=12'
                     j=0
                     do i=1,nl
                       j=j+1
                       x(i)=b(j)
                       j=j+1
                       y(i)=b(j)
                     enddo
                     c0=0.0d0
                     call renorm(nl,x,y,lang,c0,fn)
                     nrm=1
                     nbtm(1)=nl
                     ibtm(1)=lang-10
                     call linear(nrm,nbtm,ibtm,nl,x,y,etol,ymin,npmax)
                     do i=1,nl
                       f0=y(i)
                       if (f0.lt.0.0d0) y(i)=0.0d0
                     enddo
                     lang=12
                     c=1.0d0
                     call renorm(nl,x,y,lang,c,fn)
                     if (imon.gt.0) then
                       write(iou,*)' Renormalized f(u,E) at E= ',e
                       write(iou,*)' before u linearization  = ',c0
                       write(iou,*)' after  u linearization  = ',c
                       write(iou,*)' ======================'
                     endif
                     nw=2*nl
                     if (nw.le.nbmax) then
                       j=0
                       do i=1,nl
                         j=j+1
                         b(j)=x(i)
                         j=j+1
                         b(j)=y(i)
                       enddo
                     else
                       write(iou,*)' ERROR: increase size of arrays,',
     &                 ' nbmax=',nbmax,'is too small ',nw,' required'
                       stop
                     endif
                   else
                     write(iou,*)' Error: LAW=2 LANG=',lang,' is wrong'
                     stop
                   endif
                   call wrtlist(nou,mat,mf,mt,ns,c1,e,lang,l2,nw,nl,b)
                 enddo
               elseif (law.eq.7) then
c
c                law=7, converted to law=1, lang=11 or 12
c
                 open (jou, file='SIXLIN.TMP')
                 law=1
                 call wrtab1(nou,mat,mf,mt,ns,
     &                       zap,awp,lip,law,nr,nbt,ibt,np,x,y)
c
c                Copy to temporal file and explore the subsection
c
                 nsi=1
                 call readtab2(nin,c1e,c2e,l1,l2,nre,ne,nbte,ibte)
                 lang=11
                 lep=1
                 do ie=1,ne
                   call readtab2(nin,c1,e,l1,l2,nr,nmu,nbt,ibt)
                   call wrtab2(jou,mat,mf,mt,nsi,
     &               c1,e,l1,l2,nr,nmu,nbt,ibt)
                   call checklaw(nr,ibt,icod)
                   if (lang.eq.11.and.icod.ne.1) lang=12
                   do j=1,nmu
                    call readtab1(nin,c1,xmu,l1,l2,nr,nep,nbt,ibt,x,y)
                    call wrtab1(jou,mat,mf,mt,nsi,
     &                          c1,xmu,l1,l2,nr,nbt,ibt,nep,x,y)
                    call checklaw(nr,ibt,icod)
                    if (lep.eq.1.and.icod.ne.1) lep=2
                   enddo
                 enddo
                 rewind(jou)
                 write(iou,*)' LAW= 7 converted to LAW= 1 ',' LANG= ',
     &           lang,' LEP= ',lep
                 call wrtab2(nou,mat,mf,mt,ns,
     &           c1e,c2e,lang,lep,nre,ne,nbte,ibte)
c
c                prepare the union grid for secundary energy and
c                compute linearly interpolable energy-angle
c                distribution in the common grid of E' and
c                the original grid of cosine
c
                 do ie=1,ne
                   call readtab2(jou,c1,e,l1,l2,nrm,nmu0,nbtm,ibtm)
                   call readtab1(jou,c1,xmu,l1,l2,nr,nep,nbt,ibt,x,y)
                   y0(1)=xmu
                   x1(1)=xmu
                   call checklaw(nr,ibt,icod)
                   if (lep.ne.1.and.icod.ne.0) then
                     call linear(nr,nbt,ibt,nep,x,y,etol,ymin,npmax)
                   endif
                   k=0
                   nepmu(1)=nep
                   do i=1,nep
                     call inc(k,nbmax)
                     b0(k)=x(i)
                     call inc(k,nbmax)
                     b0(k)=y(i)
                   enddo
                   do i=1,nep
                     x0(i)=x(i)
                   enddo
                   nep0=nep
                   do j=2,nmu0
                     call readtab1(jou,c1,xmu,l1,l2,nr,nep,nbt,ibt,x,y)
                     y0(j)=xmu
                     x1(j)=xmu
                     call checklaw(nr,ibt,icod)
                     if (lep.ne.1.and.icod.ne.0) then
                       call linear(nr,nbt,ibt,nep,x,y,etol,ymin,npmax)
                     endif
                     nepmu(j)=nep
                     do i=1,nep
                       call inc(k,nbmax)
                       b0(k)=x(i)
                       call inc(k,nbmax)
                       b0(k)=y(i)
                     enddo
                     call union(nep0,x0,nep,x,nep1,y1,npmax)
                     nep0=nep1
                     do i=1,nep0
                       x0(i)=y1(i)
                     enddo
                   enddo
                   write(iou,*)
                   write(iou,*)' Law 7 common E'' grid at E= ',e,
     &                         ' nep= ',nep0
                   k=0
                   l=0
                   nr=1
                   ibt(1)=lep
                   do j=1,nmu0
                     nep=nepmu(j)
                     nbt(1)=nep
                     do i=1,nep
                       k=k+1
                       x(i)=b0(k)
                       k=k+1
                       y(i)=b0(k)
                     enddo
                     do i=1,nep0
                       yy=fvalue(nr,nbt,ibt,nep,x,y,x0(i))
                       if (yy.lt.0.0d0) yy=ymin
                       call inc(l,nbmax)
                       b(l)=yy
                     enddo
                   enddo
c
c                  linearize distribution for cosines on union grid
c
                   call checklaw(nrm,ibtm,icod)
                   nmu1=nmu0
                   if (lang.ne.11.and.icod.ne.0) then
                     do i=1,nep0
                       do j=1,nmu0
                         x(j)=y0(j)
                         ibj=(j-1)*nep0
                         y(j)=b(ibj+i)
                       enddo
                       nmu=nmu0
                       call linear(nrm,nbtm,ibtm,nmu,x,y,
     &                             etol,ymin,npmax)
                       call union(nmu1,x1,nmu,x,nmu2,y1,npmax)
                       nmu1=nmu2
                       do j=1,nmu1
                         x1(j)=y1(j)
                       enddo
                     enddo
                   endif
                   write(iou,*)' Law 7',
     &                         ' common cosine grid for all E'' =', nmu1
c
c                  prepare list data for law=1, lang=11 or 12
c
                   k=0
                   do i=1,nep0
                     do j=1,nmu0
                       ibj=(j-1)*nep0
                       y(j)=b(ibj+i)
                     enddo
                     c=fintval(nrm,nbtm,ibtm,nmu0,y0,y)
                     x(i)=c
                     write(iou,*)' E''=',x0(i),' integral f0 over mu=',c
                     do j=1,nmu1
                       yy=fvalue(nrm,nbtm,ibtm,nmu0,y0,y,x1(j))
                       if (yy.lt.0.0d0) yy=ymin
                       call inc(k,npmax)
                       b0(k)=yy
                     enddo
                   enddo
                   k=0
                   do i=1,nep0
                     do j=1,nmu1
                       ibj=(i-1)*nep0
                       y(j)=b0(ibj+j)
                     enddo
                     c=x(i)
                     call renorm(nmu1,x1,y,lang,c,fn)
                     if (imon.gt.0) then
                       write(iou,*)' Law 7 ',
     &                 'lin-lin renormalization factor= ',fn
                     endif
                     call inc(k,nbmax)
                     b(k)=x0(i)
                     call inc(k,nbmax)
                     b(k)=c
                     do j=1,nmu1
                       call inc(k,nbmax)
                       b(k)=x1(j)
                       call inc(k,nbmax)
                       if (c.ne.0.0d0) then
                         b(k)=0.5d0*y(j)/c
                       else
                         b(k)=0.0d0
                       endif
                     enddo
                   enddo
                   c1=0.0d0
                   nd=0
                   na=2*nmu1
                   nw=nep0*(na+2)
                   call wrtlist(nou,mat,mf,mt,ns,c1,e,nd,na,nw,nep0,b)
                 enddo
                 close (jou,status='DELETE')
               else
c
c                For MF6/LAW=0,3,4,5 & 6 data is copied
c
                 write(iou,*)' zap=',zap,' awp=',awp,' lip=',lip,
     &                       ' law=',law, ' data copied'
                 call wrtab1(nou,mat,mf,mt,ns,
     &                       zap,awp,lip,law,nr,nbt,ibt,np,x,y)
                 if (law.eq.5) then
                   call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
                   call wrtab2(nou,mat,mf,mt,ns,
     &                         c1,c2,l1,l2,nr,ne,nbt,ibt)
                   do i=1,ne
                     call readlist(nin,c1,e,l1,l2,nw,nl,b)
                     call wrtlist(nou,mat,mf,mt,ns,c1,e,l1,l2,nw,nl,b)
                   enddo
                 elseif (law.eq.6) then
                   call readcont(nin,ap,c2,l1,l2,n1,np,mat,mf,mt,nsi)
                   call wrtcont(nou,mat,mf,mt,ns,ap,c2,l1,l2,n1,np)
                 elseif (law.gt.7) then
                   write(iou,*)' Error: LAW= ',law,' not allowed on MF6'
                   stop
                 endif
               endif
c
c              end of nk subsection cycle
c
             enddo
             call wrtsend(nou,mat,mf,ns)
             call readtext(nin,line,mat,mf,mt,nsi)
             if (mt.ne.0) then
               write(iou,*)' ERROR: SEND record expected'
               write(iou,*)' MAT=',mat,' MF=',mf,' MT=',mt
               stop
             endif
            else
             write(iou,*)' Section MT=',mt,' with JP= ',jp,
     &                   ' was skipped'
             write(*,*)' Section MT=',mt,' with JP= ',jp,
     &                   ' was skipped'
             call findnextmt(nin,mf,mt)
             if (mt.lt.0) backspace(nin)
            endif
c
c           end of mf6 cycle
c
            call readcont(nin,za,awr,jp,lct,nk,n2,mat,mf,mt,nsi)
           enddo
           call wrtfend(nou,mat,ns)
           if (mf.ne.0) then
               write(iou,*)' ERROR: FEND record expected'
               write(iou,*)' MAT=',mat,' MF=',mf,' MT=',mt
               stop
           endif
         else
c
c          No MF6 data
c
           write(iou,*)' No MF6 data for this material '
         endif
c
c          Copy rest of the file to output file
c
         write(*,*)' MF7 to MF40 copied to output tape'
         write(iou,*)' MF7 to MF40 copied to output tape'
         call readtext(nin,line,mat,mf,mt,nsi)
         do while (mat.ne.0)
           if (mt.eq.0.and.mf.ne.0) then
             call wrtsend(nou,mat,mf,ns)
           elseif (mf.eq.0) then
             call wrtfend(nou,mat,ns)
           else
             call wrtext(nou,mat,mf,mt,ns,line)
           endif
           call readtext(nin,line,mat,mf,mt,nsi)
         enddo
         call wrtmend(nou,ns)
         ns=0
       endif
       call readcont(nin,za,awr,l1,l2,n1,n2,mat,mf,mt,nsi)
       call checksel(isel,mat,za,llower,lupper,nsel)
      enddo
      call wrtend(nou,ns)
      close (nin)
      close (iou)
      close (nou)
      stop
      end
C======================================================================
      subroutine checksel(isel,mat,za,matlow,matup,nsel)
c
c     Check if material is in the selection range
c         nsel=-1 end of tape
c         nsel=0 mat is not in the mat/za range
c         nsel=mat/za process material
c
      implicit real*8 (a-h, o-z)
      if (mat.eq.-1) then
        nsel=-1
      else
        if (isel.eq.0) then
          nsel=mat
        else
          nsel=int(za)
        endif
        if (nsel.lt.matlow.or.nsel.gt.matup) nsel=0
      endif
      return
      end
C======================================================================
C     Processing and linearization rutines
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
      real*8 function fvalue(nr,nbt,ibt,np,x,y,x0)
c
c      Return f(x0): value of the function f evaluated at x0.
c      Function f is given by np tabulated points (x,y) and
c      nr intervals for the interpolation law (arrays nbt & ibt)
c
        implicit real*8 (a-h, o-z)
        dimension nbt(*),ibt(*),x(*),y(*)
        if (x0.lt.x(1).or.x0.gt.x(np)) then
          fvalue=0.0d0
        else
          i=1
          do while (i.le.np.and.x(i).lt.x0)
            i=i+1
          enddo
          if (x0.eq.x(i)) then
            fvalue=y(i)
          else
            j=1
            do while (nbt(j).lt.i)
              j=j+1
            enddo
            ilaw=ibt(j)
            call terp1m(x(i-1),y(i-1),x(i),y(i),x0,y0,ilaw)
            fvalue=y0
          endif
        endif
        return
      end
C======================================================================
      subroutine terp1m(x1,y1,x2,y2,x,y,i)
c
c      interpolate one point using ENDF-6 interpolation laws
c      (borrowed and modified from NJOY)
c      (x1,y1) and (x2,y2) are the end points
c      (x,y) is the interpolated point
c      i is the interpolation law (1-5)
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
          call inc(ii,npmax)
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
              if (ilaw.eq.4) then
                xm=0.5d0*(x0+x1)
              else
                xm=dsqrt(x0*x1)
              endif
              ymlin=0.5d0*(y0+y1)
              call terp1m(x0,y0,x1,y1,xm,ymlaw,ilaw)
              dy=dabs(ymlaw-ymlin)
              ymabs=dabs(ymlaw)
              ylabs=dabs(ymlin)
              dx=dabs(xm-x0)
              h0=dabs(edelta(x0,1.0d0)-x0)
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
      subroutine renorm(np,x,y,ilaw,c,fn)
c
c     renormalize tabulated data with interpolation law=ilaw to c
c     if c=0.0d0 return the value of integral in c
c
c       lep=ilaw-(ilaw/10)*10
c              y        x
c       lep=1 constant
c       lep=2 lin      lin
c       lep=4 log      lin
c
c
      implicit real*8 (a-h, o-z)
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
          y1=y(j1)
          if (y1.ne.0.0d0.and.xr.ne.0.0d0) then
            yr=y(j)/y1
            if (yr.gt.1.0d-12) then
              sum=sum+y1*xr*(yr-1.0d0)/dlog(yr)
            else
              sum=sum+y1*xr
            endif
          endif
        enddo
      else
        write(*,*)' === Renormalization problem in renorm '
        write(*,*)'     ilaw=',ilaw,' not coded'
        stop
      endif
      if (c.eq.0.0d0) then
        c=sum
        fn=1.0d0
      else
        if (sum.eq.0.0d0) then
          fn=1.0d0
        else
          fn=c/sum
        endif
        do j=1,np
          y(j)=y(j)*fn
        enddo
      endif
      return
      end
C======================================================================
      real*8 function fintval(nr,nbt,ibt,np,x,y)
c
c     Calculate the integral value of y over x given as a TAB1 record
c
      implicit real*8 (a-h, o-z)
      dimension nbt(*),ibt(*),x(*),y(*)
      fintval=0.0d0
      i0=1
      do j=1,nr
        ilaw=ibt(j)
        npt=nbt(j)
        npx=npt-i0+1
        sum=0.0d0
        call renorm(npx,x(i0),y(i0),ilaw,sum,fn)
        fintval=fintval+sum
        i0=npt
      enddo
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
      subroutine union(np1,x1,np2,x2,np3,x3,npmax)
      implicit real*8 (a-h, o-z)
c
c     Prepare union grid x3 from x1 and x2
c
      dimension x1(*),x2(*),x3(*)
      i=1
      j=1
      k=0
      do while (i.le.np1.and.j.le.np2)
        call inc(k,npmax)
        if (x1(i).lt.x2(j)) then
          x3(k)=x1(i)
          i=i+1
        elseif (x1(i).eq.x2(j)) then
          x3(k)=x1(i)
          i=i+1
          j=j+1
        else
          x3(k)=x2(j)
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
      real*8 function bvalue(x0,n,x,y)
      implicit real*8 (a-h, o-z)
      dimension x(*), y(*)
c
c      delta function. take values at discrete x
c
      bvalue=0.0d0
      if (x0.lt.x(1).or.x0.gt.x(n)) return
      istop=0
      i=1
      do while (istop.eq.0.and.i.le.n)
        if (x(i).eq.x0) then
          bvalue=y(i)
          istop=1
        else
          i=i+1
        endif
      enddo
      return
      end
C======================================================================
      SUBROUTINE THINXS(XI,YI,N,XO,YO,M,ERR)
C-Title  : THINXS Subroutine
C-Purpose: Thin the data to within the specified tolerance
C-Version: 1996 Original code
C-V 2002/04 A.Trkov: Fix bug (missing one but last point).
C-Author : A.Trkov, ENEA, Bologna, 1996
C-Description:
C-D  Excessive data points need to be removed such that the remaining
C-D  points can be linearly interpolated to within the specified
C-D  tolerance. The formal parameters are:
C-D    N   number of points on input argument mesh
C-D   XI   input argument mesh (array of N values in momotonic order)
C-D   YI   function values corresponding to XI(i)
C-D    M   number of points in the output mesh
C-D   XO   output argument mesh (array of M values in monotonic order)
C-D   YO   interpolated function values corresponding to XO(i) (Output)
C-D  ERR   fractional tolerance on the interpolated array.
C-D
c
c      borrowed from Andrej Trkov
c
      implicit real*8 (a-h, o-z)
      DIMENSION XI(*), YI(*), XO(*), YO(*)
      FINT(X,XA,XB,YA,YB)=YA+(YB-YA)*(X-XA)/(XB-XA)
      K=3
      M=1
      I=1
      XO(M)=XI(I)
      YO(M)=YI(I)
C* Begin loop over data points
   10 IF(K.GT.N) GO TO 60
   12 L1=I+1
      L2=K-1
      DO L=L1,L2
        IF(XI(K).EQ.XI(I)) CYCLE
        Y=FINT(XI(L),XI(I),XI(K),YI(I),YI(K))
        E=ABS(Y-YI(L))
        IF(YI(L).NE.0) E=DABS(E/YI(L))
        IF(E.GT.ERR) GO TO 40
      END DO
C* Linear interpolation adequate - increase the test interval
      K=K+1
      GO TO 10
C* Add the previous point at K-1 - lin.interpolation violated
   40 M=M+1
      I=K-1
      K=I+2
      XO(M)=XI(I)
      YO(M)=YI(I)
      IF(K.LE.N) GO TO 12
C* Add the last point
   60 M=M+1
      XO(M)=XI(N)
      YO(M)=YI(N)
      RETURN
      END
C======================================================================
C      Legendre polynomials Pn(x)
C======================================================================
      subroutine legndr(x,p,np)
c     *****************************************************************
c       generate legendre polynomials at x by recursion.
c       place pl in p(l+1). (borrowed from njoy)
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
        yleg=yleg+0.5d0*(2.0d0*l-1.0d0)*a(l)*p(l)
      enddo
      return
      end
C======================================================================
      subroutine yleg2(x,a,na,yplus,yminus)
c
c      calculate yplus=y(x) and yminus=y(-x) by legendre expansion
c
      implicit real*8 (a-h,o-z)
      dimension a(*),p(65)
      call legndr(x,p,na)
      n=na+1
      yplus=0.0d0
      do l=1,n,2
        yplus=yplus+0.5d0*(2.0d0*l-1.0d0)*a(l)*p(l)
      enddo
      yminus=yplus
      if (na.gt.0) then
        do l=2,n,2
          dy=0.5d0*(2.0d0*l-1.0d0)*a(l)*p(l)
          yplus=yplus+dy
          yminus=yminus-dy
        enddo
      endif
      return
      end
C======================================================================
      subroutine setxmu(xmu,nmu,dxmax,npmax)
c
c     xmu output mu grid (mu between 0.0 and 1.0)
c     dxmax maximum cosine interval
c     nmu number of point in the output grid
c     npmax maximun dimension of array xmu
c
      implicit real*8(a-h, o-z)
      parameter (nmu0=2015)
      dimension xmu(*),x0(nmu0)
c
c       Initial unified grid x0 includes non-negative roots
c       of Pn & Pn'for n=1 to 64
c
c       if x0(i)-x0(i-1) > dxmax new points are included in between
c
      data x0/         0.00000000000d0,0.02435029266d0,0.02473672762d0,
     & 0.02512929142d0,0.02554115886d0,0.02595977230d0,0.02639966790d0,
     & 0.02684701237d0,0.02731789740d0,0.02779703529d0,0.02830230335d0,
     & 0.02881674820d0,0.02936030704d0,0.02991410980d0,0.03109833833d0,
     & 0.03173278944d0,0.03238017096d0,0.03306886477d0,0.03377219002d0,
     & 0.03452239133d0,0.03528923696d0,0.03610956813d0,0.03694894317d0,
     & 0.03784971656d0,0.03877241751d0,0.03976607080d0,0.04078514790d0,
     & 0.04188682107d0,0.04301819847d0,0.04424650945d0,0.04550982195d0,
     & 0.04688792571d0,0.04830766569d0,0.04869199548d0,0.04945218712d0,
     & 0.04986472505d0,0.05024914171d0,0.05105890671d0,0.05147184256d0,
     & 0.05190913587d0,0.05277348409d0,0.05324511049d0,0.05368251154d0,
     & 0.05460715100d0,0.05507928988d0,0.05558129096d0,0.05657275382d0,
     & 0.05711712169d0,0.05761925760d0,0.05868505430d0,0.05923009343d0,
     & 0.05981229075d0,0.06096110015d0,0.06159641178d0,0.06217877935d0,
     & 0.06342068498d0,0.06405689286d0,0.06474013799d0,0.06608692392d0,
     & 0.06683799374d0,0.06752145544d0,0.06898698016d0,0.06973927332d0,
     & 0.07055231742d0,0.07215299087d0,0.07299312179d0,0.07305454001d0,
     & 0.07386786070d0,0.07414964795d0,0.07532439550d0,0.07562325899d0,
     & 0.07652652113d0,0.07655684261d0,0.07751013794d0,0.07780933395d0,
     & 0.07912542309d0,0.07944380461d0,0.08046363021d0,0.08054593724d0,
     & 0.08152990574d0,0.08187216487d0,0.08330518682d0,0.08367040895d0,
     & 0.08477501304d0,0.08481624924d0,0.08598899710d0,0.08635451826d0,
     & 0.08797971005d0,0.08837134328d0,0.08963524465d0,0.08974909348d0,
     & 0.09096351370d0,0.09138798374d0,0.09317470156d0,0.09363106585d0,
     & 0.09501250984d0,0.09507059157d0,0.09654818818d0,0.09700469921d0,
     & 0.09726849893d0,0.09878335645d0,0.09906199255d0,0.09955531215d0,
     & 0.10037134972d0,0.10116247531d0,0.10132627352d0,0.10198460656d0,
     & 0.10286244876d0,0.10340265904d0,0.10367833389d0,0.10539987902d0,
     & 0.10569190171d0,0.10627823013d0,0.10721024249d0,0.10805494871d0,
     & 0.10814044578d0,0.10905133281d0,0.11005901340d0,0.11064502721d0,
     & 0.11099078332d0,0.11296428806d0,0.11333234989d0,0.11397258561d0,
     & 0.11504710982d0,0.11608407068d0,0.11633186888d0,0.11716780907d0,
     & 0.11833633390d0,0.11904679844d0,0.11941046900d0,0.12146281930d0,
     & 0.12169542102d0,0.12208402534d0,0.12286469261d0,0.12338111165d0,
     & 0.12411700106d0,0.12523340851d0,0.12532922362d0,0.12536665716d0,
     & 0.12658599727d0,0.12737279830d0,0.12795705948d0,0.12873610381d0,
     & 0.12920873332d0,0.12944913540d0,0.13163064152d0,0.13188486655d0,
     & 0.13239323930d0,0.13325682430d0,0.13384825060d0,0.13473482592d0,
     & 0.13615235726d0,0.13618209360d0,0.13655293285d0,0.13764520598d0,
     & 0.13855584681d0,0.13927620404d0,0.14025172418d0,0.14075314670d0,
     & 0.14105850309d0,0.14360542732d0,0.14392980951d0,0.14447196158d0,
     & 0.14556185416d0,0.14561429223d0,0.14629582875d0,0.14733228069d0,
     & 0.14787278636d0,0.14887433902d0,0.14903550861d0,0.14909859681d0,
     & 0.15024001099d0,0.15081335486d0,0.15193551564d0,0.15264424023d0,
     & 0.15278551580d0,0.15386991361d0,0.15455412068d0,0.15489059000d0,
     & 0.15516803343d0,0.15773250559d0,0.15802557795d0,0.15838534000d0,
     & 0.15913204263d0,0.16035864564d0,0.16042885854d0,0.16122235607d0,
     & 0.16251724326d0,0.16317006259d0,0.16456928213d0,0.16462194740d0,
     & 0.16527895767d0,0.16605720952d0,0.16675393024d0,0.16809117947d0,
     & 0.16899396365d0,0.16918602341d0,0.16964442042d0,0.17060675531d0,
     & 0.17134136155d0,0.17179016593d0,0.17209278708d0,0.17231064109d0,
     & 0.17501745925d0,0.17524666216d0,0.17556801478d0,0.17605106117d0,
     & 0.17685882036d0,0.17785645288d0,0.17848418150d0,0.17858118873d0,
     & 0.17960752913d0,0.18073996487d0,0.18117327477d0,0.18197702696d0,
     & 0.18343464250d0,0.18373680656d0,0.18376898173d0,0.18385549527d0,
     & 0.18557503816d0,0.18643929883d0,0.18684695184d0,0.18816582565d0,
     & 0.18924159246d0,0.18951197352d0,0.19008560142d0,0.19111886747d0,
     & 0.19219493147d0,0.19269758070d0,0.19313538225d0,0.19337823864d0,
     & 0.19361470451d0,0.19660034679d0,0.19684890401d0,0.19710611028d0,
     & 0.19757487372d0,0.19812119934d0,0.19932125339d0,0.19972915293d0,
     & 0.20037929361d0,0.20119409400d0,0.20133343222d0,0.20257045389d0,
     & 0.20290564252d0,0.20410763539d0,0.20463452925d0,0.20564748978d0,
     & 0.20623942737d0,0.20786042669d0,0.20790226416d0,0.20796713597d0,
     & 0.20929921790d0,0.20962550339d0,0.21025275077d0,0.21135228617d0,
     & 0.21191783739d0,0.21318491664d0,0.21350089232d0,0.21495624486d0,
     & 0.21535395536d0,0.21600723688d0,0.21680182880d0,0.21742364374d0,
     & 0.21760658516d0,0.21878205828d0,0.21950381217d0,0.21999202271d0,
     & 0.22034425103d0,0.22061036236d0,0.22081849750d0,0.22426358560d0,
     & 0.22448230065d0,0.22476379039d0,0.22513960563d0,0.22566669162d0,
     & 0.22645936544d0,0.22778585114d0,0.22787610026d0,0.22856678936d0,
     & 0.22946205227d0,0.23045831596d0,0.23066859656d0,0.23154355138d0,
     & 0.23238298487d0,0.23272140372d0,0.23425292221d0,0.23501148310d0,
     & 0.23539512480d0,0.23632551246d0,0.23711263411d0,0.23861918608d0,
     & 0.23928736225d0,0.23930692497d0,0.23935901433d0,0.23955170592d0,
     & 0.24115588409d0,0.24158166645d0,0.24242305505d0,0.24342181906d0,
     & 0.24386688372d0,0.24456945693d0,0.24484679325d0,0.24631512143d0,
     & 0.24685065885d0,0.24760290943d0,0.24866779279d0,0.24871376168d0,
     & 0.24928693011d0,0.25013822183d0,0.25113517861d0,0.25188622569d0,
     & 0.25200873850d0,0.25263768717d0,0.25381306417d0,0.25463692617d0,
     & 0.25542517408d0,0.25582507934d0,0.25625195419d0,0.25648752007d0,
     & 0.25675483625d0,0.26093423734d0,0.26121584067d0,0.26146545921d0,
     & 0.26192150439d0,0.26235294121d0,0.26321594372d0,0.26413568097d0,
     & 0.26468716221d0,0.26532630740d0,0.26602478360d0,0.26636265288d0,
     & 0.26701340013d0,0.26815218501d0,0.26878597402d0,0.26954315595d0,
     & 0.26978657316d0,0.26988177638d0,0.27111181080d0,0.27206162763d0,
     & 0.27266976975d0,0.27294320270d0,0.27448162118d0,0.27485381671d0,
     & 0.27584154895d0,0.27628819378d0,0.27730124489d0,0.27760909715d0,
     & 0.27870488605d0,0.27925155320d0,0.28160355078d0,0.28170880979d0,
     & 0.28172293742d0,0.28177567550d0,0.28187226662d0,0.28428151571d0,
     & 0.28486199803d0,0.28523151648d0,0.28604720149d0,0.28636517941d0,
     & 0.28736248736d0,0.28802131680d0,0.28812506853d0,0.28910887429d0,
     & 0.28939390645d0,0.29107691431d0,0.29145024436d0,0.29200483949d0,
     & 0.29249405859d0,0.29329877785d0,0.29471806998d0,0.29479527772d0,
     & 0.29575813559d0,0.29603157028d0,0.29668499534d0,0.29707009785d0,
     & 0.29817627734d0,0.29934582270d0,0.29983046890d0,0.30028760634d0,
     & 0.30106225387d0,0.30171062896d0,0.30198985651d0,0.30332751286d0,
     & 0.30423743127d0,0.30489646942d0,0.30539580403d0,0.30578720786d0,
     & 0.30610225925d0,0.30636131297d0,0.30657807946d0,0.31132287199d0,
     & 0.31151570080d0,0.31174372083d0,0.31201753212d0,0.31235246650d0,
     & 0.31277155925d0,0.31331108134d0,0.31403163787d0,0.31504267970d0,
     & 0.31609568620d0,0.31656409996d0,0.31670269369d0,0.31742358076d0,
     & 0.31829371625d0,0.31911236893d0,0.31936480943d0,0.32071558094d0,
     & 0.32093334159d0,0.32196616840d0,0.32247205451d0,0.32319500343d0,
     & 0.32425342340d0,0.32468148634d0,0.32484938284d0,0.32600294221d0,
     & 0.32651612447d0,0.32750447298d0,0.32824761338d0,0.32883742988d0,
     & 0.32929723037d0,0.33114284827d0,0.33147504697d0,0.33186860228d0,
     & 0.33312627889d0,0.33350484782d0,0.33417692014d0,0.33441085211d0,
     & 0.33550024542d0,0.33599390364d0,0.33653708285d0,0.33761790730d0,
     & 0.33839265425d0,0.33905188738d0,0.33942554197d0,0.33998104358d0,
     & 0.34193582089d0,0.34199409083d0,0.34200765360d0,0.34207248864d0,
     & 0.34214940654d0,0.34272401334d0,0.34467600914d0,0.34506880850d0,
     & 0.34576852737d0,0.34660155443d0,0.34775784689d0,0.34838758199d0,
     & 0.34875588629d0,0.34999644220d0,0.35039504491d0,0.35123176345d0,
     & 0.35142263066d0,0.35270472553d0,0.35328261286d0,0.35359103217d0,
     & 0.35557484232d0,0.35585298233d0,0.35635395157d0,0.35722015834d0,
     & 0.35752071014d0,0.35787645669d0,0.35897244048d0,0.35973251972d0,
     & 0.36117230581d0,0.36122891417d0,0.36131666551d0,0.36263185923d0,
     & 0.36307287702d0,0.36311746383d0,0.36431750042d0,0.36470514676d0,
     & 0.36596434037d0,0.36633925775d0,0.36716594408d0,0.36783149900d0,
     & 0.36811277505d0,0.36822272089d0,0.36950502264d0,0.37075818585d0,
     & 0.37164350126d0,0.37217443357d0,0.37253709592d0,0.37318489009d0,
     & 0.37370608872d0,0.37385413502d0,0.37550145786d0,0.37625151609d0,
     & 0.37722872425d0,0.37767254712d0,0.37828632473d0,0.37857935201d0,
     & 0.37900044368d0,0.37920826911d0,0.37951502715d0,0.37967005658d0,
     & 0.37990345007d0,0.38552639421d0,0.38577099704d0,0.38593390074d0,
     & 0.38625764249d0,0.38647776408d0,0.38692637647d0,0.38724016397d0,
     & 0.38790289737d0,0.38838590161d0,0.38946313758d0,0.39030103803d0,
     & 0.39141123888d0,0.39195229633d0,0.39235318371d0,0.39269596268d0,
     & 0.39341431190d0,0.39415134708d0,0.39442424896d0,0.39542385204d0,
     & 0.39687379687d0,0.39736915473d0,0.39835927776d0,0.39844627733d0,
     & 0.39953094097d0,0.39955418695d0,0.40061533828d0,0.40099189644d0,
     & 0.40227015796d0,0.40250294386d0,0.40305175512d0,0.40361303822d0,
     & 0.40451820029d0,0.40502688093d0,0.40584515138d0,0.40670050932d0,
     & 0.40679198166d0,0.40703793791d0,0.40828061129d0,0.40868648199d0,
     & 0.40972722656d0,0.40993531781d0,0.41110956696d0,0.41175116146d0,
     & 0.41195139468d0,0.41315288817d0,0.41377920437d0,0.41413398323d0,
     & 0.41436232372d0,0.41657055993d0,0.41685113203d0,0.41731103463d0,
     & 0.41820138707d0,0.41896926326d0,0.41986137603d0,0.42063805471d0,
     & 0.42072998503d0,0.42135127613d0,0.42189706075d0,0.42328988145d0,
     & 0.42372096216d0,0.42434212021d0,0.42449495897d0,0.42514331328d0,
     & 0.42636639842d0,0.42697347171d0,0.42717374158d0,0.42813754152d0,
     & 0.42920969987d0,0.42977299334d0,0.43035252204d0,0.43068379880d0,
     & 0.43339539411d0,0.43379350763d0,0.43384716943d0,0.43386406772d0,
     & 0.43392969819d0,0.43396233194d0,0.43404771720d0,0.43441503691d0,
     & 0.43715772980d0,0.43750526004d0,0.43811751132d0,0.43871727705d0,
     & 0.43987326413d0,0.44076683919d0,0.44114825175d0,0.44178115065d0,
     & 0.44292017453d0,0.44370517654d0,0.44411578328d0,0.44475374868d0,
     & 0.44503554137d0,0.44636601725d0,0.44658407310d0,0.44703376954d0,
     & 0.44721359550d0,0.44794558739d0,0.44849275104d0,0.44890178631d0,
     & 0.44980633497d0,0.45056316466d0,0.45131637321d0,0.45185001727d0,
     & 0.45247068810d0,0.45266221946d0,0.45293023223d0,0.45521081488d0,
     & 0.45557378326d0,0.45586394443d0,0.45687307561d0,0.45749915825d0,
     & 0.45801677766d0,0.45829665118d0,0.45926051231d0,0.45956515724d0,
     & 0.46070512279d0,0.46129119017d0,0.46217191207d0,0.46285066501d0,
     & 0.46457074137d0,0.46469512392d0,0.46477409494d0,0.46488816163d0,
     & 0.46650819849d0,0.46690290475d0,0.46787345771d0,0.46807961268d0,
     & 0.46884879347d0,0.46885090429d0,0.46935583799d0,0.47039088221d0,
     & 0.47058241248d0,0.47213160951d0,0.47254243619d0,0.47300273145d0,
     & 0.47352584176d0,0.47440247894d0,0.47478724799d0,0.47587422496d0,
     & 0.47602649750d0,0.47745673677d0,0.47792494981d0,0.47819378204d0,
     & 0.47872592746d0,0.47985982183d0,0.48010654519d0,0.48087895937d0,
     & 0.48171087780d0,0.48179992230d0,0.48290982109d0,0.48307580169d0,
     & 0.48425117679d0,0.48527391839d0,0.48605942189d0,0.48617194145d0,
     & 0.48696674570d0,0.48767515819d0,0.48822928568d0,0.48831053722d0,
     & 0.48888362226d0,0.48940314571d0,0.48981487519d0,0.49102411482d0,
     & 0.49197675393d0,0.49274661910d0,0.49338169867d0,0.49391453563d0,
     & 0.49436797813d0,0.49475853909d0,0.49509844947d0,0.49539696157d0,
     & 0.49566120372d0,0.49589675660d0,0.49610805165d0,0.49629865238d0,
     & 0.49647145694d0,0.50360708934d0,0.50378786656d0,0.50398771838d0,
     & 0.50420983166d0,0.50445814491d0,0.50473758386d0,0.50505439139d0,
     & 0.50541659920d0,0.50583471793d0,0.50632277324d0,0.50689990893d0,
     & 0.50759295512d0,0.50844071482d0,0.50950147785d0,0.51055403321d0,
     & 0.51086700195d0,0.51106934524d0,0.51163756581d0,0.51226728498d0,
     & 0.51269053709d0,0.51296905895d0,0.51375600678d0,0.51464463849d0,
     & 0.51524863636d0,0.51565602649d0,0.51681749885d0,0.51772881329d0,
     & 0.51816514710d0,0.51860140006d0,0.51909612921d0,0.51956431139d0,
     & 0.51974764606d0,0.52063233439d0,0.52163226288d0,0.52182366937d0,
     & 0.52316097472d0,0.52391467437d0,0.52467282046d0,0.52522504631d0,
     & 0.52553240992d0,0.52639574993d0,0.52648792454d0,0.52673574203d0,
     & 0.52788390227d0,0.52837726866d0,0.52943518887d0,0.53031177114d0,
     & 0.53068028593d0,0.53116922515d0,0.53127946402d0,0.53278662650d0,
     & 0.53312031983d0,0.53338990479d0,0.53444630965d0,0.53499286403d0,
     & 0.53533194107d0,0.53628288591d0,0.53662414814d0,0.53785999096d0,
     & 0.53832620929d0,0.53846931011d0,0.53879773272d0,0.54055156458d0,
     & 0.54061324699d0,0.54069620659d0,0.54077759250d0,0.54138539933d0,
     & 0.54278993437d0,0.54319033026d0,0.54418227174d0,0.54511065249d0,
     & 0.54542147139d0,0.54611631666d0,0.54637686630d0,0.54769736429d0,
     & 0.54820705992d0,0.54867242781d0,0.54946712510d0,0.55059854443d0,
     & 0.55063940293d0,0.55120682486d0,0.55161883589d0,0.55181747590d0,
     & 0.55303826009d0,0.55334239186d0,0.55387519688d0,0.55401932828d0,
     & 0.55429717033d0,0.55702585299d0,0.55715830451d0,0.55760517791d0,
     & 0.55787550067d0,0.55894506094d0,0.55940340949d0,0.55977083107d0,
     & 0.56004295821d0,0.56068400593d0,0.56188943929d0,0.56227890075d0,
     & 0.56324916140d0,0.56339669876d0,0.56467245319d0,0.56523532700d0,
     & 0.56544642927d0,0.56633135798d0,0.56686126418d0,0.56714664792d0,
     & 0.56732340617d0,0.56895276820d0,0.56922094161d0,0.56972047181d0,
     & 0.57062947330d0,0.57097217261d0,0.57136728774d0,0.57189564620d0,
     & 0.57270031127d0,0.57285521635d0,0.57427632365d0,0.57445602105d0,
     & 0.57529965135d0,0.57583196026d0,0.57615297233d0,0.57722472608d0,
     & 0.57735026919d0,0.57766293024d0,0.57831937726d0,0.57904113513d0,
     & 0.57965465721d0,0.57980548007d0,0.58054534475d0,0.58162500892d0,
     & 0.58215021257d0,0.58282672265d0,0.58317273803d0,0.58363613929d0,
     & 0.58731795429d0,0.58764040351d0,0.58771575724d0,0.58774459749d0,
     & 0.58775860498d0,0.58776644795d0,0.58785191841d0,0.58788294215d0,
     & 0.58794198737d0,0.58807668984d0,0.58850483432d0,0.59170018143d0,
     & 0.59177206842d0,0.59202774070d0,0.59251374538d0,0.59287769411d0,
     & 0.59359440564d0,0.59415345496d0,0.59531523208d0,0.59607889727d0,
     & 0.59628179714d0,0.59672218277d0,0.59769605436d0,0.59848414728d0,
     & 0.59862828971d0,0.59970905188d0,0.60009673125d0,0.60054530466d0,
     & 0.60080809029d0,0.60156765814d0,0.60191900571d0,0.60349072517d0,
     & 0.60403258715d0,0.60444059705d0,0.60513425964d0,0.60602464479d0,
     & 0.60625320547d0,0.60669229302d0,0.60756551530d0,0.60770292718d0,
     & 0.60798865775d0,0.60964103291d0,0.61001247602d0,0.61024234584d0,
     & 0.61115535517d0,0.61166943828d0,0.61180752430d0,0.61255388967d0,
     & 0.61318057522d0,0.61337143270d0,0.61417869996d0,0.61449625220d0,
     & 0.61538319833d0,0.61623211860d0,0.61631178520d0,0.61740636950d0,
     & 0.61787624440d0,0.61825358920d0,0.61889368960d0,0.61939434614d0,
     & 0.61960987576d0,0.62052618299d0,0.62109260841d0,0.62147734590d0,
     & 0.62175570460d0,0.62196643526d0,0.62500269199d0,0.62520997864d0,
     & 0.62548291841d0,0.62585845276d0,0.62640749128d0,0.62728529949d0,
     & 0.62767128065d0,0.62819451225d0,0.62886739678d0,0.62890813726d0,
     & 0.62976483907d0,0.63032304286d0,0.63102172708d0,0.63115952108d0,
     & 0.63222937633d0,0.63287615303d0,0.63290797195d0,0.63364612230d0,
     & 0.63373296624d0,0.63508697770d0,0.63561114735d0,0.63605368073d0,
     & 0.63630552228d0,0.63685339445d0,0.63794856689d0,0.63851917581d0,
     & 0.63854710582d0,0.63925441583d0,0.64007490127d0,0.64039310681d0,
     & 0.64234933944d0,0.64270672292d0,0.64275483242d0,0.64285804212d0,
     & 0.64293455606d0,0.64326364446d0,0.64497282849d0,0.64547127542d0,
     & 0.64588338887d0,0.64698600794d0,0.64711807708d0,0.64774374392d0,
     & 0.64809365194d0,0.64889961172d0,0.64896547125d0,0.65006505145d0,
     & 0.65022436467d0,0.65133484620d0,0.65201622828d0,0.65238870288d0,
     & 0.65317066370d0,0.65359483428d0,0.65385165543d0,0.65402263270d0,
     & 0.65465367071d0,0.65571603210d0,0.65589646569d0,0.65617321343d0,
     & 0.65665109404d0,0.65746745031d0,0.65767115922d0,0.65802823726d0,
     & 0.65889617020d0,0.65976938763d0,0.66041820261d0,0.66099731375d0,
     & 0.66120938647d0,0.66159809035d0,0.66304426693d0,0.66309844533d0,
     & 0.66324276694d0,0.66377640229d0,0.66468434775d0,0.66498774739d0,
     & 0.66594984387d0,0.66599056134d0,0.66687104088d0,0.66713880420d0,
     & 0.66737896056d0,0.66800123659d0,0.66834322118d0,0.66994535130d0,
     & 0.67031018721d0,0.67120399032d0,0.67124010526d0,0.67195668461d0,
     & 0.67258349540d0,0.67356636847d0,0.67367186450d0,0.67383550008d0,
     & 0.67487209865d0,0.67518607067d0,0.67582252811d0,0.67651012893d0,
     & 0.67674461306d0,0.67718627951d0,0.67787237963d0,0.67821453760d0,
     & 0.67918605679d0,0.67940956828d0,0.68014190423d0,0.68042975562d0,
     & 0.68126418704d0,0.68173195997d0,0.68208461269d0,0.68305446754d0,
     & 0.68345910915d0,0.68376632738d0,0.68448630913d0,0.68461281980d0,
     & 0.68523631305d0,0.68587058508d0,0.68618846908d0,0.68670150203d0,
     & 0.68729290481d0,0.68783573333d0,0.68852168077d0,0.68946797481d0,
     & 0.69004382443d0,0.69084530037d0,0.69102898063d0,0.69133557560d0,
     & 0.69168704306d0,0.69202310339d0,0.69244555119d0,0.69304181155d0,
     & 0.69340959089d0,0.69393162129d0,0.69405102606d0,0.69448726319d0,
     & 0.69611704882d0,0.69642726042d0,0.69761866136d0,0.69785049479d0,
     & 0.69875931662d0,0.69893911322d0,0.69965518373d0,0.69979868038d0,
     & 0.70037741678d0,0.70049459056d0,0.70097203119d0,0.70106951202d0,
     & 0.70147010377d0,0.70155246871d0,0.70189337938d0,0.70196388972d0,
     & 0.70225752896d0,0.70231857115d0,0.70257413177d0,0.70262749222d0,
     & 0.70285192892d0,0.71144409958d0,0.71167801054d0,0.71173311868d0,
     & 0.71200021587d0,0.71206339999d0,0.71237128303d0,0.71244445758d0,
     & 0.71280323742d0,0.71288897341d0,0.71331240270d0,0.71341423527d0,
     & 0.71392150931d0,0.71404443589d0,0.71466317903d0,0.71481450156d0,
     & 0.71558596072d0,0.71577678459d0,0.71676539864d0,0.71701347374d0,
     & 0.71832581636d0,0.71866136313d0,0.71988185017d0,0.72037817090d0,
     & 0.72048723996d0,0.72071651336d0,0.72096617734d0,0.72128146977d0,
     & 0.72166783445d0,0.72231671556d0,0.72276209975d0,0.72351509755d0,
     & 0.72367932928d0,0.72403413092d0,0.72441773136d0,0.72491847598d0,
     & 0.72553105366d0,0.72658437268d0,0.72731825519d0,0.72859405162d0,
     & 0.72869762509d0,0.72886859909d0,0.72929412345d0,0.72948917159d0,
     & 0.73012935825d0,0.73015200557d0,0.73080834474d0,0.73106617838d0,
     & 0.73176439622d0,0.73218211874d0,0.73254423081d0,0.73364933468d0,
     & 0.73418113631d0,0.73455425424d0,0.73561087801d0,0.73584617911d0,
     & 0.73665780994d0,0.73690884895d0,0.73748952528d0,0.73822714985d0,
     & 0.73843928529d0,0.73856989984d0,0.73877386510d0,0.73951373102d0,
     & 0.73970480307d0,0.74012419158d0,0.74074632932d0,0.74153118560d0,
     & 0.74154641915d0,0.74182653881d0,0.74307883398d0,0.74324592463d0,
     & 0.74369504117d0,0.74449430223d0,0.74532464832d0,0.74533709877d0,
     & 0.74614634155d0,0.74633190646d0,0.74662723545d0,0.74723049645d0,
     & 0.74760535962d0,0.74781064528d0,0.74929204741d0,0.74955224132d0,
     & 0.75006449394d0,0.75064185635d0,0.75127993569d0,0.75149420255d0,
     & 0.75234152267d0,0.75246285173d0,0.75281990726d0,0.75360812189d0,
     & 0.75389535449d0,0.75432378393d0,0.75540440836d0,0.75568590375d0,
     & 0.75572377531d0,0.75586522579d0,0.75612419401d0,0.75742883037d0,
     & 0.75767291845d0,0.75851921157d0,0.75902042271d0,0.75925926304d0,
     & 0.75998277817d0,0.76096977992d0,0.76106487663d0,0.76168038341d0,
     & 0.76211174719d0,0.76279499519d0,0.76327601117d0,0.76351968995d0,
     & 0.76417048242d0,0.76458700179d0,0.76485756743d0,0.76504499788d0,
     & 0.76505532393d0,0.76518190305d0,0.76699011936d0,0.76715903252d0,
     & 0.76740124293d0,0.76777743210d0,0.76843996348d0,0.76871641976d0,
     & 0.76916162643d0,0.76978470039d0,0.76990267419d0,0.77048859606d0,
     & 0.77071866906d0,0.77122765493d0,0.77226147925d0,0.77227289721d0,
     & 0.77257461363d0,0.77372414424d0,0.77381025229d0,0.77457668175d0,
     & 0.77459666924d0,0.77536489030d0,0.77536826095d0,0.77610689435d0,
     & 0.77638594882d0,0.77694552461d0,0.77738126299d0,0.77789730643d0,
     & 0.77830565143d0,0.77900671360d0,0.77922491535d0,0.77962859476d0,
     & 0.78151400390d0,0.78173314842d0,0.78178431259d0,0.78180388986d0,
     & 0.78193784139d0,0.78202923715d0,0.78231965924d0,0.78397235894d0,
     & 0.78434769638d0,0.78448347366d0,0.78455583290d0,0.78519295380d0,
     & 0.78557623013d0,0.78636390923d0,0.78681578113d0,0.78689035724d0,
     & 0.78769298064d0,0.78781680598d0,0.78847114505d0,0.78869373993d0,
     & 0.78971828428d0,0.79012465708d0,0.79073005708d0,0.79139232699d0,
     & 0.79177163907d0,0.79200829186d0,0.79253395260d0,0.79254171210d0,
     & 0.79279918183d0,0.79402692289d0,0.79439821707d0,0.79448379597d0,
     & 0.79526659928d0,0.79588361012d0,0.79600192608d0,0.79645920051d0,
     & 0.79666647741d0,0.79709497829d0,0.79796205326d0,0.79810174809d0,
     & 0.79847718311d0,0.79895171880d0,0.79914375417d0,0.80009728343d0,
     & 0.80016154319d0,0.80088289455d0,0.80138180446d0,0.80154134610d0,
     & 0.80157809073d0,0.80230653396d0,0.80303148693d0,0.80361508625d0,
     & 0.80370495897d0,0.80409499919d0,0.80449660231d0,0.80488840162d0,
     & 0.80564137092d0,0.80616235627d0,0.80654416761d0,0.80683596414d0,
     & 0.80706620403d0,0.80725249842d0,0.80740632791d0,0.80753549577d0,
     & 0.81059599157d0,0.81074322496d0,0.81092064774d0,0.81113857121d0,
     & 0.81141258579d0,0.81176741858d0,0.81224473178d0,0.81292048690d0,
     & 0.81326531512d0,0.81361810729d0,0.81394892761d0,0.81403478591d0,
     & 0.81453442736d0,0.81514453965d0,0.81569625122d0,0.81590629743d0,
     & 0.81641761054d0,0.81688422790d0,0.81706254273d0,0.81783062223d0,
     & 0.81818548762d0,0.81876084749d0,0.81927932164d0,0.81953752616d0,
     & 0.81991062914d0,0.82000198597d0,0.82046929856d0,0.82136811196d0,
     & 0.82158207086d0,0.82197786731d0,0.82271465654d0,0.82293422050d0,
     & 0.82307403388d0,0.82327592300d0,0.82437272203d0,0.82461223083d0,
     & 0.82518242811d0,0.82588097006d0,0.82593575929d0,0.82657913214d0,
     & 0.82674989909d0,0.82720131507d0,0.82729199217d0,0.82785298542d0,
     & 0.82823976382d0,0.82879713965d0,0.82956576238d0,0.82965109665d0,
     & 0.83022389628d0,0.83024683707d0,0.83026013821d0,0.83056933360d0,
     & 0.83057230602d0,0.83238552115d0,0.83269724980d0,0.83272120040d0,
     & 0.83337244280d0,0.83344262876d0,0.83424980965d0,0.83453543233d0,
     & 0.83528653666d0,0.83559353522d0,0.83571355432d0,0.83584716699d0,
     & 0.83603110733d0,0.83645261202d0,0.83712013990d0,0.83755273628d0,
     & 0.83790801334d0,0.83851082278d0,0.83907614180d0,0.83911697182d0,
     & 0.83944863527d0,0.83992032015d0,0.84028598327d0,0.84049457655d0,
     & 0.84062929625d0,0.84203496327d0,0.84225338094d0,0.84263589538d0,
     & 0.84316462582d0,0.84346407015d0,0.84358826162d0,0.84425298734d0,
     & 0.84440896943d0,0.84510047109d0,0.84544594279d0,0.84614051597d0,
     & 0.84617977210d0,0.84634756465d0,0.84733895218d0,0.84735371621d0,
     & 0.84809948718d0,0.84817198479d0,0.84820658341d0,0.84879139893d0,
     & 0.84914503454d0,0.84936761373d0,0.84968211984d0,0.85079577972d0,
     & 0.85115459340d0,0.85149460662d0,0.85203502193d0,0.85238260416d0,
     & 0.85246057780d0,0.85294545085d0,0.85331989872d0,0.85336336458d0,
     & 0.85367001968d0,0.85396659500d0,0.85493474484d0,0.85542976943d0,
     & 0.85567646584d0,0.85618537433d0,0.85620790802d0,0.85657643376d0,
     & 0.85718252386d0,0.85749923151d0,0.85765999530d0,0.85800965268d0,
     & 0.85900556705d0,0.85925293800d0,0.85997829318d0,0.86016247596d0,
     & 0.86071427280d0,0.86085671118d0,0.86113631159d0,0.86129055424d0,
     & 0.86140398326d0,0.86175402597d0,0.86184648236d0,0.86506336669d0,
     & 0.86563120239d0,0.86581257772d0,0.86589252257d0,0.86593463833d0,
     & 0.86595950321d0,0.86597539487d0,0.86598616285d0,0.86599379407d0,
     & 0.86599939815d0,0.86611125234d0,0.86613081320d0,0.86615790818d0,
     & 0.86619696798d0,0.86625623203d0,0.86635247601d0,0.86652432396d0,
     & 0.86687797809d0,0.86780105383d0,0.87026727702d0,0.87036109429d0,
     & 0.87072401806d0,0.87083929756d0,0.87129097153d0,0.87143601580d0,
     & 0.87174014851d0,0.87201349414d0,0.87220151169d0,0.87296578018d0,
     & 0.87321912503d0,0.87423420066d0,0.87427810075d0,0.87451992265d0,
     & 0.87463780492d0,0.87508656752d0,0.87543545407d0,0.87613646208d0,
     & 0.87620208621d0,0.87657202027d0,0.87675235827d0,0.87746154772d0,
     & 0.87802056981d0,0.87802323532d0,0.87848323721d0,0.87918634348d0,
     & 0.87923353315d0,0.87929475532d0,0.87979232242d0,0.87992980089d0,
     & 0.88023915373d0,0.88071530614d0,0.88140844557d0,0.88152383651d0,
     & 0.88226301283d0,0.88256053579d0,0.88257138016d0,0.88317878511d0,
     & 0.88345376522d0,0.88392610833d0,0.88487101721d0,0.88496397217d0,
     & 0.88504657540d0,0.88508204422d0,0.88587032851d0,0.88596797952d0,
     & 0.88612496216d0,0.88641552700d0,0.88693510431d0,0.88706259977d0,
     & 0.88737022432d0,0.88785167888d0,0.88816475372d0,0.88853423829d0,
     & 0.88914767910d0,0.88931544600d0,0.88976002995d0,0.89006229019d0,
     & 0.89027121803d0,0.89033945984d0,0.89051428261d0,0.89158269202d0,
     & 0.89185573900d0,0.89188264603d0,0.89260246650d0,0.89266571998d0,
     & 0.89303453361d0,0.89329167175d0,0.89392721256d0,0.89426591199d0,
     & 0.89433689053d0,0.89499199788d0,0.89513171174d0,0.89534904191d0,
     & 0.89613084908d0,0.89632115577d0,0.89672171225d0,0.89716711929d0,
     & 0.89718396785d0,0.89775271153d0,0.89818205788d0,0.89851031081d0,
     & 0.89920053309d0,0.89945855804d0,0.89969921820d0,0.89975799541d0,
     & 0.89988390898d0,0.90002498910d0,0.90013494068d0,0.90022257813d0,
     & 0.90029387556d0,0.90172916247d0,0.90182228628d0,0.90194132944d0,
     & 0.90209880697d0,0.90231676774d0,0.90263786198d0,0.90315590361d0,
     & 0.90354940186d0,0.90391316302d0,0.90411725637d0,0.90439611326d0,
     & 0.90479812252d0,0.90506815880d0,0.90527180074d0,0.90587913672d0,
     & 0.90606695144d0,0.90617984594d0,0.90668594476d0,0.90671346638d0,
     & 0.90748209566d0,0.90770567511d0,0.90772630278d0,0.90780967772d0,
     & 0.90848820014d0,0.90854362042d0,0.90948232068d0,0.90958565583d0,
     & 0.90972516020d0,0.90986195218d0,0.91052213708d0,0.91085683307d0,
     & 0.91087999592d0,0.91095972490d0,0.91164967852d0,0.91184993906d0,
     & 0.91223442825d0,0.91232438056d0,0.91259406065d0,0.91285426136d0,
     & 0.91307855666d0,0.91405114045d0,0.91430333969d0,0.91460092856d0,
     & 0.91494790721d0,0.91498277073d0,0.91532900829d0,0.91563302639d0,
     & 0.91592544998d0,0.91637386231d0,0.91707759100d0,0.91711730345d0,
     & 0.91740743879d0,0.91749777452d0,0.91759839922d0,0.91793817351d0,
     & 0.91842587070d0,0.91867525998d0,0.91931014427d0,0.91948612892d0,
     & 0.91953390817d0,0.91994768700d0,0.92007847618d0,0.92009933415d0,
     & 0.92042911624d0,0.92064918535d0,0.92118023295d0,0.92143554682d0,
     & 0.92178143741d0,0.92192829568d0,0.92216393672d0,0.92225921426d0,
     & 0.92242860304d0,0.92249537850d0,0.92262258138d0,0.92267196686d0,
     & 0.92504763564d0,0.92521335987d0,0.92526010000d0,0.92543379881d0,
     & 0.92549645414d0,0.92574133205d0,0.92582968285d0,0.92620004743d0,
     & 0.92633393202d0,0.92695677219d0,0.92718345873d0,0.92736092062d0,
     & 0.92772097097d0,0.92785142472d0,0.92832723948d0,0.92843488366d0,
     & 0.92850269301d0,0.92890152815d0,0.92916067588d0,0.92940914849d0,
     & 0.92956917213d0,0.93006275439d0,0.93035288025d0,0.93037832398d0,
     & 0.93075699790d0,0.93100032698d0,0.93138669071d0,0.93227305661d0,
     & 0.93229298159d0,0.93232516712d0,0.93246951420d0,0.93272696107d0,
     & 0.93281280828d0,0.93297108683d0,0.93352717102d0,0.93400143041d0,
     & 0.93409976547d0,0.93410029476d0,0.93441860024d0,0.93490607594d0,
     & 0.93498213759d0,0.93518547361d0,0.93591820858d0,0.93593449881d0,
     & 0.93597698750d0,0.93644602748d0,0.93665661894d0,0.93694271852d0,
     & 0.93712619035d0,0.93727339240d0,0.93753162614d0,0.93791463594d0,
     & 0.93818296541d0,0.93827455200d0,0.93838119768d0,0.93869437261d0,
     & 0.93892355735d0,0.93906754400d0,0.93916627612d0,0.94033014943d0,
     & 0.94047554934d0,0.94070308940d0,0.94110478095d0,0.94110898668d0,
     & 0.94134385364d0,0.94167195685d0,0.94197629696d0,0.94216239740d0,
     & 0.94236773332d0,0.94288171966d0,0.94296040139d0,0.94297457123d0,
     & 0.94349535346d0,0.94363976494d0,0.94423950912d0,0.94430302748d0,
     & 0.94457502307d0,0.94472613404d0,0.94486917021d0,0.94489927222d0,
     & 0.94514532868d0,0.94534514821d0,0.94553097516d0,0.94614274477d0,
     & 0.94636419965d0,0.94641137486d0,0.94664169100d0,0.94715906666d0,
     & 0.94720428400d0,0.94727738631d0,0.94745886804d0,0.94789305797d0,
     & 0.94827298440d0,0.94828483842d0,0.94889236345d0,0.94889630545d0,
     & 0.94910791234d0,0.94928647956d0,0.94928786293d0,0.94955965094d0,
     & 0.95067552177d0,0.95090055781d0,0.95097234326d0,0.95100396926d0,
     & 0.95102062645d0,0.95114776423d0,0.95118580461d0,0.95125371933d0,
     & 0.95139345140d0,0.95175795571d0,0.95266223579d0,0.95266755752d0,
     & 0.95297943381d0,0.95298770316d0,0.95330984664d0,0.95345210738d0,
     & 0.95346633093d0,0.95410753746d0,0.95423064891d0,0.95423300938d0,
     & 0.95425928063d0,0.95467623745d0,0.95485365867d0,0.95550542261d0,
     & 0.95572225584d0,0.95574822093d0,0.95577521232d0,0.95582394957d0,
     & 0.95628304431d0,0.95661095524d0,0.95682705844d0,0.95714015191d0,
     & 0.95728559578d0,0.95742612431d0,0.95780609308d0,0.95791681921d0,
     & 0.95826784861d0,0.95849117297d0,0.95920911587d0,0.95922536551d0,
     & 0.95925109232d0,0.95926413825d0,0.95977944976d0,0.95983182693d0,
     & 0.95990689173d0,0.95993504527d0,0.96002186497d0,0.96020815213d0,
     & 0.96028985650d0,0.96062327353d0,0.96091315351d0,0.96100879965d0,
     & 0.96130969462d0,0.96139973401d0,0.96175936534d0,0.96192743808d0,
     & 0.96237787477d0,0.96249848799d0,0.96250392509d0,0.96270764579d0,
     & 0.96316799886d0,0.96341885014d0,0.96348661301d0,0.96397192728d0,
     & 0.96398948011d0,0.96403132859d0,0.96434901733d0,0.96476225559d0,
     & 0.96509965042d0,0.96514840245d0,0.96524592650d0,0.96528387708d0,
     & 0.96528590191d0,0.96539345548d0,0.96547423650d0,0.96647608517d0,
     & 0.96654711037d0,0.96660831040d0,0.96671704356d0,0.96682290969d0,
     & 0.96701007649d0,0.96722683857d0,0.96757083027d0,0.96760620250d0,
     & 0.96762428586d0,0.96796625545d0,0.96802139185d0,0.96816023951d0,
     & 0.96861086957d0,0.96868022168d0,0.96870826253d0,0.96914655166d0,
     & 0.96934678733d0,0.96956804627d0,0.96970178877d0,0.96984580729d0,
     & 0.97006049784d0,0.97009809663d0,0.97026290146d0,0.97043761604d0,
     & 0.97059159255d0,0.97067425883d0,0.97131983489d0,0.97148223505d0,
     & 0.97160072337d0,0.97176220090d0,0.97184660317d0,0.97202769105d0,
     & 0.97232148845d0,0.97248403470d0,0.97254247122d0,0.97277258301d0,
     & 0.97286438511d0,0.97313217663d0,0.97327164536d0,0.97332682779d0,
     & 0.97349303006d0,0.97365493582d0,0.97390336802d0,0.97390652852d0,
     & 0.97397741501d0,0.97417377114d0,0.97472855597d0,0.97484632859d0,
     & 0.97488388422d0,0.97503104500d0,0.97510411400d0,0.97529469048d0,
     & 0.97581023371d0,0.97584638775d0,0.97609870933d0,0.97610555741d0,
     & 0.97615928407d0,0.97662248657d0,0.97666392146d0,0.97668632886d0,
     & 0.97678616332d0,0.97714884689d0,0.97725994998d0,0.97736181705d0,
     & 0.97751573550d0,0.97806666283d0,0.97807812450d0,0.97814668888d0,
     & 0.97822865815d0,0.97830170914d0,0.97833867356d0,0.97838544596d0,
     & 0.97861176622d0,0.97873913319d0,0.97895191065d0,0.97904722671d0,
     & 0.97934250806d0,0.97939114343d0,0.97975501469d0,0.97977453241d0,
     & 0.97992347596d0,0.98027822098d0,0.98042757396d0,0.98053235128d0,
     & 0.98054990362d0,0.98074370489d0,0.98106720175d0,0.98115183308d0,
     & 0.98128157130d0,0.98130316537d0,0.98156063425d0,0.98158141499d0,
     & 0.98167601128d0,0.98196871503d0,0.98197275606d0,0.98225594910d0,
     & 0.98254550526d0,0.98254798520d0,0.98257229660d0,0.98262638754d0,
     & 0.98273366980d0,0.98280881059d0,0.98324513530d0,0.98333625388d0,
     & 0.98344048250d0,0.98345100307d0,0.98366812328d0,0.98383143603d0,
     & 0.98412458372d0,0.98415243846d0,0.98418305472d0,0.98426628072d0,
     & 0.98438751751d0,0.98468590967d0,0.98475789591d0,0.98491541971d0,
     & 0.98503185912d0,0.98535408405d0,0.98541701345d0,0.98552715588d0,
     & 0.98561151155d0,0.98574292951d0,0.98589401696d0,0.98591599174d0,
     & 0.98628380870d0,0.98634801055d0,0.98640454264d0,0.98644619565d0,
     & 0.98645572623d0,0.98673055351d0,0.98678044968d0,0.98694703502d0,
     & 0.98702117793d0,0.98719267660d0,0.98722781641d0,0.98742063740d0,
     & 0.98758593077d0,0.98759681924d0,0.98778994493d0,0.98786894120d0,
     & 0.98793576444d0,0.98799251802d0,0.98813501918d0,0.98829371554d0,
     & 0.98858647890d0,0.98863895392d0,0.98869657765d0,0.98872741231d0,
     & 0.98907900825d0,0.98911147001d0,0.98918596321d0,0.98940093499d0,
     & 0.98944236513d0,0.98955512462d0,0.98956096373d0,0.98973945427d0,
     & 0.98978789522d0,0.98997222006d0,0.99011674523d0,0.99025153685d0,
     & 0.99030540262d0,0.99036483369d0,0.99042997119d0,0.99057547531d0,
     & 0.99072623870d0,0.99072854689d0,0.99073484376d0,0.99097298827d0,
     & 0.99101337148d0,0.99108395188d0,0.99116710970d0,0.99141370255d0,
     & 0.99156516842d0,0.99157394284d0,0.99157728834d0,0.99172550029d0,
     & 0.99195955759d0,0.99202062454d0,0.99211684435d0,0.99230024282d0,
     & 0.99231639214d0,0.99240684384d0,0.99256542223d0,0.99260893397d0,
     & 0.99264999845d0,0.99281713971d0,0.99296234891d0,0.99305629094d0,
     & 0.99305635843d0,0.99312859919d0,0.99325521099d0,0.99328369837d0,
     & 0.99346436259d0,0.99350011827d0,0.99353017227d0,0.99370624702d0,
     & 0.99375217062d0,0.99378866194d0,0.99383744364d0,0.99390272670d0,
     & 0.99403196943d0,0.99409015012d0,0.99417947547d0,0.99426126044d0,
     & 0.99429458548d0,0.99447759093d0,0.99449380945d0,0.99468191931d0,
     & 0.99476933500d0,0.99478335680d0,0.99487511702d0,0.99505065612d0,
     & 0.99505797785d0,0.99518722000d0,0.99523122608d0,0.99529792924d0,
     & 0.99539552368d0,0.99552712745d0,0.99555147660d0,0.99555696979d0,
     & 0.99569964038d0,0.99573997005d0,0.99584052512d0,0.99588570115d0,
     & 0.99593797676d0,0.99597459982d0,0.99610229632d0,0.99612249479d0,
     & 0.99617926289d0,0.99622401278d0,0.99629472187d0,0.99634011677d0,
     & 0.99644249757d0,0.99645572570d0,0.99660646054d0,0.99667944226d0,
     & 0.99674778134d0,0.99688045591d0,0.99689348407d0,0.99700517536d0,
     & 0.99708748182d0,0.99712256312d0,0.99723318271d0,0.99726386185d0,
     & 0.99733754456d0,0.99742469425d0,0.99743611187d0,0.99752930576d0,
     & 0.99757175379d0,0.99761750979d0,0.99770107390d0,0.99770656910d0,
     & 0.99778031787d0,0.99783046248d0,0.99785553444d0,0.99792699192d0,
     & 0.99794458248d0,0.99799493667d0,0.99804993054d0,0.99805959521d0,
     & 0.99812117607d0,0.99814738307d0,0.99817987150d0,0.99823585899d0,
     & 0.99823770971d0,0.99832158857d0,0.99839961899d0,0.99847233224d0,
     & 0.99854020064d0,0.99860364518d0,0.99866304213d0,0.99871872858d0,
     & 0.99877100725d0,0.99882015061d0,0.99886640442d0,0.99890999085d0,
     & 0.99895111110d0,0.99898994778d0,0.99902666687d0,0.99906141956d0,
     & 0.99909434380d0,0.99912556563d0,0.99915520041d0,0.99918335391d0,
     & 0.99925985931d0,0.99928298403d0,0.99930504174d0,1.00000000000d0/
      i=1
      j=1
      xmu(i)=x0(j)
      do j=2,nmu0
        x1=x0(j-1)
        x2=x0(j)
        dxmu=x2-x1
        if (dxmu.gt.dxmax) then
          n=dint(dxmu/dxmax)
          h=dxmu/(n+1)
          do k=1,n
            x1=x1+h
            call inc(i,npmax)
            xmu(i)=x1
          enddo
        endif
        call inc(i,npmax)
        xmu(i)=x2
      enddo
      nmu=i
      return
      end
C======================================================================
C      General routines for ENDF-6 formatted files
C======================================================================
      subroutine findnextmat(nin)
c
c     Find next material from cursor position
c
      character*66 line
      mat=10000
      do while (mat.ne.0)
        call readtext(nin,line,mat,mf,mt,ns)
      enddo
      end
C======================================================================
      subroutine findnextmt(nin,mf,mt)
c
c     find next mt on mf file for current material
c
      character*66 line
      mt1=-1
      do while (mt1.ne.0)
        call readtext(nin,line,mat1,mf1,mt1,ns1)
      enddo
      call readtext(nin,line,mat1,mf1,mt,ns1)
      if (mt.ne.0.and.mf1.eq.mf) then
        backspace nin
      else
        mt=-1
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
      subroutine wrtext(lib,mat,mf,mt,ns,line)
c
c     write TEXT record
c
      character*66 line
      ns=ns+1
      if (ns.gt.99999) ns=0
      write(lib,10)line,mat,mf,mt,ns
      return
   10 format(a66,i4,i2,i3,i5)
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
      real*8 function edelta(ffin,fdig)
      implicit real*8 (a-h,o-z)
      character*16 str16
      afdig=dabs(fdig)
      if (afdig.le.0.99999999999d0) then
        edelta=ffin
      else
        if (afdig.gt.9.0d0) fdig=sign(1.0d0,fdig)*9.0d0
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
