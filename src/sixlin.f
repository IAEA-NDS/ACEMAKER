      program sixlin
c     version 2.0
c
c     Process and linearize MF6 data
c
c     The code belongs to ACEMAKER code package to produce continuous
c     energy ACE-formatted file for Monte Carlo calculations
c
c     SIXLIN converts  laboratory angle-energy distribution (LAW=7) 
c     to continuum energy-angle distribution given by linearly 
c     interpolable tabulated data (LAW=1/LANG=12).
c     Additionally, SIXLIN linearizes secundary energy distributions
c     and yields given in MF6.
c
c     INPUT data:
c
c     Input data options should be entered on the SIXLIN.INP text file.
c
c     line 1:       isel       imon                (2i11)
c     line 2:       input  ENDF filename           (a72)
c     line 3:       output ENDF filename           (a72)
c     line 4:       ilow       iupp                (2i11)
c     line 5:        tol       ymin                (2e11.0)
c
c     where
c       isel: selection criterium      (0/1) = MAT/ZA   (Default = 0)
c       imon: Monitor printing trigger (0/1) = min/max  (Default = 0)
c       ilow: Lower material number requested
c       iupp: Upper material number requested
c            (Default ilow=iupp=0,first material on the input ENDF tape)
c        tol: Fraction tolerance for linearization/reconstruction
c             (Default = 0.01)
c       ymin: Minimum value alowable for linearization/reconstruction
c             of angle-energy distributions (1.0d-30)
c
c     Example of SIXLIN.INP:
c
c               0          0
c     \PENDF\U235.PENDF
c     \LIN\U235.SIX
c            9228       9228
c         1.0D-3     1.0D-20
c
c     Retrieve material 9228 (U-235) from tape \PENDF\U235.PENDF
c     and process/linearize MF6 file. Write results on \LIN\U235.SIX
c     tape. The minimum allowable value for angle-energy reconstruction
c     will be 1.e-20. Use a linearization tolerance of 0.1% (0.001).
c
c     Output files:
c      Output ENDF formatted file
c      SIXLIN.LST listing file
c
      implicit real*8(a-h, o-z)
      parameter (npmax=2000000, nbmax=2000000, npmumax=5001)
      parameter (etol0=0.01d0,ymin0=1.0d-30)
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
      read(in1,'(2i11)')isel,imon
      if (isel.ne.1) isel=0
      if (imon.ne.1) imon=0
      read(in1,'(a72)')fin1
      read(in1,'(a72)')fout
      read(in1,'(2i11)')llower,lupper
      if (llower.lt.0) llower=0
      if (lupper.lt.llower) lupper=llower
      read(in1,'(2e11.0)')tolmax,ymin
      etol=tolmax
      if (etol.le.0.0d0) etol=etol0
      if (ymin.le.0.0d0) ymin=ymin0
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
      write(iou,*)' Input file name  =',fin1
      write(iou,*)' Output file name =',fout
      if (llower.eq.0.and.lupper.eq.0) then
        write(iou,*)' Selection range  = first material'
      else
        write(iou,*)' Selection range  =',llower, ' to ',lupper
      endif
      write(iou,'(1x,a,1pe11.5)')' Tolerance [%]    =',etol*100.0d0
      write(iou,'(1x,a,1pe11.5)')' Minimum Y-value  =',ymin
      write(iou,*)
c
c     read tape header
c
      ns=-1
      call readtext(nin,line,mat,mf,mt,nsi)
      call wrtext(nou,mat,mf,mt,ns,line)
      call readcont(nin,za,awr,l1,l2,n1,n2,mat,mf,mt,nsi)
      if (llower.eq.0) then
        if (isel.eq.0) then
          llower=mat
        else
          llower=int(za)
        endif
      endif
      if (lupper.lt.llower) lupper=llower
      call checksel(isel,mat,za,llower,lupper,nsel)
      write(iou,*)' nsel=', nsel
      do while (nsel.ne.-1)
       if (nsel.gt.0) then
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
         line(34:66)=' (VERSION 2021-1) ***************'
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
         mf0=1
         mfn=5
         call cpmfs(nin,nou,nlst,mf0,mfn)
         write(*,*)' MF1 to MF5 copied to ouput tape'
         write(iou,*)' MF1 to MF5 copied to ouput tape'         
c
c          Process MF6 data, if any
c
         call readtext(nin,line,mat,mf,mt,nsi)
         backspace nin        
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
c              for law 1 & 7 data is processed
c
               if (law.eq.1) then
c
c                law=1
c
                 call wrtab1(nou,mat,mf,mt,ns,
     &             zap,awp,lip,law,nr,nbt,ibt,np,x,y)
                 call readtab2(nin,c1,c2,lang,lep,nr,ne,nbt,ibt)
                 lang2=lang
                 if (lep.gt.2) then
                   lep2=2
                   write(iou,*)' LANG=',lang,' LEP=',lep,' NE=',ne,
     &             ' Secundary energies will be linearized'
                 else
                   lep2=lep
                   write(iou,*)' LANG=',lang,' LEP=',lep,' NE=',ne
                 endif
                 call wrtab2(nou,mat,mf,mt,ns,c1,c2,lang,lep2,
     &             nr,ne,nbt,ibt)
                 nr=1
                 ibt(1)=lep
                 do ie=1,ne
                   call readlist(nin,c1,e,nd,na,nw,nep,b)
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
c                    saving original secondary energy grid x0(i)
c                    initializing common energy grid x(i)
c
                     j=i0+1
                     do i=1,nepc0
                       ep=b0(j)
                       x0(i)=ep
                       x(i)=ep
                       j=j+na2
                     enddo
c
c                    linearize law parameters on E' grid
c
                     do k=2,na2
                       j=i0+k
                       do i=1,nepc0
                         y(i)=b0(j)
                         j=j+na2
                       enddo
                       if (k.eq.2) then
                         f0=0.0d0
                         call renorm(nepc0,x0,y,lep,f0,fn)
                       endif
                       nbt(1)=nepc0
                       do i=1,nepc
                         b(i)=fvalue(nr,nbt,ibt,nepc0,x0,y,x(i))
                       enddo
                       nbt(1)=nepc
                       call linear(nr,nbt,ibt,nepc,x,b,etol,ymin,npmax)
                     enddo
c
c                    Prepare linearized list record
c                     1. Calculate new nep and nw
c                     2. Process discrete data
c                     3. Set new secondary energy grid for lin-lin
c                     4. Update linearly interpolable angular parameters
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
                     nbt(1)=nepc0
                     do k=2,na2
                       j=i0+k
                       do i=1,nepc0
                         y(i)=b0(j)
                         j=j+na2
                       enddo
                       j=i0+k
                       do i=1,nepc
                         b(j)=fvalue(nr,nbt,ibt,nepc0,x0,y,x(i))
                         j=j+na2
                       enddo
                     enddo
c
c                    Renormalize distributions if required
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
                       y(i)=b(j+1)
                       if (y(i).lt.0.0d0) y(i)=0.0d0
                       j=j+na2
                     enddo
                     call renorm(nepc,x,y,lep2,f0,fn)
                     j=i0+2
                     do i=1,nepc
                       b(j)=y(i)
                       j=j+na2
                     enddo                     
                     if (na.gt.0.and.lang.ne.2) then             
                       if (lang.eq.1) then
                         do i=1,nepc
                           j=i0+na2*(i-1)+2
                           do k=1,na
                             b(j+k)=b(j+k)*fn
                           enddo
                         enddo
                       else
                         nmu=na/2
                         if ((na-2*nmu).ne.0) then
                           write(iou,*)' ERROR: NA is not even'
                           stop
                         endif
                         c=0.5d0                         
                         do i=1,nepc
                           j=i0+na2*(i-1)+2
                           if (b(j).gt.0.0d0) then
                             do k=1,nmu
                               k2=2*k
                               x(k)=b(j+k2-1)
                               y(k)=b(j+k2)
                               if (y(k).lt.0.0d0) y(k)=0.0d0
                             enddo                           
                             call renorm(nmu,x,y,lang,c,fn)
                             do k=1,nmu
                               k2=2*k
                               b(j+k2)=y(k)
                             enddo
                           endif
                         enddo
                       endif
                     endif                     
                     if (imon.gt.0) then
                      write(iou,*)' ==========================='
                      write(iou,*)' Incident energy E         = ',e
                      write(iou,*)' Outgoing energies (NEP)   = ',nep
                      write(iou,*)'  - Discrete data (ND)     = ',nd
                      write(iou,*)'  - Continuous data (NEPC) = ',nepc
                      write(iou,*)' Angular parameters (NA)   = ',na
                      write(iou,*)' Integral f0(E)            = ',f0+sum
                      write(iou,*)'  - Discrete contribution  = ',sum
                      write(iou,*)'  - Continuous contribution= ',f0
                      write(iou,*)
                     endif
                   endif
                   call wrtlist(nou,mat,mf,mt,ns,c1,e,nd,na,nw,nep,b)
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
     &             c1e,c2e,lang,lep,nre,ne,nbte,ibte)
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
                   if (lep.eq.2.and.icod.ne.0) then
                     call linear(nr,nbt,ibt,nep,x,y,etol,ymin,npmax)
                   endif
                   k=0
                   nepmu(1)=nep
                   do i=1,nep
                     call inc(k,nbmax)
                     b0(k)=x(i)
                     x0(i)=x(i)
                     call inc(k,nbmax)
                     b0(k)=y(i)
                   enddo
                   nep0=nep
                   do j=2,nmu0
                     call readtab1(jou,c1,xmu,l1,l2,nr,nep,nbt,ibt,x,y)
                     y0(j)=xmu
                     x1(j)=xmu
                     call checklaw(nr,ibt,icod)
                     if (lep.eq.2.and.icod.ne.0) then
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
                       if (yy.lt.0.0d0) yy=0.0d0
                       call inc(l,nbmax)
                       b(l)=yy
                     enddo
                   enddo
c
c                  linearize distribution for cosines on union grid
c
                   nmu1=nmu0
                   call checklaw(nrm,ibtm,icod)
                   if (lang.eq.12.and.icod.ne.0) then
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
                   write(iou,*)' =================================='
                   write(iou,*)' Convertion of Law 7 to LAW 1 at E=',e
                   write(iou,*)' Common outgoing energy points nep=',
     &               nep0
                   write(iou,*)' Common outgoing cosine points nmu=', 
     &               nmu1
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
                     if (imon.gt.0) then
                       write(iou,*)' E''=',x0(i),' f0E,E'')=',c
                     endif
                     do j=1,nmu1
                       yy=fvalue(nrm,nbtm,ibtm,nmu0,y0,y,x1(j))
                       if (yy.lt.0.0d0) yy=0.0d0
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
                       write(iou,*)' Applied renormalization factor=',fn
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
                         b(k)=0.25d0
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
c                For MF6/LAW=0,2,3,4,5 & 6 data is copied
c
                 write(iou,*)' zap=',zap,' awp=',awp,' lip=',lip,
     &                       ' law=',law, ' data copied'
c
c                  MF6/LAW=0,3,4
c     
                 call wrtab1(nou,mat,mf,mt,ns,
     &                       zap,awp,lip,law,nr,nbt,ibt,np,x,y)                
                 if (law.eq.2.or.law.eq.5) then
c
c                  MF6/LAW=2,5
c                 
                   call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
                   call wrtab2(nou,mat,mf,mt,ns,
     &                         c1,c2,l1,l2,nr,ne,nbt,ibt)
                   do i=1,ne
                     call readlist(nin,c1,e,l1,l2,nw,nl,b)
                     call wrtlist(nou,mat,mf,mt,ns,c1,e,l1,l2,nw,nl,b)
                   enddo
                 elseif (law.eq.6) then
c
c                  MF6/LAW=6
c                 
                   call readcont(nin,ap,c2,l1,l2,n1,np,mat,mf,mt,nsi)
                   call wrtcont(nou,mat,mf,mt,ns,ap,c2,l1,l2,n1,np)
                 elseif (law.gt.7) then
c
c                  unknown
c                 
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
         write(*,*)' MF7 to MF28 copied to output tape'
         write(iou,*)' MF7 to MF28 copied to output tape'
         mf0=7
         mfn=28
         call cpmfs(nin,nou,nlst,mf0,mfn)
         call wrtmend(nou,ns)
         ns=0
       endif
       call findnextmat(nin)
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
        if (nsel.lt.matlow) then
          nsel=0
        elseif (nsel.gt.matup) then
          nsel=-1
        endif
      endif
      return
      end
C======================================================================
      subroutine cpmfs(nin,nou,nlst,mf0,mfn)
c
c       Copy MF files from MF=mf0 to MF=mfn for the current MAT
c      
      character*66 line          
      call readtext(nin,line,mat,mf,mt,nsi)
      if (mat.gt.0) then
        write(nlst,'(a,i2,a,i2,a,i5)')' MF=',mf0,' to MF=',mfn,
     &    ' copied for MAT=',mat      
        mat1=mat
        do while(((mf.ge.mf0.and.mf.le.mfn).or.mf.eq.0).and.mat1.eq.mat)
          if (mt.eq.0) then
            if (mf.eq.0) then
              call wrtfend(nou,mat,ns)
            else
              call wrtsend(nou,mat,mf,ns)
            endif             
          else
            call wrtext(nou,mat,mf,mt,ns,line)
          endif
          call readtext(nin,line,mat,mf,mt,nsi)
        enddo
      endif
      backspace(nin)
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
C      General routines for ENDF-6 formatted files
C======================================================================
      subroutine findnextmat(nin)
c
c     Find next material from cursor position
c
      character*66 line
      mat=10000
      do while (mat.ne.0.and.mat.ne.-1)
        read(nin,'(a66,i4,i2,i3,i5)')line,mat,mf,mt,ns
      enddo
      if (mat.eq.-1) backspace nin
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
C      The following subroutine was taken from NJOY2016 and
C      adapted/modified by D. Lopez Aldama for ACEMAKER:
C       1. subroutine terp1 (renamed as terp1m)
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
