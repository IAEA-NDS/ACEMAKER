      program gamlin
c     version: acemaker 3.0
c     
c     Process and linearize photon files (MF12,MF13 and MF14)
c
c     The code belongs to ACEMAKER code package to produce continuous
c     energy ACE-formatted file for Monte Carlo calculations
c
c     GAMLIN process photon production yield data (MF12),
c     photon production cross sections (MF13) and photon
c     angular distributions (MF14).
c
c      MF12: - Photon transition probability data (lo=2) are converted 
c              to multiplicities
c            - Multiplicities (lo=1) are linearized, if required
c      
c      MF13: - Photon production cross sections are linearized, if 
c              required
c
c      MF14: - Legendre coefficients representation is converted to
c              linearly interpolable tabular data.
c            - Tabular data are lineraized, if required
c            - MTs converted from transition probabilities to 
c              multiplicities in MF12 are checked and re-written
c              (as isotropic)
c
c     INPUT data:
c
c     Input data options should be entered on the GAMLIN.INP text file.
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
c       ilow: Lower mat/za number requested
c       iupp: Upper mat/za number requested
c            (Default ilow=iupp=0,first material on the input ENDF tape)
c        tol: Fraction tolerance for linearization/reconstruction
c             (Default = 0.01)
c       ymin: Minimum value alowable for linearization/reconstruction
c             of cross sections and angle distributions (1.0d-30)
c
c     Example of GAMLIN.INP:
c
c               0          0    
c     \PENDF\U235.PENDF
c     \LIN\U235.GAM
c            9228       9228
c         1.0e-3     1.0E-20
c
c     Retrieve material 9228 (U-235) from tape \PENDF\U235.PENDF
c     and process/linearize photon files. Write results on \LIN\U235.GAM
c     tape. The minimum allowable value for linearization/reconstruction
c     will be 1.E-20. Use a linearizaation tolerance of 0.1% (0.001).
c
c     Output files:
c      Output ENDF formatted file
c      GAMLIN.LST listing file
c
      implicit real*8 (a-h, o-z)
      parameter (nng2=320, nnxc=450, npmax=2000000)
      parameter (ytol0=0.001d0,ymin0=1.0d-30)
      dimension ngam2(nng2),mtt(nnxc),qi(nnxc),ethr(nnxc)
      character*66 line
      character*72 fin1,fout
      data ninp/1/,nin/20/,nou/30/,nlst/40/      
c
c      read input data from GAMLIN.INP
c
      open (ninp, file='GAMLIN.INP')
      isel=0
      imon=0
      llower=0
      lupper=0
      read(ninp,'(2i11)')isel,imon
      if (isel.ne.1) isel=0
      if (imon.ne.1) imon=0
      read(ninp,'(a72)')fin1
      read(ninp,'(a72)')fout
      read(ninp,'(2i11)')llower,lupper
      if (llower.lt.0) llower=0
      if (lupper.lt.llower) lupper=llower
      read(ninp,'(2e11.0)')tolmax,ymin
      ytol=tolmax
      if (ytol.le.0.0d0) ytol=ytol0
      if (ymin.le.0.0d0) ymin=ymin0
      close(ninp)
c
c     printing input data
c
      open (nlst, file='GAMLIN.LST')
      write(nlst,*)' PROGRAM GAMLIN: Process and linearize photon files'
      write(nlst,*)' =================================================='
      write(nlst,*)      
      write(nlst,*)' Input parameters'
      write(nlst,*)' ================'
      write(nlst,*)' MAT/ZA selection =',isel
      write(nlst,*)' Printing option  =',imon
      write(nlst,*)' Input file name  =',fin1
      write(nlst,*)' Output file name =',fout
      if (llower.eq.0.and.lupper.eq.0) then
        write(nlst,*)' Selection range  = first material'
      else
        write(nlst,*)' Selection range  =',llower, ' to ',lupper
      endif
      write(nlst,'(1x,a,1pe11.5)')' Tolerance [%]    =',ytol*100.0d0
      write(nlst,'(1x,a,1pe11.5)')' Minimum Y-value  =',ymin
      write(nlst,*)
c
c       Process input ENDF tape
c       
      open (nin, file=fin1)
      open (nou, file=fout)
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
      do while (nsel.ne.-1)
        if (nsel.gt.0) then
c
c         Process material
c                       
c         mf=1/mt=451
c
          call wrtcont(nou,mat,mf,mt,ns,za,awr,l1,l2,n1,n2)
          call readcont(nin,elis,c2,l1,l2,n1,n2,mat,mf,mt,nsi)
          call wrtcont(nou,mat,mf,mt,ns,elis,c2,l1,l2,n1,n2)
          call readcont(nin,awi,elim,lrel,l2,nsub,nver,mat,mf,mt,nsi)
          call wrtcont(nou,mat,mf,mt,ns,awi,elim,lrel,l2,nsub,nver)
          call readcont(nin,temp,c2,ldrv,l2,nwd,nxc,mat,mf,mt,nsi)
          nwd1=nwd+3
          call wrtcont(nou,mat,mf,mt,ns,temp,c2,ldrv,l2,nwd1,nxc)
          do i=1,nwd
            call readtext(nin,line,mat,mf,mt,nsi)
            call wrtext(nou,mat,mf,mt,ns,line)
          enddo
          line( 1:33)=' ***************** Program GAMLIN'
          line(34:66)=' (VERSION 2021-1) ***************'
          call wrtext(nou,mat,mf,mt,ns,line)
          write(line,'(a26,1pe11.4,a24)')' For All Data Greater than',
     &      ymin,' in Absolute Value'
          call wrtext(nou,mat,mf,mt,ns,line)
          write(line,'(a42,f10.7,a9)')
     &      ' Data Linearized to Within an Accuracy of ',
     &      ytol*100.d0, ' per-cent'
          call wrtext(nou,mat,mf,mt,ns,line)
          do i=1,nxc
            call readtext(nin,line,mat,mf,mt,nsi)
            call wrtext(nou,mat,mf,mt,ns,line)
          enddo
          call readtext(nin,line,mat,mf,mt,nsi)
          call wrtsend(nou,mat,mf,ns)
c
c         Copy MF1 & MF2 to output file
c
          mf0=1
          mfn=2
          call cpmfs(nin,nou,nlst,mf0,mfn)
          call findmf(nin,mat,3,icod)
          if (icod.eq.0) then
            call getmf3(nin,nou,nlst,elim,nmtt,mtt,qi,ethr,npmax)
            write(nlst,'(a,i4,a)')' ',nmtt,' sections found on MF3'
          else
            write(nlst,'(a)')' MF=3 not found'
            stop           
          endif
          mf0=4
          mfn=10
          call cpmfs(nin,nou,nlst,mf0,mfn)          
          call findmf(nin,mat,12,icod)
          if (icod.eq.0) then
            call mf12(nin,nou,nlst,elis,elim,awi,nmtt,mtt,qi,ethr,
     &                nng2,ngam2,ymin,ytol)
            write(nlst,'(a)')' MF=12 processed'
          else
            write(nlst,'(a)')' MF=12 is not used'    
          endif
          call findmf(nin,mat,13,icod)
          if (icod.eq.0) then
            call mf13(nin,nou,nlst,ymin,ytol)
            write(nlst,'(a)')' MF=13 processed'
          else
            write(nlst,'(a)')' MF=13 is not used'            
          endif
          call findmf(nin,mat,14,icod)
          if (icod.eq.0) then
            call mf14(nin,nou,nlst,nng2,ngam2,ymin,ytol)
            write(nlst,'(a)')' MF=14 processed'
          else
            write(nlst,'(a)')' MF=14 is not used' 
          endif
          mf0=15
          mfn=28
          call cpmfs(nin,nou,nlst,mf0,mfn)
          call wrtmend(nou,ns)
        endif
        call findnextmat(nin)       
        call readcont(nin,za,awr,l1,l2,n1,n2,mat,mf,mt,nsi)
        call checksel(isel,mat,za,llower,lupper,nsel)
      enddo
      call wrtend(nou,ns)
      close (nlst)
      close (nin)
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
        ns=0
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
      backspace nin
      return
      end                
C======================================================================
      subroutine getmf3(nin,nou,nlst,elim,nmtt,mtt,qi,ethr,npmax)
c
c       Read MF3 from tape nin and copy the data to tape nou
c       Store MT, QI and Threshold values in the arrays mtt, qi 
c       and ethr respectively
c       Check if cross section data if linearly interpolable and
c       extend cross section data up to elim as null if required
c      
      implicit real*8 (a-h, o-z)
      dimension mtt(*),qi(*),ethr(*)
      dimension x(npmax),y(npmax)
      dimension ibt(20),nbt(20)
      nmtt=0
      call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns1)
      do while (mf.eq.3)
        ns=0        
        call wrtcont(nou,mat,mf,mt,ns,c1,c2,l1,l2,n1,n2)
        write(*,*)' point mt=',mt
        call readtab1(nin,qm,qii,l1,lr,nr,ne1,nbt,ibt,x,y)
        call checklaw(nr,ibt,icod)
        if (icod.ne.0) then
          write(nlst,*)' non linearly interpolable data for MF3/MT=',mt
          stop
        endif
        nmtt=nmtt+1
        mtt(nmtt)=mt
        qi(nmtt)=qii
        ethri=x(ne1-1)        
        do ie=2,ne1
          if (y(ie).gt.0.0d0) then
            ethri=x(ie-1)
            exit
          endif
        enddo
        if (ethri.lt.1.0d-5) ethri=1.0d-5
        ethr(nmtt)=ethri
        if (x(ne1).lt.elim) then
          ie=ne1+1
          x(ie)=x(ne1)
          y(ie)=0.0d0
          ie=ie+1
          x(ie)=elim
          y(ie)=0.0d0
          ne1=ie
          nbt(nr)=ne1            
        endif
        call wrtab1(nou,mat,mf,mt,ns,qm,qii,l1,lr,nr,nbt,ibt,ne1,x,y)       
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns1)
        call wrtsend(nou,mat,mf,ns)
        call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,ns1)
      enddo
      call wrtfend(nou,mat,ns)     
      return
      end
C======================================================================      
      subroutine mf12(nin,nou,nlst,elis,elim,awi,nmtt,mtt,qi,ethr,
     &  nng2,ngam2,ymin,ytol)
c
c       Convert transition probability tables to yield data
c       Linearize yield data, if required
c
      implicit real*8 (a-h, o-z)
      character*66 line
      parameter (mtmax=50,lymax=5500,nbmax=2000000)
      dimension mtt(*),qi(*),ethr(*)
      dimension tp(mtmax,mtmax),tt(mtmax,mtmax),ee(mtmax)
      dimension eg(lymax),es(lymax),y(lymax),x(2),yld(2)
      dimension b(nbmax),ibt(20),nbt(20)
      dimension ngam2(*)
c
c     Initialize
c     
      do i=1,nng2
        ngam2(i)=0.0d0
      enddo
      mtlast=0
      mt0old=0
      mt0=0
      call readcont(nin,za,awr,l0,lg,nk,n2,mat,mf,mt,nsi)
      do while (mf.eq.12)
        ns=0
        if (mt.ne.460) then
          if (mt.eq.51.or.mt.eq.601.or.mt.eq.651.or.mt.eq.701.or.
     &        mt.eq.751.or.mt.eq.801.or.mt.eq.876) then
            mtlast=mt
          elseif ((mt.gt.51.and.mt.lt.91).or.
     &            (mt.gt.601.and.mt.lt.649).or.
     &            (mt.gt.651.and.mt.lt.699).or.
     &            (mt.gt.701.and.mt.lt.749).or.
     &            (mt.gt.751.and.mt.lt.799).or.
     &            (mt.gt.801.and.mt.lt.849).or.
     &            (mt.gt.876.and.mt.lt.891)) then
            if ((mt-mtlast).ne.1) then
              write(nlst,'(a,i2,a,a)')' MF12/MT=',mt-1,' may be',
     &          ' missing, discrete photon data may be incomplete'
            endif
            mtlast=mt
          endif
          if (l0.eq.2) then
c
c           Convert transition probabilities to multiplicities
c          
            if (mt.ge.51.and.mt.lt.91.and.mt0.ne.49) then
              mt0=49
            elseif (mt.ge.601.and.mt.lt.649.and.mt0.ne.599) then
              mt0=599
            elseif (mt.ge.651.and.mt.lt.699.and.mt0.ne.649) then
              mt0=649
            elseif (mt.ge.701.and.mt.lt.749.and.mt0.ne.699) then
              mt0=699
            elseif (mt.ge.751.and.mt.lt.799.and.mt0.ne.749) then
              mt0=749
            elseif (mt.ge.801.and.mt.lt.849.and.mt0.ne.799) then
              mt0=799
            elseif (mt.ge.876.and.mt.lt.891.and.mt0.ne.874) then
              mt0=874                                   
            endif        
            if (mt0.ne.mt0old) then
              mt0old=mt0
              mt1=mt0+2
              if (mt0.eq.49) then
                mtc=91
              elseif (mt0.eq.874) then
                mtc=891
              else
                mtc=mt0+50
              endif
              do i=1,mtmax
                ee(i)=0.0d0
                do j=1,mtmax
                  if (i.eq.j) then
                    tp(i,j)=1.0d0
                  else
                    tp(i,j)=0.0d0
                  endif
                  tt(i,j)=0.0d0
                enddo
              enddo              
              do j=mt1,mtc-1
                i=iposm(nmtt,mtt,j)
                if (i.eq.0) then
                  ee(j-mt0)=0.0d0
                else
                  ee(j-mt0)=elis-qi(i)
                endif
              enddo
            endif
            call readlist(nin,c1,c2,lp,l2,nw,nt,b)
            j=mt-mt0
            ee(j)=c1
            lg1=lg+1
            jm1=j-1
            do kk=1,jm1
              k=jm1-kk+1
              eek=ee(k)
              ii=0
              do i=1,nt
                esi=b(lg1*i-lg)
                if ((esi.eq.0.0d0.and.eek.eq.0.0d0).or.
     &          ((esi.ne.0.0d0).and.(abs(esi-eek).le.0.0001d0*esi)))then
                  ii=i
                  exit
                endif
              enddo
              if (ii.eq.0) then
                tp(k,j)=0.0d0
                tt(k,j)=0.0d0
              else
                p=b(lg1*ii-lg+1)
                tp(k,j)=p
                if (lg.eq.2) then
                  tt(k,j)=p*b(lg1*ii)
                else
                  tt(k,j)=p
                endif
              endif
              if (k.ne.jm1) then
                kp1=k+1
                do i=kp1,jm1
                  tp(k,j)=tp(k,j)+tt(k,i)*tp(i,j)
                enddo
              endif
            enddo
            l=0
            ysum=0.0d0
            do i=2,j
              j1=j+2-i
              ej1=ee(j1)
              do ii=1,jm1
                j2=j-ii
                ej2=ee(j2)
                yy=tt(j2,j1)*tp(j1,j)
                if (yy.ne.0.0d0) then
                  l=l+1
                  if (l.gt.lymax) then
                    write(nlst,*)' Error: too many photons from',
     &                ' transition probabilities, increase lymax=',lymax
                    stop
                  endif
                  eg(l)=ej1-ej2
                  es(l)=ej1
                  y(l)=yy
                  ysum=ysum+yy
                endif
              enddo
            enddo
            ee(1)=0.0d0
c
c           re-ordering data in descending order of gamma energies
c            
            if (l.gt.1) then
              lm1=l-1
              do i=1,lm1
                ip1=i+1
                do ii=ip1,l
                  if (eg(i).lt.eg(ii)) then
                    temp=eg(i)
                    eg(i)=eg(ii)
                    eg(ii)=temp
                    temp=y(i)
                    y(i)=y(ii)
                    y(ii)=temp
                    temp=es(i)
                    es(i)=es(ii)
                    es(ii)=temp
                  endif 
                enddo
              enddo
            endif
            mtref=49
            if (mt.ge.600.and.mt.le.849) then
              mtref=549
            elseif (mt.ge.875.and.mt.le.891) then
              mtref=574
            endif
            ngam2(mt-mtref)=l            
c
c           prepare yields in endf format for lo=1 (multiplicities)
c            
            c1=za
            c2=awr
            l1=1
            l2=0
            n1=l
            n2=0
            i=iposm(nmtt,mtt,mt)
            if (i.gt.0) then
              elow=max(ethr(i),1.0d-5)
            else
              elow=1.0d-5
            endif
            call wrtcont(nou,mat,12,mt,ns,c1,c2,l1,l2,n1,n2)
            if (l.gt.1) then
c
c             write the sum tab1 if more than 1 photon  
c            
              c1=0.0d0
              c2=0.0d0
              l1=0
              l2=0
              n1=1
              n2=2
              ibt(1)=2
              nbt(1)=2
              x(1)=elow
              yld(1)=ysum
              x(2)=elim
              yld(2)=ysum
              call wrtab1(nou,mat,12,mt,ns,c1,c2,l1,l2,n1,nbt,ibt,
     &          n2,x,yld)
            endif
c
c           Rest of tab1 records
c            
            do i=1,l
              c1=eg(i)
              c2=es(i)
              l1=0
              l2=2
              n1=1
              n2=2
              ibt(1)=2
              nbt(1)=2
              x(1)=elow
              yld(1)=y(i)
              x(2)=elim
              yld(2)=y(i)
              call wrtab1(nou,mat,12,mt,ns,c1,c2,l1,l2,n1,nbt,ibt,
     &          n2,x,yld)            
            enddo
            read(nin,*)
            call wrtsend(nou,mat,12,ns)
          else
c
c           Multiplicities (linearize if required)
c          
            call sectlin(nin,nou,za,awr,l0,nk,mat,mf,mt,ymin,ytol)
          endif     
        else
c
c         copy section MT460
c        
          call wrtcont(nou,mat,mf,mt,ns,za,awr,l0,lg,nk,n2)
          call readtext(nin,line,mat,mf,mt,nsi)
          do while (mt.ne.0)
            call wrtext(nou,mat,mf,mt,ns,line)
            call readtext(nin,line,mat,mf,mt,nsi)
          enddo
          call wrtsend(nou,mat,mf,ns)
        endif
        call readcont(nin,za,awr,l0,lg,nk,n2,mat,mf,mt,nsi)
      enddo
      call wrtfend(nou,mat,ns)      
      return
      end
C======================================================================
      subroutine sectlin(nin,nou,za,awr,l0,nk,mat,mf,mt,ymin,ytol)
      implicit real*8 (a-h, o-z)
      parameter (npmax=2000000)
      dimension x(npmax),y(npmax),xt(npmax),yt(npmax)
      dimension nbt(20),ibt(20)
      data nscr/90/,nerr/0/
      ns=0
      l2=0
      n2=0
      call wrtcont(nou,mat,mf,mt,ns,za,awr,l0,l2,nk,n2)
      call readtab1(nin,c1,c2,l1,l2,nr,np,nbt,ibt,x,y)      
      if (nk.gt.1) then
        open (nscr,file='sectlin.tmp')
        nsc=0
        nty=2
        yt(1)=x(1)
        yt(2)=x(np)
        nk1=nk
        do k=1,nk1
          call readtab1(nin,c1,c2,l1,l2,nr,np,nbt,ibt,x,y)
          call linear(nr,nbt,ibt,np,x,y,ytol,ymin,npmax)
          call checkdis(nerr,np,x,y,icod)
          nr=1
          nbt(1)=np
          ibt(1)=2
          call wrtab1(nscr,mat,mf,mt,nsc,c1,c2,l1,l2,nr,nbt,ibt,np,x,y)
          call union(nty,yt,np,x,ntx,xt,npmax)
          if (k.lt.nk1) then
            nty=ntx
            do kk=1,ntx
              yt(kk)=xt(kk)
            enddo
          endif          
        enddo
        rewind(nscr)
        do i=1,np
          yt(i)=0.0d0
        enddo
        do k=1,nk1
          call readtab1(nscr,c1,c2,l1,l2,nr,np,nbt,ibt,x,y)
          do i=1,ntx
            yt(i)=fvalin(np,x,y,xt(i))+yt(i)
          enddo
        enddo
        rewind(nscr)
        call checkdis(nerr,ntx,xt,yt,icod)
        nr=1
        nbt(1)=ntx
        ibt(1)=2
        c1=0.0d0
        c2=0.0d0
        l1=0
        l2=0
        call wrtab1(nou,mat,mf,mt,ns,c1,c2,l1,l2,nr,nbt,ibt,ntx,xt,yt)
        do k=1,nk1
          call readtab1(nscr,c1,c2,l1,l2,nr,np,nbt,ibt,x,y)
          call wrtab1(nou,mat,mf,mt,ns,c1,c2,l1,l2,nr,nbt,ibt,np,x,y)
        enddo
        close(nscr, status='DELETE')
      else
        call linear(nr,nbt,ibt,np,x,y,ytol,ymin,npmax)
        call checkdis(nerr,np,x,y,icod)
        nr=1
        nbt(1)=np
        ibt(1)=2
        call wrtab1(nou,mat,mf,mt,ns,c1,c2,l1,l2,nr,nbt,ibt,np,x,y)
      endif
      read(nin,*)
      call wrtsend(nou,mat,mf,ns)      
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
      subroutine mf13(nin,nou,nlst,ymin,ytol)
      implicit real*8 (a-h, o-z)
      call readcont(nin,za,awr,l0,l1,nk,n2,mat,mf,mt,nsi) 
      do while (mf.ne.13)
        call sectlin(nin,nou,za,awr,l0,nk,mat,mf,mt,ymin,ytol)
        call readcont(nin,za,awr,l0,l1,nk,n2,mat,mf,mt,nsi)
      enddo
      call wrtfend(nou,mat,ns)
      return
      end
C======================================================================      
      subroutine mf14(nin,nou,nlst,nng2,ngam2,ymin,ytol)
c
c      process/linearize MF14 (angular distributions for gammas) 
c      
      implicit real*8 (a-h, o-z)
      parameter (nkmax=1000,npmax=2000000,nbmax=2000000)
      dimension ngam2(*)
      dimension x0(npmax),x1(npmax)
      dimension x(npmax),y(npmax)
      dimension nbt(20),ibt(20)
      dimension b(nbmax)
      data nscr1/91/,nscr2/92/
      call readcont(nin,za,awr,li,ltt,nk,ni,mat,mf,mt,nsi)
      do while (mf.eq.14)
        if (mt.ge.50.and.mt.le.91) then
          nkk=ngam2(mt-49)
        elseif (mt.ge.600.and.mt.le.849) then
          nkk=ngam2(mt-549)
        elseif (mt.ge.875.and.mt.le.891) then
          nkk=ngam2(mt-574)
        else
          nkk=0
        endif
        if (nkk.gt.0) then
c
c         Transition probabilities were converted to multiplicities
c         Isotropic distribution will be assumed.
c
          li=1
          ltt=0
          nk=nkk
          ni=0
          call wrtcont(nou,mat,mf,mt,ns,za,awr,li,ltt,nk,ni)         
        elseif (li.eq.1) then 
c
c         full isotropic distribution (li=1)
c        
          call wrtcont(nou,mat,mf,mt,ns,za,awr,li,ltt,nk,ni)
        elseif (li.eq.0.and.(ltt.eq.1.or.ltt.eq.2)) then
c
c         Legendre coefficients representation (LI=0,LTT=1) OR
c         Tabulated angular distribution       (LI=0,LTT=2)
c          
          ltt2=2
          call wrtcont(nou,mat,mf,mt,ns,za,awr,li,ltt2,nk,ni)
          if (ni.gt.0) then
            do i=1,ni
              call readcont(nin,c1,c2,l1,l2,n1,n2,mat,mf,mt,nsi)
              call wrtcont(nou,mat,mf,mt,ns,c1,c2,l1,l2,n1,n2)
            enddo
          endif
          nkk=nk-ni
          do k=1,nkk
            open (nscr1,file='mf14scr1.tmp')
            call readtab2(nin,c1,c2,l1,l2,nr,ne,nbt,ibt)
            call wrtab2(nscr1,mat,mf,mt,nsc,c1,c2,l1,l2,nr,ne,nbt,ibt)
            if (ltt.eq.1) then
c
c             Legendre coefficients representation (LI=0,LTT=1)
c             Convert to linearly interpolable tabulated data (LTT=2)
c            
              nr=1
              ibt(1)=2
              do j=1,ne
                call readlist(nin,c1,c2,l1,l2,nl,n2,b(2))
                b(1)=1.0d0
                call leg2lin(nl,b,nmu,x,y,ytol,ymin,npmax)
                nbt(1)=nmu
                call wrtab1(nscr1,mat,mf,mt,nsc,c1,c2,l1,l2,
     &                      nr,nbt,ibt,nmu,x,y)
              enddo
            else
c
c             Tabulated angular distribution(LI=0,LTT=2)
c             Linearize data, if required
c
              do j=1,ne
                call readtab1(nin,c1,c2,l1,l2,nr,np,nbt,ibt,x,y)
                call linear(nr,nbt,ibt,np,x,y,ytol,ymin,npmax)
                call wrtab1(nscr1,mat,mf,mt,nsc,c1,c2,l1,l2,
     &                      nr,nbt,ibt,np,x,y)                
              enddo
            endif
c           
c           write the complete tab2 record. Linearize tab2 data
c           as a function of incident energy E, if required 
c              
            call lintab2(nscr1,nou,mat,mf,mt,ns,ymin,ytol,npmax)
            close (nscr1,status='DELETE')
          enddo          
        else
          write(nlst,'(a)')' Error: LI/LTT option not valid in MF14'
          write(nlst,'(a,i2,a,i2)')'        LI=',LI,' LTT=',LTT
          write(*,'(a)')' Error: LI/LTT option not valid in MF14'
          write(*,'(a,i2,a,i2)')'        LI=',LI,' LTT=',LTT
          stop   
        endif
        read(nin,*)
        call wrtsend(nou,mat,mf,ns)
        call readcont(nin,za,awr,li,ltt,nk,ni,mat,mf,mt,nsi)
      enddo
      call wrtfend(nou,mat,ns)
      return
      end
C======================================================================
C      General routines for ENDF-6 formatted files
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
      subroutine lintab2(nscr1,nou,mat,mf,mt,ns,ymin,ytol,npmax)
c
c       linearize a tab2 record
c      
      implicit real*8 (a-h, o-z)
      parameter (kmax=20)
      data nscr2/91/
      dimension nbt(20),ibt(20),nbte(20),ibte(20)
      allocatable x(:),y(:),x0(:),y0(:),x1(:),y1(:),x2(:),y2(:)
      allocatable nps(:),e(:),xs(:,:),ys(:,:)
      rewind (nscr1)
      call readtab2(nscr1,eg,es,l1e,l2e,nre,ne,nbte,ibte)
      call checklaw(nre,ibte,icod)
      if (icod.eq.0) then
        nscr2=nscr1
        nne=ne
      else
        open(nscr2,file='mf1402.tmp')
        allocate(e(npmax),x(npmax),y(npmax))
        j=1
        call readtab1(nscr1,c1,e0,l1,l2,nr,np0,nbt,ibt,x,y)
        call wrtab1(nscr2,mat,mf,mt,nsc,c1,e0,l1,l2,nr,nbt,ibt,np0,x,y)                     
        nne=1
        allocate(x0(np0),y0(np0),e(kmax),nps(kmax))
        do i=1,np0
          x0(i)=x(i)
          y0(i)=y(i)
        enddo
        do j=2,ne
          call readtab1(nscr1,c1,e1,l1,l2,nr,np1,nbt,ibt,x,y)    
          allocate(x1(np1),y1(np1))
          do i=1,np1
            x1(i)=x(i)
            y1(i)=y(i)
          enddo
          call union(np0,x0,np1,x1,np2,x,npmax)
          allocate(xs(kmax,np2),ys(kmax,np2))
          l=1
          do while (l.lt.nre.and.nbte(l).lt.j)
            l=l+1
          enddo
          ilaw=ibte(l)          
          k=0
          nostop=1
          do while (nostop.eq.1)
            em=0.5d0*(e0+e1)
            iflag=1
            do i=1,np2
              xx=x(i)
              yy0=fvalin(np0,x0,y0,xx)
              yy1=fvalin(np1,x1,y1,xx)
              yym=0.5d0*(yy0+yy1)
              call terp1m(e0,yy0,e1,yy1,em,yylaw,ilaw)
              if (abs(yylaw-yym).gt.(abs(yylaw)*ytol)) iflag=0
              y(i)=yylaw
            enddo
            if (iflag.gt.0.or.k.eq.kmax) then
              nr=1
              ibt(1)=2
              nbt(1)=np1
              call wrtab1(nscr2,mat,mf,mt,nsc,c1,e1,l1,l2,
     &                    nr,nbt,ibt,np1,x1,y1)
              nne=nne+1
              if (k.eq.0) then
                nostop=0
              else
                deallocate (x0,y0)
                np0=np1
                e0=e1
                allocate(x0(np0),y0(np0))
                do i=1,np0
                  x0(i)=x1(i)
                  y0(i)=y1(i)
                enddo
                deallocate (x1,y1)
                np1=nps(k)
                e1=e(k)
                allocate(x1(np1),y1(np1))
                do i=1,np1
                  x1(i)=xs(k,i)
                  y1(i)=ys(k,i)
                enddo
                k=k-1
              endif
            else
              k=k+1
              nps(k)=np1
              e(k)=e1
              do i=1,np1
                xs(k,i)=x1(i)
                ys(k,i)=y1(i)
              enddo
              np1=np2
              e1=em
              do i=1,np2
                x1(i)=x(i)
                y1(i)=y(i)
              enddo
            endif  
          enddo
          deallocate (x0,y0)
          np0=np1
          e0=e1
          allocate (x0(np1),y0(np1))
          do i=1,np1
            x0(i)=x1(i)
            y0(i)=y1(i)
          enddo
          deallocate (x1,y1,xs,ys)
        enddo
        rewind(nscr2)
      endif
c
c       write linearized tab2 data
c
      nre=1
      nbte(1)=nne
      ibte(1)=2
      call wrtab2(nou,mat,mf,mt,ns,eg,es,l1e,l2e,nre,nne,nbte,ibte)
      do i=1,nne
        call readtab1(nscr2,c1,c2,l1,l2,nr,np,nbt,ibt,x,y)
        call wrtab1(nou,mat,mf,mt,ns,c1,c2,l1,l2,nr,nbt,ibt,np,x,y)        
      enddo
      if (nscr2.ne.nscr1) close(nscr2,status='DELETE')
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
C     Legendre polynomials
C======================================================================
      subroutine leg2lin(na,b,np,xt,yt,epy,ymin,npmax)
      implicit real*8 (a-h, o-z)
      parameter (nmu0=5, ns=25, h0=1.0d-8)
      dimension b(*),xt(*),yt(*)
      dimension xs(ns),ys(ns),xmu0(nmu0)
c
c     Starting convertion from Legendre representation to tabular data
c     
      if (na.gt.0) then
        xmu0(1)=-1.0d0
        xmu0(2)=-0.5d0
        xmu0(3)=0.0d0
        xmu0(4)=0.5d0
        xmu0(5)=1.0d0
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
          if (yt(i).lt.0.0d0) yt(i)=1.0d-30
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
C      The following subroutines were taken from NJOY2016 and
C      adapted/modified by D. Lopez Aldama for ACEMAKER:
C       1. subroutine legndr
C       2. subroutine terp1 (renamed as terp1m)
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
c      (borrowed and modified from NJOY)
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
