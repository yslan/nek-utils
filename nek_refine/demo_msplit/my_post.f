cc-----------------------------------------------------------------------
c      subroutine set_outfld
c
cc     Check if we are going to checkpoint at this timestep
cc     and set ifoutfld accordingly
c
c      include 'SIZE'
c      include 'TOTAL'
c      include 'CTIMER'
c
c      common /rdump/ ntdump
c
c      ifoutfld = .false.
c
c      if (iostep.le.0 .and. timeio.le.0) return
c
cc      if (istep.ge.nsteps) lastep=1
c
c      if (iostep.gt.0) then
c         if(mod(istep,iostep).eq.0) ifoutfld=.true.
c      else if (timeioe.ne.0.0) then
c         if (dnekclock_sync()-etimes .ge. (ntdump + 1)*timeio) then
c            ntdump=ntdump+1
c            ifoutfld=.true.
c         endif 
c      else if (timeio.ne.0.0) then
c         if (time.ge.(ntdump + 1)*timeio) then
c            ntdump=ntdump+1
c            ifoutfld=.true.
c         endif
c      endif
c
c      if (ioinfodmp.ne.0 .or. lastep.eq.1) ifoutfld=.true. 
c
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine check_ioinfo
c
cc     Check for io request in file 'ioinfo'
c
c      include 'SIZE'
c      include 'TSTEP'
c      include 'INPUT'
c
c      parameter (lxyz=lx1*ly1*lz1)
c      parameter (lpsc9=ldimt1+9)
c      common /ctmp1/ tdump(lxyz,lpsc9)
c      real*4         tdump
c      real           tdmp(4)
c      equivalence   (tdump,tdmp)
c
c      integer maxstep
c      save    maxstep
c      data    maxstep /999999999/
c
c      character*132 fname
c      character*1   fname1(132)
c      equivalence  (fname,fname1)
c
c      ioinfodmp=0
c      if (nid.eq.0 .and. (mod(istep,10).eq.0 .or. istep.lt.200)) then
c         call blank(fname1,size(fname1))
c         len = ltrunc(path,132)
c         call chcopy(fname1,path,len)
c         call chcopy(fname1(len+1),'ioinfo',6)
c         open(unit=87,file=fname,status='old',err=88)
c         read(87,*,end=87,err=87) idummy
c         if (ioinfodmp.eq.0) ioinfodmp=idummy
c         if (idummy.ne.0) then  ! overwrite last i/o request
c            rewind(87)
c            write(87,86)
c   86       format(' 0')
c         endif
c   87    continue
c         close(unit=87)
c   88    continue
c         if (ioinfodmp.ne.0) write(6,*) 'Output:',ioinfodmp
c      endif
c
c      tdmp(1)=ioinfodmp
c      call gop(tdmp,tdmp(3),'+  ',1)
c      ioinfodmp=tdmp(1)
c      if (ioinfodmp.lt.0) maxstep=abs(ioinfodmp)
c      if (istep.ge.maxstep.or.ioinfodmp.eq.-2) lastep=1
c
c      return
c      end
cc-----------------------------------------------------------------------
      subroutine my_prepost(elist,ifdoin,prefin)

c     Store results for later postprocessing
c
c     Recent updates:
c
c     p65 now indicates the number of parallel i/o files; iff p66 >= 6
c
c     we now check whether we are going to checkpoint in set_outfld
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)

      character*3    prefin,prefix

      logical  ifdoin
      integer  elist(1)

      if (ioinfodmp.eq.-2) return

c#ifdef TIMER
      etime1=dnekclock_sync()
c#endif

      prefix = prefin
      if (prefix.eq.'his') prefix = '   '

      if (ifdoin) then
         icalld=icalld+1
         nprep=icalld

         call prepost_map(0) ! map pr and axisymm. arrays
         call my_outfld(elist,prefix)
         call prepost_map(1) ! map back axisymm. arrays

c#ifdef TIMER
         tprep=tprep+dnekclock_sync()-etime1
c#endif
      endif

      return
      end
cc-----------------------------------------------------------------------
c      subroutine prepost_map(isave) ! isave=0-->fwd, isave=1-->bkwd
c
cc     Store results for later postprocessing
c
c      include 'SIZE'
c      include 'TOTAL'
cC
cC     Work arrays and temporary arrays
cC
c      common /scruz/ vxax   (lx1,ly1,lelv)
c     $             , vyax   (lx1,ly1,lelv)
c     $             , prax   (lx2,ly2,lelv)
c     $             , yax    (lx1,ly1,lelt)
c      common /scrmg/ tax    (lx1,ly1,lelt,ldimt)
c      common /scrcg/ pm1    (lx1,ly1,lz1,lelv)
cC
cc
c      common /prepst/ pa(lx1,ly2,lz2),pb(lx1,ly1,lz2)
c      integer e
c
c      if (isave.eq.0) then ! map to GLL grid
c
c         if (ifaxis) then
c            ntotm1 = lx1*ly1*nelt
c            call copy (yax,ym1,ntotm1)
c            do 5 e=1,nelt
c               if (ifrzer(e)) then
c                  call mxm  (ym1(1,1,1,e),lx1,iatjl1,ly1,pb,ly1)
c                  call copy (ym1(1,1,1,e),pb,lx1*ly1)
c               endif
c    5       continue
c            if (ifflow) then
c               ntotm1 = lx1*ly1*nelt
c               ntotm2 = lx2*ly2*nelt
c               call copy (vxax,vx,ntotm1)
c               call copy (vyax,vy,ntotm1)
c               call copy (prax,pr,ntotm2)
c               do 10 e=1,nelt
c                  if (ifrzer(e)) then
c                     call mxm  (vx(1,1,1,e),lx1,iatjl1,ly1,pb,ly1)
c                     call copy (vx(1,1,1,e),pb,lx1*ly1)
c                     call mxm  (vy(1,1,1,e),lx1,iatjl1,ly1,pb,ly1)
c                     call copy (vy(1,1,1,e),pb,lx1*ly1)
c                     call mxm  (pr(1,1,1,e),lx2,iatjl2,ly2,pb,ly2)
c                     call copy (pr(1,1,1,e),pb,lx2*ly2)
c                  endif
c 10            continue
c            endif
c            if (ifheat) then
c               ntotm1 = lx1*ly1*nelt
c               do 15 ifldt=1,npscal+1
c                  call copy (tax(1,1,1,ifldt),t(1,1,1,1,ifldt),ntotm1)
c 15            continue
c               do 30 e=1,nelt
c                  if (ifrzer(e)) then
c                    do 25 ifldt=1,npscal+1
c                      call mxm  (t(1,1,1,e,ifldt),lx1,iatjl1,ly1,
c     $                                                  pb,ly1)
c                      call copy (t(1,1,1,e,ifldt),pb,lx1*ly1)
c 25                 continue
c                  endif
c 30            continue
c            endif
c         endif
cC        Map the pressure onto the velocity mesh
cC
c         ntott = lx1*ly1*lz1*nelt
c         ntot1 = lx1*ly1*lz1*nelt
c         nyz2  = ly2*lz2
c         nxy1  = lx1*ly1
c         nxyz  = lx1*ly1*lz1
c         nxyz2 = lx2*ly2*lz2
cC
c         
c         call rzero(pm1,ntott)
c         if (ifsplit) then
c            call copy(pm1,pr,ntot1)
c         elseif (if_full_pres) then
c            call rzero(pm1,ntot1)
c            do e=1,nelt
c               call copy(pm1(1,1,1,e),pr(1,1,1,e),nxyz2)
c            enddo
c         else
c            do 1000 e=1,nelt
c               call mxm (ixm21,lx1,pr(1,1,1,e),lx2,pa(1,1,1),nyz2)        
c               do 100 iz=1,lz2
c                  call mxm (pa(1,1,iz),lx1,iytm21,ly2,pb(1,1,iz),ly1)
c  100          continue
c               call mxm (pb(1,1,1),nxy1,iztm21,lz2,pm1(1,1,1,e),lz1)
c 1000       continue
c         endif
c
c      else       ! map back
c
c         if (ifaxis) then
c            ntot1 = lx1*ly1*nelt
c            call copy (ym1,yax,ntot1)
c            if (ifflow) then
c               ntot1 = lx1*ly1*nelt
c               ntot2 = lx2*ly2*nelt
c               call copy (vx,vxax,ntot1)
c               call copy (vy,vyax,ntot1)
c               call copy (pr,prax,ntot2)
c            endif
c            if (ifheat) then
c               ntot1 = lx1*ly1*nelt
c               do 3000 ifldt=1,npscal+1
c                  call copy (t(1,1,1,1,ifldt),tax(1,1,1,ifldt),ntot1)
c 3000          continue
c            endif
c         endif
c
c      endif
c      return
c      end
cc-----------------------------------------------------------------------
      subroutine my_outfld(elist,prefix)

c     output .fld file 

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
C
C     Work arrays and temporary arrays
C
      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)
c
c     note, this usage of CTMP1 will be less than elsewhere if NELT ~> 3.
      parameter (lxyz=lx1*ly1*lz1)
      parameter (lpsc9=ldimt1+9)
      common /ctmp1/ tdump(lxyz,lpsc9)
      real*4         tdump
      real           tdmp(4)
      equivalence   (tdump,tdmp)

      real*4         test_pattern

      character*3    prefix
      character*1    fhdfle1(132)
      character*132   fhdfle
      equivalence   (fhdfle,fhdfle1)

      character*1 excode(30)
      character*10 frmat

      common /nopenf/ nopen(99)

      common /rdump/ ntdump
      data ndumps / 0 /

      logical ifxyo_s
      integer elist(1)

      if(nio.eq.0) then 
        WRITE(6,1001) istep,time
 1001   FORMAT(/,i9,1pe12.4,' Write checkpoint')
      endif
      call nekgsync()      

      p66 = param(66)
      if (abs(p66).eq.6) then
c         call err_chk(ierr,'mfo not supported in my_post$')
         call my_mfo_outfld(elist,prefix)        
         return
      endif

      ifxyo_s = ifxyo              ! Save ifxyo

      iprefix = i_find_prefix(prefix,99)

      ierr = 0
      if (nid.eq.0) then

c       Open new file for each dump on /cfs
        nopen(iprefix)=nopen(iprefix)+1

        if (prefix.eq.'   '.and.nopen(iprefix).eq.1) ifxyo = .true. ! 1st file

        if (prefix.eq.'rst'.and.max_rst.gt.0) 
     $         nopen(iprefix) = mod1(nopen(iprefix),max_rst) ! restart

        call file2(nopen(iprefix),prefix)
        if (p66.lt.1.0) then
           open(unit=24,file=fldfle,form='formatted',status='unknown')
        else
           call byte_open (fldfle,ierr)
c          write header as character string
           call blank(fhdfle,132)
        endif
      endif
      call bcast(ifxyo,lsize)
      if(p66.ge.1.0)
     $   call err_chk(ierr,'Error opening file in outfld. Abort. $')

C     Figure out what goes in EXCODE
      CALL BLANK(EXCODE,30)
      NDUMPS=NDUMPS+1
      i=1
      if (mod(p66,1.0).eq.0.0) then !old header format
         IF(IFXYO) then
            EXCODE(1)='X'
            EXCODE(2)=' '
            EXCODE(3)='Y'
            EXCODE(4)=' '
            i = 5
            IF(IF3D) THEN
              EXCODE(i)  ='Z'
              EXCODE(i+1)=' '
              i = i + 2
            ENDIF
         ENDIF
         IF(IFVO) then
            EXCODE(i)  ='U'
            EXCODE(i+1)=' '
            i = i + 2
         ENDIF
         IF(IFPO) THEN
           EXCODE(i)='P'
           EXCODE(i+1)=' '
           i = i + 2
         ENDIF
         IF(IFTO) THEN
           EXCODE(i)='T '
           EXCODE(i+1)=' '
           i = i + 1
         ENDIF
         do iip=1,ldimt1
            if (ifpsco(iip)) then
              write(excode(iip+I)  ,'(i1)') iip
              write(excode(iip+I+1),'(a1)') ' '
              i = i + 1 
            endif
         enddo
      else
         !new header format
         IF (IFXYO) THEN
            EXCODE(i)='X'
            i = i + 1
         ENDIF
         IF (IFVO) THEN
            EXCODE(i)='U'
            i = i + 1
         ENDIF
         IF (IFPO) THEN
            EXCODE(i)='P'
            i = i + 1
         ENDIF
         IF (IFTO) THEN
            EXCODE(i)='T'
            i = i + 1
         ENDIF
         IF (LDIMT.GT.1) THEN
            NPSCALO = 0
            do k = 1,ldimt-1
              if(ifpsco(k)) NPSCALO = NPSCALO + 1
            enddo
            IF (NPSCALO.GT.0) THEN
               EXCODE(i) = 'S'
               WRITE(EXCODE(i+1),'(I1)') NPSCALO/10
               WRITE(EXCODE(i+2),'(I1)') NPSCALO-(NPSCALO/10)*10
            ENDIF
         ENDIF
      endif
     

C     Dump header
      ierr = 0
      nelgt_bak=nelgt
      nelgt=iglsum(elist,nelt)
      if (nid.eq.0) write(*,*) 'my_outpost,nelgt',nelgt,nelgt_bak
      if (nid.eq.0) call dump_header(excode,p66,ierr)
      call err_chk(ierr,'Error dumping header in outfld. Abort. $')
      nelgt=nelgt_bak

      call get_id(id)

      nxyz  = lx1*ly1*lz1

      ierr = 0
      do ieg=1,nelgt

         jnid = gllnid(ieg)
         ie   = gllel (ieg)

         ifdump=0
         if (nid.eq.0) then
            if (jnid.eq.0) then
               call fill_tmp(tdump,id,ie)
               ifdump=elist(ie)
            else
               mtype=2000+ie
               len=4*id*nxyz
               dum1=0.
               call csend (mtype,dum1,wdsize,jnid,nullpid)
               call crecv2 (mtype,tdump,len,jnid)
               call crecv2 (mtype,ifdump,4,jnid)
            endif
c           write(*,*)'ddd',ieg,jnid,ie,ifdump,ierr
            if(ifdump.gt.0.AND.ierr.eq.0) call out_tmp(id,p66,ierr)
c           if(ierr.eq.0) call out_tmp(id,p66,ierr)
         elseif (nid.eq.jnid) then
            call fill_tmp(tdump,id,ie)
            dum1=0.
            mtype=2000+ie
            len=4*id*nxyz
            call crecv2 (mtype,dum1,wdsize,node0)
            call csend (mtype,tdump,len,node0,nullpid)
            call csend (mtype,elist(ie),4,node0,nullpid)
         endif
      enddo
      call err_chk(ierr,'Error writing file in outfld. Abort. $')

      ifxyo = ifxyo_s           ! restore ifxyo

      if (nid.eq.0) call close_fld(p66,ierr)
      call err_chk(ierr,'Error closing file in outfld. Abort. $')

      return
      end
cc-----------------------------------------------------------------------
c      subroutine file2(nopen,PREFIX)
cC----------------------------------------------------------------------
cC
cC     Defines machine specific input and output file names.
cC
cC----------------------------------------------------------------------
c      include 'SIZE'
c      include 'INPUT'
c      include 'TSTEP'
c      include 'PARALLEL'
cC
c      CHARACTER*132 NAME
c      CHARACTER*1   SESS1(132),PATH1(132),NAM1(132)
c      EQUIVALENCE  (SESSION,SESS1)
c      EQUIVALENCE  (PATH,PATH1)
c      EQUIVALENCE  (NAME,NAM1)
c      CHARACTER*1  DMP(4),FLD(4),REA(4),HIS(4),SCH(4) ,ORE(4), NRE(4)
c      CHARACTER*4  DMP4  ,FLD4  ,REA4  ,HIS4  ,SCH4   ,ORE4  , NRE4
c      EQUIVALENCE (DMP,DMP4), (FLD,FLD4), (REA,REA4), (HIS,HIS4)
c     $          , (SCH,SCH4), (ORE,ORE4), (NRE,NRE4)
c      CHARACTER*1  NUMRL(0:9)
c      DATA DMP4,FLD4,REA4 /'.dmp','.fld','.rea'/
c      DATA HIS4,SCH4      /'.his','.sch'/
c      DATA ORE4,NRE4      /'.ore','.nre'/
c      DATA NUMRL          /'0','1','2','3','4','5','6','7','8','9'/
c      CHARACTER*78  STRING
cc
c      character*1    prefix(3)
cC
c      call blank(name  ,132)
c      call blank(fldfle,132)
cC
c      LS=LTRUNC(SESSION,132)
c      LPP=LTRUNC(PATH,132)
c      LSP=LS+LPP
c      l = 0
c
cc     Construct file names containing full path to host:
cc      DO 100 I=1,LPP
cc         l = l+1
cc         NAM1(l)=PATH1(I)
cc  100 CONTINUE
cC
c      if (prefix(1).ne.' '.and.prefix(2).ne.' '.and.
c     $     prefix(3).ne.' ') then
c         do i=1,3
c            l = l+1
c            NAM1(l)=prefix(i)
c         enddo
c      endif
cC
c      DO 200 I=1,LS
c         l = l+1
c         NAM1(l)=SESS1(I)
c  200 CONTINUE
cC
cC .fld file
c      DO 300 I=1,4
c         l = l+1
c         NAM1(l)=FLD(I)
c  300 CONTINUE
c      if (nopen.lt.100) then
cC        less than 100 dumps....
c         ITEN=NOPEN/10
c         l = l+1
c         NAM1(l)=NUMRL(ITEN)
c         IONE=MOD(NOPEN,10)
c         l = l+1
c         NAM1(l)=NUMRL(IONE)
c      elseif (nopen.lt.1000) then
cC        less than 1000 dumps....
c         IHUN=NOPEN/100
c         l = l+1
c         NAM1(l)=NUMRL(IHUN)
c         ITEN=MOD(NOPEN,100)/10
c         l = l+1
c         NAM1(l)=NUMRL(ITEN)
c         IONE=MOD(NOPEN,10)
c         l = l+1
c         NAM1(l)=NUMRL(IONE)
c      elseif (nopen.lt.10000) then
cC        less than 10000 dumps....
c         ITHO=NOPEN/1000
c         l = l+1
c         NAM1(l)=NUMRL(ITHO)
c         IHUN=MOD(NOPEN,1000)/100
c         l = l+1
c         NAM1(l)=NUMRL(IHUN)
c         ITEN=MOD(NOPEN,100)/10
c         l = l+1
c         NAM1(l)=NUMRL(ITEN)
c         IONE=MOD(NOPEN,10)
c         l = l+1
c         NAM1(l)=NUMRL(IONE)
c      endif
c      FLDFLE=NAME
cC
cC     Write the name of the .fld file to the logfile.
cC
c      if (nio.eq.0) then
c         call chcopy(string,fldfle,78)
c         write(6,1000) istep,time,string
c 1000    format(/,i9,1pe12.4,' OPEN: ',a78)
c      endif
c 
c      return
c      end
cc=======================================================================
c      subroutine rzero4(a,n)
c      real*4 A(1)
c      DO 100 I = 1, N
c 100     A(I ) = 0.0
c      return
c      end
cc=======================================================================
c      subroutine copyX4(a,b,n)
c      REAL*4 A(1)
c      REAL   B(1)
c      DO 100 I = 1, N
c 100     A(I) = B(I)
c      return
c      end
cc=======================================================================
c      subroutine copy4r(a,b,n)
c      real   a(1)
c      real*4 b(1)
c      do i = 1, n
c         a(i) = b(i)
c      enddo
c      return
c      end
cc=======================================================================
c      function i_find_prefix(prefix,imax)
cc
c      character*3 prefix
c      character*3 prefixes(99)
c      save        prefixes
c      data        prefixes /99*'...'/
cc
c      integer nprefix
c      save    nprefix
c      data    nprefix /0/
cc
cc     Scan existing list of prefixes for a match to "prefix"
cc
c      do i=1,nprefix
c         if (prefix.eq.prefixes(i)) then
c            i_find_prefix = i
c            return
c         endif
c      enddo
cc
cc     If we're here, we didn't find a match.. bump list and return
cc
c      nprefix                = nprefix + 1
c      prefixes(nprefix)      = prefix
c      i_find_prefix          = nprefix
cc
cc     Array bounds check on prefix list
cc
c      if (nprefix.gt.99.or.nprefix.gt.imax) then
c         write(6,*) 'Hey! nprefix too big! ABORT in i_find_prefix'
c     $      ,nprefix,imax
c         call exitt
c      endif
cc
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine dump_header(excodein,p66,ierr)
c
c      include 'SIZE'
c      include 'TOTAL'
c
c      character*30  excodein
c
c      character*30 excode
c      character*1  excode1(30)
c      equivalence (excode,excode1) 
c
c      real*4         test_pattern
c
c      character*1 fhdfle1(132)
c      character*132 fhdfle
c      equivalence (fhdfle,fhdfle1)
c
c      write(excode,'(A30)') excodein
c
c      ikstep = istep
c      do ik=1,10
c         if (ikstep.gt.9999) ikstep = ikstep/10
c      enddo
c
c      call blank(fhdfle,132)
c
cc       write(6,111)               !       print on screen
cc     $     nelgt,lx1,ly1,lz1,time,istep,excode
cc
c      if (mod(p66,1.0).eq.0.0) then !       old header format
c         if (p66.lt.1.0) then       !ASCII
c           if(nelgt.lt.10000) then
c            WRITE(24,'(4i4,1pe14.7,I5,1X,30A1,1X,A12)')
c     $           NELGT,lx1,ly1,lz1,TIME,ikstep,(EXCODE1(I),I=1,30),
c     $           'NELT,NX,NY,N'
c           else
c            WRITE(24,'(i10,3i4,1pe18.9,I9,1X,30A1,1X,A12)')
c     $           NELGT,lx1,ly1,lz1,TIME,ikstep,(EXCODE1(I),I=1,30),
c     $           'NELT,NX,NY,N'
c           endif
c         else                       !Binary
c            if (nelgt.lt.10000) then
c               WRITE(fhdfle,'(4I4,1pe14.7,I5,1X,30A1,1X,A12)')
c     $              NELGT,lx1,ly1,lz1,TIME,ikstep,(EXCODE1(I),I=1,30),
c     $              ' 4 NELT,NX,NY,N'
c            else
c               write(fhdfle,'(i10,3i4,1P1e18.9,i9,1x,30a1)')
c     $         nelgt,lx1,ly1,lz1,time,istep,(excode1(i),i=1,30)
c            endif
c            call byte_write(fhdfle,20,ierr)
c         endif
c      else                        !       new header format
c         if (p66.eq.0.1) then
c            write(24,111)
c     $           nelgt,lx1,ly1,lz1,time,istep,excode
c        else       
c             write(fhdfle,111)
c     $            nelgt,lx1,ly1,lz1,time,istep,excode
c             call byte_write(fhdfle,20,ierr)
c        endif
c 111    FORMAT(i10,1x,i2,1x,i2,1x,i2,1x,1P1e18.9,1x,i9,1x,a)
c      endif
c
c      if(ierr.ne.0) return
c
c      CDRROR=0.0
c      if (p66.LT.1.0) then       !       formatted i/o
c         WRITE(24,'(6G11.4)')(CDRROR,I=1,NELGT)   ! dummy 
c      else
cC       write byte-ordering test pattern to byte file...
c        test_pattern = 6.54321
c        call byte_write(test_pattern,1,ierr)
c      endif
c
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine fill_tmp(tdump,id,ie)
cC
c      include 'SIZE'
c      include 'TOTAL'
cc
c      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)
cC
cC     Fill work array
cC
c      parameter (lxyz=lx1*ly1*lz1)
c      parameter (lpsc9=ldimt1+9)
c      real*4 tdump(lxyz,lpsc9)
cC
c      nxyz = lx1*ly1*lz1
cc
c      ID=0
c      IF(IFXYO)then
c         ID=ID+1
c         CALL COPYx4(TDUMP(1,ID),XM1(1,1,1,IE),NXYZ)
c         ID=ID+1
c         CALL COPYx4(TDUMP(1,ID),YM1(1,1,1,IE),NXYZ)
c         IF(IF3D) then
c            ID=ID+1
c            CALL COPYx4(TDUMP(1,ID),ZM1(1,1,1,IE),NXYZ)
c         ENDIF
c      ENDIF
cc
c      IF(IFVO)then
c         IF (IE.LE.NELT) then
c            ID=ID+1
c            CALL COPYx4(TDUMP(1,ID),VX(1,1,1,IE),NXYZ)
c            ID=ID+1
c            CALL COPYx4(TDUMP(1,ID),VY(1,1,1,IE),NXYZ)
c            IF(IF3D)then
c               ID=ID+1
c               CALL COPYx4(TDUMP(1,ID),VZ(1,1,1,IE),NXYZ)
c            ENDIF
c         ELSE
c            ID=ID+1
c            CALL RZERO4(TDUMP(1,ID),NXYZ)
c            ID=ID+1
c            CALL RZERO4(TDUMP(1,ID),NXYZ)
c            IF(IF3D)then
c               ID=ID+1
c               CALL RZERO4(TDUMP(1,ID),NXYZ)
c            ENDIF
c         ENDIF
c      ENDIF
c      IF(IFPO)then
c         IF (IE.LE.NELT) then
c            ID=ID+1
c            CALL COPYx4(TDUMP(1,ID),PM1(1,1,1,IE),NXYZ)
c         ELSE
c            ID=ID+1
c            CALL RZERO4(TDUMP(1,ID),NXYZ)
c         ENDIF
c      ENDIF
c      IF(IFTO)then
c         ID=ID+1
c         CALL COPYx4(TDUMP(1,ID),T(1,1,1,IE,1),NXYZ)
c      ENDIF
cC     PASSIVE SCALARS
c      do iip=1,ldimt1
c         if (ifpsco(iip)) then
c            id=id+1
c            call copyX4(tdump(1,id),t(1,1,1,ie,iip+1),nxyz)
c        endif
c      enddo
cc
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine get_id(id) !  Count amount of data to be shipped
c
c      include 'SIZE'
c      include 'TOTAL'
c
c      id=0
c
c      if (ifxyo) id=id+ldim
c      if (ifvo)  id=id+ldim
c      if (ifpo)  id=id+1
c      if (ifto)  id=id+1
c
c      do iip=1,ldimt1
c         if (ifpsco(iip)) id=id+1     !     Passive scalars
c      enddo
c
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine close_fld(p66,ierr)
c
c      include 'SIZE'
c      include 'TOTAL'
c
c      if (nid.eq.0) then
c         if (p66.lt.1) then
c            close(unit=24)
c         else
c            call byte_close(ierr)
c         endif
c      endif
c
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine out_tmp(id,p66,ierr)
c
c      include 'SIZE'
c      include 'TOTAL'
c
c      parameter (lxyz=lx1*ly1*lz1)
c      parameter (lpsc9=ldimt1+9)
c
c      common /ctmp1/ tdump(lxyz,lpsc9)
c      real*4         tdump
c
c      character*11 frmat
c
c      nxyz = lx1*ly1*lz1
c
c      call blank(frmat,11)
c      if (id.le.9) then
c         WRITE(FRMAT,1801) ID
c 1801    FORMAT('(1p',I1,'e14.6)')
c      else
c         WRITE(FRMAT,1802) ID
c 1802    FORMAT('(1p',I2,'e14.6)')
c      endif
c
c      if (p66.lt.1.0) then
cC       formatted i/o
c        WRITE(24,FRMAT)
c     $      ((TDUMP(I,II),II=1,ID),I=1,NXYZ)
c      else
cC        C binary i/o
c         do ii=1,id
c            call byte_write(tdump(1,ii),nxyz,ierr)
c            if(ierr.ne.0) goto 101
c         enddo
c      endif
c 101  continue
c
c      return
c      end
c-----------------------------------------------------------------------
      subroutine my_mfo_outfld(elist,prefix)  ! muti-file output

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)  ! mapped pressure

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzo8
      integer elist(1)
      character*3 prefix
      logical ifxyo_s
 
      common /SCRUZ/  ur1(lxo*lxo*lxo*lelt)
     &              , ur2(lxo*lxo*lxo*lelt)
     &              , ur3(lxo*lxo*lxo*lelt)

      integer lglel_bak(lelt),lglel2(lelt)

      tiostart=dnekclock_sync()

      call io_init
      if(nid.eq.0)write(*,*)'my_mfo_outfld',param(66),ifmpiio

      ifxyo_s = ifxyo 
      ifxyo_  = ifxyo
      nout = nelt ! map2reg
c      nout = ivlsum(elist,nelt) !nelt
      nxo  = lx1
      nyo  = ly1
      nzo  = lz1
      nelgt2= iglsum(elist,nelt)
      nout2 = ivlsum(elist,nelt)

      nn = nout2
      nelB2 = igl_running_sum(nn)
      nelB2 = nelB2 - nout2

      j=0
      call izero(lglel2,lelt)
      do ie=1,nelt
        ieg=lglel(ie)
        if (elist(ie).gt.0) then
          j = j + 1
          lglel2(j)=ieg
        endif
      enddo

      if (ifreguo) then ! dump on regular (uniform) mesh
         if (nrg.gt.lxo) then
            if (nid.eq.0) write(6,*) 
     &         'WARNING: nrg too large, reset to lxo!'
            nrg = lxo
         endif
         nxo  = nrg
         nyo  = nrg
         nzo  = 1
         if(if3d) nzo = nrg
      endif
c      offs0 = iHeaderSize + 4 + isize*nelgt
      offs0 = iHeaderSize + 4 + isize*nelgt2
      if(nid.eq.0)write(*,*)'my_mfo_outfld'
     $      ,nelB,nelB2,nout,nout2,nelgt,nelgt2

      ierr=0
      if (nid.eq.pid0) then
         call mfo_open_files(prefix,ierr)         ! open files on i/o nodes
      endif
      call err_chk(ierr,'Error opening file in mfo_open_files. $')
      call bcast(ifxyo_,lsize)
      ifxyo = ifxyo_

      nelgt_bak=nelgt
      nelt_bak=nelt
      nelB_bak=nelB
      call icopy(lglel_bak,lglel,lelt)

      nelgt=nelgt2
      nelt=nout2
      nelB=nelB2
      call icopy(lglel,lglel2,lelt)

      call mfo_write_hdr                     ! create element mapping +

      nelgt=nelgt_bak
      nelt=nelt_bak
      nelB=nelB_bak
      call icopy(lglel,lglel_bak,lelt)

c     call exitti('this is wdsizo A:$',wdsizo)
                                             ! write hdr
      nxyzo8  = nxo*nyo*nzo
c      strideB = nelB * nxyzo8*wdsizo
c      stride  = nelgt* nxyzo8*wdsizo
      strideB = nelB2 * nxyzo8*wdsizo
      stride  = nelgt2* nxyzo8*wdsizo

      ioflds = 0
      ! dump all fields based on the t-mesh to avoid different
      ! topologies in the post-processor
      if (ifxyo) then
         offs = offs0 + ldim*strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifreguo) then
            call map2reg(ur1,nrg,xm1,nout)
            call map2reg(ur2,nrg,ym1,nout)
            if (if3d) call map2reg(ur3,nrg,zm1,nout)
            call my_mfo_outv(elist,ur1,ur2,ur3,nout,nxo,nyo,nzo)
         else
            call my_mfo_outv(elist,xm1,ym1,zm1,nout,nxo,nyo,nzo)
         endif
         ioflds = ioflds + ldim
      endif
c      ierr = 0
c      if (nid.eq.pid0) then
c         if(ifmpiio) then
c           call byte_close_mpi(ifh_mbyte,ierr)
c         else
c           call byte_close(ierr)
c         endif
c      endif
c      call err_chk(ierr,'Error closing file in mfo_outfld. Abort. $')   
c      return

      if (ifvo ) then
         offs = offs0 + ioflds*stride + ldim*strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifreguo) then
            call map2reg(ur1,nrg,vx,nout)
            call map2reg(ur2,nrg,vy,nout)
            if (if3d) call map2reg(ur3,nrg,vz,nout)
            call my_mfo_outv(elist,ur1,ur2,ur3,nout,nxo,nyo,nzo) 
         else
            call my_mfo_outv(elist,vx,vy,vz,nout,nxo,nyo,nzo)  ! B-field handled thru outpost
         endif
         ioflds = ioflds + ldim
      endif
      if (ifpo ) then
         offs = offs0 + ioflds*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifreguo) then
            call map2reg(ur1,nrg,pm1,nout)
            call my_mfo_outs(elist,ur1,nout,nxo,nyo,nzo)
         else
            call my_mfo_outs(elist,pm1,nout,nxo,nyo,nzo)
         endif
         ioflds = ioflds + 1
      endif
      if (ifto ) then
         offs = offs0 + ioflds*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifreguo) then
            call map2reg(ur1,nrg,t,nout)
            call my_mfo_outs(elist,ur1,nout,nxo,nyo,nzo)
         else
            call my_mfo_outs(elist,t,nout,nxo,nyo,nzo)
         endif
         ioflds = ioflds + 1
      endif
      do k=1,ldimt-1
         if(ifpsco(k)) then
           offs = offs0 + ioflds*stride + strideB
           call byte_set_view(offs,ifh_mbyte)
           if (ifreguo) then
              call map2reg(ur1,nrg,t(1,1,1,1,k+1),nout)
              call my_mfo_outs(elist,ur1,nout,nxo,nyo,nzo)
           else
              call my_mfo_outs(elist,t(1,1,1,1,k+1),nout,nxo,nyo,nzo)
           endif
           ioflds = ioflds + 1
         endif
      enddo
c      dnbyte = 1.*ioflds*nout*wdsizo*nxo*nyo*nzo
      dnbyte = 1.*ioflds*nout2*wdsizo*nxo*nyo*nzo

      if (if3d) then
         offs0   = offs0 + ioflds*stride
c         strideB = nelB *2*4   ! min/max single precision
c         stride  = nelgt*2*4
         strideB = nelB2 *2*4   ! min/max single precision
         stride  = nelgt2*2*4
         ioflds  = 0
         ! add meta data to the end of the file
         if (ifxyo) then
            offs = offs0 + ldim*strideB
            call byte_set_view(offs,ifh_mbyte)
            call my_mfo_mdatav(elist,xm1,ym1,zm1,nout)
            ioflds = ioflds + ldim
         endif
         if (ifvo ) then
            offs = offs0 + ioflds*stride + ldim*strideB
            call byte_set_view(offs,ifh_mbyte)
            call my_mfo_mdatav(elist,vx,vy,vz,nout)
            ioflds = ioflds + ldim
         endif
         if (ifpo ) then
            offs = offs0 + ioflds*stride + strideB
            call byte_set_view(offs,ifh_mbyte)
            call my_mfo_mdatas(elist,pm1,nout)
            ioflds = ioflds + 1
         endif
         if (ifto ) then
            offs = offs0 + ioflds*stride + strideB
            call byte_set_view(offs,ifh_mbyte)
            call my_mfo_mdatas(elist,t,nout)
            ioflds = ioflds + 1
         endif
         do k=1,ldimt-1
            offs = offs0 + ioflds*stride + strideB
            call byte_set_view(offs,ifh_mbyte)
            if(ifpsco(k)) call my_mfo_mdatas(elist,t(1,1,1,1,k+1),nout)
            ioflds = ioflds + 1
         enddo
c         dnbyte = dnbyte + 2.*ioflds*nout*wdsizo
         dnbyte = dnbyte + 2.*ioflds*nout2*wdsizo
      endif

      ierr = 0
      if (nid.eq.pid0) then 
         if(ifmpiio) then
           call byte_close_mpi(ifh_mbyte,ierr)
         else
           call byte_close(ierr)
         endif
      endif
      call err_chk(ierr,'Error closing file in mfo_outfld. Abort. $')

      tio = dnekclock_sync()-tiostart
      if (tio.le.0) tio=1.

      dnbyte = glsum(dnbyte,1)
c      dnbyte = dnbyte + iHeaderSize + 4. + isize*nelgt
      dnbyte = dnbyte + iHeaderSize + 4. + isize*nelgt2
      dnbyte = dnbyte/1024/1024
      if(nio.eq.0) write(6,7) istep,time,dnbyte,dnbyte/tio,
     &             nfileo
    7 format(/,i9,1pe12.4,' done :: Write checkpoint',/,
     &       30X,'file size = ',3pG12.2,'MB',/,
     &       30X,'avg data-throughput = ',0pf7.1,'MB/s',/,
     &       30X,'io-nodes = ',i5,/)

      ifxyo = ifxyo_s ! restore old value

      return
      end
c-----------------------------------------------------------------------
c      subroutine io_init ! determine which nodes will output
c      character*132 hname
c
c      include 'SIZE'
c      include 'INPUT'
c      include 'PARALLEL'
c      include 'RESTART'
c
c      ifdiro = .false.
c
c      ifmpiio = .false.
c      if(abs(param(65)).eq.1 .and. abs(param(66)).eq.6) ifmpiio=.true.
c#ifdef NOMPIIO
c      ifmpiio = .false.
c#endif
c
c      wdsizo = 4
c      if (param(63).gt.0) wdsizo = 8 ! 64-bit .fld file
c      nrg = lxo
c
c      if(ifmpiio) then
c        nfileo  = np
c        nproc_o = 1
c        fid0    = 0
c        pid0    = nid
c        pid1    = 0
c      else
c        if(param(65).lt.0) ifdiro = .true. !  p65 < 0 --> multi subdirectories
c        nfileo  = abs(param(65))
c        if(nfileo.eq.0) nfileo = 1
c        if(np.lt.nfileo) nfileo=np   
c        nproc_o = np / nfileo              !  # processors pointing to pid0
c        fid0    = nid/nproc_o              !  file id
c        pid0    = nproc_o*fid0             !  my parent i/o node
c        pid1    = min(np-1,pid0+nproc_o-1) !  range of sending procs
c      endif
c
c      ! how many elements are present up to rank nid
c      nn = nelt
c      nelB = igl_running_sum(nn)
c      nelB = nelB - nelt
c     
c      pid00 = glmin(pid0,1)
c
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine mfo_open_files(prefix,ierr) ! open files
c
c      include 'SIZE'
c      include 'INPUT'
c      include 'PARALLEL'
c      include 'RESTART'
c
c      character*1 prefix(3)
c      character*3 prefx
c
c      character*132  fname
c      character*1    fnam1(132)
c      equivalence   (fnam1,fname)
c
c      character*6  six,str
c      save         six
c      data         six / "??????" /
c
c
c      character*1 slash,dot
c      save        slash,dot
c      data        slash,dot  / '/' , '.' /
c
c      integer nopen(99,2)
c      save    nopen
c      data    nopen  / 198*0 /
c
c      call blank(fname,132)      !  zero out for byte_open()
c
c      iprefix        = i_find_prefix(prefix,99)
c      if (ifreguo) then
c         nopen(iprefix,2) = nopen(iprefix,2)+1
c         nfld             = nopen(iprefix,2)
c      else
c         nopen(iprefix,1) = nopen(iprefix,1)+1
c         nfld             = nopen(iprefix,1)
c      endif
c
c      call chcopy(prefx,prefix,3)        ! check for full-restart request
c      if (prefx.eq.'rst'.and.max_rst.gt.0) nfld = mod1(nfld,max_rst)
c
c      call restart_nfld( nfld, prefix ) ! Check for Restart option.
c      if (prefx.eq.'   '.and.nfld.eq.1) ifxyo_ = .true. ! 1st file
c
c      if(ifmpiio) then
c        rfileo = 1
c      else
c        rfileo = nfileo
c      endif
c      ndigit = log10(rfileo) + 1
c
c      lenp = ltrunc(path,132)
c      call chcopy(fnam1(1),path,lenp)    
c      k = 1 + lenp 
c 
c      if (ifdiro) then                                  !  Add directory
c         call chcopy(fnam1(k),'A',1)
c         k = k + 1
c         call chcopy(fnam1(k),six,ndigit)  ! put ???? in string
c         k = k + ndigit
c         call chcopy(fnam1(k),slash,1)
c         k = k + 1
c      endif
c
c      if (prefix(1).ne.' '.and.prefix(2).ne.' '.and.    !  Add prefix
c     $    prefix(3).ne.' ') then
c         call chcopy(fnam1(k),prefix,3)
c         k = k + 3
c      endif
c
c      len=ltrunc(session,132)                           !  Add SESSION
c      call chcopy(fnam1(k),session,len)
c      k = k+len
c     
c      if (ifreguo) then
c         len=4
c         call chcopy(fnam1(k),'_reg',len)
c         k = k+len
c      endif
c
c      call chcopy(fnam1(k),six,ndigit)                  !  Add file-id holder
c      k = k + ndigit
c
c      call chcopy(fnam1(k  ),dot,1)                     !  Add .f appendix
c      call chcopy(fnam1(k+1),'f',1)
c      k = k + 2
c
c      write(str,4) nfld                                 !  Add nfld number
c    4 format(i5.5)
c      call chcopy(fnam1(k),str,5)
c      k = k + 5
c
c      call addfid(fname,fid0)
c
c      if(ifmpiio) then
c        if(nio.eq.0)    write(6,*) '      FILE:',fname 
c        call byte_open_mpi(fname,ifh_mbyte,.false.,ierr) 
c      else
c        if(nid.eq.pid0) write(6,*) '      FILE:',fname 
c        call byte_open(fname,ierr)
c      endif
c 
c      return
c      end
cc-----------------------------------------------------------------------
c
c      subroutine restart_nfld( nfld, prefix ) 
c      include 'SIZE' ! For nio
c      character*3 prefix
cc
cc     Check for Restart option and return proper nfld value.
cc     Also, convenient spot to explain restart strategy.
cc
cc
cc     The approach is as follows:
cc
cc         Prefix rs4 would indicate 4 files in the restart cycle.
cc         
cc         This would be normal usage for velocity only, with
cc         checkpoints taking place in synch with standard io.
cc
cc         The resultant restart sequence might look like:
cc
cc         blah.fld09           Step 0
cc         rs4blah.fld01             1
cc         rs4blah.fld02             2
cc
cc         which implies that fld09 would be used as the i.c.
cc         in the restart, rs4blah.fld01 would overwrite the
cc         solution at Step 1, and rs4blah.fld02 would overwrite
cc         Step 2.   Net result is that Steps 0-2 of the restart
cc         session have solutions identical to those computed in
cc         the prior run.   (It's important that both runs use
cc         the same dt in this case.)
cc
cc
cc         Another equally possible restart sequence would be:
cc
cc
cc         blah.fld10           Step 0
cc         rs4blah.fld03             1
cc         rs4blah.fld04             2
cc
cc         Why the 3 & 4 ?   If one were to use only 1 & 2, there
cc         is a risk that the system crashes while writing, say,
cc         rs4blah.fld01, in which case the restart is compromised --
cc         very frustrating at the end of a run that has been queued
cc         for a week.  By providing a toggled sequence in pairs such as
cc
cc         (1,2),   (3,4),  (1,2), ...
cc
cc         ensures that one always has at least one complete restart
cc         sequence.   In the example above, the following files would
cc         be written, in order:
cc
cc         :
cc         :
cc         blah.fld09
cc         rs4blah.fld01
cc         rs4blah.fld02
cc         blah.fld10
cc         rs4blah.fld03
cc         rs4blah.fld04
cc         blah.fld11
cc         rs4blah.fld01       (overwriting existing rs4blah.fld01)
cc         rs4blah.fld02       (    "           "        "  .fld02)
cc         blah.fld12
cc         rs4blah.fld03       (   etc.  )
cc         rs4blah.fld04
cc         :
cc         :
cc
cc
cc         Other strategies are possible, according to taste.
cc
cc         Here is a data-intensive one:
cc
cc         MHD + double-precision restart, but single-precision std files
cc
cc         In this case, single-precision files are kept as the running
cc         file sequence (i.e., for later post-processing) but dbl-prec.
cc         is required for restart.  A total of 12 temporary restart files
cc         must be saved:  (3 for velocity, 3 for B-field) x 2 for redundancy.
cc
cc         This is expressed, using hexadecimal notation (123456789abc...),
cc         as prefix='rsc'.
cc         
cc         
c      character*16 kst
c      save         kst
c      data         kst / '0123456789abcdef' /
c      character*1  ks1(0:15),kin
c      equivalence (ks1,kst)
c
cc
cc
c      if (indx1(prefix,'rs',2).eq.1) then
c         read(prefix,3) kin
c    3    format(2x,a1)
c         do kfld=1,15
c            if (ks1(kfld).eq.kin) goto 10
c         enddo
c   10    if (kfld.eq.16) kfld=4 ! std. default
c         nfln = mod1(nfld,kfld) ! Restart A (1,2) and B (3,4)
c         if (nio.eq.0) write(6,*) nfln,nfld,kfld,' kfld'
c         nfld = nfln
c      endif
c
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine full_restart_save(iosave)
c
c      integer iosave,save_size,nfld_save
c      logical if_full_pres_tmp
c
c      include 'SIZE'
c      include 'INPUT'
c      include 'TSTEP'
c
c      if (PARAM(27).lt. 0) then
c          nfld_save=abs(PARAM(27))  ! For full restart
c      else 
c          nfld_save=3
c      endif
c      save_size=8  ! For full restart
c
c      dtmp = param(63)
c      if_full_pres_tmp = if_full_pres     
c
c      param(63) = 1 ! Enforce 64-bit output
c      if_full_pres = .true. !Preserve mesh 2 pressure
c
c      if (lastep.ne.1) call restart_save(iosave,nfld_save)
c
c      param(63) = dtmp
c      if_full_pres = if_full_pres_tmp 
c
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine restart_save(iosave,nfldi)
c
c      integer iosave,nfldi
c
c
cc     Save current fields for later restart.
cc
cc     Input arguments:
cc
cc       .iosave plays the usual triggering role, like iostep
cc
cc       .nfldi is the number of rs files to save before overwriting
cc
c
c      include 'SIZE'
c      include 'TOTAL'
c      include 'RESTART'
c
c      character*3 prefix
c
c      character*17 kst
c      save         kst
c      data         kst / '0123456789abcdefx' /
c      character*1  ks1(0:16)
c      equivalence (ks1,kst)
c
c      logical if_full_pres_tmp
c
c      iosav = iosave
c
c      if (iosav.eq.0) iosav = iostep
c      if (iosav.eq.0) return
c
c      iotest = 0
cc     if (iosav.eq.iostep) iotest = 1  ! currently spoiled because of 
cc                                      ! incompatible format of .fld
cc                                      ! and multi-file i/o;  the latter
cc                                      ! is the only form used for restart
c
c      nfld  = nfldi*2
c      nfld2 = nfld/2
c      mfld  = min(17,nfld)
c      if (ifmhd) nfld2 = nfld/4
c
c      i2 = iosav/2
c      m1 = istep+iosav-iotest
c      mt = mod(istep+iosav-iotest,iosav)
c      prefix = '   '
c
c      if (istep.gt.iosav/2  .and.
c     $   mod(istep+iosav-iotest,iosav).lt.nfld2) then ! save
c         write(prefix,'(A)') 'rs_'
cc         write(prefix,3) ks1(mfld)
cc    3    format('rs',a1)
c
c         p66 = param(66)
c         param(66) = 6
c         if (ifmhd) call outpost2(bx,by,bz,pm,t,0,prefix)  ! first B
c         call prepost (.true.,prefix)
c         param(66) = p66
c
c      endif
c
c      return
c      end
c-----------------------------------------------------------------------
      subroutine my_outpost(elist,v1,v2,v3,vp,vt,name3)

      include 'SIZE'
      include 'INPUT'

      real v1(1),v2(1),v3(1),vp(1),vt(1)
      integer elist(1)
      character*3 name3


      itmp=0
      if (ifto) itmp=1
      call my_outpost2(elist,v1,v2,v3,vp,vt,itmp,name3)

      return
      end
c-----------------------------------------------------------------------
      subroutine my_outpost2(elist,v1,v2,v3,vp,vt,nfldt,name3)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'

      parameter(ltot1=lx1*ly1*lz1*lelt)
      parameter(ltot2=lx2*ly2*lz2*lelv)
      common /outtmp/  w1(ltot1),w2(ltot1),w3(ltot1),wp(ltot2)
     &                ,wt(ltot1,ldimt)
c
      real v1(1),v2(1),v3(1),vp(1),vt(ltot1,1)
      integer elist(1)
      character*3 name3
      logical if_save(ldimt)
c
      ntot1  = lx1*ly1*lz1*nelt
      ntot1t = lx1*ly1*lz1*nelt
      ntot2  = lx2*ly2*lz2*nelt

      if(nfldt.gt.ldimt) then
        write(6,*) 'ABORT: outpost data too large (nfldt>ldimt)!'
        call exitt
      endif

c foolproof 
      do ie=1,nelt
        if(elist(ie).gt.0) then
          elist(ie)=1
        else
          elist(ie)=0
        endif
      enddo

      itmp1=iglmax(elist,nelt)
      itmp2=iglmin(elist,nelt)
      if (itmp1.eq.0.AND.itmp2.eq.0) then
         write(6,*)'null elist in my_outpost, return'
         return
      endif
      if (itmp1.lt.0.OR.itmp2.gt.1) then
         if (nid.eq.0) write(6,*)'ABORT: elist out of range',itmp1,itmp2
         call exitt
      endif

c store solution
      call copy(w1,vx,ntot1)
      call copy(w2,vy,ntot1)
      call copy(w3,vz,ntot1)
      call copy(wp,pr,ntot2)
      do i = 1,nfldt
         call copy(wt(1,i),t(1,1,1,1,i),ntot1t)
      enddo

c swap with data to dump
      call copy(vx,v1,ntot1)
      call copy(vy,v2,ntot1)
      call copy(vz,v3,ntot1)
      call copy(pr,vp,ntot2)
      do i = 1,nfldt
         call copy(t(1,1,1,1,i),vt(1,i),ntot1t)
      enddo

c dump data
      if_save(1) = ifto
      ifto = .false.
      if(nfldt.gt.0) ifto = .true. 
      do i = 1,ldimt-1
         if_save(i+1) = ifpsco(i)
         ifpsco(i) = .false.   
         if(i+1.le.nfldt) ifpsco(i) = .true.
      enddo

      call my_prepost(elist,.true.,name3)

      ifto = if_save(1)
      do i = 1,ldimt-1
         ifpsco(i) = if_save(i+1) 
      enddo

c restore solution data
      call copy(vx,w1,ntot1)
      call copy(vy,w2,ntot1)
      call copy(vz,w3,ntot1)
      call copy(pr,wp,ntot2)
      do i = 1,nfldt
         call copy(t(1,1,1,1,i),wt(1,i),ntot1t)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine my_mfo_mdatav(elist,u,v,w,nel)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(lx1*ly1*lz1,1),v(lx1*ly1*lz1,1),w(lx1*ly1*lz1,1)

      real*4 buffer(1+6*lelt)

      integer e,elist(1)

      call nekgsync() ! clear outstanding message queues.

      nxyz = lx1*ly1*lz1
      n    = 2*ldim
      len  = 4 + 4*(n*lelt)   ! recv buffer size
      leo  = 4 + 4*(n*nelt) 
      nel2 = ivlsum(elist,nel)
      leo2 = 4 + 4*(n*nel2)
      ierr = 0

      ! Am I an I/O node?
      if (nid.eq.pid0) then
         j = 1
         do e=1,nel
         if (elist(e).gt.0) then
            buffer(j+0) = vlmin(u(1,e),nxyz) 
            buffer(j+1) = vlmax(u(1,e),nxyz)
            buffer(j+2) = vlmin(v(1,e),nxyz) 
            buffer(j+3) = vlmax(v(1,e),nxyz)
            j = j + 4
            if(if3d) then
              buffer(j+0) = vlmin(w(1,e),nxyz) 
              buffer(j+1) = vlmax(w(1,e),nxyz)
              j = j + 2
            endif
         endif
         enddo

         ! write out my data
c         nout = n*nel
         nout = n*nel2
         if(ierr.eq.0) then
           if(ifmpiio) then
             call byte_write_mpi(buffer,nout,-1,ifh_mbyte,ierr)
           else
             call byte_write(buffer,nout,ierr)
           endif
         endif

         ! write out the data of my childs
         idum  = 1
         do k=pid0+1,pid1
            mtype = k
            call csend(mtype,idum,4,k,0)           ! handshake
            call crecv(mtype,buffer,len)
            inelp = buffer(1)
            nout  = n*inelp
            if(ierr.eq.0) then 
              if(ifmpiio) then 
                call byte_write_mpi(buffer(2),nout,-1,ifh_mbyte,ierr)
              else
                call byte_write(buffer(2),nout,ierr)
              endif
            endif
         enddo
      else
         j = 1
c         buffer(j) = nel
         buffer(j) = nel2
         j = j + 1
         do e=1,nel
         if (elist(e).gt.0) then
            buffer(j+0) = vlmin(u(1,e),nxyz) 
            buffer(j+1) = vlmax(u(1,e),nxyz)
            buffer(j+2) = vlmin(v(1,e),nxyz) 
            buffer(j+3) = vlmax(v(1,e),nxyz)
            j = j + 4
            if(n.eq.6) then
              buffer(j+0) = vlmin(w(1,e),nxyz) 
              buffer(j+1) = vlmax(w(1,e),nxyz)
              j = j + 2
            endif
         endif
         enddo

         ! send my data to my pararent I/O node
         mtype = nid
         call crecv(mtype,idum,4)                ! hand-shake
         call csend(mtype,buffer,leo2,pid0,0)     ! u4 :=: u8
      endif

      call err_chk(ierr,'Error writing data to .f00 in mfo_mdatav. $')

      return
      end
c-----------------------------------------------------------------------
      subroutine my_mfo_mdatas(elist,u,nel)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(lx1*ly1*lz1,1)

      real*4 buffer(1+2*lelt)

      integer e,elist(1)

      call nekgsync() ! clear outstanding message queues.

      nxyz = lx1*ly1*lz1
      n    = 2
      len  = 4 + 4*(n*lelt)    ! recv buffer size
      leo  = 4 + 4*(n*nelt)
      nel2 = ivlsum(elist,nel)
      leo2 = 4 + 4*(n*nel2)
      ierr = 0

      ! Am I an I/O node?
      if (nid.eq.pid0) then
         j = 1
         do e=1,nel
         if (elist(e).gt.0) then
            buffer(j+0) = vlmin(u(1,e),nxyz) 
            buffer(j+1) = vlmax(u(1,e),nxyz)
            j = j + 2
         endif
         enddo

         ! write out my data
c         nout = n*nel
         nout = n*nel2
         if(ierr.eq.0) then 
           if(ifmpiio) then
             call byte_write_mpi(buffer,nout,-1,ifh_mbyte,ierr)
           else
             call byte_write(buffer,nout,ierr)
           endif
         endif

         ! write out the data of my childs
         idum  = 1
         do k=pid0+1,pid1
            mtype = k
            call csend(mtype,idum,4,k,0)           ! handshake
            call crecv(mtype,buffer,len)
            inelp = buffer(1)
            nout  = n*inelp
            if(ierr.eq.0) then 
              if(ifmpiio) then
                call byte_write_mpi(buffer(2),nout,-1,ifh_mbyte,ierr)
              else
                call byte_write(buffer(2),nout,ierr)
              endif
            endif
         enddo
      else
         j = 1
c         buffer(j) = nel
         buffer(j) = nel2
         j = j + 1
         do e=1,nel
         if (elist(e).gt.0) then
            buffer(j+0) = vlmin(u(1,e),nxyz) 
            buffer(j+1) = vlmax(u(1,e),nxyz)
            j = j + 2
         endif
         enddo

         ! send my data to my pararent I/O node
         mtype = nid
         call crecv(mtype,idum,4)                ! hand-shake
         call csend(mtype,buffer,leo2,pid0,0)     ! u4 :=: u8
      endif

      call err_chk(ierr,'Error writing data to .f00 in mfo_mdatas. $')

      return
      end
c-----------------------------------------------------------------------
      subroutine my_mfo_outs(elist,u,nel,mx,my,mz)   ! output a scalar field

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(mx,my,mz,1)

      common /SCRNS/ u4(2+lxo*lxo*lxo*2*lelt)
      real*4         u4
      real*8         u8(1+lxo*lxo*lxo*1*lelt)
      equivalence    (u4,u8)

      integer e,elist(1)

      umax = glmax(u,nel*mx*my*mz)
      umin = glmin(u,nel*mx*my*mz)
      if(nid.eq.0) write(6,'(A,2g13.5)') ' min/max:', umin,umax

      call nekgsync() ! clear outstanding message queues.
      if(mx.gt.lxo .or. my.gt.lxo .or. mz.gt.lxo) then
        if(nid.eq.0) write(6,*) 'ABORT: lxo too small'
        call exitt
      endif

      nxyz = mx*my*mz
      len  = 8 + 8*(lelt*nxyz)  ! recv buffer size
      leo  = 8 + wdsizo*(nel*nxyz)
      ntot = nxyz*nel

      nel2 = ivlsum(elist,nel)
      ntot2= nxyz*nel2
      leo2 = 8 + wdsizo*(nel2*nxyz)

      idum = 1
      ierr = 0

      if (nid.eq.pid0) then

         if (wdsizo.eq.4) then             ! 32-bit output
c             call copyx4 (u4,u,ntot)
             j = 0
             do iel = 1,nel
             if (elist(iel).gt.0) then
                call copyx4   (u4(j+1),u(1,1,1,iel),nxyz)
                j = j + nxyz
             endif
             enddo
         else
c             call copy   (u8,u,ntot)
             j = 0
             do iel = 1,nel
             if (elist(iel).gt.0) then
                call copy     (u8(j+1),u(1,1,1,iel),nxyz)
                j = j + nxyz
             endif
             enddo
         endif
c         nout = wdsizo/4 * ntot
         nout = wdsizo/4 * ntot2
         if(ierr.eq.0) then 
           if(ifmpiio) then
             call byte_write_mpi(u4,nout,-1,ifh_mbyte,ierr)
           else
             call byte_write(u4,nout,ierr)          ! u4 :=: u8
           endif
         endif

         ! write out the data of my childs
         idum  = 1
         do k=pid0+1,pid1
            mtype = k
            call csend(mtype,idum,4,k,0)       ! handshake
            call crecv(mtype,u4,len)
            nout  = wdsizo/4 * nxyz * u8(1)
            if (wdsizo.eq.4.and.ierr.eq.0) then
               if(ifmpiio) then
                 call byte_write_mpi(u4(3),nout,-1,ifh_mbyte,ierr)
               else
                 call byte_write(u4(3),nout,ierr)
               endif
            elseif(ierr.eq.0) then
               if(ifmpiio) then
                 call byte_write_mpi(u8(2),nout,-1,ifh_mbyte,ierr)
               else
                 call byte_write(u8(2),nout,ierr)
               endif
            endif
         enddo

      else

         u8(1)= nel2
         if (wdsizo.eq.4) then             ! 32-bit output
c             call copyx4 (u4(3),u,ntot2)
             j = 2
             do iel = 1,nel
             if (elist(iel).gt.0) then
                 call copyx4 (u4(j+1),u(1,1,1,iel),nxyz)
                 j = j + nxyz
             endif
             enddo
         else
c             call copy   (u8(2),u,ntot2)
             j = 1
             do iel = 1,nel
             if (elist(iel).gt.0) then
                 call copy   (u8(j+1),u(1,1,1,iel),nxyz)
                 j = j + nxyz
             endif
             enddo
         endif

         mtype = nid
         call crecv(mtype,idum,4)            ! hand-shake
         call csend(mtype,u4,leo2,pid0,0)     ! u4 :=: u8

      endif

      call err_chk(ierr,'Error writing data to .f00 in mfo_outs. $')

      return
      end
c-----------------------------------------------------------------------

      subroutine my_mfo_outv(elist,u,v,w,nel,mx,my,mz)   ! output a vector field

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(mx*my*mz,1),v(mx*my*mz,1),w(mx*my*mz,1)

      common /SCRNS/ u4(2+lxo*lxo*lxo*6*lelt)
      real*4         u4
      real*8         u8(1+lxo*lxo*lxo*3*lelt)
      equivalence    (u4,u8)

      integer e,elist(1)

      umax = glmax(u,nel*mx*my*mz)
      vmax = glmax(v,nel*mx*my*mz)
      wmax = glmax(w,nel*mx*my*mz)
      umin = glmin(u,nel*mx*my*mz)
      vmin = glmin(v,nel*mx*my*mz)
      wmin = glmin(w,nel*mx*my*mz)
      if(nid.eq.0) write(6,'(A,6g13.5)') ' min/max:', 
     $             umin,umax, vmin,vmax, wmin,wmax

      call nekgsync() ! clear outstanding message queues.
      if(mx.gt.lxo .or. my.gt.lxo .or. mz.gt.lxo) then
        if(nid.eq.0) write(6,*) 'ABORT: lxo too small'
        call exitt
      endif

      nxyz = mx*my*mz
      len  = 8 + 8*(lelt*nxyz*ldim)   ! recv buffer size (u4)
      leo  = 8 + wdsizo*(nel*nxyz*ldim)
      idum = 1
      ierr = 0
      nel2 = ivlsum(elist,nel) ! Lan
      leo2 = 8 + wdsizo*(nel2*nxyz*ldim)

      if (nid.eq.pid0) then
         j = 0 
         if (wdsizo.eq.4) then             ! 32-bit output
             do iel = 1,nel
             if (elist(iel).gt.0) then
                call copyx4   (u4(j+1),u(1,iel),nxyz)
                j = j + nxyz
                call copyx4   (u4(j+1),v(1,iel),nxyz)
                j = j + nxyz
                if(if3d) then
                  call copyx4 (u4(j+1),w(1,iel),nxyz)
                  j = j + nxyz
                endif
             endif
             enddo
         else
             do iel = 1,nel
             if (elist(iel).gt.0) then
                call copy     (u8(j+1),u(1,iel),nxyz)
                j = j + nxyz
                call copy     (u8(j+1),v(1,iel),nxyz)
                j = j + nxyz
                if(if3d) then
                  call copy   (u8(j+1),w(1,iel),nxyz)
                  j = j + nxyz
                endif
             endif
             enddo
         endif
c         nout = wdsizo/4 * ldim*nel * nxyz
         nout = wdsizo/4 * ldim*nel2 * nxyz
         if(ierr.eq.0) then 
           if(ifmpiio) then
             call byte_write_mpi(u4,nout,-1,ifh_mbyte,ierr)
           else
             call byte_write(u4,nout,ierr)          ! u4 :=: u8
           endif
         endif

         ! write out the data of my childs
         do k=pid0+1,pid1
            mtype = k
            call csend(mtype,idum,4,k,0)           ! handshake
            call crecv(mtype,u4,len)
            nout  = wdsizo/4 * ldim*nxyz * u8(1)

            if (wdsizo.eq.4.and.ierr.eq.0) then
               if(ifmpiio) then
                 call byte_write_mpi(u4(3),nout,-1,ifh_mbyte,ierr)
               else
                 call byte_write(u4(3),nout,ierr)
               endif
            elseif(ierr.eq.0) then
               if(ifmpiio) then
                 call byte_write_mpi(u8(2),nout,-1,ifh_mbyte,ierr)
               else
                 call byte_write(u8(2),nout,ierr)
               endif
            endif
         enddo
      else

c         u8(1) = nel
         u8(1) = nel2
         if (wdsizo.eq.4) then             ! 32-bit output
             j = 2
             do iel = 1,nel
             if (elist(iel).gt.0) then
                call copyx4   (u4(j+1),u(1,iel),nxyz)
                j = j + nxyz
                call copyx4   (u4(j+1),v(1,iel),nxyz)
                j = j + nxyz
                if(if3d) then
                  call copyx4 (u4(j+1),w(1,iel),nxyz)
                  j = j + nxyz
                endif
             endif
             enddo
         else
             j = 1
             do iel = 1,nel
             if (elist(iel).gt.0) then
                call copy     (u8(j+1),u(1,iel),nxyz)
                j = j + nxyz
                call copy     (u8(j+1),v(1,iel),nxyz)
                j = j + nxyz
                if(if3d) then
                  call copy   (u8(j+1),w(1,iel),nxyz)
                  j = j + nxyz
                endif
             endif
             enddo
         endif

         mtype = nid
         call crecv(mtype,idum,4)            ! hand-shake
         call csend(mtype,u4,leo2,pid0,0)     ! u4 :=: u8

      endif

      call err_chk(ierr,'Error writing data to .f00 in mfo_outv. $')
      return
      end
c-----------------------------------------------------------------------
c      subroutine mfo_write_hdr          ! write hdr, byte key, els.
c
c      include 'SIZE'
c      include 'SOLN'
c      include 'INPUT'
c      include 'PARALLEL'
c      include 'RESTART'
c      include 'TSTEP'
c      real*4 test_pattern
c      common /ctmp0/ lglist(0:lelt)
c
c      character*132 hdr
c      integer*8 ioff
c      logical if_press_mesh
c
c      call nekgsync()
c      idum = 1
c
c      if(ifmpiio) then
c        nfileoo = 1   ! all data into one file
c        nelo = nelgt
c      else
c        nfileoo = nfileo
c        if(nid.eq.pid0) then                ! how many elements to dump
c          nelo = nelt
c          do j = pid0+1,pid1
c             mtype = j
c             call csend(mtype,idum,4,j,0)   ! handshake
c             call crecv(mtype,inelp,4)
c             nelo = nelo + inelp
c          enddo
c        else
c          mtype = nid
c          call crecv(mtype,idum,4)          ! hand-shake
c          call csend(mtype,nelt,4,pid0,0)   ! u4 :=: u8
c        endif 
c      endif
c
c      ierr = 0
c      if(nid.eq.pid0) then
c
c      call blank(hdr,132)              ! write header
c      call blank(rdcode1,10)
c      i = 1
c      IF (IFXYO) THEN
c         rdcode1(i)='X'
c         i = i + 1
c      ENDIF
c      IF (IFVO) THEN
c         rdcode1(i)='U'
c         i = i + 1
c      ENDIF
c      IF (IFPO) THEN
c         rdcode1(i)='P'
c         i = i + 1
c      ENDIF
c      IF (IFTO) THEN
c         rdcode1(i)='T'
c         i = i + 1
c      ENDIF
c      IF (LDIMT.GT.1) THEN
c         NPSCALO = 0
c         do k = 1,ldimt-1
c           if(ifpsco(k)) NPSCALO = NPSCALO + 1
c         enddo
c         IF (NPSCALO.GT.0) THEN
c            rdcode1(i) = 'S'
c            WRITE(rdcode1(i+1),'(I1)') NPSCALO/10
c            WRITE(rdcode1(i+2),'(I1)') NPSCALO-(NPSCALO/10)*10
c         ENDIF
c      ENDIF
c
cc     check pressure format
c      if_press_mesh = .false.
c      if (.not.ifsplit.and.if_full_pres) if_press_mesh = .true.
c 
c      write(hdr,1) wdsizo,nxo,nyo,nzo,nelo,nelgt,time,istep,fid0,nfileoo
c     $            ,(rdcode1(i),i=1,10),p0th,if_press_mesh
c    1 format('#std',1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,
c     &       1x,i9,1x,i6,1x,i6,1x,10a,1pe15.7,1x,l1)
c
c      test_pattern = 6.54321           ! write test pattern for byte swap
c
c      if(ifmpiio) then
c        ! only rank0 (pid00) will write hdr + test_pattern
c        call byte_write_mpi(hdr,iHeaderSize/4,pid00,ifh_mbyte,ierr)
c        call byte_write_mpi(test_pattern,1,pid00,ifh_mbyte,ierr)
c      else
c        call byte_write(hdr,iHeaderSize/4,ierr)
c        call byte_write(test_pattern,1,ierr)
c      endif
c
c      endif
c
c      call err_chk(ierr,'Error writing header in mfo_write_hdr. $')
c
c      ! write global element numbering for this group
c      if(nid.eq.pid0) then
c        if(ifmpiio) then
c          ioff = iHeaderSize + 4 + nelB*isize
c          call byte_set_view (ioff,ifh_mbyte)
c          call byte_write_mpi(lglel,nelt,-1,ifh_mbyte,ierr)
c        else
c          call byte_write(lglel,nelt,ierr)
c        endif
c
c        do j = pid0+1,pid1
c           mtype = j
c           call csend(mtype,idum,4,j,0)   ! handshake
c           len = 4*(lelt+1)
c           call crecv(mtype,lglist,len)
c           if(ierr.eq.0) then
c             if(ifmpiio) then
c              call byte_write_mpi(lglist(1),lglist(0),-1,ifh_mbyte,ierr)
c             else
c              call byte_write(lglist(1),lglist(0),ierr)
c             endif
c           endif
c        enddo
c      else
c        mtype = nid
c        call crecv(mtype,idum,4)          ! hand-shake
c        
c        lglist(0) = nelt
c        call icopy(lglist(1),lglel,nelt)
c
c        len = 4*(nelt+1)
c        call csend(mtype,lglist,len,pid0,0)  
c      endif 
c
c      call err_chk(ierr,'Error writing global nums in mfo_write_hdr$')
c      return
c      end
cc-----------------------------------------------------------------------
