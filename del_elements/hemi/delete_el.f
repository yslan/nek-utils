c-----------------------------------------------------------------------
      subroutine delete_elements(elist,new_bID)
c     elist
c       if elist(e)=1, we delete e-th element
c     new_bID
c       if new_bID<0, assign boundaryID = abs(new_bID) for all new boundary faces
c       if new bID>0, fill boundaryID from the opposite face of the deleted element, then 
c                     assign the rest with new_bID
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer elist(lelt),new_bID
      integer is_connect(6,lelt),bdryID_new(6,lelt,2) 

      integer f2opp(6)
      data f2opp /3,4,1,2,6,5/

      integer e,f,fopp,enew,nxyz,n,ndel

      nxyz=lx1*ly1*lz1
      n=lx1*ly1*lz1*nelv

      ! Step 0: some chks
      if (lx1.lt.3) then
        if (nio.eq.0) then
          write(*,*)'delete_el requires lx1>=3 because I am lazy and'
          write(*,*)'I use volume array to do face communication' 
        endif
        call exitt
      endif

      if (nelv.ne.nelt) then
        if (nio.eq.0) write(*,*) 'ERR del: Not support cht case for now'
        call exitt
      endif

      call del_chk_elist(elist,nelv,ndel)
      if (ndel.eq.-1) return ! elist is empty, do nothing

      call del_chk_per(elist,nelv)


      ! Step 1: get boundaryID candicates for new boundary faces
      call izero(bdryID_new,6*lelt*2)

      if (new_bID.gt.0) then ! get from neighbors first
        do e=1,nelv
        do f=1,2*ldim
          fopp = f2opp(f)  
          bdryID_new(fopp,e,1) = boundaryID(f,e)
          bdryID_new(fopp,e,2) = boundaryIDt(f,e)
        enddo
        enddo
        call exchange_boundaryID(bdryID_new)

      elseif (new_bID.lt.0) then ! user specified value
        call ifill(bdryID_new,abs(new_bID),6*lelt*2)
      endif


      ! Step 2: check if (e,f) connects to someone
      call set_is_onnect(is_connect,elist)


      ! Step 3: backup
      do e=1,nelv
        call del_copy_el_to_cache(e,e)
      enddo


      ! Step 4: delete elements
      enew=0
      do e=1,nelv
      if (elist(e).eq.0) then ! keep this element

        enew = enew + 1
        call del_copy_el_from_cache(e,enew)

        do f=1,2*ldim
        if (is_connect(f,e).eq.0) then ! no neighbor, assign new bdryID

          if (boundaryID(f,enew).eq.0) then  ! exclude bdry
            if (bdryID_new(f,e,1).ne.0) then ! if it has candicate
              boundaryID(f,enew) = bdryID_new(f,e,1)
            else
              boundaryID(f,enew) = new_bID 
            endif
          endif

          if (boundaryIDt(f,enew).eq.0) then ! exclude bdry
            if (bdryID_new(f,e,2).ne.0) then
              boundaryIDt(f,enew) = bdryID_new(f,e,2)
            else
              boundaryIDt(f,enew) = new_bID 
            endif
          endif

        endif ! neighbor
        enddo ! f

      endif ! if elist(e).eq.0
      enddo ! e

      ! Step 4: update nelg
      call update_nels(enew) ! TODO: thos won't work for DPROMAP


      ! Step 5: update vertex for con/co2
ccc      call del_set_uniq_vtx(vertex)! require special parRSB

      return
      end
c-----------------------------------------------------------------------
      subroutine del_chk_elist(elist,nel,ndel)
      implicit none
      include 'SIZE'
      integer elist(1),e,nel,ndel,iglmin,iglmax,iglsum,itmp1,itmp2

c     fool proof
      do e=1,nel
        if(elist(e).gt.0) then
          elist(e)=1
        else
          elist(e)=0
        endif
      enddo

      itmp1=iglmax(elist,nel)
      itmp2=iglmin(elist,nel)

      if (itmp1.eq.0.AND.itmp2.eq.0) then
        if (nio.eq.0) write(6,*)'del: null elist, do nothing'
        ndel=0
        return
      endif
      if (itmp1.lt.0.OR.itmp2.gt.1) then
        if(nio.eq.0)write(6,*)'ERR del: elist out of range',itmp1,itmp2
        call exitt
      endif      

      ndel = iglsum(elist,nel)
      return
      end
c-----------------------------------------------------------------------
      subroutine del_chk_per(elist,nel)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer e,f,nel,elist(lelt),nfail,iglsum

      nfail = 0

      do e=1,nel
      if (elist(e).eq.1) then

        do f=1,2*ldim       
          if (cbc(f,e,1).eq.'P  '.OR.cbc(f,e,1).eq.'p  ') then
            nfail = nfail+1   
          endif
          if (cbc(f,e,2).eq.'P  '.OR.cbc(f,e,2).eq.'p  ') then
            nfail = nfail+1
          endif
        enddo

      endif
      enddo

      nfail = iglsum(nfail,1)
      if (nfail.gt.0) then
        if (nio.eq.0) 
     $    write(*,*)'ERR del: delete e with per BC is not supported'
        call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine exchange_boundaryID(bdryID)
      implicit none
      include 'SIZE'

      integer ix,iy,iz,e,f,bdryID(6,lelt,2)
      integer idd,ix_fctr(6),iy_fctr(6),iz_fctr(6) ! get face center index
      parameter(idd=ldim-1)
      data ix_fctr /  2,lx1,  2,  1,  2,  2/
      data iy_fctr /  1,  2,ly1,  2,  2,  2/
      data iz_fctr /idd,idd,idd,idd,  1,lz1/

      real tmp1(lx1,ly1,lz1,lelt),tmp2(lx1,ly1,lz1,lelt)

      do e=1,nelv
      do f=1,2*ldim
        call facev(tmp1,e,f,real(bdryID(f,e,1)),lx1,ly1,lz1)
        call facev(tmp2,e,f,real(bdryID(f,e,2)),lx1,ly1,lz1)
      enddo
      enddo

      call dssum(tmp1,lx1,ly1,lz1)
      call dssum(tmp2,lx1,ly1,lz1)

      do e=1,nelv
      do f=1,2*ldim
        ix=ix_fctr(f)
        iy=iy_fctr(f)
        iz=iz_fctr(f)

        bdryID(f,e,1) = int(tmp1(ix,iy,iz,e)) - bdryID(f,e,1)
        bdryID(f,e,2) = int(tmp2(ix,iy,iz,e)) - bdryID(f,e,2)
      enddo
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine set_is_onnect(is_connect,elist)
      implicit none
      include 'SIZE'

      integer is_connect(6,lelt),elist(lelt)
      real tmp(lx1,ly1,lz1,lelt)

      integer ix,iy,iz,e,f,nxyz
      integer idd,ix_fctr(6),iy_fctr(6),iz_fctr(6) ! get face center index
      parameter(idd=ldim-1)
      data ix_fctr /  2,lx1,  2,  1,  2,  2/
      data iy_fctr /  1,  2,ly1,  2,  2,  2/
      data iz_fctr /idd,idd,idd,idd,  1,lz1/      

      nxyz=lx1*ly1*lz1
      call ione(is_connect,6*lelt) ! start from 1
      call rone(tmp,lx1*ly1*lz1*lelt) ! start from 1
      call dssum(tmp,lx1,ly1,lz1)     ! connected face have 2

      ! detect if (f,e) is on boundary
      do e=1,nelv

        do f=1,2*ldim
          ix=ix_fctr(f)
          iy=iy_fctr(f)
          iz=iz_fctr(f)
  
          if (tmp(ix,iy,iz,e).lt.1.5) then   ! if on bdry
            call facev(tmp,e,f,0.0,lx1,ly1,lz1)
            is_connect(f,e)=0 
          else                               ! back to 1
            call facev(tmp,e,f,1.0,lx1,ly1,lz1)
          endif  
        enddo

        if (elist(e).eq.1) call rzero(tmp(1,1,1,e),nxyz) ! this e will be deleted

      enddo

      call dssum(tmp,lx1,ly1,lz1)

      do e=1,nelv
      do f=1,2*ldim
        ix=ix_fctr(f)
        iy=iy_fctr(f)
        iz=iz_fctr(f)
        if (tmp(ix,iy,iz,e).lt.1.5) is_connect(f,e)=0
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine del_set_nel_shift(nel_shift,nel) ! copy from my refine.f
c     this helps to re-assign global element id
      implicit none
      include 'SIZE'

      integer mtype,inid,nel_tmp,nel_acc
      integer nel_shift,nel

      integer nidd,np,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidd,np,nekcomm,nekgroup,nekreal

      call nekgsync()

      nel_shift = 0
      nel_tmp = 0
      nel_acc = 0

      if (nid.eq.0) nel_acc = nel_acc + nel

      !    0: send nel_acc = sum_{i,inid-1} nel
      ! inid: send nel to nid=0
      !    0: nel_acc += nel(inid)
      if (nid.eq.0) then
        do inid=1,np-1
          mtype = inid
          call csend(mtype,nel_acc,4,inid,0)
          call crecv(mtype,nel_tmp,4)
          nel_acc = nel_acc + nel_tmp
        enddo
      else
        mtype = nid
        call crecv(mtype,nel_shift,4)
        call csend(mtype,nel,4,0,0)
      endif

      call nekgsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine update_nels(enew) ! copy from mine refine.f
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer e,enow,e1,e2,iglsum
      integer mfield,nfldt
      integer nelv_shift,enew

      nelv=enew
      nelt=nelv
      nelgv=iglsum(nelv,1)
      nelgt=iglsum(nelt,1) ! outpost needs nelgt      

      ! TODO: this looks messy
      MFIELD=2
      IF (IFFLOW) MFIELD=1
      IF (IFMVBD) MFIELD=0
c     Set up TEMPORARY value for NFIELD - NFLDT
      NFLDT = 1
      IF (IFHEAT) NFLDT = 2 + NPSCAL
      DO 1200 IFIELD=MFIELD,NFLDT
         IF (IFTMSH(IFIELD)) THEN
            NELG(IFIELD)      = NELGT
         ELSE
            NELG(IFIELD)      = NELGV
         ENDIF
 1200 CONTINUE

      call del_set_nel_shift(nelv_shift,nelv)

      do e=1,nelgv
        gllel(e) = 0
        gllnid(e) = 0
      enddo
      do enow=1,nelv
        e1 = enow
        e2 = nelv_shift + enow

        lglel(e1) = e2
        gllel(e2) = e1
        gllnid(e2) = nid
      enddo
      do e=1,nelgv ! TODO use igop?
        gllel(e) = iglsum(gllel(e),1)
        gllnid(e) = iglsum(gllnid(e),1)
      enddo

      return
      end
cccc-----------------------------------------------------------------------
ccc      subroutine del_set_uniq_vtx(vertex) ! require special parRSB
ccc      implicit none
ccc      include 'SIZE'
ccc
ccc      integer lxv,lyv,lzv
ccc      parameter(lxv=2,lyv=2,lzv=ldim-1)
ccc
ccc      integer*8 vertex((2**ldim)*lelt), vertexu((2**ldim)*lelt)
ccc      integer nv,ierr,ialgo,lb,iglmax
ccc
ccc      integer nidd,np,nekcomm,nekgroup,nekreal
ccc      common /nekmpi/ nidd,np,nekcomm,nekgroup,nekreal
ccc
ccc      ialgo = 0 ! 0=bin sort, 1=hyprecube sort
ccc      lb = 0    ! 1=load balanced
ccc
ccc      nv = lxv*lyv*lzv*nelv
ccc
ccc      call FPSORT_INT8(vertexu,nv,vertex,ialgo,lb,nekcomm,ierr)
ccc      call i8copy(vertex,vertexu,nv)
ccc
ccc      return
ccc      end
cccc-----------------------------------------------------------------------
      subroutine del_dump_mesh(s3,imid,ifbinary) ! copy from refine.f
c     imid = 0  ! No midside node defs
c     imid = 1  ! Midside defs where current curve sides don't exist
c     imid = 2  ! All nontrivial midside node defs
      implicit none
      integer ierr,imid
      logical ifbinary
      character*3 s3

ccc REMOVE extra dependencies for portability...
ccc      ! chk mesh ! my_errchk.f
ccc      call my_mesh_metric
ccc      call my_mshchk_vb(s3//'   ',.false.,ierr)
ccc      call my_cbc_chk(s3)
ccc
ccc      ! dump mesh ! my_errchk.f
ccc      if (ifbinary) then
ccc        call my_gen_re2(s3,imid)
cccc       call del_gen_co2(s3) ! need special parrsb
ccc      else
ccc        call my_gen_rea(s3,imid)
cccc       call del_gen_con(s3)
ccc      endif

      if (ifbinary) then
        call gen_re2(imid)
      else
        call gen_rea(imid)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine del_gen_co2(s3) ! from refine.f
      implicit none
      include 'SIZE'
      include 'PARALLEL'
      character*3 s3

      integer*8 vertex
      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only, TODO check set up

      integer eg,e,iprog,ii,nlen,mtype,inid,idum,ierr,nv,iglsum
      integer iwrk(1+8)

      character*132 hdr
      character*5   version
      real*4 test
      data   test  / 6.54321 /

      version = '#v001'
      nv = 2**ldim

      ierr = 0
      if (nid.eq.0) then
        ! Get filenames
        call byte_open('newco2'//s3//'.co2'//char(0),ierr)

        ! write header
        call blank(hdr,132)
        write(hdr,1) version,nelgt,nelgv,nv
        write(6,*)'del_gen_co2 ',hdr
    1   format(a5,3i12)
        call byte_write(hdr,132/4,ierr)
        call byte_write(test,1,ierr) ! write the endian discriminator
      endif
      ierr = iglsum(ierr,1)
      if (ierr.gt.0) call exitti('dump co2, open file$',ierr)

      iprog=1
      do eg=1,nelgt
        inid = gllnid(eg)

        if (nid.eq.0) then
          if (inid.eq.0) then
            e = gllel(eg)
            ii = (e-1)*nv+1
            iwrk(1) = eg
            call icopy84(iwrk(2),vertex(ii),nv)
          else ! recv
            mtype= inid
            nlen = (nv+1)*4
            call csend(mtype,idum,4,inid) ! handshake
            call crecv(mtype,iwrk,nlen)
          endif

          ! Dump
          call byte_write(iwrk,nv+1,ierr)
          if(ierr.gt.0) call exitti('invalid id$',ierr)

          ! print progress
          if (eg*10.ge.iprog*nelgt) then
             write(*,3)'  progress = ',iprog*10,'%',eg,nelgt

             iprog=iprog+1
          endif
 
        else 
          if (inid.eq.nid) then ! send to nid=0
            e = gllel(eg)
            ii = (e-1)*nv+1
            iwrk(1) = eg
            call icopy84(iwrk(2),vertex(ii),nv)
  
            mtype= inid
            nlen = (nv+1)*4
            call crecv(mtype,idum,4)
            call csend(mtype,iwrk,nlen,0)
          endif 
        endif
      enddo
      ierr = iglsum(ierr,1)
      if (ierr.gt.0) call exitti('dump co2, dump file$',ierr)

    2 format(9i12)
    3 format(a13,i3a1,2i12)

      if (nid.eq.0) call byte_close(ierr)
      ierr = iglsum(ierr,1)
      if (ierr.gt.0) call exitti('dump co2, close file$',ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine del_gen_con(s3) ! from refine.f
      implicit none
      include 'SIZE'
      include 'PARALLEL'
      character*3 s3

      integer*8 vertex
      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only, TODO check set up

      integer eg,e,iprog,ii,nlen,mtype,inid,idum,nv
      integer iwrk(1+8)

      character*5   version
      real*4 test
      data   test  / 6.54321 /

      version = '#v001'
      nv = 2**ldim

      if (nid.eq.0) then
        ! Get filenames
        open (unit=29,file='newcon'//s3//'.con',status='unknown')

        ! write header
        write(29,1) version,nelgt,nelgv,nv
        write(6,*)'del_gen_con ',version,nelgt,nelgv,nv
    1   format(a5,3i12)
      endif

      iprog=1
      do eg=1,nelgt
        inid = gllnid(eg)

        if (nid.eq.0) then

          if (inid.eq.0) then 
            e = gllel(eg)
            ii = (e-1)*nv+1
            iwrk(1) = eg
            call icopy84(iwrk(2),vertex(ii),nv)

          else ! recv
            mtype= inid
            nlen = (nv+1)*4
            call csend(mtype,idum,4,inid) ! handshake
            call crecv(mtype,iwrk,nlen)
          endif

          ! Dump 
          write(29,2) (iwrk(ii), ii=1,nv+1)

          ! print progress
          if (eg*10.ge.iprog*nelgt) then
             write(*,3)'  progress = ',iprog*10,'%',eg,nelgt
             iprog=iprog+1
          endif

        elseif (inid.eq.nid.AND.nid.ne.0) then ! send to nid=0
          e = gllel(eg)
          ii = (e-1)*nv+1
          iwrk(1) = eg
          call icopy84(iwrk(2),vertex(ii),nv)

          mtype= inid
          nlen = (nv+1)*4
          call crecv(mtype,idum,4)
          call csend(mtype,iwrk,nlen,0)
        endif
      enddo

    2 format(9i12)
    3 format(a13,i3a1,2i12)

      if (nid.eq.0) close(unit=29)

      return
      end
c-----------------------------------------------------------------------
      subroutine del_boundaryID_to_bc5
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer e,f
      do e=1,nelv
      do f=1,2*ldim
        bc(5,f,e,1) = real(boundaryID(f,e))
        bc(5,f,e,2) = real(boundaryIDt(f,e))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine del_copy_el_to_cache(ei,eo)
c     copy nek variables and save into cache
c          nek(element ei) -> cache(element eo)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer*8 vertex
      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only

      real xm1_bak(lx1,ly1,lz1,lelt)
     $   , ym1_bak(lx1,ly1,lz1,lelt)
     $   , zm1_bak(lx1,ly1,lz1,lelt)
     $   , vx_bak (lx1,ly1,lz1,lelv)
     $   , vy_bak (lx1,ly1,lz1,lelv)
     $   , vz_bak (lx1,ly1,lz1,lelv)
     $   , pr_bak (lx2,ly2,lz2,lelv)
     $   , t_bak  (lx1,ly1,lz1,lelt,ldimt)

      real bc_bak(5,6,lelt,2),curve_bak(6,12,lelt)

      common /rdel_cache/ xm1_bak,ym1_bak,zm1_bak
     $                  , vx_bak,vy_bak,vz_bak,pr_bak,t_bak
     $                  , bc_bak,curve_bak

      integer bdryID_bak(6,lelt,2)
      common /idel_cache/ bdryID_bak

      character*3 cbc_bak(6,lelt,2)
      character*1 ccurve_bak(12,lelt)
      common /cdel_cache/ cbc_bak,ccurve_bak

      integer*8 vertex_bak((2**ldim)*lelt)
      common /i8del_cache/ vertex_bak

      integer i,f,ii,jj,ifld,ei,eo

      ! volume
      do i=1,lx1*ly1*lz1
        xm1_bak(i,1,1,eo) = xm1(i,1,1,ei)
        ym1_bak(i,1,1,eo) = ym1(i,1,1,ei)
        zm1_bak(i,1,1,eo) = zm1(i,1,1,ei)

        vx_bak(i,1,1,eo) = vx(i,1,1,ei)
        vy_bak(i,1,1,eo) = vy(i,1,1,ei)
        vz_bak(i,1,1,eo) = vz(i,1,1,ei)

        do ifld=2,ldimt+1
          t_bak(i,1,1,eo,ifld-1) = t(i,1,1,ei,ifld-1)
        enddo 
      enddo

      ! volume lx2
      do i=1,lx2*ly2*lz2
        pr_bak(i,1,1,eo) = pr(i,1,1,ei)
      enddo

      ! face
      do f=1,2*ldim
        bdryID_bak(f,eo,1) = boundaryID(f,ei)
        bdryID_bak(f,eo,2) = boundaryIDt(f,ei)

        cbc_bak(f,eo,1) = cbc(f,ei,1)
        cbc_bak(f,eo,2) = cbc(f,ei,2)

        call copy(bc_bak(1,f,eo,1),bc(1,f,ei,1),5)
        call copy(bc_bak(1,f,eo,2),bc(1,f,ei,2),5)
      enddo

      ! edges
      do i=1,12
        call copy(curve_bak(1,i,eo),curve(1,i,ei),6)
        ccurve_bak(i,eo) = ccurve(i,ei)
      enddo

      ! vertex
      do i=1,2**ldim
        ii=(ei-1)*(2**ldim) + i
        jj=(eo-1)*(2**ldim) + i
        vertex_bak(jj) = vertex(ii)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine del_copy_el_from_cache(ei,eo)
c     restore nek variables from cache
c          nek(element eo) <- cache(element ei)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer*8 vertex
      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only

      real xm1_bak(lx1,ly1,lz1,lelt)
     $   , ym1_bak(lx1,ly1,lz1,lelt)
     $   , zm1_bak(lx1,ly1,lz1,lelt)
     $   , vx_bak (lx1,ly1,lz1,lelv)
     $   , vy_bak (lx1,ly1,lz1,lelv)
     $   , vz_bak (lx1,ly1,lz1,lelv)
     $   , pr_bak (lx2,ly2,lz2,lelv)
     $   , t_bak  (lx1,ly1,lz1,lelt,ldimt)

      real bc_bak(5,6,lelt,2),curve_bak(6,12,lelt)

      common /rdel_cache/ xm1_bak,ym1_bak,zm1_bak
     $                  , vx_bak,vy_bak,vz_bak,pr_bak,t_bak
     $                  , bc_bak,curve_bak

      integer bdryID_bak(6,lelt,2)
      common /idel_cache/ bdryID_bak

      character*3 cbc_bak(6,lelt,2)
      character*1 ccurve_bak(12,lelt)
      common /cdel_cache/ cbc_bak,ccurve_bak

      integer*8 vertex_bak((2**ldim)*lelt)
      common /i8del_cache/ vertex_bak

      integer i,f,ii,jj,ifld,ei,eo

      ! volume
      do i=1,lx1*ly1*lz1
        xm1(i,1,1,eo) = xm1_bak(i,1,1,ei)
        ym1(i,1,1,eo) = ym1_bak(i,1,1,ei)
        zm1(i,1,1,eo) = zm1_bak(i,1,1,ei)

        vx(i,1,1,eo) = vx_bak(i,1,1,ei)
        vy(i,1,1,eo) = vy_bak(i,1,1,ei)
        vz(i,1,1,eo) = vz_bak(i,1,1,ei)

        do ifld=2,ldimt+1
          t(i,1,1,eo,ifld-1) = t_bak(i,1,1,ei,ifld-1)
        enddo 
      enddo

      ! volume lx2
      do i=1,lx2*ly2*lz2
        pr(i,1,1,eo) = pr_bak(i,1,1,ei)
      enddo

      ! face
      do f=1,2*ldim
        boundaryID(f,eo)  = bdryID_bak(f,ei,1)
        boundaryIDt(f,eo) = bdryID_bak(f,ei,2) 

        cbc(f,eo,1) = cbc_bak(f,ei,1)
        cbc(f,eo,2) = cbc_bak(f,ei,2)

        call copy(bc(1,f,eo,1),bc_bak(1,f,ei,1),5)
        call copy(bc(1,f,eo,2),bc_bak(1,f,ei,2),5)
      enddo

      ! edges
      do i=1,12
        call copy(curve(1,i,eo),curve_bak(1,i,ei),6)
        ccurve(i,eo) = ccurve_bak(i,ei)
      enddo

      ! vertex
      do i=1,2**ldim
        ii=(ei-1)*(2**ldim) + i
        jj=(eo-1)*(2**ldim) + i
        vertex(jj) = vertex_bak(ii)
      enddo

      return   
      end
c-----------------------------------------------------------------------

