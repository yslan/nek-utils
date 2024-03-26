c     In userdat2:
c        call rfn_backup_boundary_id(1)
c        ! set bdry for simulation
c     In userchk:
c        call rfn_split
c        call rfn_dump_mesh
c     TODO:
c        cbc ifield=2
c        if co2/ma2
c-----------------------------------------------------------------------
      subroutine rfn_split(cbc3,bID,Nref,Rref,Rlim)
      implicit none
      include 'SIZE'
      include 'PARALLEL' ! nelgv
      include 'REFINE'
      character*3 cbc3
      integer Nref,ie_ref(lelt),bID
      real Rref(lref-1),Rlim(4,lref-1)

      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only
      integer*8 vertex,nvtx,i8glmax

      if (Nref.le.1) return

      ! initialize global variables, and check inputs
      call rfn_init(ie_ref,cbc3,bID,Nref,Rref,Rlim)

      ! copy coordinates to local work array (rotated)
      ! rotate s.t. face 1 (2D) or face 5 (3D) will always face the target bdry
      if (ldim.eq.2) then
        call rfn_rotate2d_in(ie_ref) 
      else
        call rfn_rotate3d_in(ie_ref) 
      endif

      ! get #uniq vtx (need vertex, vertexc)
      call rfn_set_nvtx

      ! get nel_lv_shift (need nel_lv)
      call rfn_set_nel_lv_shift

      ! Main refine
      if (ldim.eq.2) then
        call rfn_refine2d(ie_ref,Nref,Rref,Rlim)
      else
        call rfn_refine3d(ie_ref,Nref,Rref,Rlim)
      endif

      ! copy coordinates back to xm1, ym1, zm1 (rotate back)
      if (ldim.eq.2) then
        call rfn_rotate2d_out(ie_ref)
      else
        call rfn_rotate3d_out(ie_ref)
      endif

c     TODO: ambiguous future plan
c      call rfn_interp_sol ! only works when limiter 2 is turning off
c      call rfn_upd_topo   ! take care gllnid, refresh gs setup
c      call rfn_part_update! redistribute elements

      nvtx = i8glmax(vertex,lxv*lyv*lzv*nelt)
      if (nid.eq.0) then
        write(*,26)'rfn REFINE END,  E_new=',nelgv
     $            ,', vtx: old=',nvtx0,', new=',nvtx,', lyr=',nvtx_lv
        write(*,'(a)')'rfn'
      endif

   26 format(4(ai9))
      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_extrude(cbc3,bID,igeo,icht) ! only extrude 1 layer
c     requires usr_extrude_pj in usr file
c     igeo is for distinguish difference faces for extrude two surface at the same time
      implicit none
      include 'SIZE'
      include 'PARALLEL'
      include 'REFINE'

      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only
      integer*8 vertex,nvtx,i8glmax

      character*3 cbc3
      integer ie_ref(lelt),igeo,bID,icht
      real Rref(lref-1),Rlim(4,lref-1) ! un-used dummy variable

      call rfn_init(ie_ref,cbc3,bID,1,Rref,Rlim)

      if (ldim.eq.2) then
        call rfn_rotate2d_in(ie_ref)
      else
        call rfn_rotate3d_in(ie_ref)
      endif

      call rfn_set_nvtx

      call rfn_set_nel_lv_shift

      if (ldim.eq.2) then
        call rfn_extrude_lyr_2d(ie_ref,igeo,icht)
      else
        call rfn_extrude_lyr_3d(ie_ref,igeo,icht)
      endif

      if (ldim.eq.2) then
        call rfn_rotate2d_out(ie_ref)
      else
        call rfn_rotate3d_out(ie_ref)
      endif


      nvtx = i8glmax(vertex,lxv*lyv*lzv*nelt)
      if (nid.eq.0) then
        write(*,26)'rfn EXTRUDE END, E_new=',nelgv
     $            ,', vtx: old=',nvtx0,', new=',nvtx,', lyr=',nvtx_lv
        write(*,'(a)')'rfn'
      endif

   26 format(4(ai9))
      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_extrude_pj2d(xf,yf,zf,ig,iside,ie)
c     Given xf_i on face 3 (2D), project it onto xf_o for face 1 (2D)
      implicit none
      include 'SIZE'
      real xf(lx1,ly1,lz1),yf(lx1,ly1,lz1),zf(lx1,ly1,lz1)
      integer ig,iside,ie

      call usr_extrude_pj(xf(1,  1,1),yf(1,  1,1),zf(1,  1,1)
     $                   ,xf(1,ly1,1),yf(1,ly1,1),zf(1,ly1,1),ig
     $                   ,iside,ie)

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_extrude_pj3d(xf,yf,zf,ig,iside,ie)
c     Given xf_i on face 6 (3D), project it onto xf_o for face 5 (3D)
      implicit none
      include 'SIZE'
      real xf(lx1,ly1,lz1),yf(lx1,ly1,lz1),zf(lx1,ly1,lz1)
      integer ig,iside,ie
      
      call usr_extrude_pj(xf(1,1,  1),yf(1,1,  1),zf(1,1,  1)
     $                   ,xf(1,1,lz1),yf(1,1,lz1),zf(1,1,lz1),ig
     $                   ,iside,ie)

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_init(ie_ref,cbc3,bID,Nref,Rref,Rlim)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'
      integer ifld,e,f
      character*3 cbc3
      integer Nref,ie_ref(lelt),ie_test(lelt),bID
      real Rref(lref-1),Rlim(4,lref-1)

      real tol,rmn,rmx
      integer ierr,ilv,ntmp,iglmax,iglsum,ivlsum,ie,iface

c     Initialize global variables
      call rzero(xmc,lx1*ly1*lz1*lelt) 
      call rzero(ymc,lx1*ly1*lz1*lelt) 
      call rzero(zmc,lx1*ly1*lz1*lelt) 
      call i8zero(vertexc,lxv*lyv*lzv*lelt)
      call i8zero(vertexu,lxv*lzv*lelt)

      if (bID.gt.0) then
        do ifld=1,2
        do e=1,lelt
        do f=1,6
          bID_wk(f,e,ifld) = 0
        enddo
        enddo
        enddo
      else
        do ifld=1,2
        do e=1,lelt
        do f=1,6
          cbc_wk(f,e,ifld) = 'E  '
        enddo
        enddo
        enddo
      endif

c     Check inputs
      if (Nref.gt.1) then
        ! chk Nref/Rref
        if (Nref.gt.lref) call exitti('Err rfn: increase lref>=$',Nref)

        ierr=0
        do ilv=1,Nref-2
          Rref(ilv)=max(Rref(ilv),0.0) ! Rref in (0,1)
          Rref(ilv)=min(Rref(ilv),1.0)
          if (Rref(ilv).ge.Rref(ilv+1)) ierr=1
        enddo
        if (ierr.eq.1) then
          if(nio.eq.0) then
            write(*,*) 'rfn Rref = ',(Rref(ilv),ilv=1,Nref-1)
            write(*,*) 'Err rfn: Rref must be monotonic increasing'
          endif
          call exitt
        else
          if(nio.eq.0) write(*,25)'rfn Rref='
     $                            ,(Rref(ilv),ilv=1,min(Nref-1,10))
        endif


        ! chk Rlim, make sure Rlim(2,ilv) = max > Rlim(1,ilv) = min
        ! if Rlim is 0, then no limiter
        ierr=0
        tol=1.e-10 ! TODO: rel tol from min(edge)?
        do ilv=1,Nref-1
          rmn = Rlim(1,ilv)
          rmx = Rlim(2,ilv)
          if (abs(rmn).gt.tol.AND.abs(rmx).gt.tol) then ! TODO, chk monotone
            if(rmx.le.rmn) ierr = ierr + 1
            if(rmn.lt.0.0.OR.rmx.lt.0.0) ierr = ierr + 1000
          endif
          rmn = Rlim(3,ilv)
          rmx = Rlim(4,ilv)
          if (abs(rmn).gt.tol.AND.abs(rmx).gt.tol) then ! TODO, chk monotone
            if(rmx.le.rmn) ierr = ierr + 1
            if(rmn.lt.0.0.OR.rmx.lt.0.0) ierr = ierr + 1000
          endif
        enddo
        if (ierr.gt.0) call exitti('Err rfn: Rlim min>max, or <0$',ierr)

      endif

      ! tag e that needs refine, assign face id to ie_ref
      call izero(ie_ref,lelv)
      call izero(ie_test,lelv)
      if (bID.gt.0) then
        do ie=1,nelt
        do iface=1,ldim*2
          if (boundaryID(iface,ie).eq.bID) then
            ie_ref(ie) = iface
            ie_test(ie) = ie_test(ie) + 1
          endif
        enddo
        enddo
      else
        do ie=1,nelt
        do iface=1,ldim*2
          if (cbc(iface,ie,1).eq.cbc3) then
            ie_ref(ie) = iface
            ie_test(ie) = ie_test(ie) + 1
          endif
        enddo
        enddo
      endif
     
      ! chk e, and size
      ntmp = iglmax(ie_test,nelt)
      if (ntmp.gt.1) call exitti('Err rfn: multi-refine per e$',ntmp)
      ntmp = iglsum(ie_test,nelt)
      if (ntmp.eq.0) call exitti('Err rfn: no e found$',ntmp)
      ntmp = ivlsum(ie_test,nelt)*(Nref-1) ! local new elements counts
      if (ntmp+nelt.gt.lelt) 
     $  call exitti('Err rfn: lelv too small, need$',ntmp)

      ! Not supports periodic BC for now FIXME
      ierr=0
      do ie=1,nelt
      do iface=1,ldim*2
        if (cbc(iface,ie,1).eq.'P  '.or.cbc(iface,ie,1).eq.'p  ') then
          ierr=1
        endif
      enddo
      enddo
      ierr=iglsum(ierr,1)
      if (ierr.gt.0) 
     $   call exitti('Err rfn: periodic BC not supported yet$',ierr)
      


c     Fill/backup global variables
      nelv0 = nelt
      nelt0 = nelt 
      nel_lv = 0
      do e=1,nelt
        if(ie_ref(e).gt.0) nel_lv = nel_lv + 1
      enddo

      nelgv0 = nelgv
      nelgt0 = nelgt
      nelg_lv = iglsum(nel_lv,1)

      if (Nref.eq.1) then
        if (nio.eq.0) write(*,21)'rfn EXTRUDE START,   E='
     $                          ,nelgv0,', E_lyr=',nelg_lv
      else
        if (nio.eq.0) write(*,22)'rfn REFINE START,    E='
     $                          ,nelgv0,', E_lyr',nelg_lv,'Nref',Nref
      endif

   21 format(2(ai9))
   22 format(2(ai9),ai3)
   25 format(a,10g13.5)
      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_backup_boundary_id(istore)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'
      integer e,f,istore
    
      if (istore.gt.2) call exitti('cbc max slots=2$',istore)

      do e=1,nelt
      do f=1,2*ldim
        bID_bkup(f,e,istore) = boundaryID(f,e)
      enddo
      enddo

      do e=1,nelt
      do f=1,2*ldim
        cbc_bkup(f,e,istore) = cbc(f,e,1)
      enddo
      enddo

      icbc_bkup(istore) = 1

      if (nio.eq.0) write(*,24)'rfn cbc backup to slot, id=',istore
     $                        ,', status=',icbc_bkup(1),icbc_bkup(2)

   24 format(a,i4,a,2i4)

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_recover_boundary_id(istore)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'
      integer e,f,istore

      if (istore.gt.2) call exitti('cbc max slots=2$',istore)
      if (icbc_bkup(istore).eq.0) then
        if (nio.eq.0) write(*,*)'rfn cbc slot not filled', istore
     $                         ,icbc_bkup(1),icbc_bkup(2)
        return ! the slot is not used, don't overwrite
      endif
      
      do e=1,nelt
      do f=1,2*ldim
        boundaryID(f,e) = bID_bkup(f,e,istore)
      enddo
      enddo

      do e=1,nelt
      do f=1,2*ldim
        cbc(f,e,1) = cbc_bkup(f,e,istore)
      enddo
      enddo

      icbc_bkup(istore) = 2

      if (nio.eq.0) write(*,24)'rfn cbc backup to slot, id=',istore
     $                        ,', status=',icbc_bkup(1),icbc_bkup(2)

   24 format(a,i4,a,2i4)

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_set_nvtx
      implicit none
      include 'SIZE'
      include 'REFINE'

      integer*8 vertex,i8glmax
      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only

      ! get uniq vertices on bdry face
      call get_uniq_vtx_fbot(vertexu,vertexc,nel_lv)

      nvtx0   = i8glmax(vertex,lxv*lyv*lzv*nelt)
      nvtx_lv = i8glmax(vertexu,lxv*lzv*nel_lv)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_uniq_vtx_fbot(vertexu,vertexc,nel)
      implicit none
      include 'SIZE'

      integer lxv,lyv,lzv
      parameter(lxv=2,lyv=2,lzv=ldim-1)

      integer*8 vertexu(lxv,lzv,lelt),vertexc(lxv,lyv,lzv,lelt)
      integer*8 vtx(lxv*lzv,lelt)
      integer nvtx,nxyzv,nxyv,nel,e,ierr,ialgo,lb,iglmax

      integer nidd,np,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidd,np,nekcomm,nekgroup,nekreal

      ialgo = 0 ! 0=bin sort, 1=hyprecube sort
      lb = 0 ! 1=load balanced

      nxyv  = lxv*lzv
      nxyzv = lxv*lyv*lzv
      nvtx  = nel*nxyv

      call i8zero(vtx,nvtx)
      call i8zero(vertexu,nvtx)

      do e=1,nel ! only take face 1 (2D) or face 5 (3D)
        call i8copy(vtx(1,e),vertexc(1,1,1,e),nxyv)
      enddo

c      call FPSORT_INT8(vertexu,nvtx,vtx,ialgo,lb,nekcomm,ierr)
      ierr = iglmax(ierr,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_set_nel_lv_shift
      implicit none
      include 'SIZE'
      include 'REFINE'

      integer mtype,inid,nel_tmp,nel_acc

      integer nidd,np,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidd,np,nekcomm,nekgroup,nekreal

      call nekgsync()

      nel_lv_shift = 0
      nel_tmp = 0
      nel_acc = 0

      if (nid.eq.0) nel_acc = nel_acc + nel_lv

      !    0: send nel_acc = sum_{i,inid-1} nel_lv
      ! inid: send nel_lv to nid=0
      !    0: nel_acc += nel_lv(inid)
      if (nid.eq.0) then
        do inid=1,np-1
          mtype = inid
          call csend(mtype,nel_acc,4,inid,0)
          call crecv(mtype,nel_tmp,4)
          nel_acc = nel_acc + nel_tmp
        enddo
      else
        mtype = nid
        call crecv(mtype,nel_lv_shift,4)
        call csend(mtype,nel_lv,4,0,0) 
      endif

      call nekgsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_rotate2d_in(ie_ref)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'
      integer ie_ref(lelt)
      integer enow,e,f,ix,iy,ixo,iyo,ii,nxyz,iglsum,i

      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only, TODO check set up
      integer*8 vertex,i8glmax

      integer irt(4,4)
      data irt / 1,2,3,4 
     $         , 2,3,4,1
     $         , 3,4,1,2
     $         , 4,1,2,3 /

      if (lx1.ne.ly1) call exitti('Error rfn rot: lx1 ne ly1$',ly1)

      nxyz = lx1*ly1*lz1
            
c      v3__F3__v4 
c       |      |  
c     F4|      |F2
c       |______|  
c      v1  F1  v2 

      ! copy into work array + rotate
      enow = 0
      do e=1,nelt
      if (ie_ref(e).gt.0) then
        enow = enow + 1

        ie_ref_wk(enow) = ie_ref(e)
        f = ie_ref(e)

        do i=1,4
           bID_wk(i,enow,1) = boundaryID(irt(i,f),e)
           bID_wk(i,enow,2) = bID_bkup(irt(i,f),e,1)
           cbc_wk(i,enow,1) = cbc(irt(i,f),e,1)
           cbc_wk(i,enow,2) = cbc_bkup(irt(i,f),e,1)
           call copy(bc_wk(1,i,enow,1), bc(1,irt(i,f),e,1), 5)
           call copy(bc_wk(1,i,enow,2), bc(1,irt(i,f),e,2), 5)
        enddo


        if (f.eq.1) then
          ! no rotate, copy only
          call copy(xmc(1,1,1,enow),xm1(1,1,1,e),nxyz)
          call copy(ymc(1,1,1,enow),ym1(1,1,1,e),nxyz)

          call i8copy(vertexc(1,1,1,enow),vertex((e-1)*lxyzv+1),lxyzv)


        elseif (f.eq.2) then
          do iy=1,lx1
          do ix=1,lx1
            ixo = ly1-iy+1 
            iyo = ix
            xmc(ix,iy,1,enow) = xm1(ixo,iyo,1,e)
            ymc(ix,iy,1,enow) = ym1(ixo,iyo,1,e)
          enddo
          enddo

          do iy=1,lxv
          do ix=1,lyv
            ixo = lyv-iy+1
            iyo = ix
            ii = (e-1)*(2**ldim) + (iyo-1)*lxv + ixo
            vertexc(ix,iy,1,enow) = vertex(ii)
          enddo
          enddo


        elseif (f.eq.3) then
          do iy=1,lx1
          do ix=1,lx1
            ixo = lx1-ix+1 
            iyo = ly1-iy+1
            xmc(ix,iy,1,enow) = xm1(ixo,iyo,1,e)
            ymc(ix,iy,1,enow) = ym1(ixo,iyo,1,e)
          enddo
          enddo

          do iy=1,lxv
          do ix=1,lyv
            ixo = lxv-ix+1
            iyo = lyv-iy+1
            ii = (e-1)*(2**ldim) + (iyo-1)*lxv + ixo
            vertexc(ix,iy,1,enow) = vertex(ii)
          enddo
          enddo


        elseif (f.eq.4) then
          do iy=1,lx1
          do ix=1,lx1
            ixo = iy 
            iyo = lx1-ix+1
            xmc(ix,iy,1,enow) = xm1(ixo,iyo,1,e)
            ymc(ix,iy,1,enow) = ym1(ixo,iyo,1,e)
          enddo
          enddo

          do iy=1,lxv
          do ix=1,lyv
            ixo = iy
            iyo = lxv-ix+1
            ii = (e-1)*(2**ldim) + (iyo-1)*lxv + ixo
            vertexc(ix,iy,1,enow) = vertex(ii)
          enddo
          enddo

        endif ! select f
      endif ! if e
      enddo ! e

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_rotate3d_in(ie_ref)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'
      integer ie_ref(lelt)
      integer enow,e,f,ix,iy,iz,ixo,iyo,izo,ii,nxyz,iglsum,i

      common /ivrtx/ vertex ((2**ldim)*lelt) ! read only, TODO check set up
      integer*8 vertex,i8glmax

      integer irt(6,6)

      data irt / 6,2,5,4,1,3
     $         , 6,3,5,1,2,4
     $         , 6,4,5,2,3,1
     $         , 6,1,5,3,4,2
     $         , 1,2,3,4,5,6
     $         , 3,2,1,4,6,5 /

      if (lx1.ne.ly1) call exitti('Error rfn rot: lx1 ne ly1$',ly1)
      if (lx1.ne.lz1) call exitti('Error rfn rot: lx1 ne lz1$',lz1)

      nxyz = lx1*ly1*lz1

      ! copy into work array + rotate
      enow = 0
      do e=1,nelt
      if (ie_ref(e).gt.0) then
        enow = enow + 1

        ie_ref_wk(enow) = ie_ref(e)
        f = ie_ref(e)

        do i=1,6
           bID_wk(i,enow,1) = boundaryID(irt(i,f),e)
           bID_wk(i,enow,2) = bID_bkup(irt(i,f),e,1)
           cbc_wk(i,enow,1) = cbc(irt(i,f),e,1)
           cbc_wk(i,enow,2) = cbc_bkup(irt(i,f),e,1)
           call copy(bc_wk(1,i,enow,1), bc(1,irt(i,f),e,1), 5)
           call copy(bc_wk(1,i,enow,2), bc(1,irt(i,f),e,2), 5)
        enddo 

        if (f.eq.5) then
          ! no rotate, copy only
          call copy(xmc(1,1,1,enow),xm1(1,1,1,e),nxyz)
          call copy(ymc(1,1,1,enow),ym1(1,1,1,e),nxyz)
          call copy(zmc(1,1,1,enow),zm1(1,1,1,e),nxyz)

          call i8copy(vertexc(1,1,1,enow),vertex((e-1)*lxyzv+1),lxyzv)

        elseif (f.eq.6) then           !   v           
          do iz=1,lz1                  !  v5__F1__v6   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F4|  F6  |F2  
            ixo = ix                   !   |______|___r
            iyo = ly1-iy+1             !  v7  F3  v8   
            izo = lz1-iz+1             !   |s          
            xmc(ix,iy,iz,enow) = xm1(ixo,iyo,izo,e)
            ymc(ix,iy,iz,enow) = ym1(ixo,iyo,izo,e)
            zmc(ix,iy,iz,enow) = zm1(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = ix
            iyo = lyv-iy+1
            izo = lzv-iz+1
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertexc(ix,iy,iz,enow) = vertex(ii)
          enddo
          enddo
          enddo


        elseif (f.eq.1) then           !   v           
          do iz=1,lz1                  !  v1__F5__v2   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F4|  F1  |F2  
            ixo = ix                   !   |______|___r
            iyo = lz1-iz+1             !  v5  F6  v6   
            izo = ly1-iy+1             !   |t          
            xmc(ix,iy,iz,enow) = xm1(ixo,iyo,izo,e)
            ymc(ix,iy,iz,enow) = ym1(ixo,iyo,izo,e)
            zmc(ix,iy,iz,enow) = zm1(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = ix
            iyo = lzv-iz+1
            izo = lyv-iy+1
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertexc(ix,iy,iz,enow) = vertex(ii)
          enddo
          enddo
          enddo


        elseif (f.eq.2) then           !   v           
          do iz=1,lz1                  !  v2__F5__v4   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F1|  F2  |F3  
            ixo = iy                   !   |______|___s
            iyo = lz1-iz+1             !  v6  F6  v8   
            izo = lx1-ix+1             !   |t          
            xmc(ix,iy,iz,enow) = xm1(ixo,iyo,izo,e)
            ymc(ix,iy,iz,enow) = ym1(ixo,iyo,izo,e)
            zmc(ix,iy,iz,enow) = zm1(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = iy
            iyo = lzv-iz+1
            izo = lxv-ix+1
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertexc(ix,iy,iz,enow) = vertex(ii)
          enddo
          enddo
          enddo


        elseif (f.eq.3) then           !   v           
          do iz=1,lz1                  !  v4__F5__v3   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F2|  F3  |F4  
            ixo = lx1-ix+1             !r__|______|__< 
            iyo = lz1-iz+1             !  v8  F6  v7   
            izo = ly1-iy+1             !   |t          
            xmc(ix,iy,iz,enow) = xm1(ixo,iyo,izo,e)
            ymc(ix,iy,iz,enow) = ym1(ixo,iyo,izo,e)
            zmc(ix,iy,iz,enow) = zm1(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = lxv-ix+1
            iyo = lzv-iz+1
            izo = lyv-iy+1
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertexc(ix,iy,iz,enow) = vertex(ii)
          enddo
          enddo
          enddo


        elseif (f.eq.4) then           !   v           
          do iz=1,lz1                  !  v3__F5__v1   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F3|  F4  |F1  
            ixo = ly1-iy+1             !s__|______|__< 
            iyo = lz1-iz+1             !  v7  F6  v5   
            izo = ix                   !   |t          
            xmc(ix,iy,iz,enow) = xm1(ixo,iyo,izo,e)
            ymc(ix,iy,iz,enow) = ym1(ixo,iyo,izo,e)
            zmc(ix,iy,iz,enow) = zm1(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = lyv-iy+1
            iyo = lzv-iz+1
            izo = ix
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertexc(ix,iy,iz,enow) = vertex(ii)
          enddo
          enddo
          enddo

        endif ! select f
      endif ! if e
      enddo ! e

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_rotate2d_out(ie_ref)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'
      integer ie_ref(lelt)
      integer enow,e,f,ix,iy,ixo,iyo,ii,nxyz,i

      common /ivrtx/ vertex ((2**ldim)*lelt) ! write only, TODO check set up
      integer*8 vertex
      integer*8 nvtx,i8glmax

      integer irt(4,4)
      data irt / 1,2,3,4 
     $         , 2,3,4,1 
     $         , 3,4,1,2 
     $         , 4,1,2,3 /

      if (lx1.ne.ly1) call exitti('Error rfn rot: lx1 ne ly1$',ly1)

      nxyz = lx1*ly1*lz1

      ! prolong ie_ref, follows the parent element
      do e=1,nelt
      if (e.gt.nelt0) then
        enow = mod(e - nelt0 - 1 , nel_lv) + 1
        ie_ref(e) = ie_ref_wk(enow)
      endif
      enddo

      enow = 0            
      do e=1,nelt
      if (ie_ref(e).gt.0) then
        enow = enow + 1
        f = ie_ref(e)

        do i=1,4
           boundaryID(irt(i,f),e)   = bID_wk(i,enow,1)
           bID_bkup  (irt(i,f),e,1) = bID_wk(i,enow,2)
           cbc       (irt(i,f),e,1) = cbc_wk(i,enow,1)
           cbc_bkup  (irt(i,f),e,1) = cbc_wk(i,enow,2)
           call copy(bc(1,irt(i,f),e,1), bc_wk(1,i,enow,1), 5)
           call copy(bc(1,irt(i,f),e,2), bc_wk(1,i,enow,2), 5)
        enddo


        if (f.eq.1) then
          ! no rotate, copy only
          call copy(xm1(1,1,1,e),xmc(1,1,1,enow),nxyz)
          call copy(ym1(1,1,1,e),ymc(1,1,1,enow),nxyz)

          call i8copy(vertex((e-1)*lxyzv+1),vertexc(1,1,1,enow),lxyzv)
          

        elseif (f.eq.2) then
          do iy=1,lx1
          do ix=1,lx1
            ixo = ly1-iy+1 
            iyo = ix
            xm1(ixo,iyo,1,e) = xmc(ix,iy,1,enow)
            ym1(ixo,iyo,1,e) = ymc(ix,iy,1,enow)
          enddo
          enddo

          do iy=1,lxv
          do ix=1,lyv
            ixo = lyv-iy+1
            iyo = ix
            ii = (e-1)*(2**ldim) + (iyo-1)*lxv + ixo
            vertex(ii) = vertexc(ix,iy,1,enow)
          enddo
          enddo


        elseif (f.eq.3) then
          do iy=1,lx1
          do ix=1,lx1
            ixo = lx1-ix+1 
            iyo = ly1-iy+1
            xm1(ixo,iyo,1,e) = xmc(ix,iy,1,enow)
            ym1(ixo,iyo,1,e) = ymc(ix,iy,1,enow)
          enddo
          enddo

          do iy=1,lxv
          do ix=1,lyv
            ixo = lxv-ix+1
            iyo = lyv-iy+1
            ii = (e-1)*(2**ldim) + (iyo-1)*lxv + ixo
            vertex(ii) = vertexc(ix,iy,1,enow)
          enddo
          enddo


        elseif (f.eq.4) then
          do iy=1,lx1
          do ix=1,lx1
            ixo = iy 
            iyo = lx1-ix+1
            xm1(ixo,iyo,1,e) = xmc(ix,iy,1,enow)
            ym1(ixo,iyo,1,e) = ymc(ix,iy,1,enow)
          enddo
          enddo

          do iy=1,lxv
          do ix=1,lyv
            ixo = iy
            iyo = lxv-ix+1
            ii = (e-1)*(2**ldim) + (iyo-1)*lxv + ixo
            vertex(ii) = vertexc(ix,iy,1,enow)
          enddo
          enddo

        endif ! select f
      endif ! if e
      enddo ! e

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_rotate3d_out(ie_ref)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'
      integer ie_ref(lelt)
      integer enow,e,f,ix,iy,iz,ixo,iyo,izo,ii,nxyz,i

      common /ivrtx/ vertex ((2**ldim)*lelt) ! write only, TODO check set up
      integer*8 vertex
      integer*8 nvtx,i8glmax

      integer irt(6,6) 
      data irt / 6,2,5,4,1,3
     $         , 6,3,5,1,2,4
     $         , 6,4,5,2,3,1
     $         , 6,1,5,3,4,2 
     $         , 1,2,3,4,5,6
     $         , 3,2,1,4,6,5 /

      if (lx1.ne.ly1) call exitti('Error rfn rot: lx1 ne ly1$',ly1)

      nxyz = lx1*ly1*lz1

      ! prolong ie_ref, follows the parent element
      do e=1,nelt
      if (e.gt.nelt0) then
        enow = mod(e - nelt0 - 1 , nel_lv) + 1
        ie_ref(e) = ie_ref_wk(enow)
      endif
      enddo

      enow = 0
      do e=1,nelt
      if (ie_ref(e).gt.0) then
        enow = enow + 1
        f = ie_ref(e)

        do i=1,6
           boundaryID(irt(i,f),e)   = bID_wk(i,enow,1)
           bID_bkup  (irt(i,f),e,1) = bID_wk(i,enow,2)
           cbc       (irt(i,f),e,1) = cbc_wk(i,enow,1)
           cbc_bkup  (irt(i,f),e,1) = cbc_wk(i,enow,2)
           call copy(bc(1,irt(i,f),e,1), bc_wk(1,i,enow,1), 5)
           call copy(bc(1,irt(i,f),e,2), bc_wk(1,i,enow,2), 5)
        enddo


        if (f.eq.5) then
          ! no rotate, copy only
          call copy(xm1(1,1,1,e),xmc(1,1,1,enow),nxyz)
          call copy(ym1(1,1,1,e),ymc(1,1,1,enow),nxyz)
          call copy(zm1(1,1,1,e),zmc(1,1,1,enow),nxyz)

          call i8copy(vertex((e-1)*lxyzv+1),vertexc(1,1,1,enow),lxyzv)

        elseif (f.eq.6) then           !   v           
          do iz=1,lz1                  !  v5__F1__v6   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F4|  F6  |F2  
            ixo = ix                   !   |______|___r
            iyo = ly1-iy+1             !  v7  F3  v8   
            izo = lz1-iz+1             !   |s          
            xm1(ixo,iyo,izo,e) = xmc(ix,iy,iz,enow)
            ym1(ixo,iyo,izo,e) = ymc(ix,iy,iz,enow)
            zm1(ixo,iyo,izo,e) = zmc(ix,iy,iz,enow)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = ix
            iyo = lyv-iy+1
            izo = lzv-iz+1
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertex(ii) = vertexc(ix,iy,iz,enow)
          enddo
          enddo
          enddo


        elseif (f.eq.1) then           !   v           
          do iz=1,lz1                  !  v1__F5__v2   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F4|  F1  |F2  
            ixo = ix                   !   |______|___r
            iyo = lz1-iz+1             !  v5  F6  v6   
            izo = ly1-iy+1             !   |t          
            xm1(ixo,iyo,izo,e) = xmc(ix,iy,iz,enow)
            ym1(ixo,iyo,izo,e) = ymc(ix,iy,iz,enow)
            zm1(ixo,iyo,izo,e) = zmc(ix,iy,iz,enow)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = ix
            iyo = lzv-iz+1
            izo = lyv-iy+1
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertex(ii) = vertexc(ix,iy,iz,enow)
          enddo
          enddo
          enddo


        elseif (f.eq.2) then           !   v           
          do iz=1,lz1                  !  v2__F5__v4   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F1|  F2  |F3  
            ixo = iy                   !   |______|___s
            iyo = lz1-iz+1             !  v6  F6  v8   
            izo = lx1-ix+1             !   |t          
            xm1(ixo,iyo,izo,e) = xmc(ix,iy,iz,enow)
            ym1(ixo,iyo,izo,e) = ymc(ix,iy,iz,enow)
            zm1(ixo,iyo,izo,e) = zmc(ix,iy,iz,enow)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = iy
            iyo = lzv-iz+1
            izo = lxv-ix+1
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertex(ii) = vertexc(ix,iy,ix,enow)
          enddo
          enddo
          enddo


        elseif (f.eq.3) then           !   v           
          do iz=1,lz1                  !  v4__F5__v3   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F2|  F3  |F4  
            ixo = lx1-ix+1             !r__|______|__< 
            iyo = lz1-iz+1             !  v8  F6  v7   
            izo = ly1-iy+1             !   |t          
            xm1(ixo,iyo,izo,e) = xmc(ix,iy,iz,enow)
            ym1(ixo,iyo,izo,e) = ymc(ix,iy,iz,enow)
            zm1(ixo,iyo,izo,e) = zmc(ix,iy,iz,enow)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = lxv-ix+1
            iyo = lzv-iz+1
            izo = lyv-iy+1
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertex(ii) = vertexc(ix,iy,iz,enow)
          enddo
          enddo
          enddo


        elseif (f.eq.4) then           !   v           
          do iz=1,lz1                  !  v3__F5__v1   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F3|  F4  |F1  
            ixo = ly1-iy+1             !s__|______|__< 
            iyo = lz1-iz+1             !  v7  F6  v5   
            izo = ix                   !   |t          
            xm1(ixo,iyo,izo,e) = xmc(ix,iy,iz,enow)
            ym1(ixo,iyo,izo,e) = ymc(ix,iy,iz,enow)
            zm1(ixo,iyo,izo,e) = zmc(ix,iy,iz,enow)
          enddo
          enddo
          enddo

          do iz=1,lzv
          do iy=1,lyv
          do ix=1,lxv
            ixo = lyv-iy+1
            iyo = lzv-iz+1
            izo = ix
            ii = (e-1)*lxv*lyv*lzv + (izo-1)*lxv*lyv + (iyo-1)*lxv + ixo
            vertex(ii) = vertexc(ix,iy,iz,enow)
          enddo
          enddo
          enddo

        endif ! select f
      endif ! if e
      enddo ! e

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_rotate3d_fld_in(ie_ref,uo,ui)
      implicit none
      include 'SIZE'

c      include 'REFINE'
      integer ie_ref(lelt)
      integer enow,e,f,ix,iy,iz,ixo,iyo,izo,ii,nxyz

      real ui(lx1,ly1,lz1,lelt)
      real uo(lx1,ly1,lz1,lelt)

      integer irt(6,6)
      data irt / 6,2,5,4,1,3
     $         , 6,3,5,1,2,4
     $         , 6,4,5,2,3,1
     $         , 6,1,5,3,4,2
     $         , 1,2,3,4,5,6
     $         , 3,2,1,4,6,5 /      

      nxyz = lx1*ly1*lz1

      ! copy into work array + rotate
      enow = 0
      do e=1,nelv
      if (ie_ref(e).gt.0) then
        enow = enow + 1

        f = ie_ref(e)

        if (f.eq.5) then
          ! no rotate, copy only
          call copy(uo(1,1,1,enow),ui(1,1,1,e),nxyz)

        elseif (f.eq.6) then           !   v           
          do iz=1,lz1                  !  v5__F1__v6   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F4|  F6  |F2  
            ixo = ix                   !   |______|___r
            iyo = ly1-iy+1             !  v7  F3  v8   
            izo = lz1-iz+1             !   |s          
            uo(ix,iy,iz,enow) = ui(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

        elseif (f.eq.1) then           !   v           
          do iz=1,lz1                  !  v1__F5__v2   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F4|  F1  |F2  
            ixo = ix                   !   |______|___r
            iyo = lz1-iz+1             !  v5  F6  v6   
            izo = ly1-iy+1             !   |t          
            uo(ix,iy,iz,enow) = ui(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

        elseif (f.eq.2) then           !   v           
          do iz=1,lz1                  !  v2__F5__v4   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F1|  F2  |F3  
            ixo = iy                   !   |______|___s
            iyo = lz1-iz+1             !  v6  F6  v8   
            izo = lx1-ix+1             !   |t          
            uo(ix,iy,iz,enow) = ui(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

        elseif (f.eq.3) then           !   v           
          do iz=1,lz1                  !  v4__F5__v3   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F2|  F3  |F4  
            ixo = lx1-ix+1             !r__|______|__< 
            iyo = lz1-iz+1             !  v8  F6  v7   
            izo = ly1-iy+1             !   |t          
            uo(ix,iy,iz,enow) = ui(ixo,iyo,izo,e)
          enddo
          enddo
          enddo

        elseif (f.eq.4) then           !   v           
          do iz=1,lz1                  !  v3__F5__v1   
          do iy=1,ly1                  !   |      |    
          do ix=1,lx1                  ! F3|  F4  |F1  
            ixo = ly1-iy+1             !s__|______|__< 
            iyo = lz1-iz+1             !  v7  F6  v5   
            izo = ix                   !   |t          
            uo(ix,iy,iz,enow) = ui(ixo,iyo,izo,e)
          enddo
          enddo
          enddo
        endif

      endif ! if e
      enddo ! e

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_refine2d(ie_ref,Nref,Rref,Rlim)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'

      integer Nref,ie_ref(lelt)
      real Rref(lref-1),Rlim(4,lref-1)

      real xlv(lx1,ly1,lz1,lelt)
      real ylv(lx1,ly1,lz1,lelt)

      integer ilv,i,e,e1,e2,f,nxyz,enow,enew,iglsum,ivlsum,nedge
      real z0,z1
      real edge_len(2*(ldim-1),lelt) ! 2D 2 edge, 3D 4 edges
      integer mfield,nfldt

      integer*8 vtmp

      nxyz = lx1*ly1*lz1


      call get_edge_len(edge_len,xmc,ymc,zmc,nel_lv)

      enew = 0 ! num of new elements
      do ilv=1,Nref ! bottom-top

        ! Get xyz, interp ref domain from [-1,1] to [z0,z1]
        if (ilv.eq.1) then
          z0 = -1.0
          z1 = Rref(ilv)*2.0-1.0
        elseif (ilv.eq.Nref) then
          z0 = z1
          z1 = 1.0
        else 
          z0 = z1
          z1 = Rref(ilv)*2.0-1.0
        endif

        ! limiter 1: change z0,z1 via edge length, globally
        ! TODO (hard to implement): 
        !       take min/max subject to seperated connectivity
        !       for example, min/max zlv for each spheres
        if (ilv.ne.Nref) then ! top lv doesn't have new lyr
          call rfn_apply_limiter1(z1,edge_len
     $                           ,Rlim(1,ilv),Rlim(2,ilv),nel_lv)
        endif
c        if(nid.eq.0) write(*,*)'dbg',ilv,z0,z1

        ! get xyz
        do enow=1,nel_lv
          call int_gll_pt_2d_e(xlv(1,1,1,enow),z0,z1,lx1
     $                        ,xmc(1,1,1,enow),lx1)
          call int_gll_pt_2d_e(ylv(1,1,1,enow),z0,z1,lx1
     $                        ,ymc(1,1,1,enow),lx1)
        enddo

        ! assign xyz and cbc
        if (ilv.eq.Nref) then ! update xyz to the old elements
          vtmp = nvtx0 + (ilv-2)*nvtx_lv
          do enow=1,nel_lv
            call copy(xmc(1,1,1,enow),xlv(1,1,1,enow),nxyz)
            call copy(ymc(1,1,1,enow),ylv(1,1,1,enow),nxyz)

            ! bc
            if (cbc_wk(1,enow,1).eq.'P  ') then
              bc_wk(1,1,enow,1) = 0.0
              bc_wk(2,1,enow,1) = 0.0
            endif
            if (cbc_wk(1,enow,2).eq.'P  ') then
              bc_wk(1,1,enow,2) = 0.0
              bc_wk(2,1,enow,2) = 0.0
            endif

            ! bc
            call rzero(bc_wk(1,1,enow,1),5)
            call rzero(bc_wk(1,1,enow,2),5)

            ! cbc
            bID_wk(1,enow,1) = 0
            bID_wk(1,enow,2) = 0
            cbc_wk(1,enow,1) = 'E  ' ! this will be the last
            cbc_wk(1,enow,2) = 'E  ' ! this will be the last

            ! vtx, keep face 3, update face 1
            do i=1,lxv
              vertexc(i,1,1,enow) = vertexu(i,1,enow) + vtmp ! face 1
            enddo
          enddo
        else ! new elements

          vtmp = nvtx0 + (ilv-2)*nvtx_lv
          enow = 0
          do e=1,nelt0
          if (ie_ref(e).gt.0) then
            enow = enow + 1
            enew = enew + 1
            e1 = nel_lv + enew
            e2 = nelt0 + enew
            call copy(xmc(1,1,1,e1),xlv(1,1,1,enow),nxyz)
            call copy(ymc(1,1,1,e1),ylv(1,1,1,enow),nxyz)

            ! rotate index
            ie_ref(e2) = ie_ref(e)

            ! bc 
            call copy (bc_wk(1,2,e1,1),bc_wk(1,2,enow,1),5)
            call copy (bc_wk(1,4,e1,1),bc_wk(1,4,enow,1),5)
            call rzero(bc_wk(1,3,e1,1),5)
            call copy (bc_wk(1,2,e1,2),bc_wk(1,2,enow,2),5)
            call copy (bc_wk(1,4,e1,2),bc_wk(1,4,enow,2),5)
            call rzero(bc_wk(1,3,e1,2),5)
            ! TODO: take care bc for periodic

            ! cbc
            bID_wk(2,e1,1) = bID_wk(2,enow,1)
            bID_wk(4,e1,1) = bID_wk(4,enow,1)
            bID_wk(3,e1,1) = 0
            bID_wk(2,e1,2) = bID_wk(2,enow,2)
            bID_wk(4,e1,2) = bID_wk(4,enow,2)
            bID_wk(3,e1,2) = 0

            cbc_wk(2,e1,1) = cbc_wk(2,enow,1)
            cbc_wk(4,e1,1) = cbc_wk(4,enow,1)
            cbc_wk(3,e1,1) = 'E  '
            cbc_wk(2,e1,2) = cbc_wk(2,enow,2)
            cbc_wk(4,e1,2) = cbc_wk(4,enow,2)
            cbc_wk(3,e1,2) = 'E  '
            if (ilv.eq.1) then
              call copy (bc_wk(1,1,e1,1),bc_wk(1,1,enow,1),5)
              call copy (bc_wk(1,1,e1,2),bc_wk(1,1,enow,2),5)
              bID_wk(1,e1,1) = bID_wk(1,enow,1) ! this will copy before re-assign
              bID_wk(1,e1,2) = bID_wk(1,enow,2)
              cbc_wk(1,e1,1) = cbc_wk(1,enow,1) ! this will copy before re-assign
              cbc_wk(1,e1,2) = cbc_wk(1,enow,2)
            else
              call rzero(bc_wk(1,1,e1,1),5)
              call rzero(bc_wk(1,1,e1,2),5)
              bID_wk(1,e1,1) = 0
              bID_wk(1,e1,2) = 0
              cbc_wk(1,e1,1) = 'E  '
              cbc_wk(1,e1,2) = 'E  '
            endif

            ! vtx ! TODO: separate if/loop
            if (ilv.eq.1) then
              do i=1,lxv
                vertexc(i,  1,1,e1) = vertexc(i,1,1,enow)              ! face 1
                vertexc(i,lyv,1,e1) = vertexu(i,1,enow) + vtmp+nvtx_lv ! face 3
              enddo
            else
              do i=1,lxv
                vertexc(i,  1,1,e1) = vertexu(i,1,enow) + vtmp         ! face 1
                vertexc(i,lyv,1,e1) = vertexu(i,1,enow) + vtmp+nvtx_lv ! face 3
              enddo
            endif
          endif
          enddo
          nelv = nelv + nel_lv
          nelt = nelt + nel_lv
        endif

      enddo ! lv

      ! limiter 2: locally adjust hex27 to reconstruct xyz
      call rfn_apply_limiter2_2d(xmc,ymc,Rlim,Nref,nel_lv) ! This is so buggy...

      call update_nels(Nref-1)
cc#ifdef DPROCMAP ! TODO: this won't work
cc            call exitti('gllel,gllnid not set for DPROCMAP$',nid)
cc#endif
cc      ! dbg chk
cc      do e=1,nelv
cc        write(*,*)'rfn lglel',nid,e,lglel(e)
cc      enddo
cc      do e=1,nelgv
cc        write(*,*)'rfn gllel',nid,e,gllel(e),gllnid(e)
cc      enddo

      ! TODO: need on-the-fly fix_geom...

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_refine3d(ie_ref,Nref,Rref,Rlim)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'

      integer Nref,ie_ref(lelt)
      real Rref(lref-1),Rlim(4,lref-1)

      real xlv(lx1,ly1,lz1,lelt)
      real ylv(lx1,ly1,lz1,lelt)
      real zlv(lx1,ly1,lz1,lelt)

      integer ilv,i,e,e1,e2,f,nxyz,enow,enew,iglsum,ivlsum,nedge
      real z0,z1
      real edge_len(2*(ldim-1),lelt) ! 2D 2 edge, 3D 4 edges
      real edge_min,edge_max,glmin,glmax
      integer mfield,nfldt

      integer*8 vtmp

      nxyz = lx1*ly1*lz1

      write(*,*)'DDDref',Nref


      call get_edge_len(edge_len,xmc,ymc,zmc,nel_lv)

      enew = 0 ! num of new elements
      do ilv=1,Nref ! bottom-top

        ! Get xyz, interp ref domain from [-1,1] to [z0,z1]
        if (ilv.eq.1) then
          z0 = -1.0
          z1 = Rref(ilv)*2.0-1.0
        elseif (ilv.eq.Nref) then
          z0 = z1
          z1 = 1.0
        else 
          z0 = z1
          z1 = Rref(ilv)*2.0-1.0
        endif

        ! limiter 1: change z0,z1 via edge length, globally
        ! TODO (hard to implement): 
        !       take min/max subject to seperated connectivity
        !       for example, min/max zlv for each spheres
        if (ilv.ne.Nref) then ! top lv doesn't have new lyr
          call rfn_apply_limiter1(z1,edge_len
     $                           ,Rlim(1,ilv),Rlim(2,ilv),nel_lv)
        endif
c        if(nid.eq.0) write(*,*)'dbg',ilv,z0,z1

        ! get xyz
        do enow=1,nel_lv
          call int_gll_pt_3d_e(xlv(1,1,1,enow),z0,z1,lx1
     $                        ,xmc(1,1,1,enow),lx1)
          call int_gll_pt_3d_e(ylv(1,1,1,enow),z0,z1,lx1
     $                        ,ymc(1,1,1,enow),lx1)
          call int_gll_pt_3d_e(zlv(1,1,1,enow),z0,z1,lx1
     $                        ,zmc(1,1,1,enow),lx1)
        enddo

        ! assign xyz and cbc
        if (ilv.eq.Nref) then ! update xyz to the old elements
          vtmp = nvtx0 + (ilv-2)*nvtx_lv
          do enow=1,nel_lv
            call copy(xmc(1,1,1,enow),xlv(1,1,1,enow),nxyz)
            call copy(ymc(1,1,1,enow),ylv(1,1,1,enow),nxyz)
            call copy(zmc(1,1,1,enow),zlv(1,1,1,enow),nxyz)

            ! bc
            call rzero(bc_wk(1,5,enow,1),5)
            call rzero(bc_wk(1,5,enow,2),5)

            ! cbc
            bID_wk(5,enow,1) = 0
            bID_wk(5,enow,2) = 0
            cbc_wk(5,enow,1) = 'E  ' ! this will be the last
            cbc_wk(5,enow,2) = 'E  ' ! this will be the last

            ! vtx, keep face 6, update face 5
            do i=1,lxv*lzv
              vertexc(i,1,1,enow) = vertexu(i,1,enow) + vtmp ! face 5
            enddo
          enddo
        else ! new elements

          vtmp = nvtx0 + (ilv-2)*nvtx_lv
          enow = 0
          do e=1,nelt0
          if (ie_ref(e).gt.0) then
            enow = enow + 1
            enew = enew + 1
            e1 = nel_lv + enew
            e2 = nelt0 + enew
            call copy(xmc(1,1,1,e1),xlv(1,1,1,enow),nxyz)
            call copy(ymc(1,1,1,e1),ylv(1,1,1,enow),nxyz)
            call copy(zmc(1,1,1,e1),zlv(1,1,1,enow),nxyz)

            ! rotate index
            ie_ref(e2) = ie_ref(e)

            ! bc
            call copy (bc_wk(1,1,e1,1),bc_wk(1,1,enow,1),5)
            call copy (bc_wk(1,2,e1,1),bc_wk(1,2,enow,1),5)
            call copy (bc_wk(1,3,e1,1),bc_wk(1,3,enow,1),5)
            call copy (bc_wk(1,4,e1,1),bc_wk(1,4,enow,1),5)
            call rzero(bc_wk(1,6,e1,1),5)
            call copy (bc_wk(1,1,e1,2),bc_wk(1,1,enow,2),5)
            call copy (bc_wk(1,2,e1,2),bc_wk(1,2,enow,2),5)
            call copy (bc_wk(1,3,e1,2),bc_wk(1,3,enow,2),5)
            call copy (bc_wk(1,4,e1,2),bc_wk(1,4,enow,2),5)
            call rzero(bc_wk(1,6,e1,2),5)

            ! cbc
            bID_wk(1,e1,1) = bID_wk(1,enow,1)
            bID_wk(2,e1,1) = bID_wk(2,enow,1)
            bID_wk(3,e1,1) = bID_wk(3,enow,1)
            bID_wk(4,e1,1) = bID_wk(4,enow,1)
            bID_wk(6,e1,1) = 0
            bID_wk(1,e1,2) = bID_wk(1,enow,2)
            bID_wk(2,e1,2) = bID_wk(2,enow,2)
            bID_wk(3,e1,2) = bID_wk(3,enow,2)
            bID_wk(4,e1,2) = bID_wk(4,enow,2)
            bID_wk(6,e1,2) = 0

            cbc_wk(1,e1,1) = cbc_wk(1,enow,1)
            cbc_wk(2,e1,1) = cbc_wk(2,enow,1)
            cbc_wk(3,e1,1) = cbc_wk(3,enow,1)
            cbc_wk(4,e1,1) = cbc_wk(4,enow,1)
            cbc_wk(6,e1,1) = 'E  '
            cbc_wk(1,e1,2) = cbc_wk(1,enow,2)
            cbc_wk(2,e1,2) = cbc_wk(2,enow,2)
            cbc_wk(3,e1,2) = cbc_wk(3,enow,2)
            cbc_wk(4,e1,2) = cbc_wk(4,enow,2)
            cbc_wk(6,e1,2) = 'E  '
            if (ilv.eq.1) then
              call copy (bc_wk(1,5,e1,1),bc_wk(1,5,enow,1),5)
              call copy (bc_wk(1,5,e1,2),bc_wk(1,5,enow,2),5)
              bID_wk(5,e1,1) = bID_wk(5,enow,1) ! this will copy before re-assign
              bID_wk(5,e1,2) = bID_wk(5,enow,2)
              cbc_wk(5,e1,1) = cbc_wk(5,enow,1) ! this will copy before re-assign
              cbc_wk(5,e1,2) = cbc_wk(5,enow,2)
            else
              call rzero(bc_wk(1,5,e1,1),5)
              call rzero(bc_wk(1,5,e1,2),5)
              bID_wk(5,e1,1) = 0
              bID_wk(5,e1,2) = 0
              cbc_wk(5,e1,1) = 'E  '
              cbc_wk(5,e1,2) = 'E  '
            endif

            ! vtx ! TODO: separate if/loop
            if (ilv.eq.1) then
              do i=1,lxv*lzv
                vertexc(i,1,  1,e1) = vertexc(i,1,1,enow)              ! face 5
                vertexc(i,1,lzv,e1) = vertexu(i,1,enow) + vtmp+nvtx_lv ! face 6
              enddo
            else
              do i=1,lxv*lzv
                vertexc(i,1,  1,e1) = vertexu(i,1,enow) + vtmp         ! face 5
                vertexc(i,1,lzv,e1) = vertexu(i,1,enow) + vtmp+nvtx_lv ! face 6
              enddo
            endif
          endif
          enddo
          nelv = nelv + nel_lv
          nelt = nelt + nel_lv
        endif

      enddo ! lv

      enow=1
      do ilv=1,Nref
        call get_edge_len(edge_len,xmc(1,1,1,enow)
     $                   ,ymc(1,1,1,enow),zmc(1,1,1,enow),nel_lv)
        enow=enow+nel_lv

        edge_min = glmin(edge_len,2*(ldim-1)*nel_lv)
        edge_max = glmax(edge_len,2*(ldim-1)*nel_lv) 
        if(nio.eq.0) write(*,*)'rfn elen bfr lim2',ilv,edge_min,edge_max
      enddo
      ! limiter 2: locally adjust hex27 to reconstruct xyz
      call rfn_apply_limiter2_3d(xmc,ymc,zmc,Rlim,Nref,nel_lv) ! This is so buggy...

      call update_nels(Nref-1)
cc#ifdef DPROCMAP ! TODO: this won't work
cc            call exitti('gllel,gllnid not set for DPROCMAP$',nid)
cc#endif
cc      ! dbg chk
cc      do e=1,nelv
cc        write(*,*)'rfn lglel',nid,e,lglel(e)
cc      enddo
cc      do e=1,nelgv
cc        write(*,*)'rfn gllel',nid,e,gllel(e),gllnid(e)
cc      enddo

      ! TODO: need on-the-fly fix_geom...

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_extrude_lyr_2d(ie_ref,igeo,icht)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'

      integer ie_ref(lelt),igeo,icht

      real xlv(lx1,ly1,lz1,lelt)
      real ylv(lx1,ly1,lz1,lelt)
      real zlv(lx1,ly1,lz1,lelt)

      integer ilv,i,iside,e,e1,e2,nxyz,enow,enew,iglsum,ivlsum
      integer mfield,nfldt

      integer j

      nxyz = lx1*ly1*lz1
      enew = 0 ! num of new elements

      call rzero(xlv,nxyz*nelt)
      call rzero(ylv,nxyz*nelt)
      call rzero(zlv,nxyz*nelt)

      ! get xyz
c     do enow=1,nel_lv
      enow = 0
      do e=1,nelt0
      iside = ie_ref(e)
      if (iside.gt.0) then
        enow = enow + 1
        call copy(xlv(1,ly1,1,enow),xmc(1,1,1,enow),lx1*lz1) ! E_old face 1 to E_new face 3
        call copy(ylv(1,ly1,1,enow),ymc(1,1,1,enow),lx1*lz1) ! E_old face 1 to E_new face 3

        call rfn_extrude_pj2d(xlv(1,1,1,enow),ylv(1,1,1,enow)
     $                       ,zlv(1,1,1,enow),igeo,iside,e) ! E_new face 3 to E_new face 1
        call int_recover2d_from_f13(xlv(1,1,1,enow),lx1) ! gen xyz from face 1 and face 3
        call int_recover2d_from_f13(ylv(1,1,1,enow),lx1)
      endif
      enddo

      ! gen new elements
      enow = 0
      do e=1,nelt0
      if (ie_ref(e).gt.0) then
        enow = enow + 1
        enew = enew + 1
        e1 = nel_lv + enew
        e2 = nelt0 + enew
        call copy(xmc(1,1,1,e1),xlv(1,1,1,enow),nxyz)
        call copy(ymc(1,1,1,e1),ylv(1,1,1,enow),nxyz)

        ! rotate index
        ie_ref(e2) = ie_ref(e)

        ! cbc
        if (icht.eq.0) then
        bID_wk(2,e1,1) = bID_wk(2,enow,1)
        bID_wk(4,e1,1) = bID_wk(4,enow,1)
        bID_wk(3,e1,1) = 0
        bID_wk(1,e1,1) = bID_wk(1,enow,1) ! this will copy before re-assign
        endif

        bID_wk(2,e1,2) = bID_wk(2,enow,2)
        bID_wk(4,e1,2) = bID_wk(4,enow,2)
        bID_wk(3,e1,2) = 0
        bID_wk(1,e1,2) = bID_wk(1,enow,2)

        if (icht.eq.0) then
        cbc_wk(2,e1,1) = cbc_wk(2,enow,1)
        cbc_wk(4,e1,1) = cbc_wk(4,enow,1)
        cbc_wk(3,e1,1) = 'E  '
        cbc_wk(1,e1,1) = cbc_wk(1,enow,1) ! this will copy before re-assign
        endif

        cbc_wk(2,e1,2) = cbc_wk(2,enow,2)
        cbc_wk(4,e1,2) = cbc_wk(4,enow,2)
        cbc_wk(3,e1,2) = 'E  '
        cbc_wk(1,e1,2) = cbc_wk(1,enow,2)

        ! vtx
        do i=1,lxv
          vertexc(i,  1,1,e1) = vertexu(i,1,enow) + nvtx0   ! face 1
          vertexc(i,lyv,1,e1) = vertexc(i,1,1,enow)         ! face 3
        enddo
      endif
      enddo
      if (icht.eq.0) nelv = nelv + nel_lv
      nelt = nelt + nel_lv

      ! update old elements
      if (icht.eq.0) then
      do enow=1,nel_lv
        ! cbc
        bID_wk(1,enow,1) = 0
        cbc_wk(1,enow,1) = 'E  ' ! reassign
      enddo
      endif

      do enow=1,nel_lv
        ! cbc
        bID_wk(1,enow,2) = 0
        cbc_wk(1,enow,2) = 'E  ' 
      enddo

      call update_nels(1)

cc      ! dbg chk
cc      do e=1,nelv
cc        write(*,*)'rfn lglel',nid,e,lglel(e)
cc      enddo
cc      do e=1,nelgv
cc        write(*,*)'rfn gllel',nid,e,gllel(e),gllnid(e)
cc      enddo

      ! TODO: need on-the-fly fix_geom...



      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_extrude_lyr_3d(ie_ref,igeo,icht)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'

      integer ie_ref(lelt),igeo,icht

      real xlv(lx1,ly1,lz1,lelt)
      real ylv(lx1,ly1,lz1,lelt)
      real zlv(lx1,ly1,lz1,lelt)

      integer ilv,i,e,iside,e1,e2,nxyz,enow,enew,iglsum,ivlsum
      integer mfield,nfldt

      integer j

      nxyz = lx1*ly1*lz1
      enew = 0 ! num of new elements

      call rzero(xlv,nxyz*nelt)
      call rzero(ylv,nxyz*nelt)
      call rzero(zlv,nxyz*nelt)

      ! get xyz
c      do enow=1,nel_lv
      enow = 0
      do e=1,nelt0
      iside = ie_ref(e)
      if (iside.gt.0) then
        call copy(xlv(1,1,lz1,enow),xmc(1,1,1,enow),lx1*lz1) ! E_old face 5 to E_new face 6
        call copy(ylv(1,1,lz1,enow),ymc(1,1,1,enow),lx1*lz1) ! E_old face 5 to E_new face 6
        call copy(zlv(1,1,lz1,enow),zmc(1,1,1,enow),lx1*lz1) ! E_old face 5 to E_new face 6

        call rfn_extrude_pj3d(xlv(1,1,1,enow),ylv(1,1,1,enow)
     $                       ,zlv(1,1,1,enow),igeo,iside,e)  ! E_new face 6 to E_new face 5
        call int_recover3d_from_f56(xlv(1,1,1,enow),lx1) ! gen xyz from face 5 and face 6
        call int_recover3d_from_f56(ylv(1,1,1,enow),lx1)
        call int_recover3d_from_f56(zlv(1,1,1,enow),lx1)
      endif
      enddo

      ! gen new elements
      enow = 0
      do e=1,nelt0
      if (ie_ref(e).gt.0) then
        enow = enow + 1
        enew = enew + 1
        e1 = nel_lv + enew
        e2 = nelt0 + enew
        call copy(xmc(1,1,1,e1),xlv(1,1,1,enow),nxyz)
        call copy(ymc(1,1,1,e1),ylv(1,1,1,enow),nxyz)
        call copy(zmc(1,1,1,e1),zlv(1,1,1,enow),nxyz)

        ! rotate index
        ie_ref(e2) = ie_ref(e)

        ! cbc
        if (icht.eq.0) then
        bID_wk(1,e1,1) = bID_wk(1,enow,1)
        bID_wk(2,e1,1) = bID_wk(2,enow,1)
        bID_wk(3,e1,1) = bID_wk(3,enow,1)
        bID_wk(4,e1,1) = bID_wk(4,enow,1)
        bID_wk(6,e1,1) = 0
        bID_wk(5,e1,1) = bID_wk(5,enow,1) ! this will copy before re-assign
        endif

        bID_wk(1,e1,2) = bID_wk(1,enow,2)
        bID_wk(2,e1,2) = bID_wk(2,enow,2)
        bID_wk(3,e1,2) = bID_wk(3,enow,2)
        bID_wk(4,e1,2) = bID_wk(4,enow,2)
        bID_wk(6,e1,2) = 0
        bID_wk(5,e1,2) = bID_wk(5,enow,2)

        if (icht.eq.0) then
        cbc_wk(1,e1,1) = cbc_wk(1,enow,1)
        cbc_wk(2,e1,1) = cbc_wk(2,enow,1)
        cbc_wk(3,e1,1) = cbc_wk(3,enow,1)
        cbc_wk(4,e1,1) = cbc_wk(4,enow,1)
        cbc_wk(6,e1,1) = 'E  '
        cbc_wk(5,e1,1) = cbc_wk(5,enow,1) ! this will copy before re-assign
        endif

        cbc_wk(1,e1,2) = cbc_wk(1,enow,2)
        cbc_wk(2,e1,2) = cbc_wk(2,enow,2)
        cbc_wk(3,e1,2) = cbc_wk(3,enow,2)
        cbc_wk(4,e1,2) = cbc_wk(4,enow,2)
        cbc_wk(6,e1,2) = 'E  '
        cbc_wk(5,e1,2) = cbc_wk(5,enow,2)

        ! vtx
        do i=1,lxv
          vertexc(i,1,  1,e1) = vertexu(i,1,enow) + nvtx0   ! face 5
          vertexc(i,1,lz1,e1) = vertexc(i,1,1,enow)         ! face 6
        enddo
      endif
      enddo
      if (icht.eq.0) nelv = nelv + nel_lv
      nelt = nelt + nel_lv

      ! update old elements
      if (icht.eq.0) then
      do enow=1,nel_lv
        ! cbc
        bID_wk(5,enow,1) = 0
        cbc_wk(5,enow,1) = 'E  ' ! reassign
      enddo
      endif

      do enow=1,nel_lv
        ! cbc
        bID_wk(5,enow,2) = 0
        cbc_wk(5,enow,2) = 'E  ' 
      enddo

      call update_nels(1)

      return
      end
c-----------------------------------------------------------------------
      subroutine update_nels(Nlv_new)
c     Input: nelv
c       nelv0,nelgv0,nel_lv,nel_lv_shift (from REFINE)
c     Update:
c       nelt, nelgv, nelgt, nelg, lglel, gllel, gllnid
c     NOT support DPROCMAP, but I don't have check on this
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'

      integer e,enow,e1,e2,iglsum
      integer mfield,nfldt
      integer ilv,Nlv_new

      integer ivlmin,ivlmax ! TODO: dbg

      ! Input: nelv

c      nelt = nelv ! TODO cht
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

      DO 100 IFIELD=MFIELD,nfldt+(LDIMT-1 - NPSCAL)
         IF (IFTMSH(IFIELD)) THEN
             NELFLD(IFIELD) = NELT
         ELSE
             NELFLD(IFIELD) = NELV
         ENDIF
 100  CONTINUE


      ! set lglel,gllel,gllnid for new e
      do e=nelgv0+1,nelgv
        gllel(e) = 0
        gllnid(e) = 0
      enddo
      do ilv=1,Nlv_new
      do enow=1,nel_lv
        e1 = nelt0 + nel_lv*(ilv-1) + enow
        e2 = nelgt0 + nelg_lv*(ilv-1) + nel_lv_shift + enow

        lglel(e1) = e2
        gllel(e2) = e1
        gllnid(e2) = nid
      enddo
      enddo
      do e=nelgv0+1,nelgv ! TODO use igop?
        gllel(e) = iglsum(gllel(e),1)
        gllnid(e) = iglsum(gllnid(e),1)
      enddo

c      write(*,*)'LL',nelv0,nelgv0,nel_lv_shift,nel_lv,Nlv_new
c      write(*,*)'NN',nid,nelv,nelt,nelgv,nelgt,nelg(1),nelg(2)
c      write(*,*)'GG',nid,ivlmin(gllel,nelgv),ivlmax(gllel,nelgv)


      return
      end
c-----------------------------------------------------------------------
      subroutine test_rotation
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'

      real vxtmp(lx1,ly1,lz1,lelv)
     $   , vytmp(lx1,ly1,lz1,lelv)
     $   , vztmp(lx1,ly1,lz1,lelv)
     $   , prtmp(lx2,ly2,lz2,lelv)
     $   , ttmp (lx1,ly1,lz1,lelt,ldimt)
      real xm1_bak(lx1,ly1,lz1,lelt)
      real ym1_bak(lx1,ly1,lz1,lelt)
      real zm1_bak(lx1,ly1,lz1,lelt)

      integer ie_ref(lelt),e,f,ix,iy,iz
      character s1

      ifxyo = .true.
      call opcopy(xm1_bak,ym1_bak,zm1_bak,xm1,ym1,zm1)


      do e=1,nelt
        do iz=1,lz1
        do iy=1,ly1
        do ix=1,lx1
          vxtmp(ix,iy,iz,e) = real(ix)
          vytmp(ix,iy,iz,e) = real(iy)
          vztmp(ix,iy,iz,e) = real(iz)
        enddo
        enddo
        enddo

        call cfill(prtmp(1,1,1,e),real(lglel(e)),lx1*ly1*lz1)

        do f=1,2*ldim
          call facev(ttmp,e,f,real(f),lx1,ly1,lz1)
        enddo
      enddo

      do f=1,2*ldim
        write(s1,'(I1)')f

        call izero(ie_ref,lelt)
        call ifill(ie_ref,f,nelt)
        call opcopy(xm1,ym1,zm1,xm1_bak,ym1_bak,zm1_bak)
  
        if (ldim.eq.2) then
          call rfn_rotate2d_in(ie_ref)
        else
          call rfn_rotate3d_in(ie_ref)
        endif
  
        call opcopy(xm1,ym1,zm1,xmc,ymc,zmc)
        call outpost(vxtmp,vytmp,vztmp,prtmp,ttmp,'rt'//s1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine test_sort
      implicit none
      integer*8 ivtx(10),ind(10)
      integer nlen,ierr,ialgo,lb,i

      integer nidd,np,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidd,np,nekcomm,nekgroup,nekreal

      if(nidd.eq.0) write(*,*)'vtx test'

      if(nidd.eq.0)then ! np=2
        ivtx(1) = 4
        ivtx(2) = 3
        nlen = 2
      else
        ivtx(1) = 5
        ivtx(2) = 7
        ivtx(3) = 3
        ivtx(4) = 1
        nlen = 4
      endif

      call nekgsync()
      write(*,*)'vtx in',nidd,(ivtx(i),i=1,nlen)
      call nekgsync()

      ialgo = 0 ! 0=bin sort, 1=hyprecube sort
      lb = 0 ! 1=load balanced

c      call FPSORT_INT8(ind,nlen,ivtx,ialgo,lb,nekcomm,ierr) ! FIXME
cc    call fparrsb_fpsort_int8(ind,nlen,ivtx,ialgo,lb,nekcomm,ierr)
      call nekgsync()

      write(*,*)'vtx out',nidd,(ind(i),i=1,nlen)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_edge_len(edge_len,xmc,ymc,zmc,nel) ! TODO: edge mode, face mode
c     Compte the edges' length for the linearlized elements
      implicit none
      include 'SIZE'
      real edge_len(2*(ldim-1),lelt),tmn,tmx,glmin,glmax
      real xmc(lx1,ly1,lz1,lelt)
      real ymc(lx1,ly1,lz1,lelt)
      real zmc(lx1,ly1,lz1,lelt)
      integer e,nel

      if (ldim.eq.2) then

        call rzero(edge_len,2*lelt)

        ! face 1 - face 3 TODO: use dist2d/dist3d
        do e=1,nel
          edge_len(1,e) = sqrt( (xmc(  1,1,1,e)-xmc(  1,ly1,1,e))**2
     $                         +(ymc(  1,1,1,e)-ymc(  1,ly1,1,e))**2 )

          edge_len(2,e) = sqrt( (xmc(lx1,1,1,e)-xmc(lx1,ly1,1,e))**2
     $                         +(ymc(lx1,1,1,e)-ymc(lx1,ly1,1,e))**2 )
        enddo 

        tmn=glmin(edge_len,2*nel)
        tmx=glmax(edge_len,2*nel)
        if (nio.eq.0) write(*,23)'rfn edges',tmn,tmx

      else

        call rzero(edge_len,4*lelt)

        ! face 5 - face 6
        do e=1,nel
          edge_len(1,e)=sqrt( (xmc(  1,  1,1,e)-xmc(  1,  1,lz1,e))**2
     $                       +(ymc(  1,  1,1,e)-ymc(  1,  1,lz1,e))**2 
     $                       +(zmc(  1,  1,1,e)-zmc(  1,  1,lz1,e))**2 )
          edge_len(2,e)=sqrt( (xmc(lx1,  1,1,e)-xmc(lx1,  1,lz1,e))**2
     $                       +(ymc(lx1,  1,1,e)-ymc(lx1,  1,lz1,e))**2 
     $                       +(zmc(lx1,  1,1,e)-zmc(lx1,  1,lz1,e))**2 )
          edge_len(3,e)=sqrt( (xmc(  1,ly1,1,e)-xmc(  1,ly1,lz1,e))**2
     $                       +(ymc(  1,ly1,1,e)-ymc(  1,ly1,lz1,e))**2 
     $                       +(zmc(  1,ly1,1,e)-zmc(  1,ly1,lz1,e))**2 )
          edge_len(4,e)=sqrt( (xmc(lx1,ly1,1,e)-xmc(lx1,ly1,lz1,e))**2
     $                       +(ymc(lx1,ly1,1,e)-ymc(lx1,ly1,lz1,e))**2 
     $                       +(zmc(lx1,ly1,1,e)-zmc(lx1,ly1,lz1,e))**2 )
        enddo

        tmn=glmin(edge_len,4*nel)
        tmx=glmax(edge_len,4*nel)
        if (nio.eq.0) write(*,23)'rfn edges',tmn,tmx

      endif

   23 format(a,2g13.5)


      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_apply_limiter1(zlv,edge_len,rmin,rmax,nel)
      implicit none
      include 'SIZE'

      real edge_len(2*(ldim-1),lelt) ! 2D 2 edge, 3D 4 edges
      real zlv,rmin,rmax

      integer nedge,nel,e
      real tol,vlmin,vlmax,glmin,glmax,rmn,rmx,zlv0

      nedge = (ldim-1)*2
      tol = 1.e-10

      zlv0 = zlv

      ! min control
      if (abs(rmin).gt.tol) then
        do e=1,nel
          rmn = vlmin(edge_len(1,e),nedge)
          zlv = max(zlv, rmin/rmn*2.0-1.0)
        enddo
        zlv = glmax(zlv,1)
        if(nio.eq.0) write(*,23)'rfn lim1 min',zlv0,zlv,rmin
      endif
 
      ! max control
      if (abs(rmax).gt.tol) then
        do e=1,nel
          rmx = vlmax(edge_len(1,e),nedge)
          zlv = min(zlv, rmax/rmx*2.0-1.0)
        enddo
        zlv = glmin(zlv,1)
        if(nio.eq.0) write(*,23)'rfn lim1 max',zlv0,zlv,rmax
      endif

   23 format(a,3g13.5)

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_apply_limiter2_2d(xmc,ymc,Rlim,Nref,nel_lv)
      implicit none
      include 'SIZE'

      real Rlim(4,1)
      real xmc(lx1,ly1,lz1,lelt)
      real ymc(lx1,ly1,lz1,lelt)
      real xhex27(3,3,lelt)
      real yhex27(3,3,lelt)

      real xface_bak(lx1,lelt)
      real yface_bak(lx1,lelt)
      real xface9(3,lelt)
      real yface9(3,lelt)

      real rmin,rmax,tol,dist2d,scale,dis0,dis1,dis2,dis3,dis4,dis5
     $   , xxtop,yytop,xxmid,yymid,xxref0,yyref0,xxbot,yybot
      integer Nref,nel_lv,ilv,ncnt,ncnt_lv,iglsum,e,enow,i,enext,itmp
      integer e_upd(lelt)
      logical ifbot,iftop

      tol = 1.e-10
      ncnt = 0

      ! chk if lim2 is on
      itmp = 0
      do ilv=1,Nref-1
        rmin = Rlim(3,ilv)
        rmax = Rlim(4,ilv)
        if (abs(rmin).gt.tol.OR.abs(rmax).gt.tol) itmp = itmp + 1
      enddo
      itmp = iglsum(itmp,1)
      if(itmp.eq.0) return
      if(nio.eq.0) write(*,28)'rfn lim2 start',itmp

      ! copy original face 1 to hex27
      do e=1,nel_lv
        enow = e + nel_lv ! get first lyr on bdry
        call int_gll_2d_f1(xface9(1,e),3,xmc(1,1,1,enow),lx1)
        call int_gll_2d_f1(yface9(1,e),3,ymc(1,1,1,enow),lx1)
      enddo

      call izero(e_upd,lelt)
      do ilv=1,Nref ! top layer doesn't have new layers
        ncnt_lv = 0

        if (ilv.ne.Nref) then
          rmin = Rlim(3,ilv)
          rmax = Rlim(4,ilv)
        else
          rmin = 0.0
          rmax = 0.0
        endif
  
        do e=1,nel_lv
          enow = ilv*nel_lv + e
          if (ilv.eq.Nref) enow = e
        
          call copy(xface_bak(1,e),xmc(1,1,1,enow),lx1)! backup face 1
          call copy(yface_bak(1,e),ymc(1,1,1,enow),lx1) 
  
          ! get hex27 
          call int_gll_2d_e(xhex27(1,1,enow),3,xmc(1,1,1,enow),lx1)
          call int_gll_2d_e(yhex27(1,1,enow),3,ymc(1,1,1,enow),lx1)
  
          do i=1,3

            xxtop = xhex27(i,3,enow) ! face 3 of current e
            yytop = yhex27(i,3,enow)
            xxmid = xhex27(i,2,enow) ! mid lv of current e
            yymid = yhex27(i,2,enow)
            xxbot = xhex27(i,1,enow) ! face 1 of current e
            yybot = yhex27(i,1,enow)

            if (e_upd(e).gt.0) then ! FIXME force linear ?
              xxmid = (xxtop + xxbot)*0.5
              yymid = (yytop + yybot)*0.5
            endif

            xxref0 = xface9(i,e) ! face 1 of original mesh
            yyref0 = yface9(i,e) 
  
            dis0 = dist2d(xxtop,yytop,xxref0,yyref0) ! lyr to bdry
            dis1 = dist2d(xxtop,yytop,xxbot,yybot) ! lyr to face1
 
            ! scale top lyr 
            scale = 1.0
            if (abs(rmin).gt.tol.AND.dis0.lt.rmin) scale = rmin/dis0
            if (abs(rmax).gt.tol.AND.dis0.gt.rmax) scale = rmax/dis0
            if (abs(1.0-scale).gt.tol) e_upd(e) = 1 ! bottom will also update top lyr
            xxtop = (xxtop - xxref0)*scale + xxref0
            yytop = (yytop - yyref0)*scale + yyref0
 
            ! scale mid lyr 
            dis2 = dist2d(xxtop,yytop,xxbot,yybot) ! lyr to face1, new
            scale = dis2/dis1
            xxmid = (xxmid-xxbot)*scale + xxbot
            yymid = (yymid-yybot)*scale + yybot

            ! chk inverted e
            dis3 = dist2d(xxtop,yytop,xxref0,yyref0) ! top to bdry
            dis4 = dist2d(xxmid,yymid,xxref0,yyref0) ! mid to bdry
            dis5 = dist2d(xxbot,yybot,xxref0,yyref0) ! bot to bdry
            if (dis4.lt.dis5) then ! mid lv is inverted, force linear
              xxmid = (xxtop + xxbot)*0.5
              yymid = (yytop + yybot)*0.5
            endif
            if (dis3.lt.dis4) then
              write(*,*)'rfn lim2, Err: invert e',ilv,e,i
     $                 , dis3,dis4,dis5
              call exitt
            endif
         
            xhex27(i,3,enow) = xxtop
            yhex27(i,3,enow) = yytop
            xhex27(i,2,enow) = xxmid
            yhex27(i,2,enow) = yymid
          enddo

          if (e_upd(e).gt.0) then ! this e is updated
            ncnt_lv = ncnt_lv + 1
            ncnt = ncnt + 1

            ifbot = .false.
            iftop = .false.
            if(ilv.eq.1) ifbot=.true.
            if(ilv.eq.Nref) iftop=.true.

c           call copy(xmc(1,1,1,enow),xface_bak(1,e),lx1)! backup face 1
c           call copy(ymc(1,1,1,enow),yface_bak(1,e),lx1) 

            ! hex27 to lx1, with frozen face 1
            call int_hex27_recover2d_f1(xmc(1,1,1,enow),lx1
     $                                 ,xhex27(1,1,enow),3,ifbot,iftop)
            call int_hex27_recover2d_f1(ymc(1,1,1,enow),lx1
     $                                 ,yhex27(1,1,enow),3,ifbot,iftop)


            call copy(xmc(1,1,1,enow),xface_bak(1,e),lx1)! backup face 1
            call copy(ymc(1,1,1,enow),yface_bak(1,e),lx1)  ! TODO: use fix_geom than this

            ! copy face3 to next layer's face 1
            if(ilv.ne.Nref) then
              enext = enow + nel_lv
              if (ilv.eq.Nref-1) enext = e
              call copy(xmc(1,1,1,enext),xmc(1,ly1,1,enow),lx1)
              call copy(ymc(1,1,1,enext),ymc(1,ly1,1,enow),lx1)
            endif

          endif ! update e
        enddo ! e

        ncnt_lv = iglsum(ncnt_lv,1)
        if(nio.eq.0) write(*,27)'rfn lim2 lv=',ilv,ncnt_lv

      enddo ! lv

      ncnt = iglsum(ncnt,1)
      if(nio.eq.0) write(*,28)'rfn lim2 all=',ncnt

   27 format(a,2i3)
   28 format(a,i3)
      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_apply_limiter2_3d(xmc,ymc,zmc,Rlim,Nref,nel_lv)
      implicit none
      include 'SIZE'

      real Rlim(4,1)
      real xmc(lx1,ly1,lz1,lelt)
      real ymc(lx1,ly1,lz1,lelt)
      real zmc(lx1,ly1,lz1,lelt)
      real xhex27(9,3,lelt), xface9(9,lelt), xface_bak(lx1,lelt)
      real yhex27(9,3,lelt), yface9(9,lelt), yface_bak(lx1,lelt)
      real zhex27(9,3,lelt), zface9(9,lelt), zface_bak(lx1,lelt)

      real rmin,rmax,tol
     $   , dist3d,scale,dis0,dis1,dis2,dis3,dis4,dis5
     $   , xxtop,yytop,zztop,xxmid,yymid,zzmid
     $   , xxref0,yyref0,zzref0,xxbot,yybot,zzbot
      integer Nref,nel_lv
      integer ilv,ncnt,ncnt_lv,iglsum,e,enow,i,enext,itmp,e_upd(lelt)
      logical ifbot,iftop

      tol = 1.e-10
      ncnt = 0

      ! chk if lim2 is on
      itmp = 0
      write(*,*)'DDDD',Nref
      do ilv=1,Nref-1
        rmin = Rlim(3,ilv)
        rmax = Rlim(4,ilv)
        if (abs(rmin).gt.tol.OR.abs(rmax).gt.tol) itmp = itmp + 1
        write(*,*)'DDDD',nid,ilv,rmin,rmax,itmp
      enddo
      itmp = iglsum(itmp,1)
      if(itmp.eq.0) return
      if(nio.eq.0) write(*,28)'rfn lim2 start',itmp

      ! copy original face 1 to hex27
      do e=1,nel_lv
        enow = e + nel_lv ! get first lyr on bdry ! TODO
        call int_gll_3d_f5(xface9(1,e),3,xmc(1,1,1,enow),lx1)
        call int_gll_3d_f5(yface9(1,e),3,ymc(1,1,1,enow),lx1)
        call int_gll_3d_f5(zface9(1,e),3,zmc(1,1,1,enow),lx1)
      enddo

      call izero(e_upd,lelt)
      do ilv=1,Nref ! top layer doesn't have new layers
        ncnt_lv = 0

        if (ilv.ne.Nref) then
          rmin = Rlim(3,ilv)
          rmax = Rlim(4,ilv)
        else
          rmin = 0.0
          rmax = 0.0
        endif
  
        do e=1,nel_lv
          enow = ilv*nel_lv + e
          if (ilv.eq.Nref) enow = e
        
          call copy(xface_bak(1,e),xmc(1,1,1,enow),lx1*lz1)! backup face 5
          call copy(yface_bak(1,e),ymc(1,1,1,enow),lx1*lz1) 
          call copy(zface_bak(1,e),zmc(1,1,1,enow),lx1*lz1) 
  
          ! get hex27  ! TODO
          call int_gll_3d_e(xhex27(1,1,enow),3,xmc(1,1,1,enow),lx1)
          call int_gll_3d_e(yhex27(1,1,enow),3,ymc(1,1,1,enow),lx1)
          call int_gll_3d_e(zhex27(1,1,enow),3,zmc(1,1,1,enow),lx1)
  
          do i=1,9

            xxtop = xhex27(i,3,enow) ! face 6 of current e
            yytop = yhex27(i,3,enow)
            zztop = zhex27(i,3,enow)
            xxmid = xhex27(i,2,enow) ! mid lv of current e
            yymid = yhex27(i,2,enow)
            zzmid = zhex27(i,2,enow)
            xxbot = xhex27(i,1,enow) ! face 5 of current e
            yybot = yhex27(i,1,enow)
            zzbot = zhex27(i,1,enow)

            if (e_upd(e).gt.0) then ! FIXME force linear ?
              xxmid = (xxtop + xxbot)*0.5
              yymid = (yytop + yybot)*0.5
              zzmid = (zztop + zzbot)*0.5
            endif

            xxref0 = xface9(i,e) ! face 1 of original mesh
            yyref0 = yface9(i,e) 
            zzref0 = zface9(i,e) 
  
            dis0 = dist3d(xxtop,yytop,zztop,xxref0,yyref0,zzref0) ! lyr to bdry
            dis1 = dist3d(xxtop,yytop,zztop,xxbot,yybot,zzbot) ! lyr to face1
 
            ! scale top lyr 
            scale = 1.0
            if (abs(rmin).gt.tol.AND.dis0.lt.rmin) scale = rmin/dis0
            if (abs(rmax).gt.tol.AND.dis0.gt.rmax) scale = rmax/dis0
            if (abs(1.0-scale).gt.tol) e_upd(e) = 1 ! bottom will also update top lyr
            xxtop = (xxtop - xxref0)*scale + xxref0
            yytop = (yytop - yyref0)*scale + yyref0
            zztop = (zztop - zzref0)*scale + zzref0
 
            ! scale mid lyr 
            dis2 = dist3d(xxtop,yytop,zztop,xxbot,yybot,zzbot) ! lyr to face1, new
            scale = dis2/dis1
            xxmid = (xxmid-xxbot)*scale + xxbot
            yymid = (yymid-yybot)*scale + yybot
            zzmid = (zzmid-zzbot)*scale + zzbot

            ! chk inverted e
            dis3 = dist3d(xxtop,yytop,zztop,xxref0,yyref0,zzref0) ! top to bdry
            dis4 = dist3d(xxmid,yymid,zzmid,xxref0,yyref0,zzref0) ! mid to bdry
            dis5 = dist3d(xxbot,yybot,zzbot,xxref0,yyref0,zzref0) ! bot to bdry
            if (dis4.lt.dis5) then ! mid lv is inverted, force linear
              xxmid = (xxtop + xxbot)*0.5
              yymid = (yytop + yybot)*0.5
              zzmid = (zztop + zzbot)*0.5
              dis4 = dist3d(xxmid,yymid,zzmid,xxref0,yyref0,zzref0) ! mid to bdry
            endif
            if (dis3.lt.dis4) then
              write(*,*)'rfn lim2, Err: invert e',ilv,nid,e,i
     $                 , dis3,dis4,dis5
              call exitt
            endif
         
            xhex27(i,3,enow) = xxtop
            yhex27(i,3,enow) = yytop
            zhex27(i,3,enow) = zztop
            xhex27(i,2,enow) = xxmid
            yhex27(i,2,enow) = yymid
            zhex27(i,2,enow) = zzmid
          enddo

          if (e_upd(e).gt.0) then ! this e is updated
            ncnt_lv = ncnt_lv + 1
            ncnt = ncnt + 1

            ifbot = .false.
            iftop = .false.
            if(ilv.eq.1) ifbot=.true.
            if(ilv.eq.Nref) iftop=.true.

c           call copy(xmc(1,1,1,enow),xface_bak(1,e),lx1)! backup face 1
c           call copy(ymc(1,1,1,enow),yface_bak(1,e),lx1) 
c           call copy(zmc(1,1,1,enow),zface_bak(1,e),lx1) 

            ! hex27 to lx1, with frozen face 5 ! TODO
            call int_hex27_recover3d_f5(xmc(1,1,1,enow),lx1
     $                                 ,xhex27(1,1,enow),3,ifbot,iftop)
            call int_hex27_recover3d_f5(ymc(1,1,1,enow),lx1
     $                                 ,yhex27(1,1,enow),3,ifbot,iftop)
            call int_hex27_recover3d_f5(zmc(1,1,1,enow),lx1
     $                                 ,zhex27(1,1,enow),3,ifbot,iftop)


            call copy(xmc(1,1,1,enow),xface_bak(1,e),lx1*lz1)! backup face 5
            call copy(ymc(1,1,1,enow),yface_bak(1,e),lx1*lz1)  ! TODO: use fix_geom than this
            call copy(zmc(1,1,1,enow),zface_bak(1,e),lx1*lz1)

            ! copy face6 to next layer's face 5
            if(ilv.ne.Nref) then
              enext = enow + nel_lv
              if (ilv.eq.Nref-1) enext = e
              call copy(xmc(1,1,1,enext),xmc(1,1,lz1,enow),lx1*lz1)
              call copy(ymc(1,1,1,enext),ymc(1,1,lz1,enow),lx1*lz1)
              call copy(zmc(1,1,1,enext),zmc(1,1,lz1,enow),lx1*lz1)
            endif

          endif ! update e
        enddo ! e

        ncnt_lv = iglsum(ncnt_lv,1)
        if(nio.eq.0) write(*,27)'rfn lim2 lv=',ilv,ncnt_lv

      enddo ! lv

      ncnt = iglsum(ncnt,1)
      if(nio.eq.0) write(*,28)'rfn lim2 all=',ncnt

   27 format(a,i3i9)
   28 format(a,i9)
      return
      end
c-----------------------------------------------------------------------
      subroutine int_gll_pt_2d_e(uf,z0,z1,n,uc,m) 
c     [z0,z1] in [-1,1]
c     interpolate uc(m,m) -> uf(m,n) 
c          from 1D-GLL(m) -> 1D-GLL rescale to [z0,z1] along s-axis

      parameter (l=50)
      real uf(m,n),uc(m,m),j,jt,z0,z1,wtmp(l)
      common /my_cmap2d_gll/ j(l*l),jt(l*l),w(l*l),z(l)

      integer mo,no,to1,to2,itmp1,itmp2
      save    mo,no,to1,to2
      data    mo,no,to1,to2 / 0,0,0,0 /

      itmp1 = int(z0*1.e4)
      itmp2 = int(z1*1.e4)
      if (m.ne.mo .or. n.ne.no .or. itmp1.ne.to1 .or. itmp2.ne.to2) then
          if (m.gt.l) call exitti('int_gll_pt_2d_e memory 1$',m)
          if (n.gt.l) call exitti('int_gll_pt_2d_e memory 2$',n)

          call zwgll (z,wtmp,m)
          call zwgll (w,wtmp,n)

          do i=1,n
            w(i) = (w(i)+1.0)*0.5*(z1-z0) + z0
          enddo

          call gen_int_gz(j,jt,w,n,z,m)
          
          to1=itmp1
          to2=itmp2
          mo=m
          no=n
      endif

c     call mxm(j,n,uc,m,w ,m)
c     call mxm(w,n,jt,m,uf,n)

      call mxm(uc,m,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_gll_pt_3d_e(uf,z0,z1,n,uc,m) 
c     [z0,z1] in [-1,1]
c     interpolate uc(m,m,m) -> uf(m,m,n) 
c          from 1D-GLL(m) -> 1D-GLL rescale to [z0,z1] along s-axis

      parameter (l=50)
      real uf(m,m,n),uc(m,m,m),j,jt,z0,z1,wtmp(l)
      common /my_cmap3d_gll/ j(l*l),jt(l*l),w(l*l),z(l)

      integer mo,no,to1,to2,itmp1,itmp2
      save    mo,no,to1,to2
      data    mo,no,to1,to2 / 0,0,0,0 /

      itmp1 = int(z0*1.e4)
      itmp2 = int(z1*1.e4)
      if (m.ne.mo .or. n.ne.no .or. itmp1.ne.to1 .or. itmp2.ne.to2) then
          if (m.gt.l) call exitti('int_gll_pt_3d_e memory 1$',m)
          if (n.gt.l) call exitti('int_gll_pt_3d_e memory 2$',n)

          call zwgll (z,wtmp,m)
          call zwgll (w,wtmp,n)

          do i=1,n
            w(i) = (w(i)+1.0)*0.5*(z1-z0) + z0
          enddo

          call gen_int_gz(j,jt,w,n,z,m)
          
          to1=itmp1
          to2=itmp2
          mo=m
          no=n
      endif

      mm = m*m
      mn = m*n
      nn = n*n

      call mxm(uc,mm,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_gll_2d_e(uf,n,uc,m)
c     interpolate uc(m,m) -> uf(n,n) 
c          from 1D-GLL(m) -> 1D-GLL(n)

      parameter (l=50)
      real uf(n,n),uc(m,m),j,jt,wtmp(l)
      common /my_cmap2d_e/ j(l*l),jt(l*l),w(l*l),z(l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.gt.l) call exitti('int_gll_pt_2d_e memory 1$',m)
      if (n.gt.l) call exitti('int_gll_pt_2d_e memory 2$',n)

      if (m.ne.mo .or. n.ne.no) then
          call zwgll (z,wtmp,m)
          call zwgll (w,wtmp,n)

          call gen_int_gz(j,jt,w,n,z,m)

          mo=m
          no=n
      endif

      call mxm(j,n,uc,m,w ,m)
      call mxm(w,n,jt,m,uf,n)

c     call mxm(j,n,uc,m,uf,m)
c     call mxm(uc,m,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_gll_3d_e(uf,n,uc,m)
c     interpolate uc(m,m) -> uf(n,n) 
c          from 1D-GLL(m) -> 1D-GLL(n)

      parameter (l=50)
      real uf(n,n,n),uc(m,m,m),j,jt,wtmp(l),v
      common /my_cmap3d_e/ j(l*l),jt(l*l),w(l*l*l),z(l),v(l*l*l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no) then
          if (m.gt.l) call exitti('int_gll_3d_e memory 1$',m)
          if (n.gt.l) call exitti('int_gll_3d_e memory 2$',n)

          call zwgll (z,wtmp,m)
          call zwgll (w,wtmp,n)

          call gen_int_gz(j,jt,w,n,z,m)

          mo=m
          no=n
      endif

      mm = m*m
      mn = m*n
      nn = n*n

      call mxm(j,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxm(v(iv),n,jt,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxm(w,nn,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_gll_2d_f1(uf,n,uc,m)
c     interpolate uc(m,1) -> uf(n,1) 
c          from 1D-GLL(m) -> 1D-GLL(n)

      parameter (l=50)
      real uf(n,1),uc(m,1),j,jt,wtmp(l)
      common /my_cmap2d_f/ j(l*l),jt(l*l),w(l*l),z(l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no) then
          if (m.gt.l) call exitti('int_gll_2d_f1 memory 1$',m)
          if (n.gt.l) call exitti('int_gll_2d_f1 memory 2$',n)

          call zwgll (z,wtmp,m)
          call zwgll (w,wtmp,n)

          call gen_int_gz(j,jt,w,n,z,m)

          mo=m
          no=n
      endif

c     call mxm(j,n,uc,m,w ,m)
c     call mxm(w,n,jt,m,uf,n)

c     mxm(a,n1,b,n2,c,n3)
c     Compute matrix-matrix product C = A*B
c     for contiguously packed matrices A,B, and C.
c     real a(n1,n2),b(n2,n3),c(n1,n3)

c      call mxm(j,n,uc,m,uf,n)

      call mxm(j,n,uc,m,uf,1)

c      write(*,*)'dbg9 0',m,mo,n,no
c     call mxm(uc,m,jt,m,uf,n)
c      write(*,*)'dbg9 1',(uc(i,1),i=1,m)
c      write(*,*)'dbg9 2',(jt(i),i=1,m*n)
c      write(*,*)'dbg9 3',(uf(i,1),i=1,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_gll_3d_f5(uf,n,uc,m)
c     interpolate uc(m,1) -> uf(n,1) 
c          from 1D-GLL(m) -> 1D-GLL(n)

      parameter (l=50)
      real uf(n,n,1),uc(m,m,1),j,jt,wtmp(l),v
      common /my_cmap3d_f/ j(l*l),jt(l*l),w(l*l),z(l),v(l*l*l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no) then
          if (m.gt.l) call exitti('int_gll_3d_f5 memory 1$',m)
          if (n.gt.l) call exitti('int_gll_3d_f5 memory 2$',n)

          call zwgll (z,wtmp,m)
          call zwgll (w,wtmp,n)

          call gen_int_gz(j,jt,w,n,z,m)

          mo=m
          no=n
      endif

      call mxm(j,n,uc,m,v ,m)
      call mxm(v,n,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_hex27_recover2d_f1(uf,n,uc,m,ifface1,ifface3)
c     interpolate uc(m,m) -> uf(n,n) 
c          from 1D-GLL(m) -> 1D-GLL(n)
c          ifface=T, then freeze face 1 from uf

      parameter (l=50)
      real uf(n,n),uc(m,m),j,jt,ud(n,n),uf0(n,n),wtmp(l),zg(l),uc0(m,m)
      common /my_cmap2d_f2/ j(l*l),jt(l*l),w(l*l),z(l)

      real w1(n*n*n),w2(n*n*n)
      logical ifface1,ifface3

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no) then
          if (m.gt.l) call exitti('int_hex27_recover2d_f1 memory 1$',m)
          if (n.gt.l) call exitti('int_hex27_recover2d_f1 memory 2$',n)

          call zwgll (z,wtmp,m)
          call zwgll (w,wtmp,n)

          call gen_int_gz(j,jt,w,n,z,m)

          mo=m
          no=n
      endif


      call copy(uf0,uf,n*n) ! backup uf
      call copy(uc0,uc,m*m) ! backup uc

      call mxm(j,n,uc,m,w ,m)
      call mxm(w,n,jt,m,uf,n) ! interp uc to uf

c      if(.not.ifface) return


      call zwgll(zg,wtmp,n)

      if(ifface1)then ! FIXME some how thih works... idk why... 
        do k1=1,3 ! gorden-hall mapping
          call rzero(ud,n*n)
          call sub3(ud,uf,uf0,n*n) ! monotone to monotone ! FIXME
          if(ifface1) call rzero(ud(1,1),n) ! freeze face 1, uf = uf0 = bottom
          if(ifface3) call rzero(ud(1,n),n) ! freeze face 1, uf = uf0 = bottom

          call gh_face_extend(ud,zg,n,k1,w1,w2)

          call add3(uf,uf0,ud,n*n)
          if(ifface1) call copy(uf(1,1),uf0(1,1),n)
          if(ifface3) call copy(uf(1,n),uf0(1,n),n)
        enddo
      else
        return
        do k1=1,3 ! gorden-hall mapping
          call rzero(ud,n*n)
          call copy(ud,uf,n*n)
          call copy(uf(1,1),uf0(1,1),n)
   
          call gh_face_extend(ud,zg,n,k1,w1,w2)

          call copy(uf,ud,n*n)
          call copy(uf(1,1),uf0(1,1),n)
        enddo

c        do k1=1,3 ! gorden-hall mapping
c          call rzero(ud,n*n)
c          call sub3(ud,uf,uf0,n*n) ! monotone to monotone ! FIXME
c          call rzero(ud(1,1),n) ! freeze face 1, uf = uf0 = bottom
c
c          call gh_face_extend(ud,zg,n,k1,w1,w2)
c
c          call add3(uf,uf0,ud,n*n)
c          call copy(uf(1,1),uf0(1,1),n)
c        enddo
      
      endif

c      write(*,*)'dd 0',(uf0(i,1),i=1,n*n)
c      write(*,*)'dd 1',(uf(i,1),i=1,n*n)
c      write(*,*)'dd 2',(ud(i,1),i=1,n*n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_hex27_recover3d_f5(uf,n,uc,m,ifface5,ifface6)
c     interpolate uc(m,m) -> uf(n,n) 
c          from 1D-GLL(m) -> 1D-GLL(n)
c          ifface=T, then freeze face 5 from uf

      parameter (l=50)
      real uf(n,n,n),uc(m,m,m),j,jt
     $   , ud(n,n,n),uf0(n,n,n),wtmp(l),zg(l),uc0(m,m)
      common /my_cmap3d_f2/ j(l*l),jt(l*l),w(l*l*l),v(l*l*l),z(l)

      real w1(n*n*n),w2(n*n*n)
      logical ifface5,ifface6

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no) then
          if (m.gt.l) call exitti('int_hex27_recover3d_f5 memory 1$',m)
          if (n.gt.l) call exitti('int_hex27_recover3d_f5 memory 2$',n)

          call zwgll (z,wtmp,m)
          call zwgll (w,wtmp,n)

          call gen_int_gz(j,jt,w,n,z,m)

          mo=m
          no=n
      endif

      mm = m*m
      mn = m*n
      nn = n*n

      call copy(uf0,uf,n*n*n) ! backup uf
      call copy(uc0,uc,m*m*m) ! backup uc

      ! interp uc to uf
      call mxm(j,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxm(v(iv),n,jt,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxm(w,nn,jt,m,uf,n)


c      if(.not.ifface) return

      call zwgll(zg,wtmp,n)

      if(ifface5)then ! FIXME some how thih works... idk why... 
        do k1=1,3 ! gorden-hall mapping
          call rzero(ud,n*n*n)
          call sub3(ud,uf,uf0,n*n*n) ! monotone to monotone ! FIXME
          if(ifface5) call rzero(ud(1,1,1),n*n) ! freeze face 1, uf = uf0 = bottom
          if(ifface6) call rzero(ud(1,1,n),n*n) ! freeze face 1, uf = uf0 = bottom

          call gh_face_extend(ud,zg,n,k1,w1,w2)

          call add3(uf,uf0,ud,n*n*n)
          if(ifface5) call copy(uf(1,1,1),uf0(1,1,1),n*n)
          if(ifface6) call copy(uf(1,1,n),uf0(1,1,n),n*n)
        enddo
      else
        return
        do k1=1,3 ! gorden-hall mapping
          call rzero(ud,n*n*n)
          call copy(ud,uf,n*n*n)
          call copy(uf(1,1,1),uf0(1,1,1),n*n)
   
          call gh_face_extend(ud,zg,n,k1,w1,w2)

          call copy(uf,ud,n*n*n)
          call copy(uf(1,1,1),uf0(1,1,1),n*n)
        enddo

c        do k1=1,3 ! gorden-hall mapping
c          call rzero(ud,n*n)
c          call sub3(ud,uf,uf0,n*n) ! monotone to monotone ! FIXME
c          call rzero(ud(1,1),n) ! freeze face 1, uf = uf0 = bottom
c
c          call gh_face_extend(ud,zg,n,k1,w1,w2)
c
c          call add3(uf,uf0,ud,n*n)
c          call copy(uf(1,1),uf0(1,1),n)
c        enddo
      
      endif

c      write(*,*)'dd 0',(uf0(i,1),i=1,n*n)
c      write(*,*)'dd 1',(uf(i,1),i=1,n*n)
c      write(*,*)'dd 2',(ud(i,1),i=1,n*n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_recover2d_from_f13(u,n)
      implicit none
      real u(n,n),u0(n,n),ud(n,n),uc(2,2)

      integer l,n,k1
      parameter(l=50)
      real wtmp(l),zg(l),w1(n*n*n),w2(n*n*n),j,jt,w,zn,z2
      common /my_cmap2d_f13/ j(l*l),jt(l*l),w(l*l),zn(l),z2(l)

      integer icalld
      save icalld
      data icalld /0/

      call copy(u0,u,n*n)

      ! Step 1: regen linear mesh from vertices
      if (icalld.eq.0) then
         call zwgll (z2,wtmp,2)
         call zwgll (zn,wtmp,n)

         call gen_int_gz(j,jt,zn,n,z2,2) ! from 2 to n

         icalld = 1
      endif

      ! Get vertix 
      uc(1,1) = u(1,1)
      uc(1,2) = u(1,n)
      uc(2,1) = u(n,1)
      uc(2,2) = u(n,n)

      call mxm(j,n,uc,2,w,2)
      call mxm(w,n,jt,2,u,n) ! interp uc to u

      ! Step 2: fix curved sides from face 1 and face 3
      do k1=1,3 ! gorden-hall mapping
        call rzero(ud,n*n)
        call sub3(ud(1,1),u0(1,1),u(1,1),n) ! face 1
        call sub3(ud(1,n),u0(1,n),u(1,n),n) ! face 3

        call gh_face_extend(ud,zn,n,k1,w1,w2)

        call add2(u,ud,n*n)
      enddo

      call copy(u(1,1),u0(1,1),n)
      call copy(u(1,n),u0(1,n),n)

      return
      end
c-----------------------------------------------------------------------
      subroutine int_recover3d_from_f56(u,n)
      implicit none
      real u(n,n,n),u0(n,n,n),ud(n,n,n),uc(2,2,2)

      integer l,l3,n,k,iv,iw,nn,mn,mm,m
      parameter(l=50,l3=l*l*l)
      real wtmp(l),zg(l),w1(n*n*n),w2(n*n*n),j,jt,v,w,zn,z2
      common /my_cmap3d_f56/ j(l*l),jt(l*l),w(l3),v(l3),zn(l),z2(l)

      integer icalld
      save icalld
      data icalld /0/

      call copy(u0,u,n*n*n)

      ! Step 1: regen linear mesh from vertices
      if (icalld.eq.0) then
         call zwgll (z2,wtmp,2)
         call zwgll (zn,wtmp,n)

         call gen_int_gz(j,jt,zn,n,z2,2) ! from 2 to n

         icalld = 1
      endif

      ! Get vertix, uc=Hex8
      uc(1,1,1) = u(1,1,1)
      uc(1,2,1) = u(1,n,1)
      uc(2,1,1) = u(n,1,1)
      uc(2,2,1) = u(n,n,1)
      uc(1,1,2) = u(1,1,n)
      uc(1,2,2) = u(1,n,n)
      uc(2,1,2) = u(n,1,n)
      uc(2,2,2) = u(n,n,n)

c      call mxm(j,n,uc,2,w,2)
c      call mxm(w,n,jt,2,u,n) ! interp uc to u

      ! interp uc to u
      m  = 2
      mm = m*m
      mn = m*n
      nn = n*n
      call mxm(j,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxm(v(iv),n,jt,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxm(w,nn,jt,m,u,n)

      ! Step 2: fix curved sides from face 5 and face 6
      do k=1,3 ! gorden-hall mapping
        call rzero(ud,n*n*n)
        call sub3(ud(1,1,1),u0(1,1,1),u(1,1,1),n*n) ! face 5
        call sub3(ud(1,1,n),u0(1,1,n),u(1,1,n),n*n) ! face 6

        call gh_face_extend(ud,zn,n,k,w1,w2)

        call add2(u,ud,n*n*n)
      enddo

      call copy(u(1,1,1),u0(1,1,1),n*n)
      call copy(u(1,1,n),u0(1,1,n),n*n)

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_chk_idx(s3,cbc3,bID)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      character*3 s3,cbc3
      integer e,f,ix,iy,iz,n,bID

      real vxtmp(lx1,ly1,lz1,lelv)
     $   , vytmp(lx1,ly1,lz1,lelv)
     $   , vztmp(lx1,ly1,lz1,lelv)
     $   , prtmp(lx2,ly2,lz2,lelv)
     $   , ttmp (lx1,ly1,lz1,lelt,ldimt)

      n=lx1*ly1*lz1*nelt

      call rzero(vxtmp,n)
      call rzero(vytmp,n)
      call rzero(vztmp,n)
      call rzero(prtmp,n)
      call rzero(ttmp ,n)

c     debug plot
      do e=1,nelt
        do iz=1,lz1
        do iy=1,ly1
        do ix=1,lx1
          vxtmp(ix,iy,iz,e) = real(ix)
          vytmp(ix,iy,iz,e) = real(iy)
          vztmp(ix,iy,iz,e) = real(iz)
        enddo
        enddo
        enddo

        do f=1,ldim*2
          if (bID.gt.0) then
            if (boundaryID(f,e).eq.bID) then
              call facev(prtmp,e,f,1.0,lx1,ly1,lz1)
            endif
          else
            if (cbc(f,e,1).eq.cbc3) then
              call facev(prtmp,e,f,1.0,lx1,ly1,lz1)
            endif
          endif

          call facev(ttmp,e,f,real(f),lx1,ly1,lz1)
        enddo
      enddo

      ifxyo = .true.
      call outpost(vxtmp,vytmp,vztmp,prtmp,ttmp,s3)
      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_dump_mesh(s3,ifbinary) ! TODO: remove my_errchk??
      implicit none
      integer ierr
      logical ifbinary
      character*3 s3

      ! swap cbc
      call rfn_backup_boundary_id(2)   ! backup current cbc
      call rfn_recover_boundary_id(1)  ! get original cbc from mesh

      ! chk mesh ! my_errchk.f
      call my_mesh_metric
      call my_mshchk_vb(s3//'   ',.false.,ierr)
      call my_cbc_chk(s3)

      ! dump mesh
      if (ifbinary) then
        call my_gen_re2(s3,2) 
c        call my_gen_re2(s3,1) 
c        call my_gen_re2(s3,0) 
        call rfn_gen_co2(s3)
      else
        call my_gen_rea(s3,1)
        call rfn_gen_con(s3)
      endif


      call rfn_recover_boundary_id(2)  ! recover current cbc

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_gen_co2(s3)
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
        write(6,*)'rfn_gen_co2 ',hdr
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
      subroutine rfn_gen_con(s3)
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
        write(6,*)'rfn_gen_con ',version,nelgt,nelgv,nv
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
