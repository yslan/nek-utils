c     In userdat2:
c        call rfn_backup_boundary_id(1)
c        ! set bdry for simulation
c     In userchk:
c        call rfn_split
c        call rfn_dump_mesh
c     TODO:
c        cbc ifield=2
c        if co2/ma2
      include 'rfn_rotate.f' 
      include 'rfn_interpolate.f'
      include 'rfn_ssplit.f'
      include 'my_errchk.f' ! mesh chk
      include 'my_post.f'   ! selected outpost
c-----------------------------------------------------------------------
      subroutine rfn_ssplit(cbc3,bID,Ncut) ! surface split
      implicit none
      include 'SIZE'
      include 'PARALLEL' ! nelgv
      include 'REFINE'
      character*3 cbc3
      integer Nref,ie_ref(lelt),bID,Ncut,Nseg,Nnew_p_e
      real Rref(lref-1),Rlim(4,lref-1)

      Nseg = Ncut + 1 ! # segment per direction
      Nnew_p_e = Nseg*Nseg-1 ! each element will create m^2-1 new elements
   
      if (ldim.ne.3) call exitti('Err: ssplit only works in 3D$',1)
      if (Ncut.eq.0) return

      call rfn_init
      call rfn_tag_e(ie_ref,cbc3,bID,1) ! tag bdry elements
      call rfn_trace_layers(ie_ref)     ! propagate the tag
      call rfn_bak_nelg(ie_ref,3,Ncut)  ! bak nel, set nel_lv
      call rfn_chk_Ncut(Ncut)           ! also chk size (lelg,lelv)
      call rfn_set_nel_lv_shift         ! prepare new element id offset

      ! copy coordinates to local work array (rotated)
      ! rotate s.t. face 1 (2D) or face 5 (3D) will always face the target bdry
      if (ldim.eq.2) then
        call rfn_rotate2d_in(ie_ref)
      else
        call rfn_rotate3d_in(ie_ref)
      endif

      ! split
      call rfn_ssplit_work(ie_ref,Ncut)

      call update_nels(Nnew_p_e)

      ! copy coordinates back to xm1, ym1, zm1 (rotate back)
      if (ldim.eq.2) then
        call rfn_rotate2d_out(ie_ref)
      else
        call rfn_rotate3d_out(ie_ref)
      endif

      if (nid.eq.0) then
        write(*,26)'rfn SSPLIT END,  E_new=',nelgv
        write(*,'(a)')'rfn'
      endif

   26 format((ai9))
      return
      end
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
      call rfn_init                          ! zero global arrays
      call rfn_chk_Rref(Nref,Rref,Rlim)      ! chk refine inputs
      call rfn_tag_e(ie_ref,cbc3,bID,Nref)   ! tag elements
      call rfn_bak_nelg(ie_ref,1,Nref)       ! save nelgv/nelgt 

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

      call rfn_init
      call rfn_tag_e(ie_ref,cbc3,bID,1)
      call rfn_bak_nelg(ie_ref,2,1)

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
      subroutine rfn_init
      implicit none
      include 'SIZE'
      include 'REFINE'
      integer ifld,e,f

c     Initialize global variables
      call rzero(xmc,lx1*ly1*lz1*lelt) 
      call rzero(ymc,lx1*ly1*lz1*lelt) 
      call rzero(zmc,lx1*ly1*lz1*lelt) 
      call i8zero(vertexc,lxv*lyv*lzv*lelt)
      call i8zero(vertexu,lxv*lzv*lelt)

      do ifld=1,2
        do e=1,lelt
        do f=1,6
          bID_wk(f,e,ifld) = 0
          cbc_wk(f,e,ifld) = 'E  '
        enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_tag_e(ie_ref,cbc3,bID,Nref)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      character*3 cbc3
      integer Nref,ie_ref(lelt),ie_test(lelt),bID
      integer ifld,e,f,iglmax,iglsum,ivlsum,ierr,ntmp

      ! tag e that needs refine, assign face id to ie_ref
      call izero(ie_ref,lelv)
      call izero(ie_test,lelv)
      if (bID.gt.0) then
        do e=1,nelt
        do f=1,ldim*2
          if (boundaryID(f,e).eq.bID) then
            ie_ref(e) = f
            ie_test(e) = ie_test(e) + 1
          endif
        enddo
        enddo
      else
        do e=1,nelt
        do f=1,ldim*2
          if (cbc(f,e,1).eq.cbc3) then
            ie_ref(e) = f
            ie_test(e) = ie_test(e) + 1
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
      do e=1,nelt
      do f=1,ldim*2
        if (cbc(f,e,1).eq.'P  '.or.cbc(f,e,1).eq.'p  ') then
          ierr=1
        endif
      enddo
      enddo
      ierr=iglsum(ierr,1)
      if (ierr.gt.0) 
     $   call exitti('Err rfn: periodic BC not supported yet$',ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine update_face_tag(face_tag,ie_ref,istatus)
      implicit none
      include 'SIZE'
      include 'PARALLEL' ! lglel
      integer face_tag(2*ldim,lelt),ie_ref(lelt)
      integer e,f,fo,fa,istatus,i,ftmp0,ftmp1,ftmp2,ntmp1,ntmp2,iglmax
      
      ! opposite dace id
      integer fopp(6), fadj(4,6) ! opposite face id, adjacent face id
      data    fopp / 3,4,1,2,6,5 /
      data    fadj / 2,4,5,6 
     $             , 1,3,5,6 
     $             , 2,4,5,6 
     $             , 1,3,5,6 
     $             , 1,2,3,4 
     $             , 1,2,3,4 /

      ! from ie_ref to face_tag
      do e=1,nelv
        if (ie_ref(e).gt.0) then
          f = ie_ref(e)
          fo = fopp(f)
          face_tag(f,e) = 1
          face_tag(fo,e) = 1
          do i=1,4   
            fa = fadj(i,f)
            face_tag(fa,e) = 2
          enddo
        endif
      enddo

      ! from face_tag to face_tag via opposite
      do e=1,nelv
      do f=1,2*ldim
        fo = fopp(f)
        if (face_tag(f,e).eq.1) then
          face_tag(fo,e) = 1
        elseif (face_tag(f,e).eq.2) then
          face_tag(fo,e) = 2
        endif
      enddo
      enddo

      ! fill adj face tag + update ie_ref
      istatus = 0
      do e=1,nelv

        ntmp1=0
        ntmp2=0
        ftmp0=0
        ftmp1=0
        ftmp2=0
        do f=1,2*ldim
          if (face_tag(f,e).eq.1) then
            ntmp1 = ntmp1 + 1
            ftmp1 = f
          elseif (face_tag(f,e).eq.1) then
            ntmp2 = ntmp2 + 1
            ftmp2 = f
          else ! face_tag = 0
            ftmp0 = f
          endif
        enddo
        if     (ntmp1.eq.2) then
          if (ie_ref(e).eq.0) ie_ref(e) = ftmp1
          face_tag(ftmp1,e) = 1
          face_tag(fopp(ftmp1),e) = 1
          do i=1,4
            fa = fadj(i,ftmp1) 
            face_tag(fa,e) = 2
          enddo
        elseif (ntmp2.eq.4.AND.ntmp1.eq.0) then
          ie_ref(e) = ftmp0
          face_tag(ftmp0,e) = 1
          face_tag(fopp(ftmp0),e) = 1
          do i=1,4
            fa = fadj(i,ftmp0) 
            face_tag(fa,e) = 2
          enddo
        elseif (ntmp2.eq.2.AND.ntmp1.eq.0) then ! ok for now, err at the end
          istatus = 1
        elseif (ntmp2.eq.0.AND.ntmp1.eq.0) then 
          ! do nothing
        else
          write(*,*)'Err: face types',lglel(e),ntmp1,ntmp2
          istatus = 2
        endif
        
      enddo
      istatus = iglmax(istatus,1)
      if (istatus.eq.2) call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_trace_layers(ie_ref)
      implicit none
      include 'SIZE'
      include 'PARALLEL' ! nelgv
      integer ie_ref(lelt),e,f,fo,fa,n,i,ix,iy,iz
      integer iter,max_iter,ierr,istatus,ntag0,ntag1,ntag,ndif,iglsum
      
      integer face_tag(2*ldim,lelt), face_tag0(2*ldim,lelt)
      real tmpo(lx1,ly1,lz1,lelt)
      real tmpa(lx1,ly1,lz1,lelt)
      
      ! opposite dace id
      integer fopp(6), fadj(4,6) ! opposite face id, adjacent face id
      data    fopp / 3,4,1,2,6,5 /
      data    fadj / 2,4,5,6
     $             , 1,3,5,6 
     $             , 2,4,5,6 
     $             , 1,3,5,6 
     $             , 1,2,3,4 
     $             , 1,2,3,4 /

      ! (one of the) center point of face
      integer ix_fctr(6), iy_fctr(6), iz_fctr(6), idd
      parameter(idd=ldim-1)
      data ix_fctr /  2,lx1,  2,  1,  2,  2/
      data iy_fctr /  1,  2,ly1,  2,  2,  2/
      data iz_fctr /idd,idd,idd,idd,  1,lz1/

      ! TODO: use gs dg face
      if (lx1.lt.3) call exitti('WARN: trace layer needs lx1>=3 $',1)

      n = lx1*ly1*lz1*nelv
      max_iter = nelgt ! worse case senario
      ierr = 0
      call izero(face_tag0,2*ldim*lelt)
      call izero(face_tag,2*ldim*lelt)

      call update_face_tag(face_tag,ie_ref,istatus)
      call icopy(face_tag0,face_tag0,2*ldim*lelt)

      ntag = 0
      do e=1,nelv
        if (ie_ref(e).ne.0) ntag=ntag+1
      enddo
      ntag=iglsum(ntag,1)
      ntag0=ntag
      if (nio.eq.0) write(*,26) 0,ntag,ntag

      do iter=1,max_iter

        call rzero(tmpo,n)
        call rzero(tmpa,n)
        do e=1,nelv
        do f=1,2*ldim
          if (face_tag(f,e).eq.1) call facev(tmpo,e,f,1.0,lx1,ly1,lz1)
          if (face_tag(f,e).eq.2) call facev(tmpa,e,f,1.0,lx1,ly1,lz1)
        enddo
        enddo

        call dssum(tmpo,lx1,ly1,lz1)
        call dssum(tmpa,lx1,ly1,lz1)

        do e=1,nelv
          do f=1,2*ldim
            ix = ix_fctr(f)
            iy = iy_fctr(f)
            iz = iz_fctr(f)
            fo = fopp(f)
           
            ! propagate via opposite face 
            if (tmpo(ix,iy,iz,e).gt.0.0) then
              if (ie_ref(e).eq.0) ie_ref(e) = f
              if (ie_ref(e).eq.f.OR.ie_ref(e).eq.fo) then
                face_tag(f,e) = 1 
                face_tag(fo,e) = 1 
                do i=1,4
                  fa = fadj(i,f)
                  face_tag(fa,e) = 2
                enddo
              else
                ierr = ierr + 1
              endif
            endif

            ! propagate via adjacent face 
            if (tmpa(ix,iy,iz,e).gt.0.0) then
              if (ie_ref(e).ne.f.AND.ie_ref(e).ne.fo) then
                face_tag(f,e) = 2
                face_tag(fo,e) = 2
              else
                ierr = ierr + 1
              endif
            endif
          enddo
        enddo
        ierr = iglsum(ierr,1)
        if (ierr.gt.0) goto 90

        call update_face_tag(face_tag,ie_ref,istatus)

        ! check conv
        ntag1 = ntag
        ntag = 0
        ndif = 0
        do e=1,nelv
          if (ie_ref(e).gt.0) ntag = ntag + 1
          do f=1,2,ldim
            if (face_tag(f,e).ne.face_tag0(f,e)) ndif=ndif+1
          enddo
        enddo
        ntag = iglsum(ntag,1)
        ndif = iglsum(ndif,1)
        if (nio.eq.0) write(*,26) iter,ntag,ndif
        if (ndif.eq.0.AND.ntag.eq.ntag1) goto 90
        if (ntag.eq.nelgv) goto 90 ! all e is assigned
        call icopy(face_tag0,face_tag0,2*ldim*lelt)

      enddo

   90 continue ! err or conv
      if (ierr.gt.0) then
        if (nio.eq.0) 
     $    write(*,*) 'ERROR: cross path detected in trace_layers'
     $             , iter,ierr
        call exitt
      endif

      call update_face_tag(face_tag,ie_ref,istatus)
      if (istatus.eq.1) then
        call exitti('ERROR: type 2,0 is not supported$',1)
      endif

      if (nio.eq.0) then
        write(*,27)'rfn trace lyr Niter',iter
     $            ,', tag0=',ntag0,' tag1=',ntag
      endif

   26 format('rfn trace lyr iter=',3i9)
   27 format(3(ai9))
      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_bak_nelg(ie_ref,imode,Nref)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'REFINE'
      integer imode,ie_ref(lelt),e,Nref,iglsum
      
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

      if (imode.eq.1) then
        if (nio.eq.0) write(*,22)'rfn REFINE START,    E='
     $                          ,nelgv0,', E_lyr',nelg_lv,'Nref',Nref
      elseif (imode.eq.2) then
        if (nio.eq.0) write(*,21)'rfn EXTRUDE START,   E='
     $                          ,nelgv0,', E_lyr=',nelg_lv
      elseif (imode.eq.3) then 
        if (nio.eq.0) write(*,21)'rfn SSPLIT START,   E='
     $                          ,nelgv0,', E_lyr=',nelg_lv,'Ncut',Nref
      endif

   21 format(2(ai9))
   22 format(2(ai9),ai3)

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_chk_Rref(Nref,Rref,Rlim)
      implicit none
      include 'SIZE'
      include 'REFINE'

      integer Nref,ierr,ilv
      real Rref(lref-1),Rlim(4,lref-1)
      real tol,rmn,rmx

      if (Nref.gt.1) return

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
     $                          ,(Rref(ilv),ilv=1,min(Nref-1,10))
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

   25 format(a,10g13.5)
      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_chk_Ncut(Ncut)
      implicit none
      include 'SIZE'
      include 'PARALLEL' ! nelgv
      include 'REFINE'
      integer Ncut,nelv_est,nelgv_est

      if (Ncut.lt.0) call exitti('Err rfn: invalid Ncut$',Ncut)
      if (Ncut.gt.lcut) call exitti('Err rfn: Ncut>lcut$',Ncut)
      
      nelv_est = nelv + nel_lv*( (Ncut+1)*(Ncut+1) - 1 )
      nelgv_est = nelgv + nelg_lv*( (Ncut+1)*(Ncut+1) - 1 )
      if (nelgv_est.gt.lelg.OR.nelv_est.gt.lelv)
     $  call exitti('Err rfn: need larger lelg/lelv$',Ncut)

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
