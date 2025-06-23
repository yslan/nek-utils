c-----------------------------------------------------------------------
      subroutine mirror_mesh(dir)
      include 'SIZE'
      include 'GEOM'

      character dir
 
c     Mirror coordinates
      if(dir.eq.'x'.or.dir.eq.'X')
     &  call cmult(xm1,-1.0,lx1*ly1*lz1*nelt)
      if(dir.eq.'y'.or.dir.eq.'Y')
     &  call cmult(ym1,-1.0,lx1*ly1*lz1*nelt)
      if(dir.eq.'z'.or.dir.eq.'Z')
     &  call cmult(zm1,-1.0,lx1*ly1*lz1*nelt)

c     fix left-handedness
      call flip_elements(.true.)

c     connectivity needs to be regenerated
      call gen_rea(2)
      call exitt !
       
      return
      end 
c-----------------------------------------------------------------------
c     Yu-Hsiang's code below
c-----------------------------------------------------------------------
      subroutine flip_elements(if_con_in_bc)
c     This subroutine flip the elements from lhs to rhs
c     Note that this is not the general fix for rhs_check.
c     We simply mirrors the elements here
c     call this subroutine from usrdat2
c     
c     if_con_in_bc = T, will also fix the connectivity stored in bc(1) and bc(2)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer elist(lelt),nei_flist(6,lelt),ierr
      integer ifld,e,f,f_old,iz,iz_old,iv,iv_old,nxyz
      logical if_per_in_bc

      integer nei_flist_bak(2*ldim)
      real xbak(lx1,ly1,lz1), ybak(lx1,ly1,lz1), zbak(lx1,ly1,lz1)
      real bc_bak(5,2*ldim)
      character*3 cbc_bak(2*ldim)

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer*8 vertex,vertex_bak(8)

      integer flist(6),vlist(8)
      data flist /1,2,3,4,6,5/     ! map from new face to old face
      data vlist /5,6,7,8,1,2,3,4/ ! map from new vtx to old vtx

      nxyz=lx1*ly1*lz1

      if (ldim.ne.3) call exitti('this only works for 3d$',1)
      call my_mshchk('before')

      ! chk rhs for each elements
      call izero(elist,lelt)
      call my_chkrhs(ierr,elist) 

      ! see if neighbor is flip
      call get_neighbor_elist(nei_flist,elist)

      ! flip elements
      do e=1,nelt

        if (elist(e).gt.0) then ! this element fails rhs check

          ! flip coordinates through "t" axis
          call copy(xbak,xm1(1,1,1,e),nxyz)
          call copy(ybak,ym1(1,1,1,e),nxyz)
          call copy(zbak,zm1(1,1,1,e),nxyz)
          do iz=1,lz1
            iz_old = lz1-iz+1
            call copy(xm1(1,1,iz,e),xbak(1,1,iz_old),lx1*ly1)
            call copy(ym1(1,1,iz,e),ybak(1,1,iz_old),lx1*ly1)
            call copy(zm1(1,1,iz,e),zbak(1,1,iz_old),lx1*ly1)
          enddo

          ! swap face 5-6
          do ifld=1,nfield
            do f=1,2*ldim
              cbc_bak(f) = cbc(f,e,ifld)
              call copy(bc_bak(1,f),bc(1,f,e,ifld),5)
              nei_flist_bak(f) = nei_flist(f,e)
            enddo

            do f=1,2*ldim
              f_old = flist(f)
              cbc(f,e,ifld) = cbc_bak(f_old)
              call copy(bc(1,f,e,ifld),bc_bak(1,f_old),5)
              nei_flist(f,e) = nei_flist_bak(f_old)
            enddo
          enddo

          ! swap vtx ! TODO: This is not tested
          do iv=1,2**ldim
            vertex_bak(iv) = vertex(iv,e)
          enddo
          do iv=1,2**ldim
            iv_old = vlist(iv)
            vertex(iv,e) = vertex_bak(iv_old)
          enddo

        endif

      enddo

      if (if_con_in_bc) then      
        do e=1,nelt
        do f=1,2*ldim
c          if (CBC(f,e,1).eq.'P  '.OR.CBC(f,e,1).eq.'p  ') then ! only need this info for PER
          if (nei_flist(f,e).eq.1) then ! your neighbor is fliped.
            do ifld=1,nfield
              f_old = bc(2,f,e,ifld)
              bc(2,f,e,ifld) = flist(f_old) ! permute: fnew = flist(fold) and fold = flist(fnew)
            enddo
          endif
c          endif
        enddo
        enddo
      endif

      call my_mshchk('after ')

      return
      end
c-----------------------------------------------------------------------
      subroutine get_neighbor_elist(nei_flist,elist)
c     Given elist that assigning each element an integer, this subroutine 
c     copy neighbor elements' values. So at element e and face f, the 
c     neighbor will have value nei_flist(f,e)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer nei_flist(6,lelt),elist(lelt),n,nxyz
      real tmp(lx1,ly1,lz1,lelt)

      ! face center
      integer ix_fctr(6),iy_fctr(6),iz_fctr(6),idd,ix,iy,iz,e,f
      parameter(idd=ldim-1)
      data ix_fctr /  2,lx1,  2,  1,  2,  2/
      data iy_fctr /  1,  2,ly1,  2,  2,  2/
      data iz_fctr /idd,idd,idd,idd,  1,lz1/

      if (lx1.lt.3) call exitti('This need lx1>=3$',1) ! because we are lazy

      n = lx1*ly1*lz1*nelt
      nxyz = lx1*ly1*lz1

      call izero(nei_flist,6*lelt)
      call rzero(tmp,n)

      ! face gather-scatter
      do e=1,nel
        call cfill(tmp(1,e), real(elist(e)), nxyz)
      enddo
      call dssum(tmp,lx1,ly1,lz1)

      do e=1,nel
      do f=1,2*ldim
         ix = ix_fctr(f)
         iy = iy_fctr(f)
         iz = iz_fctr(f)
         nei_flist(f,e) = int(tmp(ix,iy,iz,e)) - elist(e)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine my_mshchk(s6)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      character*6 s6
      integer ierr,ierr1,ierr2,ierr3
      integer e,nel,nxyz,iwk(lelt)
      integer inei,nnei,nplt0,nplt1,iglsum
      real tmp(lx1,ly1,lz1,lelt),dummy(lx1,ly1,lz1,lelt), vlmax
      logical ifdump

      integer err_elist
      common /my_errchki/ err_elist(lelt)

      nel=nelv
      nxyz=lx1*ly1*lz1
      nnei=3

      ierr=0
      call izero(err_elist,lelt)

      ! rhs
      call izero(iwk,lelt)
      call my_chkrhs(ierr2,iwk)
      if(ierr2.gt.0) then
        ierr=max(ierr,ierr2)
        call my_iadd2(err_elist,iwk,nel)
      endif

      ! jac
      call izero(iwk,lelt)
      call my_chkjac(ierr3,iwk)
      if(ierr3.gt.0) then
        ierr=max(ierr,ierr3)
        call my_iadd2(err_elist,iwk,nel)
      endif

      ! summary
      if(ierr.eq.0) then
        if(nid.eq.0)write(*,*)'PASSED mshchk ',s6
      else
        if(nid.eq.0)write(*,*)'FAILED mshchk ',s6
     $                       ,'  rhs fail',ierr2,'  neg jac ',ierr3
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine my_chkjac(ierr,elist)
c     Check negative Jacobian
c     Note, Nek only requires Jac in the same sign of V(1)
c           At here, we want Jac>0 for all vertices
      implicit none
      include 'SIZE'
      integer e,nel,nxyz,ierr,iglsum,elist(lelt)
      real Jac(lx1*ly1*lz1,lelt),vlmin

      nel=nelv
      nxyz=lx1*ly1*lz1

      ierr=0
      call izero(elist,nel)

      call chk_comp_Jac0(Jac)

      do e=1,nel
        if(vlmin(Jac(1,e),nxyz).lt.0.0) then
          elist(e)=1
          ierr=ierr+1
        endif
      enddo
      ierr=iglsum(ierr,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_comp_Jac0(Jac)
      implicit none
      include 'SIZE'
      real Jac (lx1*ly1*lz1*lelt)

      if (ldim.eq.3) then
        call chk_comp_Jac0_3D(Jac)
      else
        call chk_comp_Jac0_2D(Jac)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_comp_Jac0_2D(Jac)
      implicit none
      include 'SIZE'
      real Jac (lx1*ly1*lz1*lelt), J(lx1*ly1*lz1*lelt,ldim*ldim)
      real dmy (lx1*ly1*lz1*lelt,5)
      integer i,nel,ntot

      nel=nelv
      ntot=lx1*ly1*lz1*nel

      call xyzrst(J(1,1),J(1,2),dmy(1,1)
     $           ,J(1,3),J(1,4),dmy(1,2)
     $           ,dmy(1,3),dmy(1,4),dmy(1,5),.false.)               ! mesh 1 ! FIXME 2D

      do i=1,ntot
        Jac(i) = J(i,1)*J(i,4) - J(i,3)*J(i,2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_comp_Jac0_3D(Jac)
      implicit none
      include 'SIZE'
      integer nel,ntot,i
      real Jac(lx1*ly1*lz1*lelt),J(lx1*ly1*lz1*lelt,ldim*ldim)

      nel=nelv
      ntot=lx1*ly1*lz1*nel

      call xyzrst(J(1,1),J(1,2),J(1,3),J(1,4),J(1,5),J(1,6)
     $           ,J(1,7),J(1,8),J(1,9),.false.)               ! mesh 1

      do i=1,ntot
        Jac(i) = J(i,1)*J(i,5)*J(i,9)
     $         + J(i,7)*J(i,2)*J(i,6)
     $         + J(i,4)*J(i,8)*J(i,3)
     $         - J(i,1)*J(i,8)*J(i,6)
     $         - J(i,4)*J(i,2)*J(i,9)
     $         - J(i,7)*J(i,5)*J(i,3)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine my_chkrhs(ierr,elist)
      implicit none
      include 'SIZE'
      integer ierr,elist(lelt)

      if (ldim.eq.3) then
        call my_chkrhs3D(ierr,elist)
      else
        call my_chkrhs2D(ierr,elist)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine xm1toxc  ! xc: vertices
      implicit none
      include 'SIZE' ! nelv
      include 'GEOM' ! xm1
      include 'INPUT'! xc
      integer ie,nel

      nel=nelv
      do ie=1,nel
        xc(1,ie) = XM1(1  ,1  ,1  ,ie)
        xc(2,ie) = XM1(lx1,1  ,1  ,ie)
        xc(3,ie) = XM1(lx1,ly1,1  ,ie)
        xc(4,ie) = XM1(1  ,ly1,1  ,ie)
        xc(5,ie) = XM1(1  ,1  ,lz1,ie)
        xc(6,ie) = XM1(lx1,1  ,lz1,ie)
        xc(7,ie) = XM1(lx1,ly1,lz1,ie)
        xc(8,ie) = XM1(1  ,ly1,lz1,ie)
        ! y
        yc(1,ie) = YM1(1  ,1  ,1  ,ie)
        yc(2,ie) = YM1(lx1,1  ,1  ,ie)
        yc(3,ie) = YM1(lx1,ly1,1  ,ie)
        yc(4,ie) = YM1(1  ,ly1,1  ,ie)
        yc(5,ie) = YM1(1  ,1  ,lz1,ie)
        yc(6,ie) = YM1(lx1,1  ,lz1,ie)
        yc(7,ie) = YM1(lx1,ly1,lz1,ie)
        yc(8,ie) = YM1(1  ,ly1,lz1,ie)
        ! z
        if(ldim.eq.3) then
          zc(1,ie) = ZM1(1  ,1  ,1  ,ie)
          zc(2,ie) = ZM1(lx1,1  ,1  ,ie)
          zc(3,ie) = ZM1(lx1,ly1,1  ,ie)
          zc(4,ie) = ZM1(1  ,ly1,1  ,ie)
          zc(5,ie) = ZM1(1  ,1  ,lz1,ie)
          zc(6,ie) = ZM1(lx1,1  ,lz1,ie)
          zc(7,ie) = ZM1(lx1,ly1,lz1,ie)
          zc(8,ie) = ZM1(1  ,ly1,lz1,ie)
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine my_chkrhs2D(ierr,elist)
      implicit none
      include 'SIZE' ! nelv
      include 'GEOM' ! xm1
      include 'INPUT'! xc
      include 'SCRCT'! xyz

      integer ierr,ie,nel,j,ivtx,iglsum,elist(lelt)
      real CRSS2D,c1,c2,c3,c4
      integer indx(4) ! local usage, todo: include TOPOL?
      data indx /1,2,4,3/

      nel=nelv

      ierr=0
      call izero(elist,nel)
      call xm1toxc
      call rzero(xyz,24*nel)

      do ie=1,nel

        do j=1,4
          ivtx = indx(j)
          xyz(1,ivtx,ie) = xc(j,ie)
          xyz(2,ivtx,ie) = yc(j,ie)
        enddo

        C1=CRSS2D(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,1,IE))
        C2=CRSS2D(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,2,IE))
        C3=CRSS2D(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,3,IE))
        C4=CRSS2D(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,4,IE))

        if (C1.LE.0.0.OR.C2.LE.0.0.OR.
     $      C3.LE.0.0.OR.C4.LE.0.0 ) then
          elist(ie)=1
          ierr=ierr+1
        endif

      enddo

      ierr=iglsum(ierr,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine my_chkrhs3D(ierr,elist)
c     Check right-handness
      implicit none
      include 'SIZE' ! nelv
      include 'GEOM' ! xm1
      include 'INPUT'! xc
      include 'SCRCT'! xyz

      integer ierr,ie,nel,j,ivtx,iglsum,elist(lelt)
      real VOLUM0,v1,v2,v3,v4,v5,v6,v7,v8
      integer indx(8) ! local usage, todo: include TOPOL?
      data indx /1,2,4,3,5,6,8,7/

      nel=nelv

      ierr=0
      call izero(elist,nel)
      call xm1toxc
      call rzero(xyz,24*nel)

      do ie=1,nel

        do j=1,8
          ivtx = indx(j)
          xyz(1,ivtx,ie) = xc(j,ie)
          xyz(2,ivtx,ie) = yc(j,ie)
          xyz(3,ivtx,ie) = zc(j,ie)
        enddo

        V1= VOLUM0(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,5,IE),XYZ(1,1,IE))
        V2= VOLUM0(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,6,IE),XYZ(1,2,IE))
        V3= VOLUM0(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,7,IE),XYZ(1,3,IE))
        V4= VOLUM0(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,8,IE),XYZ(1,4,IE))
        V5=-VOLUM0(XYZ(1,6,IE),XYZ(1,7,IE),XYZ(1,1,IE),XYZ(1,5,IE))
        V6=-VOLUM0(XYZ(1,8,IE),XYZ(1,5,IE),XYZ(1,2,IE),XYZ(1,6,IE))
        V7=-VOLUM0(XYZ(1,5,IE),XYZ(1,8,IE),XYZ(1,3,IE),XYZ(1,7,IE))
        V8=-VOLUM0(XYZ(1,7,IE),XYZ(1,6,IE),XYZ(1,4,IE),XYZ(1,8,IE))

        if (V1.LE.0.0.OR.V2.LE.0.0.OR.
     $      V3.LE.0.0.OR.V4.LE.0.0.OR.
     $      V5.LE.0.0.OR.V6.LE.0.0.OR.
     $      V7.LE.0.0.OR.V8.LE.0.0    ) then
          elist(ie)=1
          ierr=ierr+1
        endif

      enddo

      ierr=iglsum(ierr,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine my_iadd2(a,b,n)
      integer a(1),b(1)

      do i=1,n
        a(i)=a(i)+b(i)
      enddo

      return
      end
