c     dependencies
c        my_outpost.f
c           my_outpost
c-----------------------------------------------------------------------
      subroutine my_domain_chk
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real x,y,rad
      real glmin,glmax,rmn,rmx,xmn,xmx,ymn,ymx,zmn,zmx
      integer i,n

      n=lx1*ly1*lz1*nelv

      rmn=1.e22
      rmx=0.0

      ! cylinder
      do i=1,n

        x=xm1(i,1,1,1)
        y=ym1(i,1,1,1)

        rad=sqrt(x*x+y*y)
        rmx=max(rmx,rad)
        rmn=min(rmn,rad)

      enddo
      rmx=glmax(rmx,1)
      rmn=glmin(rmn,1)

      ! plane
      call domain_size(xmn,xmx,ymn,ymx,zmn,zmx)

      if (nio.eq.0) then
        write(6,18) rmn,xmn,ymn,zmn
        write(6,19) rmx,xmx,ymx,zmx
   18   format('new rxyz min  ',1p4e16.8)
   19   format('new rxyz max  ',1p4e16.8)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine my_mesh_metric

      call geom_reset(2)
      call xm1toxc
      call mesh_metrics ! read jacm, xc

      return
      end
c-----------------------------------------------------------------------
      subroutine my_mshchk(s6,ifdump,ierr)
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

      ! find neighbors of bad elements
      if(ifdump.AND.ierr.gt.0) then 

        do e=1,nel          
          if(err_elist(e).gt.0) err_elist(e)=1
        enddo
        nplt0=iglsum(err_elist,nel)

        call rzero(dummy,nel*nxyz) 
        do e=1,nel
          if(err_elist(e).gt.0) then
            call cfill(dummy(1,1,1,e),lglel(e)*1.0,nxyz)
          endif
        enddo

        do inei=1,nnei
          call rzero(tmp,nxyz*nel)
          do e=1,nel
            if(err_elist(e).gt.0) call rone(tmp(1,1,1,e),nxyz)
          enddo

          call dssum(tmp,lx1,ly1,lz1) ! use dssum to extend neighbor
          do e=1,nel          
            if(vlmax(tmp(1,1,1,e),nxyz).gt.0) err_elist(e)=1
          enddo

          do e=1,nel          
            if(err_elist(e).gt.0) err_elist(e)=1
          enddo
          nplt1=iglsum(err_elist,nel)
          if(nid.eq.0)write(*,*)'err exi',inei,nplt1
        enddo
        if(nid.eq.0)write(*,*)'err dmp',nnei,nplt0,nplt1

        ifxyo=.true.
        call my_outpost(err_elist,dummy,vy,vz,pr,t,'err')

      endif

   50 format(a,i3,i6,i6,i12,1p3e16.8)
      return
      end
c-----------------------------------------------------------------------
      subroutine my_mshchk_vb(s6,ifdump,ierr)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      character*6 s6
      integer ierr,ierr1,ierr2,ierr3
      integer e,nel,nxyz,iwk(lelt),ic
      integer inei,nnei,nplt0,nplt1,iglsum
      real tmp(lx1,ly1,lz1,lelt), dummy(lx1,ly1,lz1,lelt), vlmax
      logical ifdump

      integer err_elist
      common /my_errchki/ err_elist(lelt)

      nel=nelv
      nxyz=lx1*ly1*lz1

      nnei=3 ! # neighbor for dumping files
      ierr=0
      call izero(err_elist,lelt)

      ! rhs
      call izero(iwk,lelt)
      call my_chkrhs(ierr2,iwk)
      if(ierr2.gt.0) then
        ierr=max(ierr,ierr2)
        call my_iadd2(err_elist,iwk,nel)
        do e=1,nel
          if(iwk(e).gt.0) then
            write(*,51)'Eerr rhs',istep,nid,lglel(e),int(bc(5,5,e,1))
     $                       ,xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e)
     $                       ,' ',(cbc(ic,e,1),ic=1,6)
          endif
        enddo
      endif

      ! jac
      call izero(iwk,lelt)
      call my_chkjac(ierr3,iwk)
      if(ierr3.gt.0) then
        ierr=max(ierr,ierr3)
        call my_iadd2(err_elist,iwk,nel)
        do e=1,nel
          if(iwk(e).gt.0) then
            write(*,51)'Eerr jac',istep,nid,lglel(e),int(bc(5,5,e,1))
     $                       ,xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e)
     $                       ,' ',(cbc(ic,e,1),ic=1,6)
          endif
        enddo
      endif

      ! summary
      if(ierr.eq.0) then
        if(nid.eq.0)write(*,*)'PASSED mshchk ',s6
      else
        if(nid.eq.0)write(*,*)'FAILED mshchk ',s6
     $                       ,'  rhs fail',ierr2,'  neg jac ',ierr3
      endif

      ! find neighbors of bad elements
      if(ifdump.AND.ierr.gt.0) then 

        do e=1,nel
          if(err_elist(e).gt.0) err_elist(e)=1
        enddo
        nplt0=iglsum(err_elist,nel)

        call rzero(dummy,nel*nxyz) 
        do e=1,nel
          if(err_elist(e).gt.0) then
            call cfill(dummy(1,1,1,e),lglel(e)*1.0,nxyz)
          endif
        enddo

        do inei=1,nnei
          call rzero(tmp,nxyz*nel)
          do e=1,nel
            if(err_elist(e).gt.0) call rone(tmp(1,1,1,e),nxyz)
          enddo

          call dssum(tmp,lx1,ly1,lz1) 
          do e=1,nel
            if(vlmax(tmp(1,1,1,e),nxyz).gt.0) err_elist(e)=1
          enddo

          do e=1,nel
            if(err_elist(e).gt.0) err_elist(e)=1
          enddo
          nplt1=iglsum(err_elist,nel)
          if(nid.eq.0)write(*,*)'err exi',inei,nplt1
        enddo
        if(nid.eq.0)write(*,*)'err dmp',nnei,nplt0,nplt1

        ifxyo=.true.
        call my_outpost(err_elist,vx,vy,dummy,pr,t,'err')

      endif

   50 format(a,i3,i6,i6,i12,1p3e16.8)
   51 format(a,i3,i6,i12,i3,1p3e16.8,a1,6a3)

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
      subroutine kirk_jac_hex8(h8,det)
      ! evaluate Jacobian of hex8 element
      include 'SIZE'
      real h8(3,8),det(8)
      real ma(3,3)
      logical ifnegjac_subhex

      if(ldim.eq.2)return
      ! pack ma for node1
      call vec_minus(ma(1,1),h8(1,2),h8(1,1))
      call vec_minus(ma(1,2),h8(1,4),h8(1,1))
      call vec_minus(ma(1,3),h8(1,5),h8(1,1))
      call determinant_3by3(ma,det(1))

      ! pack ma for node2
      call vec_minus(ma(1,1),h8(1,2),h8(1,1))
      call vec_minus(ma(1,2),h8(1,3),h8(1,2))
      call vec_minus(ma(1,3),h8(1,6),h8(1,2))
      call determinant_3by3(ma,det(2))

      ! pack ma for node3
      call vec_minus(ma(1,1),h8(1,3),h8(1,4))
      call vec_minus(ma(1,2),h8(1,3),h8(1,2))
      call vec_minus(ma(1,3),h8(1,7),h8(1,3))
      call determinant_3by3(ma,det(3))

      ! pack ma for node4
      call vec_minus(ma(1,1),h8(1,3),h8(1,4))
      call vec_minus(ma(1,2),h8(1,4),h8(1,1))
      call vec_minus(ma(1,3),h8(1,8),h8(1,4))
      call determinant_3by3(ma,det(4))

      ! pack ma for node5
      call vec_minus(ma(1,1),h8(1,6),h8(1,5))
      call vec_minus(ma(1,2),h8(1,8),h8(1,5))
      call vec_minus(ma(1,3),h8(1,5),h8(1,1))
      call determinant_3by3(ma,det(5))

      ! pack ma for node6
      call vec_minus(ma(1,1),h8(1,6),h8(1,5))
      call vec_minus(ma(1,2),h8(1,7),h8(1,6))
      call vec_minus(ma(1,3),h8(1,6),h8(1,2))
      call determinant_3by3(ma,det(6))
              
      ! pack ma for node7
      call vec_minus(ma(1,1),h8(1,7),h8(1,8))
      call vec_minus(ma(1,2),h8(1,7),h8(1,6))
      call vec_minus(ma(1,3),h8(1,7),h8(1,3))
      call determinant_3by3(ma,det(7))
 
      ! pack ma for node8
      call vec_minus(ma(1,1),h8(1,7),h8(1,8))
      call vec_minus(ma(1,2),h8(1,8),h8(1,5))
      call vec_minus(ma(1,3),h8(1,8),h8(1,4))
      call determinant_3by3(ma,det(8))
                 
      return
      end
!--------------------------------------------   
      subroutine vec_minus(c,a,b)
      real a(3),b(3),c(3)
      c(1)=(a(1)-b(1))*0.5 
      c(2)=(a(2)-b(2))*0.5
      c(3)=(a(3)-b(3))*0.5
      return
      end
!--------------------------------------------   
      subroutine determinant_3by3(ma,jac)
      real ma(9),jac
        jac =    ma(1)*ma(5)*ma(9)
     $         + ma(7)*ma(2)*ma(6)
     $         + ma(4)*ma(8)*ma(3)
     $         - ma(1)*ma(8)*ma(6)
     $         - ma(4)*ma(2)*ma(9)
     $         - ma(7)*ma(5)*ma(3)
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_kirk_comp_Jac_hex8 ! a litter help for testing Kirk's implementation
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer e,n
      real h8(3,8),jtmp(8),err,glamax
      real Jac1(lx1,ly1,lz1,lelt)
      real Jac2(lx1*ly1*lz1*lelt)
      real Jdif(lx1*ly1*lz1*lelt)

      ! use exodus ordering

      do e=1,nelv
        h8(1,1) = xm1(1  ,1  ,1  ,e)
        h8(1,2) = xm1(lx1,1  ,1  ,e)
        h8(1,3) = xm1(lx1,lx1,1  ,e) 
        h8(1,4) = xm1(1  ,lx1,1  ,e)
        h8(1,5) = xm1(1  ,1  ,lx1,e)
        h8(1,6) = xm1(lx1,1  ,lx1,e)
        h8(1,7) = xm1(lx1,lx1,lx1,e)
        h8(1,8) = xm1(1  ,lx1,lx1,e)
        h8(2,1) = ym1(1  ,1  ,1  ,e)
        h8(2,2) = ym1(lx1,1  ,1  ,e)
        h8(2,3) = ym1(lx1,lx1,1  ,e) 
        h8(2,4) = ym1(1  ,lx1,1  ,e)
        h8(2,5) = ym1(1  ,1  ,lx1,e)
        h8(2,6) = ym1(lx1,1  ,lx1,e)
        h8(2,7) = ym1(lx1,lx1,lx1,e)
        h8(2,8) = ym1(1  ,lx1,lx1,e)
        h8(3,1) = zm1(1  ,1  ,1  ,e)
        h8(3,2) = zm1(lx1,1  ,1  ,e)
        h8(3,3) = zm1(lx1,lx1,1  ,e) 
        h8(3,4) = zm1(1  ,lx1,1  ,e)
        h8(3,5) = zm1(1  ,1  ,lx1,e)
        h8(3,6) = zm1(lx1,1  ,lx1,e)
        h8(3,7) = zm1(lx1,lx1,lx1,e)
        h8(3,8) = zm1(1  ,lx1,lx1,e)
        call kirk_jac_hex8(h8,jtmp)
        jac1(1  ,1  ,1  ,e) = jtmp(1)
        jac1(lx1,1  ,1  ,e) = jtmp(2)
        jac1(lx1,lx1,1  ,e) = jtmp(3)
        jac1(1  ,lx1,1  ,e) = jtmp(4)
        jac1(1  ,1  ,lx1,e) = jtmp(5)
        jac1(lx1,1  ,lx1,e) = jtmp(6)
        jac1(lx1,lx1,lx1,e) = jtmp(7)
        jac1(1  ,lx1,lx1,e) = jtmp(8)
      enddo

      call chk_comp_Jac0(jac2)

      n=lx1*ly1*lz1*nelv
      call sub3(jdif,jac1,jac2,n)
      err=glamax(jdif,n)
      if(nid.eq.0)write(*,*)'err',err

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
c     This is not called, we leave here for reference/comparison/debugging
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
      subroutine xm1toxyz 
      implicit none
      include 'SIZE' ! nelv
      include 'GEOM' ! xm1
      include 'INPUT'! xc
      include 'SCRCT'! xyz

      integer ierr,ie,nel,j,ivtx
      integer indx(8) ! local usage, todo: include TOPOL?
      data indx /1,2,4,3,5,6,8,7/

      nel=nelv

      call xm1toxc
      call rzero(xyz,24*nel)

      if(ldim.eq.2) then
        do ie=1,nel
          do j=1,4
            ivtx = indx(j)
            xyz(1,ivtx,ie) = xc(j,ie)
            xyz(2,ivtx,ie) = yc(j,ie)
          enddo
        enddo
      else
        do ie=1,nel
          do j=1,8
            ivtx = indx(j)
            xyz(1,ivtx,ie) = xc(j,ie)
            xyz(2,ivtx,ie) = yc(j,ie)
            xyz(3,ivtx,ie) = zc(j,ie)
          enddo
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine my_iadd2(i1,i2,n)
      DIMENSION I1(1),I2(1)
C
      DO 10 I=1,N
         I1(I)=I1(I)+I2(I)
   10 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine my_prt_lglel
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SCRCT'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      logical ifverbm
C     Output the processor-element map:
      ifverbm=.false.
      if (loglevel .gt. 2) ifverbm=.true.

      loglevel_bak=loglevel
      loglevel=3

      ifverbm=.true.

      if(ifverbm) then
        idum = 1
        if(nid.eq.0) then
           N8 = min(8,nelt)
           write(6 ,1310) node-1,(lglel(ie),ie=1,n8)
           if (NELT.GT.8) write(6 ,1315) (lglel(ie),ie=9,NELT)
           if (NELT.GT.8) write(6 ,*)' '
           DO inid=1,NP-1
              mtype = inid
              call csend(mtype,idum,4,inid,0)         ! handshake
              call crecv(mtype,inelt,4)               ! nelt of other cpus
              N8 = min(8,inelt)
           ENDDO
 1310      FORMAT(' RANK',I6,' IEG',8I12)
 1315      FORMAT('     ',6X,'    ',8I12)
        else
           mtype = nid
           call crecv(mtype,idum,4)                ! hand-shake
           call csend(mtype,nelt,4,0,0)            ! nelt
           if (loglevel .gt. 2) then
              N8 = min(8,nelt)
              write(6 ,1310) node-1,(lglel(ie),ie=1,n8)
              if (NELT.GT.8) write(6 ,1315) (lglel(ie),ie=9,NELT)
              if (NELT.GT.8) write(6 ,*)' '
           endif
        endif
      endif

      loglevel=loglevel_bak

      return
      end
c----------------------------------------------------------------------
      subroutine my_dump_mesh(s3,imode)
      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'INPUT'

      character*3 s3
      integer imid,ierr

      if (imode.gt.0) then
        call fix_geom ! reset jac for metric
        call my_mshchk('bfrdmp',.true.,ierr) ! last chk before dump

        call domain_size(xmn,xmx,ymn,ymx,zmn,zmx)
        if (nio.eq.0) then
          write(6,18) xmn,ymn,zmn
          write(6,19) xmx,ymx,zmx
   18     format('new xyz min  ',5g13.5)
   19     format('new xyz max  ',5g13.5)
        endif
        call my_mesh_metric

        ifxyo = .true.
        call outpost(vx,vy,vz,pr,t,s3)
      else
        call my_domain_chk
      endif

      if (abs(imode).eq.1) then
        imid=1
c       call gen_rea(imid)
c       if (nid.eq.0) write(6,*)'done genrea'
c        call my_gen_rea(s3,imid)
c        if (nid.eq.0) write(6,*)'done genrea:','newrea'//s3//'.out'

        call my_gen_re2(s3,imid)
c       call gen_re2(imid) ! FIXME, bug for large E (in bc)
        if (nid.eq.0) write(6,*)'done genre2:','newre2'//s3//'.re2'

      elseif (abs(imode).eq.2) then
        imid=0
        call my_gen_re2_xyz
        if (nid.eq.0) write(6,*)'done genre2 (xyz only)'

c       call gen_re2(imid) ! FIXME, bug for large E (in bc) ! TODO for dbg
c       if (nid.eq.0) write(6,*)'done genre2'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine my_gen_re2(s3,imid)  ! Generate and output essential parts of .rea
                                ! And re2
                                ! Clobbers ccurve()
                                ! byte read is float size..
                                ! 4 wdsize
      include 'SIZE'
      include 'TOTAL'
      character*3 s3

      character*80 hdr
      real*4 test
      data   test  / 6.54321 /
      integer ierr


c     imid = 0  ! No midside node defs
c     imid = 1  ! Midside defs where current curve sides don't exist
c     imid = 2  ! All nontrivial midside node defs

      ierr = 0
      if (nid.eq.0) then
         call byte_open('newre2'//s3//'.re2'//char(0),ierr)
c        call byte_open('newre2.re2' // char(0), ierr)
         call blank(hdr,80)
         if(wdsize.eq.8) then
            write(hdr,112) nelgt,ldim,nelgv
         else
            write(hdr,111) nelgt,ldim,nelgv
         endif
  111    format('#v001',i9,i3,i9,' hdr')
  112    format('#v002',i9,i3,i9,' hdr')
         if(ierr.eq.0) call byte_write(hdr,20,ierr)
         if(ierr.eq.0) call byte_write(test,1,ierr) !write endian discriminator
      endif
      call err_chk(ierr,'Error opening  in gen_re2$')

      call gen_re2_xyz
      call gen_re2_curve(imid)  ! Clobbers ccurve()

      do ifld=1,nfield
         call gen_re2_bc   (ifld)
      enddo

      if (nid.eq.0) call byte_close(ierr)
      call err_chk(ierr,'Error closing in gen_re2$')

      return
      end
c-----------------------------------------------------------------------
      subroutine my_gen_rea(s3,imid)  
                                ! Clobbers ccurve()
      include 'SIZE'
      include 'TOTAL'
      character*3 s3

c     imid = 0  ! No midside node defs
c     imid = 1  ! Midside defs where current curve sides don't exist
c     imid = 2  ! All nontrivial midside node defs

      if (nid.eq.0) 
     $   open(unit=10,file='newrea'//s3//'.out',status='unknown')! clobbers existing file

      call gen_rea_xyz

      call gen_rea_curve(imid)  ! Clobbers ccurve()

      if (nid.eq.0) write(10,*)' ***** BOUNDARY CONDITIONS *****'
      do ifld=1,nfield
         call gen_rea_bc   (ifld)
      enddo

      if (nid.eq.0) close(10)

      return
      end
c-----------------------------------------------------------------------
      subroutine my_gen_re2_xyz ! <case>_xyz.re2
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer imid
      character*80 hdr
      real*4 test
      data   test  / 6.54321 /
      integer ierr

      character*32 ffname
      character*132 datfle
      character*1   datfle1(132)
      equivalence  (datfle,datfle1)
      integer lfname,ltrunc
      logical ifdat

      ierr=0
      if(nid.eq.0) then
         lfname = ltrunc(reafle,132) - 4
         call blank (datfle,132)
         call chcopy(datfle,reafle,lfname)
         call chcopy(datfle1(lfname+1),'_xyz.re2',8)
         write(*,*)'my_gen_re2: ',wdsize,datfle

         call byte_open(datfle, ierr)

         call blank(hdr,80)
         if(wdsize.eq.8) then
            write(hdr,112) nelgt,ldim,nelgv
         else
            write(hdr,111) nelgt,ldim,nelgv
         endif
  111    format('#v001',i9,i3,i9,' hdr')
  112    format('#v002',i9,i3,i9,' hdr')
         if(ierr.eq.0) call byte_write(hdr,20,ierr)
         if(ierr.eq.0) call byte_write(test,1,ierr) !write endian discriminator
      endif

      call gen_re2_xyz

      if (nid.eq.0) call byte_close(ierr)
      call err_chk(ierr,'Error in my_gen_re2$')

      return
      end
c-----------------------------------------------------------------------
      subroutine my_cbc_chk(s3)
c     This print types of current BC, counting as well
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer ncbc,lcbc,i,j,e,f,inid,mtype,idummy,ncbc_tmp
      parameter (lcbc=10) ! max # unique CBC
      integer cbc_count(lcbc),cbc_count_tmp(lcbc),ivlsum
      character*3 cbc_list(lcbc),cbc_list_tmp(lcbc),cbc3,s3
      logical cbc_in_list

      call izero(cbc_count,lcbc)

      ncbc = 1            ! # unique CBC
      cbc_list(1) = 'E  ' ! internal BC (we also merge empty CBC into this)

      do e=1,nelv
      do f=1,2*ldim
        cbc3 = cbc(f,e,1)
        cbc_in_list = .false.

        if (cbc3.eq.'E  '.or.cbc3.eq.'   ') then
          cbc_count(1) = cbc_count(1) + 1
          cbc_in_list = .true.
        else
          do i=2,ncbc ! go throught the registered CBC
            if (cbc3.eq.cbc_list(i)) then 
              cbc_count(i) = cbc_count(i) + 1
              cbc_in_list = .true.
            endif
          enddo
        endif

        if (.not.cbc_in_list) then ! new CBC detected, expand list
          ncbc = ncbc + 1
          if (ncbc.gt.lcbc)
     $      call exitti('BCchk: Increase lcbc$',ncbc)

          cbc_list(ncbc) = cbc3
          cbc_count(ncbc) = 1
        endif
      enddo
      enddo

      ! dbg-print
c      do i=1,ncbc
c        write(*,*) 'dbg',s3,'BC: ',nid,cbc_list(i),cbc_count(i)
c      enddo

      call nekgsync()

      ! All reduce to nid=0
      if (nid.eq.0) then
        do inid=1,np-1 ! get data from all other proc
          mtype = inid
          call csend(mtype,idummy,4,inid,0) ! handshake

          call crecv(mtype,ncbc_tmp,4)
          call crecv(mtype,cbc_list_tmp,3*ncbc_tmp,0,0)
          call crecv(mtype,cbc_count_tmp,4*ncbc_tmp,0,0)
         
          cbc_count(1) = cbc_count(1) + cbc_count_tmp(1) 
          do j=2,ncbc_tmp ! go through pending list
            cbc_in_list = .false.
            do i=2,ncbc ! search in the existed list
              if (cbc_list(i).eq.cbc_list_tmp(j)) then
                cbc_count(i) = cbc_count(i) + cbc_count_tmp(j)
                cbc_in_list = .true.
              endif
            enddo
            if (.not.cbc_in_list) then ! not found -> new CBC
              ncbc = ncbc + 1
              if (ncbc.gt.lcbc)
     $          call exitti('BCchk: Increase lcbc$',ncbc)

              cbc_list(ncbc) = cbc_list_tmp(j)
              cbc_count(ncbc) = cbc_count_tmp(j)
            endif
          enddo
        enddo
      else
        mtype = nid
        call crecv(mtype,idummy,4) ! ! handshake

        call csend(mtype,ncbc,4,0,0)
        call csend(mtype,cbc_list,3*ncbc,0,0)
        call csend(mtype,cbc_count,4*ncbc,0,0)
      endif

      call nekgsync()

      ! print
      if (nid.eq.0) then
        do i=1,ncbc
          write(*,*) s3,'BC: ',cbc_list(i),cbc_count(i)
        enddo

        idummy = ivlsum(cbc_count,ncbc)
        write(*,*) s3,'BC: ','sum',idummy
      endif

      return
      end
c-----------------------------------------------------------------------
