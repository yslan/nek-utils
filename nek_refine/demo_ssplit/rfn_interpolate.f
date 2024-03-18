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

