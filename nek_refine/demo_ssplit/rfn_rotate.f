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

