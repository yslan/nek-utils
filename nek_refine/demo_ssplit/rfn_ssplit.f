      subroutine rfn_ssplit_work(ie_ref,Ncut)
c     On a given surface, do surface split that a quad face is splitted into 
c     (Ncut+1)^2 elements. The vertical direction remains the same.
      implicit none
      include 'SIZE'
      include 'REFINE'
      integer ie_ref(lelt),Ncut,Nseg,iseg,jseg,Nnew_p_e,nxyz
      integer e,enow,enew,e1,e2,f

      integer ix,iy,iz ! debug

      real jmat,jtmat
      common /my_ssplit_int/ jmat(lx1*lx1,lseg),jtmat(lx1*lx1,lseg)

      real xbak(lx1,ly1,lz1)
      real ybak(lx1,ly1,lz1)
      real zbak(lx1,ly1,lz1)
      real xnew(lx1,ly1,lz1)
      real ynew(lx1,ly1,lz1)
      real znew(lx1,ly1,lz1)

      Nseg = Ncut + 1 ! # segment per direction
      Nnew_p_e = Nseg*Nseg-1 ! each element will create m^2-1 new elements
      nxyz = lx1*ly1*lz1

      call set_int_mat(jmat,jtmat,lx1,Nseg)

      enow = 0
      enew = 0
      do e=1,nelt0
      if (ie_ref(e).gt.0) then
        enow = enow + 1

        call copy(xbak,xmc(1,1,1,enow),nxyz)
        call copy(ybak,ymc(1,1,1,enow),nxyz)
        call copy(zbak,zmc(1,1,1,enow),nxyz)

        do jseg=1,Nseg
        do iseg=1,Nseg

          call int_ssplit_elem(xnew,xbak,lx1,iseg,jseg)
          call int_ssplit_elem(ynew,ybak,lx1,iseg,jseg)
          call int_ssplit_elem(znew,zbak,lx1,iseg,jseg)

          if (iseg.eq.Nseg.AND.jseg.eq.Nseg) then ! copy to old elem

            call copy(xmc(1,1,1,enow),xnew,nxyz)
            call copy(ymc(1,1,1,enow),ynew,nxyz)
            call copy(zmc(1,1,1,enow),znew,nxyz)

            ! bc 
            if (Nseg.ne.1) then ! this should always be true
              call rzero(bc_wk(1,1,enow,1),5)
              call rzero(bc_wk(1,1,enow,2),5)
              call rzero(bc_wk(1,4,enow,1),5)
              call rzero(bc_wk(1,4,enow,2),5)
              bID_wk(1,e1,1) = 0
              bID_wk(1,e1,2) = 0
              bID_wk(4,e1,1) = 0
              bID_wk(4,e1,2) = 0
              cbc_wk(1,e1,1) = 'E  '
              cbc_wk(1,e1,2) = 'E  '
              cbc_wk(4,e1,1) = 'E  '
              cbc_wk(4,e1,2) = 'E  '
            endif

            ! vtx ! TODO
            
          else                                    ! create new elem

            enew = enew + 1
            e1 = nel_lv + enew ! local (work)
            e2 = nelt0 + enew ! local (mpi rank)

            call copy(xmc(1,1,1,e1),xnew,nxyz)
            call copy(ymc(1,1,1,e1),ynew,nxyz)
            call copy(zmc(1,1,1,e1),znew,nxyz)

            ! rotate index
            ie_ref(e2) = ie_ref(e)

            ! bc 
            do f=1,2*ldim
              call rzero(bc_wk(1,f,e1,1),5)
              call rzero(bc_wk(1,f,e1,2),5)
              bID_wk(f,e1,1) = 0
              bID_wk(f,e1,2) = 0
              cbc_wk(f,e1,1) = 'E  '
              cbc_wk(f,e1,2) = 'E  '
            enddo
            if (iseg.eq.1)    call rfn_copy_cbc_face(e1,4,enow,4)
            if (iseg.eq.Nseg) call rfn_copy_cbc_face(e1,2,enow,2)
            if (jseg.eq.1)    call rfn_copy_cbc_face(e1,1,enow,1)
            if (jseg.eq.Nseg) call rfn_copy_cbc_face(e1,3,enow,3)
            call rfn_copy_cbc_face(e1,5,enow,5)
            call rfn_copy_cbc_face(e1,6,enow,6)

            ! vtx ! TODO

          endif

        enddo
        enddo

        nelv = nelv + Nnew_p_e
        nelt = nelt + Nnew_p_e

      endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_int_mat(jmat,jtmat,n,Nseg)
      implicit none
      integer Nseg,iseg,i,n
      real jmat(n*n,1),jtmat(n*n,1),wtmp(n),z(n),w(n)
      real z0,z1,zdel

      call zwgll (z,wtmp,n)
      zdel = 2.0 / real(Nseg)

      do iseg=1,Nseg

        z0 = -1.0 + (iseg-1)*zdel
        z1 = z0 + zdel

        do i=1,n
          w(i) = (z(i)+1.0)*0.5*(z1-z0) + z0
        enddo

        call gen_int_gz(jmat(1,iseg),jtmat(1,iseg),w,n,z,n)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine int_ssplit_elem(uout,uin,n,i,j)
      implicit none
      include 'SIZE'
      include 'REFINE'
      real uout(1),uin(1)
      integer n,i,j,nn,l,iv,iw,k

      real jmat,jtmat
      common /my_ssplit_int/ jmat(lx1*lx1,lseg),jtmat(lx1*lx1,lseg)

      real w,v
      parameter(l=50)
      common /my_ssplit_int_wk/ w(l*l*l),v(l*l*l)

      nn = n*n

      call mxm(jmat(1,i),n,uin,n,v ,nn)
      iv=1
      iw=1
      do k=1,n
         call mxm(v(iv),n,jtmat(1,j),n,uout(iw),n)
         iv = iv+nn
         iw = iw+nn
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rfn_copy_cbc_face(e1,f1,e2,f2)
c     copy cbc_wk from (e2,f2) to (e1,f1)
      implicit none
      include 'SIZE'
      include 'REFINE'
      integer e1,e2,f1,f2,islot

      do islot=1,2
        call copy(bc_wk(1,f1,e1,islot),bc_wk(1,f2,e2,islot),5)
        bID_wk(f1,e1,islot) = bID_wk(f2,e2,islot)
        cbc_wk(f1,e1,islot) = cbc_wk(f2,e2,islot)
      enddo

      return
      end
