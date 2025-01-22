c-----------------------------------------------------------------------
      subroutine chk_bdry(ifld,nel,bdryid,s3)
      implicit none
      include 'SIZE'
      include 'TSTEP' ! ifield
c      include 'TOTAL'

      integer ifld,nel,bdryid(2*ldim,lelt)
      character*3 s3

      real tmp(lx1,ly1,lz1,lelt)
      integer maxbdry, ibc, nbc, nbcface0, nbcface1
      integer  e, f, ix, iy, iz
      parameter(maxbdry=99)
      integer nbdrys(maxbdry),iwork(maxbdry)
      integer iglsum,iglmax

      integer ix_fctr(6), iy_fctr(6), iz_fctr(6), idd
      parameter(idd=ldim-1)
      data ix_fctr /  2,lx1,  2,  1,  2,  2/
      data iy_fctr /  1,  2,ly1,  2,  2,  2/
      data iz_fctr /idd,idd,idd,idd,  1,lz1/

      ifield=ifld

      call rone(tmp,lx1*ly1*lz1*lelt)
      call izero(nbdrys,maxbdry)
      call dssum(tmp,lx1,ly1,lz1)

      nbcface0 = 0
      nbcface1 = 0

      nbc = 0
      do e=1,nel
      do f=1,2*ldim

        ibc = bdryid(f,e)
        if (ibc.gt.0) then
          nbdrys(ibc) = nbdrys(ibc) + 1
          nbcface0 = nbcface0 + 1
          nbc = max(nbc,ibc)
        endif

        ix = ix_fctr(f)
        iy = iy_fctr(f)
        iz = iz_fctr(f)
        if (tmp(ix,iy,iz,e).lt.1.5) nbcface1 = nbcface1 + 1

      enddo
      enddo
      nbcface0 = iglsum(nbcface0,1)
      nbcface1 = iglsum(nbcface1,1)

      nbc = iglmax(nbc,1)
      call igop(nbdrys,iwork,'+  ',nbc)

      ! if any one of nbdrys(ibc) is zero, there is an error in connecivity or bc
      ! Also, nbcface0 shoud be the same as nbcface1
      if (nio.eq.0) then
        do ibc=1,nbc
          write(*,*)'chk bdryidt',s3,ibc,nbdrys(ibc)
        enddo
        write(*,*)'chk bdryidt',s3,' total',nbc,nbcface0,nbcface1
      endif

      return
      end
