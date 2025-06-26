c-----------------------------------------------------------------------
      subroutine rotate_z_to_x
      ! the mesh was extruded along z, here we rotate it to x
      implicit none
      include 'SIZE'
      include 'GEOM'

      real xo(lx1,ly1,lz1,lelt)
      real yo(lx1,ly1,lz1,lelt)
      real zo(lx1,ly1,lz1,lelt)
      integer n

      n = lx1*ly1*lz1*nelt
      call copy(xo,xm1,n)
      call copy(yo,ym1,n)
      call copy(zo,zm1,n)

      call copy(zm1,xo,n)
      call chsign(zm1,n)
      call copy(xm1,zo,n)

      call fix_geom

      return
      end
