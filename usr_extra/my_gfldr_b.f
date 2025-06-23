#ifndef NOMPIIO

      subroutine my_gfldr_b(nel_block,sourcefld)
c
c     generic field file reader
c     reads sourcefld and interpolates all avaiable fields
c     onto current mesh
c
c     memory requirement: 
c     nelgs*nxs**ldim < np*(4*lelt*lx1**ldim)
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'GFLDR'

      character sourcefld*(*)

      common /scrcg/  pm1(lx1*ly1*lz1,lelv)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      character*1   hdr(iHeaderSize)

      integer nel_block
      integer*8 dtmp8

      logical if_byte_swap_test
      real*4 bytetest

      etime_t = dnekclock_sync()
      if(nio.eq.0) write(6,*) 'call my_gfldr_b '
     $                      , nel_block,trim(sourcefld) 

      ! open source field file
      ierr = 0
      if(nid.eq.0) then
        open (90,file=sourcefld,status='old',err=100)
        close(90)
        goto 101
 100    ierr = 1
 101  endif
      call err_chk(ierr,' Cannot open source fld file!$')
      call byte_open_mpi(sourcefld,fldh_gfldr,.true.,ierr)
      if(nio.eq.0) write(6,*) 'my_gfldr_b file opened'

      ! read and parse header
      call byte_read_mpi(hdr,iHeaderSize/4,0,fldh_gfldr,ierr)
      call byte_read_mpi(bytetest,1,0,fldh_gfldr,ierr)

      call mfi_parse_hdr(hdr,ierr)
      call err_chk(ierr,' Invalid header!$')
      ifbswp = if_byte_swap_test(bytetest,ierr)
      call err_chk(ierr,' Invalid endian tag!$')
      if(nio.eq.0) write(6,*) 'my_gfldr_b header is read'

      nelgs   = nelgr
      nxs     = nxr
      nys     = nyr
      nzs     = nzr
      if(nzs.gt.1) then 
        ldims = 3
      else
        ldims = 2
      endif
      if (ifgtim) time = timer

      ! distribute elements across all ranks
      nels = nelgs/np
      do i = 0,mod(nelgs,np)-1
         if(i.eq.nid) nels = nels + 1
      enddo
      nxyzs      = nxs*nys*nzs
      dtmp8      = nels
      ntots_b    = dtmp8*nxyzs*wdsizr
      rankoff_b  = igl_running_sum(nels) - dtmp8
      rankoff_b  = rankoff_b*nxyzs*wdsizr  
      dtmp8      = nelgs
      nSizeFld_b = dtmp8*nxyzs*wdsizr
      noff0_b    = iHeaderSize + iSize + iSize*dtmp8

      ! do some checks
      if(ldims.ne.ldim) 
     $ call exitti('ldim of source does not match target!$',0)
      if(ntots_b/wdsize .gt. ltots) then
        dtmp8 = nelgs
        lelt_req = dtmp8*nxs*nys*nzs / (np*ltots/lelt)
        lelt_req = lelt_req + 1
        if(nio.eq.0) write(6,*)
     $   'ABORT: buffer too small, increase lelt > ', lelt_req
        call exitt
      endif
      if(nio.eq.0) write(6,*) 'my_gfldr_b nelgs',nelgs,nels

      ifldpos = 0
      if(ifgetxr) then
        ! read source mesh coordinates
        call gfldr_getxyz(xm1s,ym1s,zm1s)
        ifldpos = ldim
      else
        call exitti('source does not contain a mesh!$',0)
      endif

      if(if_full_pres) then
        call exitti('no support for if_full_pres!$',0)
      endif
      if(nio.eq.0) write(6,*) 'my_gfldr_b mesh is read'

      ! initialize interpolation tool using source mesh
      nxf   = 2*nxs
      nyf   = 2*nys
      nzf   = 2*nzs
      nhash = nels*nxs*nys*nzs 
      nmax  = 128

      call fgslib_findpts_setup(inth_gfldr,nekcomm,np,ldim,
     &                          xm1s,ym1s,zm1s,nxs,nys,nzs,
     &                          nels,nxf,nyf,nzf,bb_t,
     &                          nhash,nhash,nmax,tol)
      if(nio.eq.0) write(6,*) 'my_gfldr_b findpts_setup done'

      if(ifgetur) then 
        if(.not.ifgfldr.or.ifgetu) then !skip if this is a restart call and the scalar isn't requested
          if(nio.eq.0) write(6,*) 'my_gfldr_b reading vel'
          call my_gfldr_loop(vx,vy,vz,nelv,nel_block,ldim,ifldpos+1)
        endif
        ifldpos = ifldpos + ldim
      endif
      if(ifgetpr) then
        if(.not.ifgfldr.or.ifgetp) then !skip if this is a restart call and the scalar isn't requested
          if(nio.eq.0) write(6,*) 'my_gfldr_b reading pr'
          call my_gfldr_loop(pm1,dum,dum,nelv,nel_block,1,ifldpos+1)
          if (ifaxis) call axis_interp_ic(pm1) ! move to the end
          call map_pm1_to_pr(pm1,1)
        endif
        ifldpos = ifldpos + 1
      endif
      if(ifgettr .and. ifheat) then
        if(.not.ifgfldr.or.ifgett) then !skip if this is a restart call and the scalar isn't requested
          if(nio.eq.0) write(6,*) 'my_gfldr_b reading temp'
          call my_gfldr_loop(t,dum,dum,nelt,nel_block,1,ifldpos+1)
        endif
        ifldpos = ifldpos + 1
      endif
      do i = 1,ldimt-1
        if(ifgtpsr(i)) then
          if(.not.ifgfldr.or.ifgtps(i)) then !skip if this is a restart call and the scalar isn't requested
            if(nio.eq.0) write(6,*) 'my_gfldr_b reading scalar',i
            ntot = nx1*ny1*nz1*mm
            call my_gfldr_loop(t(1,1,1,1,i+1),dum,dum
     $                        ,nelt,nel_block,1,ifldpos+1)
          endif
          ifldpos = ifldpos + 1
        endif
      enddo

      if(nio.eq.0) write(6,*) 'my_gfldr_b findpts getfld done'
      call byte_close_mpi(fldh_gfldr,ierr)
      etime_t = dnekclock_sync() - etime_t
      call fgslib_findpts_free(inth_gfldr)
      if(nio.eq.0) write(6,'(A,1(1g9.2),A)')
     &                   ' done :: my_gfldr_b  ', etime_t, ' sec'

  51  format(a,i3,'/',i3,' kk=',i9,' mm=',i9)
      return
      end
cc-----------------------------------------------------------------------
      subroutine my_gfldr_loop(out1,out2,out3
     $                        ,nel_in,nel_block,nldim,ifldpos)
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'GFLDR'

      integer*8 i8glsum,nfail,nfail_sum

      real out1(*)
      real out2(*)
      real out3(*)
      integer kk, mm, nel_in, nel_block, ipass, npass
      character*3 s3

      call my_gfldr_getfld_buf(nldim,ifldpos) ! read data

      npass = 1 + nel_in / nel_block
      kk = 1
      do ipass = 1,npass ! interp data block by block
        mm = nel_in - kk + 1
        mm = min(mm,nel_block)

        if(nio.eq.0) write(6,51)'my_gfldr loop',ipass,npass,kk,mm

        ntot = nx1*ny1*nz1*mm

        ! locate points (iel,iproc,r,s,t)
        nfail = 0
        toldist = 5e-6
        if(wdsizr.eq.8) toldist = 5e-14

        call fgslib_findpts(inth_gfldr,
     &                      grcode,1,
     &                      gproc,1,
     &                      gelid,1,
     &                      grst,ldim,
     &                      gdist,1,
     &                      xm1(1,1,1,kk),1,
     &                      ym1(1,1,1,kk),1,
     &                      zm1(1,1,1,kk),1,ntot)
        if(nio.eq.0) write(6,*) 'my_gfldr loop findpts done',ipass

        do i=1,ntot
           if(grcode(i).eq.1 .and. sqrt(gdist(i)).gt.toldist)
     &       nfail = nfail + 1
           if(grcode(i).eq.2) nfail = nfail + 1
        enddo

        nfail_sum = i8glsum(nfail,1)
        if(nfail_sum.gt.0) then
          if(nio.eq.0) write(6,*)
     &      ' WARNING: Unable to find all mesh points in source fld ',
     &      nfail_sum
        endif

        ! read source fields and interpolate
        call my_gfldr_getfld_intp(out1,out2,out3,kk,ntot,nldim)

        kk = kk + mm
      enddo ! ipass

  51  format(a,i3,'/',i3,' kk=',i9,' mm=',i9)
      return
      end
cc-----------------------------------------------------------------------
c      subroutine gfldr_getxyz(xout,yout,zout)
c
c      include 'SIZE'
c      include 'GFLDR'
c      include 'RESTART'
c
c      real xout(*)
c      real yout(*)
c      real zout(*)
c
c      integer*8 ioff_b
c
c 
c      ioff_b = noff0_b + ldim*rankoff_b
c      call byte_set_view(ioff_b,fldh_gfldr)
c
c      nread = ldim*ntots_b/4
c      call byte_read_mpi(bufr,nread,-1,fldh_gfldr,ierr)
c      if(ifbswp) then
c        if(wdsizr.eq.4) call byte_reverse (bufr,nread,ierr)
c        if(wdsizr.eq.8) call byte_reverse8(bufr,nread,ierr)
c      endif
c
c      call gfldr_buf2vi (xout,1,bufr,ldim,wdsizr,nels,nxyzs)
c      call gfldr_buf2vi (yout,2,bufr,ldim,wdsizr,nels,nxyzs)
c      if(ldim.eq.3)
c     $ call gfldr_buf2vi(zout,3,bufr,ldim,wdsizr,nels,nxyzs)
c
c      return
c      end
c-----------------------------------------------------------------------
      subroutine my_gfldr_getfld_buf(nldim,ifldpos)

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'RESTART'

      integer*8 ioff_b

      logical ifpts

      integer icalld
      save    icalld
      data    icalld /0/

      ifpts = .false.
      if(icalld.eq.0) then
        ifpts = .true. ! find points
        icalld = 1
      endif

      ! read field data from source fld file
      ioff_b = noff0_b + (ifldpos-1)*nSizeFld_b
      ioff_b = ioff_b  + nldim*rankoff_b
      call byte_set_view(ioff_b,fldh_gfldr)
      nread = nldim*ntots_b/4
      call byte_read_mpi(bufr,nread,-1,fldh_gfldr,ierr)
      if(ifbswp) then
        if(wdsizr.eq.4) call byte_reverse (bufr,nread,ierr)
        if(wdsizr.eq.8) call byte_reverse8(bufr,nread,ierr)
      endif

      return
      end

      subroutine my_gfldr_getfld_intp(out1,out2,out3,e,nout,nldim)

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'RESTART'

      integer e, i
      real out1(*)
      real out2(*)
      real out3(*)

      i = (e-1)*lx1*ly1*lz1 + 1

      ! interpolate onto current mesh
      call gfldr_buf2vi  (buffld,1,bufr,nldim,wdsizr,nels,nxyzs)
      call gfldr_intp    (out1(i),nout,buffld,ifpts)
      if(nldim.eq.1) return

      call gfldr_buf2vi  (buffld,2,bufr,nldim,wdsizr,nels,nxyzs)
      call gfldr_intp    (out2(i),nout,buffld,.false.)
      if(nldim.eq.2) return

      if(nldim.eq.3) then
        call gfldr_buf2vi(buffld,3,bufr,nldim,wdsizr,nels,nxyzs)
        call gfldr_intp  (out3(i),nout,buffld,.false.)
      endif

      return
      end
c-----------------------------------------------------------------------
c      subroutine gfldr_buf2vi(vi,index,buf,ldim,wds,nel,nxyz)
c
c      real    vi(*)
c      real*4  buf(*)
c      integer wds
c
c
c      do iel = 1,nel
c         j = (iel-1)*nxyz
c         k = (iel-1)*ldim*nxyz
c
c         if(index.eq.2) k = k+nxyz
c         if(index.eq.3) k = k+2*nxyz
c
c         if(wds.eq.4) call copy4r(vi(j+1),buf(k+1)  ,nxyz)
c         if(wds.eq.8) call copy  (vi(j+1),buf(2*k+1),nxyz)
c      enddo
c
c      return
c      end
cc-----------------------------------------------------------------------
c      subroutine gfldr_intp(fieldout,nout,fieldin,iffpts)
c
c      include 'SIZE'
c      include 'RESTART'
c      include 'GEOM'
c      include 'GFLDR'
c
c      real    fieldout(nout)
c      real    fieldin (*)
c
c      ! evaluate inut field at given points
c      npt = nout
c      call fgslib_findpts_eval(inth_gfldr,
c     &                         fieldout,1,
c     &                         grcode,1,
c     &                         gproc,1,
c     &                         gelid,1,
c     &                         grst,ldim,npt,
c     &                         fieldin)
c
c      return
c      end
cc-----------------------------------------------------------------------
#else
      subroutine my_gfldr_b(sourcefld)

      character sourcefld*(*)

      call exitti("MPIIO needed for my_gfldr_b!$",0)

      return
      end
#endif
