Several modified gfldr.

This is used to kickstart the 2Bn case with a reduced elements count (5x)

- `my_gfldr_v.f` that only fills the fluid domain.

  Usage:
  ```
  #include <my_gfldr_v.f>
  
        subroutine userchk
        include 'SIZE'
        include 'TOTAL'
  
        ! fill fluid domain with gfldr
        if (istep.eq.0) call my_gfldr_v('pink0.run')
  
        return
        end
  ```


- `my_gfldr_mesh.f` reads mesh from another file (`msh.fld`) in case the source field has no mesh.

  Usage:
  ```
  #include <my_gfldr_mesh.f>
  
        subroutine userchk
        include 'SIZE'
        include 'TOTAL'
  
        if (istep.eq.0) call my_gfldr_mesh('msh.fld ', 'b.fld ')

        return
        end
  ```

- `my_gfldr_b.f` query meshes in batches. To avoid crystal router `MAX_INT` error.

  Warning: this will reads the file multiple files to loop over query batches
  ```
  findpts_setup
  
  do ibatch=1,nbatch
  
    findpts      ! reduced size
    read a field ! extra reads with full nels
    findpts_eval ! reduce size
  
  enddo
  ```

  Usage:
  ```
  #include <my_gfldr_b.f>

        subroutine userchk
        include 'SIZE'
        include 'TOTAL'

        ! load 8192 elements at a batch, loop over batches
        if (istep.eq.0) call my_gfldr_b('r.fld ', 8192)

        return
        end
  ```
