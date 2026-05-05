## Extra files 

I will (try to) collect my files here. The rule is, if I use it more than three times (and if I remember), I will save it here.
These will be isolated tiny functions that I think it won't go into the main repo.  

Please create issues for questions or bugs. Happy debugging. 

- `my_cbc_chk.f`:   
  __Usage__: quick debugging, print all distinct BC types and count the number.

   You can call the subroutine `my_cbc_chk(s3)` anywhere in userdat2, userdat3, or userchk with different tags like
   ```
   call my_cbc_chk('aaa')
   call my_cbc_chk('bbb')
   ```
   It will print something like this with ifield, type, CBC, #faces
   ```
   aaa BC:   0 PR   E         768
   aaa BC:   1 VEL  E         608
   aaa BC:   1 VEL  P          96
   aaa BC:   1 VEL  W          64
   aaa BC:   2 S00  E         608
   aaa BC:   2 S00  P          96
   aaa BC:   2 S00  t          64
   aaa BC:   3 MHD  E         608
   aaa BC:   3 MHD  P          96
   aaa BC:   3 MHD  W          64
   ```
   This subroutine is intensively used as my personal debugging routine. Works for many cases up to 72000 ranks.
   Nek repo uses this https://github.com/Nek5000/Nek5000/pull/787 to print BCs which will give you UNKNOWN when the BC is not in the table.  This one is more versatile. 


- `chk_con_via_bdry.f`:
  __Usage__: This check the connectivity by counting the boundaryID. 
  ```
  call chk_bdry(1,nelv,boundaryID,'vvv')
  call chk_bdry(2,nelt,boundaryIDt,'ttt')
  ```
  It will first print out how many faces per bcid, then compute the total number of faces that is assigned to a BC.
  That number must match with the faces with face center has value 1 after `dssum(1)`.
  Otherwise, either the connectivity is leaking or some BC is not assigned properly.
  ```
  chk bdryidtvvv            1          324
  chk bdryidtvvv            2          324
  chk bdryidtvvv            3     10935342
  chk bdryidtvvv total            3     10935990     10935990
  chk bdryidtttt            1          324
  chk bdryidtttt            2          444
  chk bdryidtttt            3     10935342
  chk bdryidtttt total            3     10936110     10936110
  ```
  

- `flip_elements/`: flip from lhs to rhs. for wired wrap with mirror direction

- `change_time.py:

   This change the time in header from atime to timestep*1.0 in order to bypass Paraview's struggle of reading nonincreasing time / timestep.
   ```
   python3 ./change_time.py avgeddy.nek5000
   ```

- gfldr variants     
  [gfldr subfolder](./gfldr/)


- rea2vtk.py: Convert rea/re2 into vtk for mesh inspection.    
  (`rea2vtk_mpi.py`): MPI version from re2 to pvtu

  ```
  python3 rea2vtk.py input.rea out.vtk
  mpiexec -np 2 python3 rea2vtk_mpi.py input.rea out.pvtu
  ```
  - It supports 2D and 3D, rea and re2 (detected from the input extension).
  - Depends on `pyvista`
  - re2 reader modifed from kth groups' pymech [here](https://github.com/eX-Mech/pymech/blob/main/src/pymech/neksuite/mesh.py)
  - has tested with a E=2M mesh

  TODO: I don't like the current ordering of xyz ndarray    
  TODO: Process data blocks by blocks to reduce memory requirement for 10x larger case.
  TODO: check re2 ver
  TODO: add examples
  TODO: rea reader is serial in rea2vtk_mpi.py

  quick test, on my laptop, Intel Core i9-12900H (12th-gen Alder Lake, 6 P-cores + 8 E-cores)
  ```
  read header ...       hdr: 62132 3 62132

  #rank    read        write       read eff    write eff
  1        1.0180e-02  4.0950e-01  1.000       1.000
  2        7.1593e-03  2.1050e-01  0.711       0.973
  4        9.4292e-03  1.2135e-01  0.270       0.844
  8        2.0745e-02  8.2203e-02  0.061       0.623
  ```

- `rotate_z_to_x.f`:
  Often, we build mesh via n2to3 which exctudes in z. However, many engineer problems set streamwise direction in x.
  Call this in usrdat2 to rotate coordinates from z+ to x+ to have peace in mind in postprocessing.

