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
