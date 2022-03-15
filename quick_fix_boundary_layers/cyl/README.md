### Cylinder in a box

To construct/alter the mesh:
  - adjust the box files (inside `mk_mesh/`)    
    `ringbox1.box`: concentric circles   
    `ringbox2.box`: circle in a 2d box
  - `./mkmsh` or `./mkmsh3d` (define the path to Nek5000's tool inside)

### Interface:

called from `userdat2`

```
      nlyr = 3   ! numbers of the layers from bdry
      itype = 1  ! for future extension

      ! axis of the cylinder, assume along z axis
      ginfo(1) = 0.0 ! x coord
      ginfo(2) = 0.0 ! y coord

      call fix_bdry('W  ',itype,ginfo,nlyr) ! boundary layers on 'W  ' bdry
```


