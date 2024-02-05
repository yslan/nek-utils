# Refinement and Extrusion for Nek Meshes

This tool can be used to refined the boundary layers and generate new layers of elements via extrusion
We use the target boundary (CBC or bdryID) to tag elements, then it increases the number of elements inward (refinement) or outward (extrusion).

Currently, we only support splitting in one direction. 

## Usage

1. Include `refine.f`, `my_errchk.f`, `my_post.f` into the usr file         
   Main algorithm is inside `refine.f`
   The `my_errchk.f` is used to do some check and dumping meshes.       
   The `my_post.f` is needed by errchk which could dump elements near bad location.
   You also need to carry `REFINE` for all working variables.   

2. (usrdat2) backup boundary conditions in original re2
   `call rfn_backup_boundary_id(1)`

3. (usrdat2) Set up dummy CBC or boundaryID     
   Also zero out igroup to avoid un-initialized numbers  

4. Add the dummy subroutine `usr_extrude_pj`       
   This is how we define geometry for extrusion, see example `demo_extrude/`     
   For refinement, leave it black.

5. (usrchk) call the interface functions     
   See examples to set up all the variables.      

6. (usrchk) Dump new mesh `rfn_dump_mesh`

7. call exitt, use the new mesh at a new run    
   Nek won't run on the newly generated mesh due to the lack of connectivity.      
   If someone finish this, maybe we can do mesh refinement at runtime, or even adaptive mesh refinement.

Extra    
a. Make sure `lelg` and `lelt` are large enough for the new mesh.       
   Especially, if the tagged elements are all in the same processor, this will create unbalance need for large a larger `lelt`.
b. `lref` inside `REFINE` has maximium number of layers. Increase this if needed


## Cases

`demo` are 2D cases with step-by-step guide. `pb` case is one of the pebble bed case.

- `demo_refine/`

- `demo_extrude/`

- `pb167s1t/`




## Note \& Features

Features:
- Multiple layers with limiter
- Run in parallel
- Support 2D/3D
- Watertight (requred a special version of parRSB for parallel sorting, FIXME)      
  input co2/ma2, dump co2 for new meshes
- Extrusion onto defined geometry
- Reusability: multiple refinement + extrusion combination (not well-tested)

Notes:
- User can use this tool to generate non-conformal mesh.     
  It's recommanded to only refine the boundary layers that already has a skin layer.

- Tested with my pebble bed mesh up to 98M elements
- Treatment for periodic BC is not implemented but I expect it work as long as tagged elements don't have periodic BC (not tested).

Future Work:
- Add more refinement, ex. octant (1 hex to 8 hex), directional, etc.
- Suppor Nek function such as dist2, `fix_geom`, etc to enhance further capability.
- On-the-fly simulation (note: unbalanced workload or redistribute? Interpolation)      
  If under-resolve, refine mesh and continue the simulations.  


