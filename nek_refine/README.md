# Refinement and Extrusion for Nek Meshes

This tool can be used to refined the boundary layers and generate new layers of elements via extrusion
We use the target boundary (CBC or bdryID) to tag elements, then it increases the number of elements inward (refinement) or outward (extrusion).

Currently, the code need to be refactored. However, examples should work individually.

upd 031924, rotation other than face5-6 seems inconsistent 


## API

We support several mode of refinement (TODO, rename / refactor).

Mesh = vertical $\otimes$ horizontal.        

1. `rfn_split`: boundary layer split    
   This is used to refine the boundary layers. See `demo_refine/`.

   `Nref=2` will split the tagged element into 2 elements.
   `Nref=3` will split the tagged element into 3 elements.

2. `rfn_extrude`: extrude layer (TODO, example)       
   `rfn_split` split the element inside the computational domain.   
   `rfn_extrude` extend the computational domain with surface extrusion.

    Given a surface, this will project the sirface onto another user-defined surface with single layer.     
    Then, one can apply another `rfn_split` to get multiple layers.     
    This is useful to generate inflow or outflow region based on a core mesh.    

3. `rfn_ssplit`: surface split (3d only)     
   The previous `rfn_split` is a vertical split on a given boundary surface.        
   This split is a horizontal split that refine the 2D surface mesh.       

   `Ncut=1` will split the tagged element into 4 elements.     
   `Ncut=2` will split the tagged element into 9 elements.  

Mesh = arbitrary.    

4. `rfn_msplit`: (volumetric) multiple split (3d only for now)    
   The "octant split" turn every elements into 8 smaller elements.   
   To generalize it, 

   `Ncut=1` will split the tagged element into 8 elements.        
   `Ncut=2` will split the tagged element into 27 elements.     


5. `rfn_asplit` (TODO: not implemented) (volumetric) adaptive split     
   Adaptive split based on a given density.        


TODO:       
table of APIs, images      
table capability, reuse, con/co2, 



## Usage (general workflow)

1. Include all `*.f` files into the usr file (The code is not finalized. See each example)      
   The `my_errchk.f` is used to do some check and dumping meshes.       
   The `my_post.f` is needed by errchk which could dump elements near bad location.    
   You also need to carry `REFINE` for all working variables.   

2. (usrdat2) backup boundary conditions in original re2     
   `call rfn_backup_boundary_id(1)`

3. (usrdat2) Set up dummy CBC or boundaryID     
   Also zero out igroup to avoid un-initialized numbers  

4. Add the dummy subroutine `usr_extrude_pj`       
   This is how we define geometry for extrusion, see example `demo_extrude/(TODO)`        
   For other splits, leave it blank.      

5. (usrchk) call the interface functions     
   See examples to set up all the variables.      

6. (usrchk) Dump new mesh `rfn_dump_mesh`

7. call exitt, use the new mesh at a new run    
   Nek won't continue the run with a newly updated mesh due to the lack of connectivity.      
   If someone finish this, maybe we can do mesh refinement at runtime, or even adaptive mesh refinement.    

Extra notes:   
a. Make sure `lelg` and `lelt` are large enough for the new mesh.       
   Especially, if the tagged elements are all in one processor, this will create unbalance need for a larger `lelt`.     
b. `lref` (or `lcut`) inside `REFINE` has maximium number of layers. Increase this if needed.


## Cases

`demo` are cases with step-by-step guide. `pb` case is one of the pebble bed case.

- `demo_refine/`
- `demo_extrude/ (TODO)`   
- `demo_ssplit/`
- `demo_msplit/`

- `pb167s1t/`



## Note \& Features

Features:
- Favor high order to preserve curvature.    
- Multiple layers with limiter
- Run in parallel
- Support 2D/3D
- Watertight (requred a special version of parRSB for parallel sorting, FIXME)      
  input co2/ma2, dump co2 for new meshes     
  Please ignore newly generated co2 file for now unless you have `FPSORT_INT8` in refine.f
- Extrusion onto defined geometry
- Reusability: multiple refinement + extrusion combination (not well-tested)

Notes:
- User can use this tool to generate non-conformal mesh.     
  It's recommended to only refine the boundary layers that already has a skin layer.      
- (`rfn_split`) Tested with my pebble bed mesh up to 98M elements.
- Treatment for periodic BC is not implemented but I expect it work as long as tagged elements don't have periodic BC (not tested).

Future Work:
- refactor, clean up...
- Add more refinement, ex. octant (1 hex to 8 hex), directional, etc.
- Suppor Nek function such as dist2, `fix_geom`, etc to enhance further capability.
- On-the-fly simulation (note: unbalanced workload or redistribute? Interpolation)      
  If under-resolve, refine mesh and continue the simulations.  


