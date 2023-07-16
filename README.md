# nek-utils   
A place to put my codes for Nek5000/NekRS.   
All tools are under a separate directory or a git-repository with README.md description and with proper example(s).   

This is still under-construction. 
Most codes are based on existed Nek5000 code. 
I will put proper reference and lisence later.

Some codes are experimental. Feel free to report issues and bugs.   

## Overview

- `usr_extra/`
   - __Main__: This store those isolated `*.f` files that I usually use is to do some checks or post-processing
   - __Usage__: Include the file in .usr file. See the README.md for individual subroutines.
   - __TODO__: I just start collecting them together. Will remember to copy files here when I see them.

- `outpost-select/`
   - __Main__: A modified nek `outpost` allowing user to only dump the user-selected elements to reduce the size of the output file and the visualization cost.   
   - __Usage__: independent file included from a nek user file.   


- `nek_connectivity/`
   - __Main__: Extra nek tools for the connectivity file `.con/.co2` such as `n2to3` and binary-ascii translation.   
   - __Usage__: [repo] Modified nek5000 tools.   
   - __TODO__: Merge into `n2to3`   


- `SphereMesh/`
   - __Main__: Mesh generator for pebble bed geometry.   
   - __Usage__: [private repo] A workflow contains MATLAB code and nek mesh smoother.   
   - __TODO__: WIP, to be released in future  


- `Tet2Hex_forNek/`
   - __Main__: [repo] Mesh generator for arbitrary input points via tet-to-hex. 
   - __Usage__: MATLAB code that generates nek meshes.
   - __TODO__: Add new features as a small subset of SphereMesh      

- `quick_fix_boundary_layers/`
  - __Main__:  fix the boundary layers for cylinder geometry   
  - __Usage__: independent file included from a nek user file. 
  - __TODO__: clean up, add general geometry
 

- `del_elements/`
  - __Main__: delete selected elements from mesh in nek
  - __Usage__: independent file included from a nek user file. 
  - __TODO__: add feature like BC tracing and dump co2 files

## TODO:
- mesh smoother   
- mesh refine / extrusion
- multi-box reader
- hpts example(s)
- web based doc
- gallery and showcases
- SEM codes (both MATLAB and Python)
- external storage for large files

