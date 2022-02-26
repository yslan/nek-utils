# nek-utils   
A place to put my codes for Nek5000/NekRS.   
All tools are under a separate directory or a git-repository with README.md description and with proper example(s).   

This is still under-construction. 
Most codes are based on existed Nek5000 code. 
I will put proper reference and lisence later.

Some codes are experimental. Feel free to report issues and bugs.   

## Overview

- `outpost-select/`
   - __Main__: A modified nek `outpost` allowing user to only dump the user-selected elements to reduce the size of the output file and the visualization cost.   
   - __Usage__: independent file included from a nek user file.   


- `nek_connectivity/`
   - __Main__: Extra nek tools for the connectivity file `.con/.co2` such as `n2to3` and binary-ascii translation.   
   - __Usage__: [repo] Modified nek5000 tools.   
   - __TODO__: Merge into `n2to3`   


- `SphereMesh/`
   - __Main__: Mesh generator for pebble bed geometry.   
   - __Usage__: [repo] A workflow contains MATLAB code and nek mesh smoother.   
   - __TODO__: WIP, to be released in future  


- `Tet2Hex_forNek/`
   - __Main__: [repo] Mesh generator for arbitrary input points via tet-to-hex. 
   - __Usage__: MATLAB code that generates nek meshes.
   - __TODO__: Add new features as a small subset of SphereMesh      


## TODO:
- mesh smoother   
- mesh refine / extrusion
- multi-box reader
- hpts example(s)
- web based doc
- gallery and showcases
- SEM codes (both MATLAB and Python)
- external storage for large files

