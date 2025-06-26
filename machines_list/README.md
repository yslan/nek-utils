# Machines List

My collection of the toys with the configuration/instruction/doc for running Nek5000/NekRS, including modules list, submit scripts, and shortcuts.

Disclaimer:       
- Using scripts from a stranger is dangerous. **Use them at your own risk!**   
- This is, and will always be, under construction.    
  Those scripts might be out-dated. However, it should still work with minimal adjustment unless the machine has gone through a major upgrade.      
  Please feel free to contact me for any issues or if you find new changes.    

- This is simply a note for special machines, NOT a general guide to build NekRS.    
  Please don't simply copy whatever in the code blocks to the terminal.      

Credits to all developers and users who built, modified, refined, debugged, fixed, and used those scripts. 


### Naming rules (if present...)

- prefix
  `nrs`: for NekRS      
  `nek`: for Nek5000    

- postfix   
  `buildonly`: used to submit the build-only job         
  `ngpu`: instead of running a full node, this enables user to run the target #gpu     
  `v23`: version v23    
  `nn`: neknek
  `pf`: profiling

nrsqsub_<machine>    
nekqsub_<machine>


### Common env flags (if present)

```
PROJ_ID
QUEUE
NEKRS_HOME

NEKRS_BUILD_ONLY
NEKRS_SKIP_BUILD_ONLY
NEKRS_CIMODE
NEKRS_DFLOAT_FP32
```


### Directories

- OLCF      
  Summit          
  Frontier        

- ALCF      
  Theta        
  ThetaGPU        
  Polaris      
  Sunspot (NDA)         
  Aurora

- LCRC      
  bebop     
  swing     

- JLSE      
  nurburg

- ANL/GCE

- NERSC     
  perlmutter         

- NCSA   
  delta        

- UIUC   
  Campus Cluster     
  monza     


      


# old alias

This is not my first time trying to put the scripts together.     
In case I haven't uploaded the script, you might find some old documentations through these links.

[www on OLCF](https://users.nccs.gov/~ylan/machines_list/)        
[private git](https://github.com/misunmin/ceed-yuhsiang/tree/master/2022_Summer/scales_crusher/scripts)

