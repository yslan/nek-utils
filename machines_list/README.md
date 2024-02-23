# Machines list

This is, and will always be, under construction.

My collection of the toys with the configuration/instruction/doc for running Nek5000/NekRS, including scipts, modules and shortcuts.

Disclaimer:       
- Using scripts from a stranger is dangerous. Use them at your own risk.   
- Those scripts might be out-dated. In general, it should still work with minimal adjustment.      
  Please feel free contact me for any issues or if you find new changes.    


Credits to all developers and users who built/modified/fixed/refined/debugged/used those scripts. 


### Naming rules (if present...)

- prefix
  `nrs`: for NekRS      
  `nek`: for Nek5000    

- postfix   
  `buildonly`: used to submit the build-only job         
  `ngpu`: instead of running a full node, this enables user to run the target #gpu     
  `v23`: version v23    

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

- NERSC     
  perlmutter         

- NCSA   
  delta        

- UIUC   
  Campus Cluster     
  monza     

- JLSE      
  nurburg
      


# old alias

This is not my first time trying to put the scripts together.     
In case I haven't uploaded the script, you might find some old documentations through these links.

[www on OLCF](https://users.nccs.gov/~ylan/machines_list/)
[private git](https://github.com/misunmin/ceed-yuhsiang/tree/master/2022_Summer/scales_crusher/scripts)

