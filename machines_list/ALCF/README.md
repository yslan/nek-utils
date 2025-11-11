# ALCF

Useful links      
- ![ALCF doc](https://docs.alcf.anl.gov/aurora/getting-started-on-aurora/)      
- ![ALCF status](https://www.alcf.anl.gov/support-center/machine-status)         
- ![ALCF announcments](https://www.alcf.anl.gov/support-center/facility-updates)    


## Aurora (NekRS)

11/11/25: This is outdated.`See nekrs_HPCsupport/Aurora/nrsqsub` script in v25.

### Aurora (Nek5000)

scripts `n5kqsub_polaris`

- env
  ```
  module load cmake
  module list
  ```

- makenek
  ```
  # source path 
  NEK_SOURCE_ROOT="/lus/flare/EnergyApps/ylan/src/Nek5000_repo"
  
  FC="mpif77"
  CC="mpicc"
  
  PPLIST="PARRSB HYPRE"
  
  FFLAGS="-mcmodel=large"
  CFLAGS="-mcmodel=large"
  ```



---

### Polaris (NekRS)

11/11/25: This is outdated.`See nekrs_HPCsupport/Polaris/nrsqsub` script in v25.


### Polaris (Nek5000)

scripts `n5kqsub_polaris`

- env
  ```
  module use /soft/modulefiles
  module load PrgEnv-gnu
  module load spack-pe-base cmake
  module list
  ```

- makenek
  ```
  # source path 
  NEK_SOURCE_ROOT="/lus/grand/projects/CSC249ADCD04/ylan/src/Nek5000_repo_polaris"
  
  FC="ftn"
  CC="cc"
  
  PPLIST="PARRSB DPROCMAP HYPRE"
  
  FFLAGS="-mcmodel=large"
  CFLAGS="-mcmodel=large"
  ```





