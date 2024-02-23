# ALCF

WIP

Useful links      
- ![ALCF doc](https://docs.alcf.anl.gov/aurora/getting-started-on-aurora/)      
- ![ALCF status](https://www.alcf.anl.gov/support-center/machine-status)         
- ![ALCF announcments](https://www.alcf.anl.gov/support-center/facility-updates)    

## Aurora

script `nrsqsub_aurora`

- version: v24.0.1
- last update: 02/21/24
- notable changes (diff from source code)
  - Use `PrgEnv-amd` for a better compatibility of HIP libraries        
    For example, `HYPRE_GPU`, `Ascent_HIP`
  - Fix s.bin, don't remove cache everytime
  - adjust module load
  - Add `NEKRS_DFLOAT_FP32` to use `nekrs-fp32`
  - Add `NEKRS_BUILD_ONLY` to submit a single node job for precompilation

- env
  ```
  module use /soft/modulefiles
  module load spack-pe-oneapi cmake
  module load oneapi/eng-compiler/2023.05.15.006
  module unload mpich/51.2/icc-all-pmix-gpu
  module load mpich/52.2
  
  module list

  #Currently Loaded Modules:
  #  1) gcc/11.2.0             6) spack-pe-oneapi/0.5-rc1
  #  2) libfabric/1.15.2.0     7) cmake/3.26.4-gcc-testing
  #  3) cray-pals/1.3.3        8) intel_compute_runtime/release/agama-devel-647
  #  4) cray-libpals/1.3.3     9) oneapi/eng-compiler/2023.05.15.006
  #  5) spack-pe-gcc/0.5-rc1  10) mpich/52.2/icc-all-pmix-gpu
  ```

- TODO: 

### Polaris

