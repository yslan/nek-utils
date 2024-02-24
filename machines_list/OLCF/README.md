# OLCF

Useful links      
- ![OLCF doc](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html)      
- ![OLCF status](https://www.olcf.ornl.gov/for-users/center-status/)         
- ![OLCF announcments](https://www.olcf.ornl.gov/center-announcement/)    

## Frontier

script `nrsqsub_frontier`

- version: v24.0.1
- last update: 02/21/24
- notable changes (diff from source code)
  - Use `PrgEnv-amd` for a better compatibility of HIP libraries        
    For example, `HYPRE_GPU`, `Ascent_HIP`
  - Fix module load
  - Add `NEKRS_DFLOAT_FP32` to use `nekrs-fp32`
  - Add `NEKRS_BUILD_ONLY` to submit a single node job for precompilation
  - Add `NEKRS_CIMODE` mode to pass `--cimode <cimode>`
  - WIP: detect neknek's `.sess` file     

- env
  ```
  module reset

  module load PrgEnv-amd
  module load craype-accel-amd-gfx90a
  module load cray-mpich
  module load rocm
  module load cmake
  module unload cray-libsci
  module unload darshan-runtime
  
  module list

  #Currently Loaded Modules:
  #  1) craype-x86-trento                       9) cray-dsmml/0.2.2
  #  2) libfabric/1.15.2.0                     10) PrgEnv-amd/8.3.3
  #  3) craype-network-ofi                     11) hsi/default
  #  4) perftools-base/22.12.0                 12) DefApps/default
  #  5) xpmem/2.6.2-2.5_2.22__gd067c3f.shasta  13) craype-accel-amd-gfx90a
  #  6) cray-pmi/6.1.8                         14) cray-mpich/8.1.23
  #  7) amd/5.3.0                              15) rocm/5.3.0
  #  8) craype/2.7.19                          16) cmake/3.23.2
  ```

- TODO: 
  - add precompialtion for neknek


## Summit (PLUS)

![OS update](https://docs.olcf.ornl.gov/software/software-news.html)

scripts `nrsqsub_summit`

- version: v23 or v23.1.1 (7237c89d)
- last update: 02/24/24
- notable changes (diff from source code)
  - Update module list due to OS udpate. Now, it's version specific.
  - Add `NEKRS_DFLOAT_FP32` to use `nekrs-fp32`
  - Add `NEKRS_BUILD_ONLY` to submit a single node job for precompilation
  - Add `QUEUE` to easily switch between `batch` and `debug`
  - Add time format check
  - keep tmp $SFILE to keep track last submission. (one can also use `bjobs -l`/`bjobs -d`)

- env
  ```
  module load DefApps-2023
  module load gcc/9.1.0 cmake/3.23.1 cuda/11.0.3
  module unload darshan-runtime

  module list

  #Currently Loaded Modules:
  #  1) lsf-tools/2.0   5) gcc/9.1.0                        9) nsight-systems/2021.3.1.54
  #  2) hsi/5.0.2.p5    6) spectrum-mpi/10.4.0.3-20210112  10) cuda/11.0.3
  #  3) xalt/1.2.1      7) cmake/3.23.1
  #  4) DefApps-2023    8) nsight-compute/2021.2.1
  ```

- note
  - IBM mpi only (officially) supports up to cuda 11.0.3 
    ![MPI Aware CUDA](https://docs.olcf.ornl.gov/systems/summit_user_guide.html#unsupported-cuda-versions-do-not-work-with-gpu-aware-mpi)

  - cuda 11 only support at most gcc 9

