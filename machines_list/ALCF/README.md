# ALCF

Useful links      
- ![ALCF doc](https://docs.alcf.anl.gov/aurora/getting-started-on-aurora/)      
- ![ALCF status](https://www.alcf.anl.gov/support-center/machine-status)         
- ![ALCF announcments](https://www.alcf.anl.gov/support-center/facility-updates)    


## Aurora

script `nrsqsub_aurora`

- version: v24.0.1
- last update: 02/21/24
- notable changes (diff from source code)
  - massively clean up
  - Fix s.bin, don't remove cache everytime
  - adjust module load
  - Add `NEKRS_DFLOAT_FP32` to use `nekrs-fp32`
  - Add `NEKRS_BUILD_ONLY` to submit a single node job for precompilation

- env
  ```
  module use /soft/modulefiles
  module load spack-pe-oneapi cmake
  module load oneapi/eng-compiler/2023.05.15.006
  module list

  #Currently Loaded Modules:
  #  1) gcc/11.2.0                    6) spack-pe-gcc/0.5-rc1
  #  2) mpich/51.2/icc-all-pmix-gpu   7) spack-pe-oneapi/0.5-rc1
  #  3) libfabric/1.15.2.0            8) cmake/3.26.4-gcc-testing
  #  4) cray-pals/1.3.3               9) intel_compute_runtime/release/agama-devel-647
  #  5) cray-libpals/1.3.3           10) oneapi/eng-compiler/2023.05.15.006
  ```

- Compile
  ```
  CC=mpicc CXX=mpic++ FC=mpif77 ./build.sh -DENABLE_HYPRE_GPU=off
  ```

- note
  - mpich/52.2 fails to run 2 nodes.

- TODO: 

---

### Polaris

scripts `nrsqsub_polaris`

- version: v23.1 (repo/next: commit 11/02/23)
- last update: 03/09/24
- notable changes    

  - Add `NEKRS_DFLOAT_FP32` and `NEKRS_BUILD_ONLY`
  - Add `FI_CXI_RX_MATCH_MODE=hybrid`
  - `LD_PRELOAD=/opt/cray/pe/gcc/11.2.0/snos/lib64/libstdc++.so.6`
  - print nodelist
  

- env    
  ```
  module use /soft/modulefiles
  module use /opt/cray/pe/lmod/modulefiles/mix_compilers
  
  module load libfabric
  module load cpe-cuda
  module load PrgEnv-gnu
  #module load gcc/11.2.0 # manually load this if you find PrgEnv-gnu load gcc-12
  module load nvhpc-mixed
  module load cmake
  
  module list
  export MPICH_GPU_SUPPORT_ENABLED=1

  # Currently Loaded Modules:
  #   1) craype-x86-rome          7) cpe-cuda/23.03     13) cray-libpals/1.2.11
  #   2) craype-network-ofi       8) craype/2.7.20      14) PrgEnv-gnu/8.3.3
  #   3) craype-accel-nvidia80    9) cray-dsmml/0.2.2   15) nvhpc-mixed/22.11
  #   4) libfabric/1.15.2.0      10) cray-mpich/8.1.25  16) cmake/3.23.2
  #   5) gcc/11.2.0              11) cray-pmi/6.1.10
  #   6) perftools-base/23.03.0  12) cray-pals/1.2.11
  ```

- Compile
  ```
  CC=cc CXX=CC FC=ftn ./nrsconfig \
   -DCMAKE_INSTALL_PREFIX=<your_install_path> 

  # optional: 
  # -DENABLE_AMGX=on -DENABLE_HYPRE_GPU=on
  ```

- note:      
  Since ss11 upgrade, default compile was changed to gcc-12 which messes up with cuda version.  
  Recently, they switch back to gcc-11.2 so the module list above works.
  However, there is still a known issues that libc++ is not under the right path.   
  Currently, it can be bypassed by
  ```
  LD_PRELOAD=/opt/cray/pe/gcc/11.2.0/snos/lib64/libstdc++.so.6
  ``` 
  See update at here: ![Polaris known issue](https://docs.alcf.anl.gov/polaris/known-issues/)










