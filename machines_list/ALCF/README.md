# ALCF

Useful links      
- ![ALCF doc](https://docs.alcf.anl.gov/aurora/getting-started-on-aurora/)      
- ![ALCF status](https://www.alcf.anl.gov/support-center/machine-status)         
- ![ALCF announcments](https://www.alcf.anl.gov/support-center/facility-updates)    


## Aurora

02/25/24: This is outdated. See nrsqsub_aurora script in v24.

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

- (TODO) version: v23.1 (repo/next: commit 11/02/23)
- version: v23 (repo/next: commit 11/02/23)
- last update: 04/25/24
- notable changes    

  - Add `NEKRS_DFLOAT_FP32` and `NEKRS_BUILD_ONLY`
  - Add `FI_CXI_RX_MATCH_MODE=hybrid`
  - `LD_PRELOAD=/opt/cray/pe/gcc/12.2.0/snos/lib64/libstdc++.so.6`
  - print nodelist
  

- env    
  ```
  module use /soft/modulefiles
  module use /opt/cray/pe/lmod/modulefiles/mix_compilers
  
  module use /soft/modulefiles
  module use /opt/cray/pe/lmod/modulefiles/mix_compilers
  module load libfabric
  module load PrgEnv-gnu
  module load nvhpc-mixed
  module load craype-x86-milan craype-accel-nvidia80
  module load spack-pe-base cmake
  
  module list
  export MPICH_GPU_SUPPORT_ENABLED=1

  #Currently Loaded Modules:
  #  1) craype-network-ofi       6) craype/2.7.30      11) cray-libpals/1.3.4      16) curl/8.4.0-2ztev25
  #  2) perftools-base/23.12.0   7) cray-dsmml/0.2.2   12) PrgEnv-gnu/8.5.0        17) cmake/3.27.7
  #  3) darshan/3.4.4            8) cray-mpich/8.1.28  13) nvhpc-mixed/23.9        18) craype-x86-milan
  #  4) libfabric/1.15.2.0       9) cray-pmi/6.1.13    14) spack-pe-base/0.6.1     19) craype-accel-nvidia80
  #  5) gcc-native/12.3         10) cray-pals/1.3.4    15) nghttp2/1.57.0-ciat5hu
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
  LD_PRELOAD=/opt/cray/pe/gcc/12.2.0/snos/lib64/libstdc++.so.6
  ``` 
  See update at here: ![Polaris known issue](https://docs.alcf.anl.gov/polaris/known-issues/)

- note hypre
  HYPRE in v23 doesn't support cuda 12, Turn `-DENABLE_HYPRE_GPU=off` (and `DENABLE_AMGX=off`) for now





