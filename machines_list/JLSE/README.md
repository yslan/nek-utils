# JLSE

## Nurburg

2 Intel Skylake Xeon 22c 2.1Ghz CPU   
4 Tesla V100 fully connnected by NVlink   


script `nrsmpi/nrsbmpi` works

- version: v26.0
- last update: 06/17/26

- env   
  ```
  module purge
  module use /soft/modulefiles
  module load gcc/14.1.0
  module load cuda/12.9.1
  module load openmpi/4.1.1-gcc
  module load cmake/3.28.3
  module list

  ulimit -s unlimited
  ```
  CUDA 13.x dropped offline compilation/library support for Volta. So we drop cuda to 12.9.1

- config
  ```
  CC=mpicc CXX=mpic++ FC=mpif77 ./build.sh \
     -DCMAKE_INSTALL_PREFIX=$NEKRS_HOME \
     -DENABLE_CPPTRACE=yes

  ```
