## Delta (gpu) + cray-mpich

Useful links      
- ![NCSA doc](https://docs.ncsa.illinois.edu/systems/delta/en/latest/index.html)
- ![Delta SS11 update](https://wiki.ncsa.illinois.edu/display/DSC/Delta+Network+Upgrade)

### NekRS v23.1.1

script `nrsqsub_delta.cray_mpich`

- last update: 02/22/24
- version: v23.1.1 (repo/next) (7237c89d) + local changes      
  There cray-mpich wrapper links `libmpi_gtl_cuda`, but it doesn't automatically link `-lcuda -lcudart` (stored in LDFLAGS). `CC` can't compile any cxx code without linking cuda enven though it's not using it.      
  To fix that, we pass the flags to cmake and make sure nekrs's compiler tests can pick i up during the test. This, however, cnnot fix the HYPRE\_GPU error which needs more time to dig in.

- git diff (`cray_v23.1.1.patch`)    
  ```
  $ git diff config/utils.cmake
  diff --git a/config/utils.cmake b/config/utils.cmake
  index 6d5d32b1..d4d232d2 100644
  --- a/config/utils.cmake
  +++ b/config/utils.cmake
  @@ -2,8 +2,13 @@ function (__MPI_find_compiler LANG QUERY_FLAG OUTPUT_VARIABLE)
     separate_arguments(_MPI_COMPILER_WRAPPER_OPTIONS NATIVE_COMMAND "${QUERY_FLAG}")
     set(DUMMYSRC "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx")
     file(WRITE ${DUMMYSRC} "int main() { return 0; }\n")
  +
  +  string(REPLACE " " ";" TMP_LIST $ENV{LDFLAGS})
  +  set (EXECUTE_COMMAND ${MPI_${LANG}_COMPILER} ${TMP_LIST} ${_MPI_COMPILER_WRAPPER_OPTIONS} ${DUMMYSRC})
  +
     execute_process(
  -    COMMAND ${MPI_${LANG}_COMPILER} ${_MPI_COMPILER_WRAPPER_OPTIONS} ${DUMMYSRC}
  +    COMMAND ${EXECUTE_COMMAND}
  +    COMMAND_ECHO     STDOUT
       OUTPUT_VARIABLE  WRAPPER_OUTPUT OUTPUT_STRIP_TRAILING_WHITESPACE
       ERROR_VARIABLE   WRAPPER_ERR ERROR_STRIP_TRAILING_WHITESPACE
       RESULT_VARIABLE  WRAPPER_RETURN)
  ```

- env
  ```
  module unload openmpi gcc
  module load PrgEnv-gnu craype-x86-milan craype-accel-ncsa cuda
  module load cmake

  module list

  #Currently Loaded Modules:
  #  1) cue-login-env/1.0    9) cray-libsci/23.09.1.1
  #  2) slurm-env/0.1       10) PrgEnv-gnu/8.4.0
  #  3) default-s11         11) craype-x86-milan
  #  4) gcc-native/11.2     12) craype-accel-ncsa/1.0
  #  5) libfabric/1.15.2.0  13) cuda/11.8.0
  #  6) craype-network-ofi  14) gcc-runtime/8.5.0
  #  7) craype/2.7.23       15) cmake/3.27.9
  #  8) cray-mpich/8.1.27
  ```

- config
  ```
  : ${NEKRS_INSTALL_DIR:=/scratch/bcla/ylan/.local/nekrs_next022024}
  CC=cc CXX=CC FC=ftn ./nrsconfig \
     -DCMAKE_INSTALL_PREFIX="${NEKRS_INSTALL_DIR}" \
     -DENABLE_HYPRE_GPU=on \
     -DENABLE_HIP=off \
     -DENABLE_DPCPP=off
  ```

### NekRS v24

### NekRS + Ascent



