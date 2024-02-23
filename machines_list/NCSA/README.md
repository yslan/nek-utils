# NCSA

Useful links      
- ![NCSA doc](https://docs.ncsa.illinois.edu/systems/delta/en/latest/index.html)
- ![Delta SS11 update](https://wiki.ncsa.illinois.edu/display/DSC/Delta+Network+Upgrade)

We switch to gnu + openmpi now. For cray-mpich, see ![README_cray](./README_cray.md)       
There is still an ongoing s11 performance related issues.   

## Delta (gpu)

- storage:
  ```
  # project storate (slow)
  /project/<project_id>  

  # project scratch (recommended)
  /project/<project_id>
  ```

### NekRS v23.1.1

script `nrsqsub_delta.openmpi`

- last update: 02/23/24
- version: v23.1.1 (repo/next)
- env
  ```
  module load cmake
  module load gcc-runtime/11.4.0
  module load openmpi/4.1.5+cuda

  module list

  #Currently Loaded Modules:
  #  1) gcc/11.4.0    3) cue-login-env/1.0   5) default-s11    7) gcc-runtime/11.4.0
  #  2) cuda/11.8.0   4) slurm-env/0.1       6) cmake/3.27.9   8) openmpi/4.1.5+cuda
  ```

- config
  ```
  # : ${NEKRS_INSTALL_DIR:=/scratch/bcla/ylan/.local/nekrs_next022024}
  CC=mpicc CXX=mpic++ FC=mpif77 ./nrsconfig \
     -DCMAKE_INSTALL_PREFIX="${NEKRS_INSTALL_DIR}" \
     -DENABLE_HYPRE_GPU=on \
     -DENABLE_HIP=off \
     -DENABLE_DPCPP=off
  ```

### NekRS v24

### NekRS + Ascent



