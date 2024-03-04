#!/bin/bash

#--------------------------------------
# Default parameters
: ${PROJ_ID:="CSC249ADCD04_CNDA"}
: ${QUEUE:="EarlyAppAccess"}
: ${NEKRS_HOME:="/lus/gecko/projects/CSC249ADCD04_CNDA/ylan/.local/nekrs_stefan_v24_022224"}
: ${NEKRS_SKIP_BUILD_ONLY:=0}
: ${NEKRS_BUILD_ONLY:=0}
: ${NEKRS_DFLOAT_FP32:=0}
: ${NEKRS_GPU_MPI:=1}
: ${NEKRS_BACKEND:="dpcpp"}
: ${OCCA_DPCPP_COMPILER_FLAGS:="-w -O3 -fsycl -gline-tables-only -fsycl-targets=intel_gpu_pvc -Xsycl-target-backend '-options -ze-intel-enable-auto-large-GRF-mode'"}

: ${GPUS_PER_NODE:=12}
: ${ONEAPI_MODULE:="oneapi/eng-compiler/2023.05.15.006"}
#: ${MPICH_MODULE:="mpich/52.2"}
#--------------------------------------
# Validate input arguments
if [ $# -ne 3 ]; then
  echo "usage: [PROJ_ID] [QUEUE] $0 <casename> <number of compute nodes> <hh:mm:ss>"
  exit 0
fi

if [ -z "$PROJ_ID" ]; then
  echo "ERROR: PROJ_ID is empty"
  exit 1
fi

if [ -z "$QUEUE" ]; then
  echo "ERROR: QUEUE is empty"
  exit 1
fi

bin=${NEKRS_HOME}/bin/nekrs
if [ $NEKRS_DFLOAT_FP32 -eq 1 ]; then
  bin=${NEKRS_HOME}/bin/nekrs-fp32
fi
case=$1
nodes=$2
cores_per_numa=16
let nn=$nodes*$GPUS_PER_NODE
let ntasks=nn
time=$3

qnodes=$nodes
jobname="nekRS_"$case
if [ $NEKRS_BUILD_ONLY -eq 1 ]; then
  NEKRS_SKIP_BUILD_ONLY=0
  qnodes=1
  jobname="nekRS_build_"$case
fi

time_fmt=`echo $time|tr ":" " "|awk '{print NF}'`
if [ "$time_fmt" -ne "3" ]; then
  echo "Warning: time is not in the format <hh:mm:ss>"
  echo $time
  exit 1
fi

if [ ! -f $bin ]; then
  echo "Cannot find" $bin
  exit 1
fi

if [ ! -f $case.par ]; then
  echo "Cannot find" $case.par
  exit 1
fi

if [ ! -f $case.udf ]; then
  echo "Cannot find" $case.udf
  exit 1
fi

if [ ! -f $case.re2 ]; then
  echo "Cannot find" $case.re2
  exit 1
fi

#--------------------------------------
# Generate the submission script
SFILE=s.bin
echo "#!/bin/bash" > $SFILE
echo "#PBS -A $PROJ_ID" >>$SFILE
echo "#PBS -N $jobname" >>$SFILE
echo "#PBS -l walltime=$time" >>$SFILE
echo "#PBS -l select=$qnodes" >>$SFILE
echo "#PBS -l place=scatter" >>$SFILE
echo "#PBS -k doe" >>$SFILE
echo "#PBS -j oe" >>$SFILE

# job to "run" from your submission directory
echo "cd \$PBS_O_WORKDIR" >> $SFILE

echo "echo Jobid: \$PBS_JOBID" >>$SFILE
echo "echo Running on host \`hostname\`" >>$SFILE
echo "echo Running on nodes \`cat \$PBS_NODEFILE\`" >>$SFILE

echo "module use /soft/modulefiles" >> $SFILE
echo "module load spack-pe-oneapi cmake" >> $SFILE
echo "module load $ONEAPI_MODULE" >> $SFILE
#echo "module unload mpich/51.2/icc-all-pmix-gpu" >> $SFILE
#echo "module load $MPICH_MODULE" >> $SFILE

echo "module list" >> $SFILE
echo "type -a mpiexec" >> $SFILE
echo "which mpicc" >> $SFILE

echo "export NEKRS_HOME=$NEKRS_HOME" >>$SFILE
echo "export NEKRS_GPU_MPI=$NEKRS_GPU_MPI" >>$SFILE
echo "export MPICH_GPU_SUPPORT_ENABLED=$NEKRS_GPU_MPI" >> $SFILE

# https://github.com/Nek5000/Nek5000/issues/759
echo "export FI_CXI_RX_MATCH_MODE=hybrid" >> $SFILE 

echo "export TZ='/usr/share/zoneinfo/US/Central'" >> $SFILE
echo "export LANG=en" >> $SFILE
echo "export OMP_PROC_BIND=spread" >> $SFILE
echo "export OMP_NUM_THREADS=1" >> $SFILE
echo "unset OMP_PLACES" >> $SFILE

echo "export PARRSB_VERBOSE_LEVEL=1" >> $SFILE
echo "export PARRSB_LEVELS=3" >> $SFILE

echo "export OCCA_DPCPP_COMPILER_FLAGS=\"$OCCA_DPCPP_COMPILER_FLAGS\"" >> $SFILE
echo "export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1" >> $SFILE

#if [ -d .cache ]; then
#  rm -r .cache
#fi

if [ $NEKRS_SKIP_BUILD_ONLY -eq 0 ]; then
echo "# precompile"  >>$SFILE
echo "date"  >>$SFILE
echo "mpiexec -n 1 $bin --setup ${case} --backend ${NEKRS_BACKEND} --device-id 0 --build-only $nn" >>$SFILE
echo "date"  >>$SFILE
fi

CMD="/soft/tools/mpi_wrapper_utils/gpu_tile_compact.sh"

if [ $NEKRS_BUILD_ONLY -eq 0 ]; then
echo "mpiexec -n $nn -ppn $GPUS_PER_NODE -d $cores_per_numa --cpu-bind depth $CMD $bin --setup ${case} --backend ${NEKRS_BACKEND} --device-id 0" >>$SFILE
fi

qsub -q $QUEUE $SFILE

# clean-up
# rm -rf $SFILE 
