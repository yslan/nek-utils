#!/bin/bash

: ${PROJ_ID:=""} # xxxx-delta-gpu where "xxxx" is yorur project id
: ${QUEUE:="gpuA100x4"}
: ${NEKRS_HOME:="$HOME/.local/nekrs"}

: ${NEKRS_DFLOAT_FP32:=0}
: ${NEKRS_BUILD_ONLY:=0}
: ${NEKRS_SKIP_BUILD_ONLY:=0}

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
gpu_per_node=4
cores_per_socket=16
let nn=$nodes*$gpu_per_node
let ntasks=nn
time=$3
backend=CUDA

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

if [ ! -f $case.oudf ]; then
  echo "Cannot find" $case.oudf
  exit 1
fi

if [ ! -f $case.re2 ]; then
  echo "Cannot find" $case.re2
  exit 1
fi


## romio setup
#export ROMIO_HINTS="$(pwd)/.romio_hint"
#if [ ! -f "$ROMIO_HINTS" ]; then
#  echo "romio_no_indep_rw true"   >$ROMIO_HINTS
#  echo "romio_cb_write enable"   >>$ROMIO_HINTS
#  echo "romio_ds_write enable"   >>$ROMIO_HINTS
#  echo "romio_cb_read enable"    >>$ROMIO_HINTS
#  echo "romio_ds_read enable"    >>$ROMIO_HINTS
#  echo "cb_buffer_size 16777216" >>$ROMIO_HINTS
#  echo "cb_config_list *:1"      >>$ROMIO_HINTS
#fi


# sbatch
SFILE=s.bin
echo "#!/bin/bash" > $SFILE
echo "#SBATCH --account=$PROJ_ID" >>$SFILE
echo "#SBATCH --job-name=$jobname" >>$SFILE
echo "#SBATCH -o %x-%j.out" >>$SFILE
echo "#SBATCH --time=$time" >>$SFILE
echo "#SBATCH --nodes=$qnodes" >>$SFILE
echo "#SBATCH --partition=$QUEUE" >>$SFILE
echo "#SBATCH --constraint=\"scratch\"" >>$SFILE
echo "#SBATCH --exclusive" >>$SFILE
echo "#SBATCH --mem=208G" >>$SFILE
echo "#SBATCH --ntasks-per-node=$gpu_per_node" >>$SFILE
echo "#SBATCH --cpus-per-task=16" >>$SFILE
echo "#SBATCH --gpu-bind=closest" >> $SFILE
echo "#SBATCH --gpus-per-node=$gpu_per_node" >> $SFILE

#echo "export SLURM_CPU_BIND=\"cores\"" >> $SFILE
#echo "export CRAY_ACCEL_TARGET=nvidia80" >>$SFILE
#echo "export MPICH_GPU_SUPPORT_ENABLED=1" >>$SFILE
echo "export OMP_NUM_THREADS=1" >>$SFILE

echo "ulimit -s unlimited" >>$SFILE
echo "export NEKRS_HOME=$NEKRS_HOME" >>$SFILE
echo "export NEKRS_GPU_MPI=1" >>$SFILE

echo "module load cmake" >>$SFILE
echo "module load gcc-runtime/11.4.0" >>$SFILE
echo "module load openmpi/4.1.5+cuda" >>$SFILE
echo "module list" >>$SFILE
echo "nvidia-smi" >>$SFILE
echo "nvcc --version" >>$SFILE
echo "cmake --version" >>$SFILE
echo "which mpicc" >>$SFILE
echo "ldd $bin" >>$SFILE

#echo "export ROMIO_HINTS=$ROMIO_HINTS" >>$SFILE

echo "# Workaround for https://github.com/Nek5000/Nek5000/issues/759" >> $SFILE
echo "export FI_OFI_RXM_RX_SIZE=32768" >> $SFILE

if [ $NEKRS_SKIP_BUILD_ONLY -eq 0 ]; then
echo "LANG=EN TZ=America/Chicago date" >>$SFILE
echo "srun --ntasks=1 --nodes=1 --gpus-per-node=1 --cpus-per-task=16 $bin --setup ${case} --backend ${backend} --device-id 0 --build-only $nn" >>$SFILE
echo "LANG=EN TZ=America/Chicago date" >>$SFILE
echo "if [ \$? -ne 0 ]; then" >> $SFILE
echo "  exit" >> $SFILE
echo "fi" >> $SFILE
fi

if [ $NEKRS_BUILD_ONLY -eq 0 ]; then
echo "srun --ntasks=$nn --nodes=$nodes --ntasks-per-node=$gpu_per_node --cpus-per-task=16 --gpus-per-node=$gpu_per_node $bin --setup ${case} --backend ${backend}"  >>$SFILE
fi

sbatch $SFILE

# clean-up
#rm -rf $SFILE $ROMIO_HINTS
