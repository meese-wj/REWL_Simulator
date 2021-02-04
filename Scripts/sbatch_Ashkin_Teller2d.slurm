#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH -p small
#SBATCH -o RFBM_K-0.25_h-0.out
#SBATCH -J K-0.25_h-0

JOB=$SLURM_ARRAY_JOB_ID
TASKID=0
hsize=0
Ksize=0.25
Lsize=40
NPROCS=160
TIMELIMIT="30m"
KILLAFTER="31m"
MODEL="Ashkin_Teller2d_RFBM"
COUPLINGS="K-${Ksize}_h-${hsize}"
SEARCHSTR="constexpr size_t SYSTEM_SIZE_L = "


SRCDIR=/home/fernand7/meese022/REWL_Simulator
BLDDIR="build_${MODEL}_${COUPLINGS}"
LOGDIR="${SRCDIR}/logs"
PARAMFILE="src/Hamiltonians/${MODEL}/${MODEL,,}_parameters.cxx"
REPLSTR="constexpr size_t SYSTEM_SIZE_L = ${Lsize};"

cd $SRCDIR; echo; pwd; echo

cd $BLDDIR; pwd

module purge
module load cmake/3.16.2 ompi gcc/9.2.0
module list

echo -e "\n\nRunning code...\n\n"

timeout -k $KILLAFTER $TIMELIMIT time mpiexec -n $NPROCS --use-hwthread-cpus ./bin/REWL_Simulator $TASKID > "${LOGDIR}/${MODEL}_L-${Lsize}_${COUPLINGS}-${TASKID}.rewl"

