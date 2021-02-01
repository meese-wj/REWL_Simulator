#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=33:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH -p small
#SBATCH -o RFBM_K-0.1_h-0.out
#SBATCH -J K-0.1_h-0

JOB=$SLURM_ARRAY_JOB_ID
TASKID=$SLURM_ARRAY_TASK_ID
hsize=0
Ksize=0.1
Lsize=40
NPROCS=160
NJOBS=64
TIMELIMIT="30m"
KILLAFTER="31m"
MODEL="Ashkin_Teller2d_RFBM"
COUPLINGS="K-${Ksize}_h-${hsize}"
SEARCHSTR="constexpr size_t SYSTEM_SIZE_L = "


SRCDIR=/home/fernand7/meese022/REWL_Simulator
#BLDDIR="build_${MODEL}_L-${Lsize}"
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

for ((jobid = 0; jobid<${NJOBS}; jobid++));
do
	TASKID=$jobid
	echo -e "\nStarting ${TASKID}...\n\n"
	timeout -k $KILLAFTER $TIMELIMIT time mpiexec -n $NPROCS --use-hwthread-cpus ./bin/REWL_Simulator $TASKID > "${LOGDIR}/${MODEL}_L-${Lsize}_${COUPLINGS}-${TASKID}.rewl"
	if [ $? == 124 ] || [ $? == 137 ]
	then
		echo -e "\n${TASKID} timed out.\n"
	fi
	echo -e "\nEnd of ${TASKID}.\n\n"
done

echo -e "\nEnd of all $NJOBS runs.\n"
