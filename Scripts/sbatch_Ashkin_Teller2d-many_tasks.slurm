#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=small
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=2gb
#SBATCH --exclusive
#SBATCH --time=56:25:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --account=fernand7
#SBATCH -o K-0_h-0.5.out
#SBATCH --job-name=K-0_h-0.5

JOB=$SLURM_JOB_ID
hsize=0.5
Ksize=0
Lsize=40
NPROCS=$SLURM_NTASKS
JOBSTART=96
JOBEND=256
NJOBS=160
NSLOTS=$SLURM_NTASKS_PER_NODE
TIMELIMIT="20m"
KILLAFTER="21m"
MODEL="Ashkin_Teller2d_RFBM"
COUPLINGS="L-${Lsize}_K-${Ksize}_h-${hsize}"
SEARCHSTR="constexpr size_t SYSTEM_SIZE_L = "


SRCDIR=/home/fernand7/meese022/REWL_Simulator
#BLDDIR="build_${MODEL}_L-${Lsize}"
BLDDIR="build_${MODEL}_${COUPLINGS}"
LOGDIR="${SRCDIR}/logs"
PARAMFILE="src/Hamiltonians/${MODEL}/${MODEL,,}_parameters.cxx"
REPLSTR="constexpr size_t SYSTEM_SIZE_L = ${Lsize};"

cd $SRCDIR; echo; pwd; echo

MYHOSTS="${SRCDIR}/${JOB}.hostfile"
HOSTNAMES=$(scontrol show hostnames)
echo $HOSTNAMES
while IFS= read -r HOSTS;
do
	for name in $HOSTS;
	do
		echo "host ${name}"
		echo -e "${name} slots=${NSLOTS}" >> $MYHOSTS
	done
done <<< $HOSTNAMES 

echo -e "\nHostfile generated:\n"
cat $MYHOSTS

cd $BLDDIR; pwd

module purge
module load cmake/3.16.2 ompi/4.1.1/gnu-9.2.0 gcc/9.2.0
module list

echo -e "\n\nRunning code...\n\n"

for ((jobid=$JOBSTART; jobid<${JOBEND}; jobid++));
do
	TASKID=$jobid
	echo -e "\nStarting ${TASKID}...\n\n"
	timeout -k $KILLAFTER $TIMELIMIT time mpiexec -np $NPROCS --hostfile $MYHOSTS ./bin/REWL_Simulator $TASKID > "${LOGDIR}/${MODEL}_L-${Lsize}_${COUPLINGS}-${TASKID}.rewl"
	if [ $? == 124 ] || [ $? == 137 ]
	then
		echo -e "\n${TASKID} timed out.\n"
	fi
	echo -e "\nEnd of ${TASKID}.\n\n"
done

echo -e "\nEnd of all $NJOBS runs.\n"
