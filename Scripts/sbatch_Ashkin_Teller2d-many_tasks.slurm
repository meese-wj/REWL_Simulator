#!/bin/bash
#SBATCH --nodes=24
#SBATCH --partition=large
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=576
#SBATCH --mem-per-cpu=2gb
#SBATCH --exclusive
#SBATCH --time=13:25:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --account=fernand7
#SBATCH -o L-120_K-0.5.out
#SBATCH --job-name=L-120_K-0.5

JOB=$SLURM_JOB_ID
hsize=0.
Ksize=0.5
Lsize=120
NPROCS=$SLURM_NTASKS
JOBSTART=0
JOBEND=32
NJOBS=32
NSLOTS=$SLURM_NTASKS_PER_NODE
TIMELIMIT="30m"
KILLAFTER="31m"
MODEL="Ashkin_Teller2d"
COUPLINGS="L-${Lsize}_K-${Ksize}"
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
module load cmake/3.16.2 ompi gcc/9.2.0
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
