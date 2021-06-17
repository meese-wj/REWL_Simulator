#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=small
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=2gb
#SBATCH --exclusive
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --account=fernand7
#SBATCH -o RFIM_h-0.6.out
#SBATCH --job-name=h-0.6

JOB=$SLURM_JOB_ID
hsize=0.6
Lsize=40
NPROCS=$SLURM_NTASKS
JOBSTART=128
JOBEND=256 
NJOBS=128
NSLOTS=$SLURM_NTASKS_PER_NODE
TIMELIMIT="25m"
KILLAFTER="26m"
MODEL="Ising2d_RFIM"
COUPLINGS="L-${Lsize}_h-${hsize}"

SRCDIR=/home/fernand7/meese022/REWL_Simulator
BLDDIR="build_${MODEL}_${COUPLINGS}"
LOGDIR="${SRCDIR}/logs"

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

for ((jobid = $JOBSTART; jobid<${JOBEND}; jobid++));
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

