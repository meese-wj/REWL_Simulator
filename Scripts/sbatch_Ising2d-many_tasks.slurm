#!/bin/bash
#SBATCH --nodes=16
#SBATCH --partition=large
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=256
#SBATCH --mem-per-cpu=2gb
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --account=fernand7
#SBATCH -o RFIM_h-0.9.out
#SBATCH --job-name=h-0.9

JOB=$SLURM_JOB_ID
hsize=0.9
Lsize=64
NPROCS=$SLURM_NTASKS
NJOBS=192
NSLOTS=$SLURM_NTASKS_PER_NODE
TIMELIMIT="20m"
KILLAFTER="21m"
MODEL="Ising2d_RFIM"
COUPLINGS="h-${hsize}"

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
module load cmake/3.16.2 ompi gcc/9.2.0
module list

echo -e "\n\nRunning code...\n\n"

for ((jobid = 0; jobid<${NJOBS}; jobid++));
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

