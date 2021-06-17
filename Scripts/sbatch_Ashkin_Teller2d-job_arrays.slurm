#!/bin/bash
#SBATCH --nodes=24
#SBATCH --partition=large
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=576
#SBATCH --mem-per-cpu=2gb
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH --account=fernand7
#SBATCH -o L-144_K-0.5.%A_%a.out
#SBATCH --job-name=L-144_K-0.5-%a

# Job array script
JOB=$SLURM_JOB_ID
TASKID=$SLURM_ARRAY_TASK_ID
hsize=0.
Ksize=0.5
Lsize=144
NPROCS=$SLURM_NTASKS
NSLOTS=$SLURM_NTASKS_PER_NODE
MODEL="Ashkin_Teller2d"
COUPLINGS="L-${Lsize}_K-${Ksize}"

SRCDIR=/home/fernand7/meese022/REWL_Simulator
#BLDDIR="build_${MODEL}_L-${Lsize}"
BLDDIR="build_${MODEL}_${COUPLINGS}"
LOGDIR="${SRCDIR}/logs"
PARAMFILE="src/Hamiltonians/${MODEL}/${MODEL,,}_parameters.cxx"
REPLSTR="constexpr size_t SYSTEM_SIZE_L = ${Lsize};"

cd $SRCDIR; echo; pwd; echo

MYHOSTS="${SRCDIR}/${JOB}-${TASKID}.hostfile"
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

echo -e "\nStarting ${TASKID}...\n\n"
timeout -k $KILLAFTER $TIMELIMIT time mpiexec -np $NPROCS --hostfile $MYHOSTS ./bin/REWL_Simulator $TASKID > "${LOGDIR}/${MODEL}_${COUPLINGS}-${TASKID}.rewl"

