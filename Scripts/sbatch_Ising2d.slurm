#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH -p large
#SBATCH -o %j_L-64.out
#SBATCH -J Ising2d_L-64

JOB=$SLURM_JOBID
Lsize=64
NPROCS=256
MODEL="Ising2d"
SEARCHSTR="constexpr size_t SYSTEM_SIZE_L = "


SRCDIR=/home/fernand7/meese022/REWL_Simulator
BLDDIR="build_${MODEL}_L-${Lsize}"
LOGDIR="${SRCDIR}/logs"
PARAMFILE="src/Hamiltonians/${MODEL}/${MODEL,,}_parameters.cxx"
REPLSTR="constexpr size_t SYSTEM_SIZE_L = ${Lsize};"

cd $SRCDIR; echo; pwd; echo

cd $BLDDIR; pwd

module purge
module load cmake/3.16.2 ompi/4.1.1/gnu-9.2.0 gcc/9.2.0
module list

echo -e "\n\nRunning code...\n\n"

time mpiexec -n $NPROCS --use-hwthread-cpus ./bin/REWL_Simulator > "${LOGDIR}/${MODEL}_L-${Lsize}_${JOB}.rewl"
