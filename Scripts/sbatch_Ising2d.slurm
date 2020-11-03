#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=0:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meese022@umn.edu
#SBATCH -p small
#SBATCH -o %j.out

Lsize=24
MODEL="Ising2d"
SEARCHSTR="constexpr size_t SYSTEM_SIZE_L = "


SRCDIR=/home/fernand7/meese022/REWL_Simulator
BLDDIR="build_${MODEL}_L-${Lsize}"
LOGDIR="${SRCDIR}/logs"
PARAMFILE="src/Hamiltonians/${MODEL}/${MODEL,,}_parameters.cxx"
REPLSTR="constexpr size_t SYSTEM_SIZE_L = ${Lsize};"

cd $SRCDIR; echo; pwd; echo

if [ ! -d $BLDDIR ]
then
	mkdir $BLDDIR
else
	rm -rf "${BLDDIR}"/*
fi

# Replace the SEARCHSTR with the REPLSTR
sed -i "/${SEARCHSTR}/c\\${REPLSTR}" $PARAMFILE

if [ $? != 0 ]
then
	echo -e "\nFailure to change System Size. Exiting.\n"
	exit 1
fi

echo -e "\nModified the file. Now building and recompiling...\n"
cd $BLDDIR; pwd

module purge
module load cmake/3.16.2 ompi gcc/9.2.0
module list

cmake .. -DISING2D=ON

if [ $? == 0 ]
then
	echo -e "\n\nCompiling code\n\n"

	make -j 
else
	echo -e "\nFailure in CMake.\n"
fi

if [ $? == 0 ]
then
	echo -e "\n\nRunning code...\n\n"

	time mpiexec -n 24 --use-hwthread-cpus ./bin/REWL_Simulator > "${LOGDIR}/${MODEL}_L-${Lsize}.rewl"
else
	echo -e "\nFailure in make\n"
fi

