#!/bin/bash

Lsize=96
MODEL="Ising2d"
SEARCHSTR="constexpr size_t SYSTEM_SIZE_L = "


SRCDIR=/home/fernand7/meese022/REWL_Simulator
BLDDIR="build_${MODEL}_L-${Lsize}"
LOGDIR="${SRCDIR}/logs"
PARAMFILE="src/Hamiltonians/${MODEL}/${MODEL,,}_parameters.cxx"
REPLSTR="constexpr size_t SYSTEM_SIZE_L = ${Lsize};"

cd $SRCDIR; echo; pwd; echo

sed -i "/${SEARCHSTR}/c\\${REPLSTR}" $PARAMFILE

if [ ! -d $BLDDIR ]
then
	mkdir $BLDDIR
else
	rm -rf "${BLDDIR}"/*
fi


# Replace the SEARCHSTR with the REPLSTR
# sed -i "/${SEARCHSTR}/c\\${REPLSTR}" $PARAMFILE

if [ $? != 0 ]
then
	echo -e "\nFailure to change System Size. Exiting.\n"
	exit 1
fi

#echo -e "\nModified the file. Now building and recompiling...\n"
cd $BLDDIR; pwd

module purge
module load cmake/3.16.2 ompi gcc/9.2.0
module list

cmake .. "-D${MODEL^^}"=ON

if [ $? == 0 ]
then
	echo -e "\n\nCompiling code\n\n"

	make -j 
else
	echo -e "\nFailure in CMake.\n"
fi

