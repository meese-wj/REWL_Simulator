#!/bin/bash

source /home/fernand7/meese022/.bashrc

SRC_DIR="/home/fernand7/meese022/REWL_Simulator"
BUILD_DIR="${SRC_DIR}/build_ashkin_teller2d"

cd $BUILD_DIR
echo
pwd
echo
ls
echo
echo "Building and running code."
pwd
export CC=`which gcc`
export CXX=`which g++`
echo $CC
echo $CXX
cmake .. -DISING2D=OFF -DASHKIN_TELLER2D=ON

if [ $? != 0 ]
then
	echo -e "\nCMake failed. Exiting now.\n\n"
	exit 1
fi

# *******************************
# Set MPI variables from CMake
DELIM="="
CACHE_LOC="CMakeCache.txt"

MPIEXEC=$( grep "MPIEXEC_EXECUTABLE:" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_NUMPROC_FLAG=$( grep "MPIEXEC_NUMPROC_FLAG:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_NUMPROC_FLAG="-np"
MPIEXEC_MAX_NUMPROCS=$( grep "MPIEXEC_MAX_NUMPROCS:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
#MPIEXEC_MAX_NUMPROCS=8
MPIEXEC_PREFLAGS=$( grep "MPIEXEC_PREFLAGS:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_PREFLAGS="${MPIEXEC_PREFLAGS} --use-hwthread-cpus --oversubscribe --hostfile ${SRC_DIR}/hostfile"
MPIEXEC_POSTFLAGS=$( grep "MPIEXEC_POSTFLAGS:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
# *******************************

make -j
if [[ $? == 0 ]]
then
    echo
    echo
    echo "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} \
          ./bin/REWL_Simulator ${MPIEXEC_POSTFLAGS}"
    time ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} \
         ./bin/REWL_Simulator ${MPIEXEC_POSTFLAGS}
else
    echo
    echo "Error in compilation. Exiting."
fi
echo

