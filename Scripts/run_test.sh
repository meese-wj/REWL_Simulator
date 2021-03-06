#!/bin/bash

model="Ising2d"
RANDFIELD="RFIM"

cd /home/joe/Linux_Code_Dev/REWL_Simulator/build
echo
pwd
echo
ls
echo
echo "Building and running code."
cmake .. "-D${model^^}=ON" "-D${RANDFIELD^^}=OFF" "-DSIMULATED_ANNEALING=OFF" "-DAT_DENSITIES=OFF" "-DJOB_ARRAYS=OFF" "-DONE_OVER_T_ALGORITHM=OFF"

# *******************************
# Set MPI variables from CMake
DELIM="="
CACHE_LOC="CMakeCache.txt"

MPIEXEC=$( grep "MPIEXEC_EXECUTABLE:" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_NUMPROC_FLAG=$( grep "MPIEXEC_NUMPROC_FLAG:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_MAX_NUMPROCS=$( grep "MPIEXEC_MAX_NUMPROCS:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_MAX_NUMPROCS=8
MPIEXEC_PREFLAGS=$( grep "MPIEXEC_PREFLAGS:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_PREFLAGS="${MPIEXEC_PREFLAGS} --use-hwthread-cpus"
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

