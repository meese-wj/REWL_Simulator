#!/bin/bash

model="Ashkin_Teller2d"
RANDFIELD="RFAT_Baxter"

cd /home/joe/Linux_Code_Dev/REWL_Simulator/build
echo
pwd
echo
ls
echo
echo "Building and running code."
rm -rf /home/joe/Linux_Code_Dev/REWL_Simulator/build/*
cmake .. -DCMAKE_BUILD_TYPE=Debug "-D${model^^}=ON" "-D${RANDFIELD^^}=OFF" "-DSIMULATED_ANNEALING=OFF"

# *******************************
# Set MPI variables from CMake
DELIM="="
CACHE_LOC="CMakeCache.txt"

MPIEXEC=$( grep "MPIEXEC_EXECUTABLE:" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_NUMPROC_FLAG=$( grep "MPIEXEC_NUMPROC_FLAG:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_MAX_NUMPROCS=$( grep "MPIEXEC_MAX_NUMPROCS:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_MAX_NUMPROCS=4
MPIEXEC_PREFLAGS=$( grep "MPIEXEC_PREFLAGS:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
MPIEXEC_POSTFLAGS=$( grep "MPIEXEC_POSTFLAGS:STRING" $CACHE_LOC | cut -d $DELIM -f2 )
# *******************************

$MPIEXEC --version

make -j
if [[ $? == 0 ]]
then
    echo
    echo
    rm nc.vg*
    echo "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS}" \
         "valgrind --verbose --leak-check=yes --show-reachable=yes --log-file=nc.vg%p ./bin/REWL_Simulator ${MPIEXEC_POSTFLAGS}"
    time ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} \
         valgrind --verbose --leak-check=yes --show-reachable=yes --track-origins=yes --log-file=nc.vg%p ./bin/REWL_Simulator ${MPIEXEC_POSTFLAGS}
else
    echo
    echo "Error in compilation. Exiting."
fi
echo

