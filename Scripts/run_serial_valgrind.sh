#!/bin/bash

cd /home/joe/Linux_Code_Dev/REWL_Simulator/build_serial
echo
pwd
echo
ls
echo
echo "Building and running code."
rm -rf /home/joe/Linux_Code_Dev/REWL_Simulator/build_serial/*
cmake .. -DMPI_ON=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo

make -j
if [[ $? == 0 ]]
then
    echo
    echo
    time valgrind --verbose --leak-check=yes ./bin/REWL_Simulator
else
    echo
    echo "Error in compilation. Exiting."
fi
echo

