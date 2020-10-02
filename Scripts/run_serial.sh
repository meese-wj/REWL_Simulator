#!/bin/bash

cd /home/joe/Linux_Code_Dev/REWL_Simulator/build_serial
echo
pwd
echo
ls
echo
echo "Building and running code."
cmake .. -DMPI_ON=OFF

make -j
if [[ $? == 0 ]]
then
    echo
    echo
    time ./bin/REWL_Simulator
else
    echo
    echo "Error in compilation. Exiting."
fi
echo

