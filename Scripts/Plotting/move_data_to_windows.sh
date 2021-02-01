#!/bin/bash

while getopts ":m:d:c:v:" arg;
do
    case $arg in
        m) MODEL=$OPTARG;;
        d) DATE=$OPTARG;;
        c) COUPLING=$OPTARG;;
        v) VALUE=$OPTARG;;
    esac
done

# Check if the model name is nonempty
if [ -z "$MODEL" ]
then
    echo -e "\nA model type is required."
    echo -e "Run code with options: -m 'MODEL' -d 'DATE' -c 'COUPLING' -v 'VALUE'\n"
    exit 1
fi

# Check if the date is nonempty
if [ -z "$DATE" ]
then
    echo -e "\nA date is required."
    echo -e "Run code with options: -m 'MODEL' -d 'DATE' -c 'COUPLING' -v 'VALUE'\n"
    exit 1
fi

# Check if the coupling string  is nonempty
if [ -z "$COUPLING" ]
then
    echo -e "\nA coupling is required."
    echo -e "Run code with options: -m 'MODEL' -d 'DATE' -c 'COUPLING' -v 'VALUE'\n"
    exit 1
fi

# Check if the coupling value  is nonempty
if [ -z "$VALUE" ]
then
    echo -e "\nA coupling value is required."
    echo -e "Run code with options: -m 'MODEL' -d 'DATE' -c 'COUPLING' -v 'VALUE'\n"
    exit 1
fi

# ********************************************************************************************
# Now check that the directories are fine

# Move into Simulation_Data
SIMDATA=/home/joe/Linux_Code_Dev/REWL_Simulator/Simulation_Data
cd $SIMDATA

# Check if the model is defined
if [ ! -d "$MODEL" ]
then
    ls
    echo -e "\n${MODEL} is not a directory\n"
    exit 1
fi

# Move into $MODEL
cd $MODEL

# Check if the date is defined
if [ ! -d "$DATE" ]
then
    ls
    echo -e "\n${DATE} is not a directory\n"
    exit 1
fi

# Move into $DATE
cd $DATE

# Check if the Figures folder is present
if [ ! -d "Figures_${COUPLING}-${VALUE}" ]
then
    echo -e "\nNo Figures_${COUPLING}-${VALUE} to copy to Windows\n"
    exit 1
fi

# ********************************************************************************************
# Now copy over to Windows

WINDIR=/mnt/c/Users/meese/Documents/UMN_Research_Projects/Wang_Landau
PLOTDIR="Data_Plots"

# Build Directories
if [ ! -d "${WINDIR}/${MODEL}" ]
then
    echo -e "\nMaking new directory in $WINDIR for $MODEL\n"
    mkdir "${WINDIR}/${MODEL}"
fi

if [ ! -d "${WINDIR}/${MODEL}/${PLOTDIR}" ]
then
    echo -e "\nMaking new directory in $WINDIR/${MODEL} for $PLOTDIR\n"
    mkdir "${WINDIR}/${MODEL}/${PLOTDIR}"
fi

if [ ! -d "${WINDIR}/${MODEL}/${PLOTDIR}/${DATE}" ]
then
    echo -e "\nMaking new directory in ${WINDIR}/${MODEL}/${PLOTDIR} for $DATE\n"
    mkdir "${WINDIR}/${MODEL}/${PLOTDIR}/${DATE}"
fi

if [ ! -d "${WINDIR}/${MODEL}/${PLOTDIR}/${DATE}/${COUPLING}_${VALUE}" ]
then
    echo -e "\nMaking a new directory in ${WINDIR}/${MODEL}/${PLOTDIR}/${DATE} for ${COUPLING}_${VALUE}\n"
    mkdir "${WINDIR}/${MODEL}/${PLOTDIR}/${DATE}/${COUPLING}_${VALUE}"
fi

if [ $? == 0 ]
then
    echo -e "\nMoving data from\n${SIMDATA}/${MODEL}/${DATE}/Figures_${COUPLING}-${VALUE} to\n${WINDIR}/${MODEL}/${PLOTDIR}/${DATE}/${COUPLING}_${VALUE}"
    cp -r "${SIMDATA}/${MODEL}/${DATE}/Figures_${COUPLING}-${VALUE}"/* "${WINDIR}/${MODEL}/${PLOTDIR}/${DATE}/${COUPLING}_${VALUE}"
fi

