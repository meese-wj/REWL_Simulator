#!/bin/bash

while getopts ":m:d:c:v:T:" arg;
do
    case $arg in
        m) MODEL=$OPTARG;;
        d) DATE=$OPTARG;;
        c) COUPLING=$OPTARG;;
        v) VALUE=$OPTARG;;
        T) TC_VALUE=$OPTARG;;  # Optional Tc for vertical line
    esac
done

# Check if the model name is nonempty
if [ -z "$MODEL" ]
then
    echo -e "\nA model type is required."
    echo -e "Run code with options: -m 'MODEL' -d 'DATE' -c 'COUPLING' -v 'VALUE' [OPTIONAL] -T 'TC_VALUE'\n"
    exit 1
fi

# Check if the date is nonempty
if [ -z "$DATE" ]
then
    echo -e "\nA date is required."
    echo -e "Run code with options: -m 'MODEL' -d 'DATE' -c 'COUPLING' -v 'VALUE' [OPTIONAL] -T 'TC_VALUE'\n"
    exit 1
fi

# Check if the coupling string  is nonempty
if [ -z "$COUPLING" ]
then
    echo -e "\nA coupling is required."
    echo -e "Run code with options: -m 'MODEL' -d 'DATE' -c 'COUPLING' -v 'VALUE' [OPTIONAL] -T 'TC_VALUE'\n"
    exit 1
fi

# Check if the coupling value  is nonempty
if [ -z "$VALUE" ]
then
    echo -e "\nA coupling value is required."
    echo -e "Run code with options: -m 'MODEL' -d 'DATE' -c 'COUPLING' -v 'VALUE' [OPTIONAL] -T 'TC_VALUE'\n"
    exit 1
fi

# ********************************************************************************************
# Now check that the directories are fine

# Move into Simulation_Data
cd /home/joe/Linux_Code_Dev/REWL_Simulator/Simulation_Data

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

# ********************************************************************************************
# Now plot Ashkin-Teller data

PLOTTER=/home/joe/Linux_Code_Dev/REWL_Simulator/Scripts/Plotting/parse_and_plot_observables.py
SCALING=/home/joe/Linux_Code_Dev/REWL_Simulator/Scripts/Plotting/parse_and_fss_observable.py
SPLITTING=/home/joe/Linux_Code_Dev/REWL_Simulator/Scripts/Plotting/rfbm_temperature_splitting.py

echo -e "\n\nPlotting microcanonical observables"
python3 $PLOTTER $MODEL "microcanonical_observables" "Intensive Observable" "$COUPLING" "$VALUE" --Tc "$TC_VALUE"

echo -e "\n\nPlotting self-averaged observables"
python3 $PLOTTER $MODEL "self_averaged_observables" "Intensive Observable" "$COUPLING" "$VALUE" --Tc "$TC_VALUE"

echo -e "\n\nPlotting nonlinear observables"
python3 $PLOTTER $MODEL "nonlinear_observables" "Intensive Nonlinear Observable" "$COUPLING" "$VALUE" --Tc "$TC_VALUE"

python3 $SPLITTING "$COUPLING" "$VALUE" "$TC_VALUE"

echo -e "\nPlotting complete.\n"
