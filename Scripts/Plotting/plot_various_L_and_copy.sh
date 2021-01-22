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

echo -e "\nCoupling value $COUPLING\n"

./Scripts/Plotting/plot_model_sizes_for_fixed_coupling.sh -m "$MODEL" -d "$DATE" -c "$COUPLING" -v "$VALUE" -T "$TC_VALUE"

if [ $? == 0 ]
then
    ./Scripts/Plotting/move_data_to_windows.sh -m "$MODEL" -d "$DATE" -c "$COUPLING" -v "$VALUE"
fi


