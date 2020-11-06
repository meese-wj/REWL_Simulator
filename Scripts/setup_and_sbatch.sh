#!/bin/bash

# This needs to come from the sbatch script
Lold=312
MODEL="Ising2d"
NPROCS=48
NTASKS_PER_NODE=24
TASKSTR="#SBATCH --ntasks-per-node="

MODELSTR="MODEL="
SETUP="Scripts/msi_setup_Lsize.sh"
SBATCH="Scripts/sbatch_${MODEL}.slurm"

for L in {24..96..24}
do	
	echo -e "\nSetting up L = ${L} with Lold = ${Lold}...\n"

	sed -i "/Lsize=/c\\Lsize=${L}" $SETUP
	sed -i "/${MODELSTR}/c\\${MODELSTR}\"${MODEL}\"" $SETUP
	
	if [ $? == 0 ]
	then	
		./${SETUP}

		sed -i "s/${Lold}/${L}/g" $SBATCH
		sed -i "/NPROCS=/c\\NPROCS=${NPROCS}" $SBATCH
		sed -i "/${TASKSTR}/c\\${TASKSTR}${NTASKS_PER_NODE}" $SBATCH

		sbatch ${SBATCH}
	else
		echo -e "\nL = ${L} failed.\n"
	fi

	Lold=$L
done
