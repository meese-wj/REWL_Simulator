"""
    Use this script to calculate the average
    side size in the square parts of the
    density of states. This is to be compared
    to the clean Ising model.
"""
import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg')  # THIS IS REQUIRED FOR WSL2
import matplotlib.pyplot as plt
from multiple_couplings import *
import ashkin_teller_density_plotter as atdp

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("micro_dir", help = "The directory of the microcanonical data.", type = str)
    parser.add_argument("sifter_value", help = "The value of the sifting coupling (usually system size).", type = str )
    parser.add_argument("coupling_symbol", help = "The constant value(s) to parse through. May be a space separated list.", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling(s). May be a space separated list.", type = str)

    return parser.parse_args()

def make_temperatures():
    """
        Define the temperature array here.
    """
    Tmin = 2.0
    Tmax = 2.8
    NTemps = 2
    temps = np.linspace( Tmin, Tmax, NTemps )
    return temps

def return_slice( density_plot, axis, off_axis_idx, density_params ):
    """
        Get the slice of the density plot
        along the the off_axis_idx
    """
    num_bins = density_params["axis_bins"]
    result = []
    if ( axis == "sigma" ):
        result = density_plot[ num_bins // 2, : ]
    else:
        result = density_plot[ :, num_bins // 2 ]

    return result

def main():

    args = setup_args()

    density_dir = atdp.correct_directory( args.micro_dir ) + "Density_Plots"
    if not os.path.isdir( density_dir ):
        print("\n%s does not exist. Exiting.\n" % density_dir)
        return 1

    sifter_coupling, coupling_tuples = parse_couplings( args.coupling_symbol, args.coupling_value )

    temp_array = make_temperatures()

    # Gather all job data
    jobnames, all_job_bins, all_job_density, density_params = atdp.collect_all_jobs( density_dir, coupling_tuples )

    # Compute all thermodynamic averages
    print("Computing all thermodynamic densities...")
    all_thermo_densities = atdp.collect_all_thermo_densities(  all_job_density, jobnames, temp_array, args.micro_dir, sifter_coupling, args.sifter_value, coupling_tuples )
    print("Densities computed.")

    # Average the thermodynamic densities
    # over all jobs
    print("\nComputing job-averages...")
    final_density_averages = atdp.compute_multiple_job_average( all_thermo_densities )
    print("Job-averages computed.")

    print(final_density_averages)
    print(final_density_averages.shape)

    print(return_slice( final_density_averages[0,:,:], "sigma", density_params["axis_bins"] // 2, density_params ))
    print(return_slice( final_density_averages[1,:,:], "tau", density_params["axis_bins"] // 2, density_params ))


    return 0

if __name__ == "__main__":
    main()
