"""
    Use this script to plot many
    density plot pngs.
"""

import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg')  # THIS IS REQUIRED FOR WSL2
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from multiple_couplings import *

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("micro_dir", help = "The directory of the microcanonical data.", type = str)
    parser.add_argument("sifter_value", help = "The value of the sifting coupling (usually system size).", type = str )
    parser.add_argument("coupling_symbol", help = "The constant value(s) to parse through. May be a space separated list.", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling(s). May be a space separated list.", type = str)
    parser.add_argument("peak_T", default = None, help = "Temperature to center the temperature array around.", type = str)

    return parser.parse_args()

def make_temperatures( peak_T ):
    """
        Define the temperature array here.
    """
    temps = np.array([ 0.250, 0.500, 0.750, 0.800, 0.850, 0.900,
                       0.925, 0.950, 0.975, 1.000, 1.025, 1.050,
                       1.075, 1.100, 1.150, 1.200, 1.250, 1.500,
                       1.750, 2.000 ])
    #temps = np.linspace( 2.0, 2.6, 60 )
    #temps = np.linspace( 2.0, 3.6, 60 )
    #temps = np.linspace( 0.1, 1.9, 38 )
    temps = np.array([1.609])
    return temps

def get_file_bin( file_string, key_string = "bin-" ):
    """
        Strip the file string for the bin number.
    """
    start = file_string.find(key_string) + len(key_string)
    end = file_string[start:].find("_")
    bin_value = int(file_string[start:start+end])
    return bin_value

def print_dictionary( dictionary ):
    """
        Print the parameter dictionary.
    """
    print("\nParameters:\n")
    for key in dictionary:
        print("\t%s =" % key, dictionary[key] )
    print("\n")
    return

def parse_density_header( file_string, comment = "#" ):
    """
        Read in the file header from the
        density plots and extract the
        parameter dictionary.
    """
    data_file = open(file_string, "r")
    data_lines = data_file.readlines()
    line = 0
    header_lines = []
    while data_lines[line][0] == comment:
        header_lines.append( data_lines[line] )
        line += 1

    axis_bins, total_bins, axis_min, axis_max = 0, 0, 0, 0
    for line in header_lines:
        line_list = line.split()
        if "axis" in line_list:
            if "bins" in line_list:
                axis_bins = int(line_list[-1])
            elif "maximum" in line_list:
                axis_max = float(line_list[-1])
            elif "minimum" in line_list:
                axis_min = float(line_list[-1])

    total_bins = axis_bins ** 2
    binwidth = (axis_max - axis_min) / axis_bins
    density_params = { "axis_bins" : axis_bins, "total_bins" : total_bins, "axis_min" : axis_min,
                       "axis_max" : axis_max, "binwidth" : binwidth }
    return density_params

def get_data( file_string, comment = "#" ):
    """
        Wrapper around loadtxt.
    """
    return np.loadtxt( file_string, delimiter = "  ", dtype=float, comments=comment )

def correct_directory( data_dir ):
    """
        Add a / symbol to the data_dir
        if it's not present.
    """
    filepath = data_dir
    if data_dir[-1] != "/":
        filepath += "/"
    return filepath

def identify_jobs( data_dir, coupling_tuples ):
    """
        Look through the data directory and
        construct a dictionary of relevant
        files.

        The keys are the jobids and
        the values are a list of all density
        bin files for that specific job.
    """
    job_marker = 'JOBID-'
    all_jobs = {}
    for f in os.listdir( data_dir ):
        file_name = correct_directory( data_dir ) + f
        boolean = "density_plots_bin" in file_name
        boolean = boolean and couplings_in_file(coupling_tuples, file_name)
        boolean = boolean and job_marker in file_name
        if boolean:
            start = file_name.find( job_marker ) + len(job_marker)
            jobid = file_name[ start : start + file_name[start:].find("_") ]
            if jobid not in all_jobs:
                all_jobs[jobid] = [f]
            else:
                all_jobs[jobid].append( f )
    return all_jobs

def yes_extract_file_data( file_name, coupling_tuples ):
    """
        Determine if the file meets the
        coupling tuple standard to be
        read.
    """
    boolean = "density_plots_bin" in file_name
    boolean = boolean and couplings_in_file(coupling_tuples, file_name)
    return boolean

def collect_single_job( data_dir, coupling_tuples, jobid = 'None', job_list = None, print_params = False ):
    """
        Go into the data_dir and then read
        in all density plots per bin for a
        single job.

        If jobid == 'None', this code will
        not care about the JOBID in the
        file names.
    """
    density_data_folder = correct_directory( data_dir )
    bin_values = []
    density_data = []
    data_tuples = []
    density_params = {}
    if jobid == 'None':
        for f in os.listdir(density_data_folder):
            filepath = density_data_folder + f
            if yes_extract_file_data( filepath, coupling_tuples ):
                data_tuples.append( (get_file_bin( filepath ), get_data( filepath )) )
                if len(data_tuples) == 1:
                    density_params = parse_density_header( filepath )
                    if print_params:
                        print_dictionary( density_params )
    else:
        # Now load all bins for a single specified jobid
        for jobfile in job_list:
            filepath = density_data_folder + jobfile
            data_tuples.append( (get_file_bin( filepath ), get_data( filepath )) )
            if len(data_tuples) == 1:
                density_params = parse_density_header( filepath )
                if print_params:
                    print_dictionary( density_params )

    data_tuples.sort(key = lambda tup: tup[0])
    for tup in data_tuples:
        bin_values.append(tup[0])
        density_data.append(tup[1])

    bin_values = np.array(bin_values)
    density_data = np.array(density_data)

    return bin_values, density_data, density_params

# ****************************************************************************************
# Canonical Statistics
def energy_weights( logDoS, energy, temperature ):
    """
        Calculate the energy weights and
        partition function for a given
        temperature.

        These weights are scaled by the
        maximum exponent.
    """
    exp_arg = logDoS - energy / temperature
    exp_arg -= np.max( exp_arg )
    weights = np.exp( exp_arg )
    return weights, np.sum(weights)

def probability_per_bin( logDoS, energy, temperature ):
    """
        Return the probability per bin.
    """
    weights, Zfunc = energy_weights( logDoS, energy, temperature )
    return weights / Zfunc

def canonical_average( observable, logDoS, energy, temperature ):
    """
        Rescale the the density plot so
        that its maximum is at one.
    """
    probs = probability_per_bin( logDoS, energy, temperature )
    output = probs[0] * observable[0]
    for endx in range( 1, observable.shape[0], 1 ):
        output += probs[endx] * observable[endx]
    # Normalize by the dividing out the maximum
    return output/np.max(output)
# ****************************************************************************************

def read_microcanonical_data( micro_file, Lsize, comment = "#" ):
    """
        Read in the microcanonical data. This returns
        the EXTENSIVE energy and logDoS data as an
        array.
    """
    data = np.loadtxt( micro_file, dtype = float, delimiter = "  ", comments = comment )
    energy_logDoS = Lsize ** 2 * data[:,0:2]
    return energy_logDoS

def load_job_logDoS_data( data_dir, sifter_coupling, sifter_value, coupling_tuples, jobname = 'None' ):
    """
        Search through the microcanonical
        directory and find the proper
        data file to read in.
    """
    energy_logDoS = []
    for f in os.listdir( data_dir ):
        filepath = correct_directory(data_dir) + f
        extract = "microcanonical" in filepath
        extract = extract and couplings_in_file(coupling_tuples, filepath)
        if jobname != 'None':
            extract = extract and ('JOBID-%s_' % jobname in filepath)
        if extract:
            Lsize = 0
            if sifter_coupling == 'L':
                Lsize = int(sifter_value)
            else:
                for tup in coupling_tuples:
                    if tup[0] == 'L':
                        Lsize = int(tup[1])
            energy_logDoS = read_microcanonical_data( filepath, Lsize )
    return energy_logDoS

def collect_all_jobs( density_dir, coupling_tuples ):
    """
        Search through the density directory
        and gather the density data per bin
        for each job.

        The data structures returned are the
        jobname dictionary, one ndarray,
        one job density dictionary, and
        then the parameter dictionary.

        The job density output needed to be
        changed to a dictionary since the jobs
        may have different numbers of energy
        bins and so they cannot be broadcast
        together.

        all_job_bins =
            ( num_jobs, num_bins, energy or logdos )

        all_job_density.shape =
            ( num_jobs, num_bins, axis_bins, axis_bins )
    """
    jobnames = identify_jobs( density_dir, coupling_tuples )
    all_job_bins = []
    all_job_density = {}
    density_params = {}
    print_params = True
    if len(jobnames) > 0:
        print("\n%d jobs found. Loading them now.\n" % len(jobnames))
        # Now iterate through the jobnames dictionary
        # to peruse all bin files per jobid
        for jobid in jobnames:
            bin_data, dens_data, temp_dict = collect_single_job( density_dir, coupling_tuples, jobid, jobnames[jobid], print_params = print_params )
            if print_params:
                print_params = False
                density_params = temp_dict
                #all_job_density = np.zeros( ( len(jobnames), dens_data.shape[0], dens_data.shape[1], dens_data.shape[2] ) )
            all_job_bins = bin_data
            all_job_density[int(jobid)] = dens_data

    else:
        print("\nNo job data found. Assuming you then mean to consider a single job...\n")
        all_job_bins = [0]
        all_job_density = {0:''}
        all_job_energy_logDoS = [0]
        all_job_bins[0], all_job_density[0], density_params = collect_single_job( density_dir, coupling_tuples, print_params = print_params )

        all_job_bins = np.array(all_job_bins)
        #all_job_density = np.array(all_job_density)

    return jobnames, all_job_bins, all_job_density, density_params

def compute_single_job_average( density_data, temp_array, micro_dir, sifter_coupling, sifter_values, coupling_tuples, jobname = 'None' ):
    """
        Compute the thermodynamic densities for
        a given temperature array for a single
        job.

        This function reads in the proper
        microcanonical data for the relevant
        job.

        The return type is an ndarray of the
        density plots.

        thermo_densities.shape =
            ( num_temps, axis_bins, axis_bins )
    """
    energy_logDoS = load_job_logDoS_data( micro_dir, sifter_coupling, sifter_values, coupling_tuples, jobname )
    axis_size = density_data.shape[1]
    thermo_densities = np.zeros( (len(temp_array), axis_size, axis_size) )
    for Tdx in range(len(temp_array)):
        thermo_densities[Tdx,:,:] = canonical_average( density_data, energy_logDoS[:,1], energy_logDoS[:,0], temp_array[Tdx] )
    return thermo_densities

def collect_all_thermo_densities(  all_job_density, jobnames, temp_array, micro_dir, sifter_coupling, sifter_value, coupling_tuples ):
    """
        From all the job densities, calculate
        the thermodynamic average for each
        job. Return all averages for each
        temperature.

        all_thermo_densities.shape =
            ( num_jobs, num_temps, axis_bins, axis_bins )
    """
    all_thermo_densities = []
    njobs = len( all_job_density )
    job_index = 0
    for jobid in jobnames:
        job = 'None'
        if njobs > 1:
            job = jobid
        temp_densities = compute_single_job_average( all_job_density[ int( jobid ) ], temp_array, micro_dir, sifter_coupling, sifter_value, coupling_tuples, jobname = job )

        if len(all_thermo_densities) == 0:
            all_thermo_densities = np.zeros( ( njobs, len(temp_array), temp_densities.shape[1], temp_densities.shape[2] ) )

        all_thermo_densities[ job_index, :, :, : ] = temp_densities
        job_index += 1

    return all_thermo_densities


def compute_multiple_job_average( all_thermo_densities ):
    """
        Average the density plots over the
        jobs for fixed temperature.
        Return an ndarry of the following
        shape:

        final_density_averages.shape =
            ( num_temps, axis_bins, axis_bins )
    """
    njobs = all_thermo_densities.shape[0]
    ntemps = all_thermo_densities.shape[1]
    final_density_averages = np.zeros( (all_thermo_densities.shape[1], all_thermo_densities.shape[2], all_thermo_densities.shape[3]) )
    for Tdx in range(ntemps):
        for jdx in range(njobs):
            final_density_averages[Tdx,:,:] += all_thermo_densities[jdx,Tdx,:,:]
        final_density_averages[Tdx,:,:] *= 1./njobs
        # Normalize by the maximum
        final_density_averages[Tdx,:,:] *= 1./np.max(final_density_averages[Tdx,:,:])
    return final_density_averages

def plot_final_thermo_densities( plot_dir, density_params, sifter_coupling, sifter_value, coupling_tuples, temp_array, final_density_averages, interpolation = 'None' ):
    """
        Plot the 2d density plots for
        each temperature.
    """
    plt.style.use("dark_background")
    axis_min = density_params["axis_min"]
    axis_max = density_params["axis_max"]

    for Tdx in range( len(temp_array) ):
        fig, ax = plt.subplots(1,1, figsize = (6,5))
        image = ax.imshow( final_density_averages[Tdx], extent = [axis_min, axis_max, axis_min, axis_max],
                           cmap = "magma", vmin = 0., vmax = 1, interpolation=interpolation)
        cbar = fig.colorbar(image, ax=ax)
    #     cbar.ax.set_ylabel("Normalized Occurences", rotation=-90)
        ax.set_xlabel(r"$\frac{1}{N} \sum_i \sigma_i$",fontsize=14)
        ax.set_ylabel(r"$\frac{1}{N} \sum_i \tau_i$",fontsize=14)
        ax.set_title(r"$T = %.3f$, $%s = %s$, %s" % (temp_array[Tdx], sifter_coupling, sifter_value, latex_couplings(coupling_tuples)[1:]))

        plot_name = "density T %.3f %s %s %s" % ( temp_array[Tdx], sifter_coupling, sifter_value, get_coupling_string(coupling_tuples) )
        if interpolation != 'None':
            plot_name = interpolation + " " + plot_name

        plt.tight_layout()
        plt.savefig( correct_directory( plot_dir ) + plot_name + ".png", facecolor="black", edgecolor="none", dpi=300 )
        plt.close()


def main():

    args = setup_args()

    density_dir = correct_directory( args.micro_dir ) + "Density_Plots"
    if not os.path.isdir( density_dir ):
        print("\n%s does not exist. Exiting.\n" % density_dir)
        return 1

    sifter_coupling, coupling_tuples = parse_couplings( args.coupling_symbol, args.coupling_value )

    temp_array = make_temperatures( args.peak_T )

    # Gather all job data
    jobnames, all_job_bins, all_job_density, density_params = collect_all_jobs( density_dir, coupling_tuples )

    # Compute all thermodynamic averages
    print("Computing all thermodynamic densities...")
    all_thermo_densities = collect_all_thermo_densities(  all_job_density, jobnames, temp_array, args.micro_dir, sifter_coupling, args.sifter_value, coupling_tuples )
    print("Densities computed.")

    # Average the thermodynamic densities
    # over all jobs
    print("\nComputing job-averages...")
    final_density_averages = compute_multiple_job_average( all_thermo_densities )
    print("Job-averages computed.")

    # Plot the job averaged densities
    print("\nPlotting densities...")
    plot_dir = ''
    if len(jobnames) > 0:
        plot_dir = correct_directory( args.micro_dir ) + "../Density_Figures/"
    else:
        plot_dir = correct_directory( args.micro_dir ) + "Density_Figures/"

    if not os.path.isdir( plot_dir ):
        print("%s does not exist. Creating directory now." % plot_dir)
        os.mkdir( plot_dir )

    plot_dir += get_coupling_string( coupling_tuples )
    if not os.path.isdir( plot_dir ):
        print("%s does not exist. Creating directory now." % plot_dir)
        os.mkdir( plot_dir )

    plot_final_thermo_densities( plot_dir, density_params, sifter_coupling, args.sifter_value, coupling_tuples, temp_array, final_density_averages, interpolation = 'None' )
    plot_final_thermo_densities( plot_dir, density_params, sifter_coupling, args.sifter_value, coupling_tuples, temp_array, final_density_averages, interpolation = 'gaussian' )
    print("Densities plotted in", plot_dir)

    return 0

if __name__ == "__main__":

    main()
