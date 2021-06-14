"""
    This script will use NumPy to re-calculate
    the canonical observables for many independent
    jobs or many disorder configurations.
"""

import sys
sys.path.insert(1, "Scripts/Plotting")

import argparse
from parse_file_header import collect_labels
from multiple_couplings import *
import numpy as np
import os
from pathlib import Path

job_id_string = "JOBID"
filename_base = "self_averaged_observables"
microname_base = "microcanonical_observables"
nonlinear_base = "nonlinear_observables"
micro_dilutor = 0.1 # fraction by which to dilute the self-averaged observables for the microcanonical ones
file_types = [ filename_base, nonlinear_base ]
observable_markers = [ "Intensive Observable Names by Column", "Intensive Nonlinear Observable Names by Column" ]

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help = "Path to Job Arrays.", type = str)
    parser.add_argument("Lsize", help = "System size to average.", type = str)
    parser.add_argument("coupling_symbol", help = "The constant value(s) to parse through. May be a space separated list.", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling(s). May be a space separated list.", type = str)

    return parser.parse_args()

def extract_file_name( file_string ):

    start = file_string.find(job_id_string + "-")
    end = start + file_string[start:].find("_")

    front = file_string[:start]
    back = file_string[end+1:]

    return front + back

def get_file_header( file_string, comment = "#" ):

    header_lines = []
    data_file = open(file_string, "r")
    data_lines = data_file.readlines()

    line = 0
    while data_lines[line][0] == comment:
        header_lines.append( data_lines[line] )
        line += 1

    return header_lines

def shift_array( arr, shift, idx=0 ):
    return arr + (shift - arr[idx])

def correct_entropy( entropy_arr, corrector, idx=0 ):
    value = np.copy(entropy_arr[idx])
    return shift_array( entropy_arr, corrector, idx ), value

def correct_free_energy( free_energy, entropy_value, corrector, temperature ):
    return free_energy - ( corrector - entropy_value ) * temperature

def read_in_data( input_path, Lsize, coupling_tuples, file_type, observable_marker, comment = "#" ):

    file_header = []
    extracted_file = ""
    labels = []
    data_tuples = []
    L_string = "L-" + ("%d" % int(Lsize))
    key_string = get_coupling_string(coupling_tuples, isfloat = True)

    for fl in os.listdir( input_path ):
        if not os.path.isdir( fl ) and ( file_type in fl and L_string in fl and key_string in fl ):

            if len(labels) == 0:
                extracted_file = extract_file_name( fl )
                file_header = get_file_header( input_path + fl )
                labels = collect_labels( input_path + fl, observable_marker, comment )

            ID = find_string_value(job_id_string, fl)

            data = np.loadtxt( input_path + fl, delimiter = "  ", dtype = "float64", comments = comment)

            # Define the entropy corrector for the Ising model
            # at T = 0. Change this for different models
            entropy_corrector = np.log(2) / int(Lsize)**2
            if 'Ashkin_Teller' in input_path:
                entropy_corrector = np.log(4) / int(Lsize)**4

            if file_type == filename_base:
                # Only look for entropy and free energy in the self-averaged observables
                free_energy_idx, entropy_idx = labels.index('Free Energy'), labels.index('Entropy')
                old_entropy_T0, data[:, entropy_idx] = correct_entropy( data[:,entropy_idx], entropy_corrector )
                data[:, free_energy_idx] = correct_free_energy( data[:,free_energy_idx], old_entropy_T0, entropy_corrector, data[:, labels.index('Temperature')] )

            data_tuples.append( (ID, data) )


    data_tuples.sort(key = lambda tup: int(tup[0]))

    return extracted_file, file_header, labels, data_tuples

# Calculate the specific heat from intensive energy
# and energy2.
def specific_heat( energy2, energy, size, temp ):

    return ( energy2 - size * energy * energy ) / ( temp ** 2. )

# Calculate error for specific heat from
# intensive energy
def specific_heat_error( en2_err, en_err, energy, size, temp ):

    return np.sqrt( en2_err ** 2. + (2 * size * energy * en_err) ** 2. ) / (temp ** 2.)

# Calculate the average of the final data
def average_job_data( Lsize, labels, data_tuples ):

    num_jobs = len(data_tuples)

    print("Observable Labels:", labels, "\n")

    final_data = np.zeros(data_tuples[0][1].shape)

    # Average Data along columns
    for col in range(0, len(labels)):
        for job in range(0, num_jobs):
            final_data[:,col] += data_tuples[job][1][:,col]
        final_data[:,col] /= num_jobs

    # Recompute the specific heat
    # TODO: This will break for higher dimensions
    # TODO: get rid of this
    """
    Nfloat = float(Lsize) ** 2.
    en_idx = labels.index("Energy")
    en2_idx = labels.index("Energy2")
    cv_idx = labels.index("Specific Heat")

    #final_data[:, cv_idx] = specific_heat( final_data[:, en2_idx], final_data[:, en_idx], Nfloat, final_data[:, 0] )
    """

    return final_data

# Calculate the standard error of the
# final data
def stderr_job_data( Lsize, labels, data_tuples, final_averages ):

    num_jobs = len(data_tuples)

    divisor = num_jobs
    if num_jobs == 1:
        divisor *= num_jobs
    else:
        divisor *= (num_jobs - 1)

    final_stderr = np.zeros(final_averages.shape)

    # Get Standard Error along columns
    for col in range(0, len(labels)):
        for job in range(0, num_jobs):
            final_stderr[:,col] += ( data_tuples[job][1][:,col] - final_averages[:,col] ) ** 2.
        final_stderr[:,col] = np.sqrt( final_stderr[:,col] / divisor )

    # Recompute the specific heat error
    # TODO: This will break for higher dimensions
    # TODO: Get rid of this
    """Nfloat = float(Lsize) ** 2.
    en_idx = labels.index("Energy")
    en2_idx = labels.index("Energy2")
    cv_idx = labels.index("Specific Heat")

    #final_stderr[:, cv_idx] = specific_heat_error( final_stderr[:, en2_idx], final_stderr[:, en_idx], final_averages[:, en_idx], Nfloat, final_averages[:, 0] )
    """

    return final_stderr

# Get the microcanonical data
def microcanonical_observable_data( labels, final_averages, final_stderr ):

    num_rows = final_averages.shape[0]
    num_cols = final_averages.shape[1] -2 # Don't include the specific heat or temperature

    #micro_rows = int( np.floor( (1 - micro_dilutor) * num_rows )  )
    micro_rows=20000
    micro_increment = num_rows // micro_rows
    micro_avgs = np.zeros((micro_rows, num_cols))
    micro_stderr = np.zeros(micro_avgs.shape)

    energy_idx = labels.index("Energy")
    entropy_idx = labels.index("Entropy")
    cv_idx = labels.index("Specific Heat")
    obs1_idx = cv_idx + 1

    micro_labels = [ "1: Energy", "2: logDoS" ]
    micro_avgs[:,0] = final_averages[::micro_increment, energy_idx]
    micro_stderr[:,0] = final_stderr[::micro_increment, energy_idx]
    micro_avgs[:,1] = final_averages[::micro_increment, entropy_idx]
    micro_stderr[:,1] = final_stderr[::micro_increment, entropy_idx]

    for idx in range(obs1_idx, len(labels)):
        micro_idx = idx - obs1_idx + 2
        label_idx = micro_idx + 1
        micro_avgs[:,micro_idx] = final_averages[::micro_increment, idx]
        micro_stderr[:,micro_idx] = final_stderr[::micro_increment, idx]
        micro_labels.append("%d: %s" % (label_idx, labels[idx]))

    return micro_labels, micro_avgs, micro_stderr

# Get the microcanonical observable header
def microcanonical_header( micro_labels, header_lines, observable_marker = observable_markers[0] ):

    full_line_obs_marker = "# " + observable_marker + "\n"
    micro_header_string = ""
    for ldx in range(0, header_lines.index(full_line_obs_marker) + 1):
        micro_header_string += header_lines[ldx].rstrip()
        if ldx != header_lines.index(full_line_obs_marker):
            micro_header_string += "\n"

    for ldx in range(0, len(micro_labels)):
        micro_header_string += "\n#    " + micro_labels[ldx]

    micro_header_string += "\n#"

    return micro_header_string

# Write data out to the proper file path
def write_out_data( output_file, input_path, num_jobs, labels, header_lines, final_averages, final_stderr, file_type ):

    path = Path(input_path)
    output_path = path.parent

    output = str(output_path) + "/" + output_file

    header_lines.insert(1, "#\n# Number of independent jobs = %d\n" % num_jobs)

    header = ""
    for ldx in range(0, len(header_lines)):
        header += header_lines[ldx].rstrip()
        if ldx != len(header_lines)-1:
            header += "\n"

    np.savetxt(output + ".job_mean", final_averages, delimiter="  ", newline="\n", header=header, comments="")
    np.savetxt(output + ".job_stderr", final_stderr, delimiter="  ", newline="\n", header=header, comments="")

    if file_type == filename_base:
        # Write out the microcanonical data from the self-averaged data
        micro_labels, micro_avgs, micro_stderr = microcanonical_observable_data( labels, final_averages, final_stderr )
        micro_header = microcanonical_header( micro_labels, header_lines )
        micro_output = str(output_path) + "/" + microname_base + "-" + output_file[ output_file.find("REWL_") :  ]
        np.savetxt(micro_output + ".job_mean", micro_avgs, delimiter="  ", newline="\n", header=micro_header, comments="")
        np.savetxt(micro_output + ".job_stderr", micro_stderr, delimiter="  ", newline="\n", header=micro_header, comments="")

    return None

def main():

    cl_args = setup_args()

    Lvalue, coupling_tuples = parse_couplings( cl_args.coupling_symbol, cl_args.coupling_value )

    # Exit the program if the coupling tuples are null
    if len(coupling_tuples) == 0:
        return

    print("="*70, "\nPost-Simulation Inter-Job Statistics Calculator\n" + "="*70 + "\n")
    for type_idx in range(len(file_types)):

        print("\nNow analyzing %s data for L = %s\n" % (file_types[type_idx], cl_args.Lsize))

        output_file, header_lines, labels, data_tuples = read_in_data( cl_args.input_path, cl_args.Lsize, coupling_tuples, file_types[type_idx], observable_markers[type_idx] )

        if len(data_tuples) == 0:
            print("\nNo data to average.\n")
            return 1
        else:
            print("\nTotal Files: %d Jobs Found.\n" % len(data_tuples))

        final_averages, final_stderr = 0,0
        final_averages = average_job_data( cl_args.Lsize, labels, data_tuples )

        final_stderr = stderr_job_data( cl_args.Lsize, labels, data_tuples, final_averages )

        write_out_data( output_file, cl_args.input_path, len(data_tuples), labels, header_lines, final_averages, final_stderr, file_types[type_idx] )

        if type_idx != len(file_types) - 1:
            print("*"*70)

    print("="*70, "\n")

    return 0


if __name__ == "__main__":
    main()
