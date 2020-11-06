"""
    This script will use NumPy to re-calculate
    the canonical observables for many independent
    jobs or many disorder configurations.
"""

import sys
sys.path.insert(1, "Scripts/Plotting")

import argparse
from parse_file_header import collect_labels
import numpy as np
import os
from pathlib import Path

job_id_string = "JOBID"
filename_base = "self_averaged_observables"
observable_marker = "Intensive Observable Names by Column"

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help = "Path to Job Arrays.", type = str)
    parser.add_argument("Lsize", help = "System size to average.", type = str)
    parser.add_argument("coupling_symbol", help = "The constant value to parse through", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling", type = str)

    return parser.parse_args()

def find_string_value( string_type, file_string ):

    start = file_string.find(string_type + "-") + len(string_type + "-")
    end = start + file_string[start:].find("_")

    return file_string[start:end]

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

def read_in_data( input_path, Lsize, coupling_symbol, coupling_value, comment = "#" ):

    file_header = []
    extracted_file = ""
    labels = []
    data_tuples = []
    L_string = "L-" + ("%d" % int(Lsize))
    key_string = coupling_symbol + "-" + ("%.6f" % float(coupling_value))

    print(os.listdir(input_path))

    for fl in os.listdir( input_path ):
        if not os.path.isdir( fl ) and ( filename_base in fl and L_string in fl and key_string in fl ):

            if len(labels) == 0:
                extracted_file = extract_file_name( fl )
                file_header = get_file_header( input_path + fl )
                labels = collect_labels( input_path + fl, observable_marker, comment )

            ID = find_string_value(job_id_string, fl)

            data = np.loadtxt( input_path + fl, delimiter = "  ", dtype = float, comments = comment)

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

    print(labels)

    final_data = np.zeros(data_tuples[0][1].shape)

    # Average Data along columns
    for col in range(0, len(labels)):
        for job in range(0, num_jobs):
            final_data[:,col] += data_tuples[job][1][:,col]
        final_data[:,col] /= num_jobs

    # Recompute the specific heat
    # TODO: This will break for higher dimensions
    Nfloat = float(Lsize) ** 2.
    en_idx = labels.index("Energy")
    en2_idx = labels.index("Energy2")
    cv_idx = labels.index("Specific Heat")

    final_data[:, cv_idx] = specific_heat( final_data[:, en2_idx], final_data[:, en_idx], Nfloat, final_data[:, 0] )

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
    Nfloat = float(Lsize) ** 2.
    en_idx = labels.index("Energy")
    en2_idx = labels.index("Energy2")
    cv_idx = labels.index("Specific Heat")

    final_stderr[:, cv_idx] = specific_heat_error( final_stderr[:, en2_idx], final_stderr[:, en_idx], final_averages[:, en_idx], Nfloat, final_averages[:, 0] )

    return final_stderr

# Write data out to the proper file path
def write_out_data( output_file, input_path, num_jobs, header_lines, final_averages, final_stderr ):

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

    return None

def main():

    cl_args = setup_args()

    output_file, header_lines, labels, data_tuples = read_in_data( cl_args.input_path, cl_args.Lsize, cl_args.coupling_symbol, cl_args.coupling_value )

    if len(data_tuples) == 0:
        print("\nNo data to average.\n")
        return 1

    final_averages, final_stderr = 0,0
    final_averages = average_job_data( cl_args.Lsize, labels, data_tuples )

    final_stderr = stderr_job_data( cl_args.Lsize, labels, data_tuples, final_averages )

    write_out_data( output_file, cl_args.input_path, len(data_tuples), header_lines, final_averages, final_stderr )

    return 0


if __name__ == "__main__":
    main()