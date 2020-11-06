"""
    This script will recalculate the
    nonlinear observables from a job-array
    average canonical observables.
"""

import sys
sys.path.insert(1, "Scripts/Plotting")
from parse_file_header import collect_labels
import argparse
import numpy as np
from pathlib import Path

# Calculate the intensive susceptibility
# from an intensive order parameter
def susceptibility( m, m2, size, temp ):
    def avg( m, m2, size, temp ):
        return ( m2 - size * m ** 2 ) / temp
    def stderr( m_err, m2_err, m, size, temp ):
        return np.sqrt( m2_err ** 2 + (2 * size * m * m_err) ** 2 ) / temp

# Calculate the Binder cumulant from an
# intensive one-component order parameter
def scalar_binder_cumulant( ):
    def avg(m2, m4, size):
        return 1. - m4 / (3. * size * m2 ** 2)
    def stderr( m2_err, m4_err, m2, m4, size ):
        first_term = ( m4_err / (3 * size * m2 ** 2) ) ** 2
        second_term = ( 2 * m4 * m2_err / (3 * size * m2 ** 3) ) ** 2
        return np.sqrt( first_term + second_term )

# Calculate the Binder cumulant from an
# intensive two-component order parameter
# Calculate the Binder cumulant from an
# intensive one-component order parameter
def two_component_binder_cumulant( ):
    def avg(m2, m4, size):
        return 2. * ( 1. - m4 / (2. * size * m2 ** 2) )
    def stderr( m2_err, m4_err, m2, m4, size ):
        first_term = ( m4_err / (size * m2 ** 2) ) ** 2
        second_term = ( 2 * m4 * m2_err / (size * m2 ** 3) ) ** 2
        return np.sqrt( first_term + second_term )


def two_component_binder_cumulant( m2, m4, size ):
    return 2 * ( 1 - m4 / (2 * size * m2 ** 2) )

def setup_args():
    parser = argparse.ArgumentParser()
    parser = argparse.add_argument("model", help="Model type.", type=str)
    parser = argparse.add_argument("averages", help="Path to job-averaged canonical observable averages.", type=str)
    parser = argparse.add_argument("stderr", help="Path to job-averaged canonical observable standard error.", type=str)
    parser = argparse.add_argument("Lsize", help="Linear system size to consider.", type=str)
    parser.add_argument("coupling_symbol", help = "The constant value to parse through", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling", type = str)

    return parser.parse_args()

def read_in_data( average_file, stderr_file, comment = "#" ):

    labels = collect_labels( average_file, "Intensive Observable Names by Column", comment )

    header_lines = []
    data_file = open(average_file, "r")
    data_lines = data_file.readlines()

    line = 0
    cutoff = 0
    while data_lines[line][0] == comment:
        header_lines.append( data_lines[line] )
        if "Intensive Observable Names by Column" in data_lines[line]:
            cutoff == line
        line += 1

    avg_data = np.loadtxt( average_file, delimiter = "  ", comments = comment )
    stderr_data = np.loadtxt( stderr_file, delimiter = "  ", comments = comment )

    return labels, header_lines[:cutoff], avg_data, stderr_data

def calculate_Ising_observables( labels, Lsize, avg_data, stderr_data, dim = 2 ):

    Nfloat = float(Lsize) ** dim
    T_idx = labels.index("Temperature")
    m_idx = labels.index("Magnetization")
    m2_idx = labels.index("Magnetization2")
    m4_idx = labels.index("Magnetization4")

    nonlinear_labels = [ "1: Temperature", "2: Susceptibility", "3: Binder Cumulant" ]

    num_obs = 2

    nonlinear_avgs = np.zeros((avg_data.shape[0], 1 + num_obs))
    nonlinear_stderr = np.zeros(nonlinear_avgs.shape)

    # Set the temperature
    nonlinear_avgs[:,T_idx] = avg_data[:,T_idx]
    nonlinear_stderr[:,T_idx] = stderr_data[:,T_idx]

    # Calculate the susceptibility and error
    nonlinear_avgs[:,1] = susceptibility.avg( avg_data[:, m_idx], avg_data[:, m2_idx], Nfloat, avg_data[:, T_idx] )
    nonlinear_stderr[:,1] = susceptibility.stderr( stderr_data[:, m_idx], stderr_data[:, m2_idx], avg_data[:, m_idx], Nfloat, avg_data[:, T_idx] )

    # Calculate the Binder cumulant and error
    nonlinear_avgs[:,2] = scalar_binder_cumulant.avg( avg_data[:, m2_idx], avg_data[:, m4_idx], Nfloat )
    nonlinear_stderr[:,2] = scalar_binder_cumulant.stderr( stderr_data[:, m2_idx], stderr_data[:, m4_idx], avg_data[:, m2_idx], avg_data[:, m4_idx], Nfloat )

    return nonlinear_labels, nonlinear_avgs, nonlinear_stderr

# TODO: Calculate Ashkin Teller observables

def write_output_data( input_file, nonlinear_labels, header_lines, nonlinear_avgs, nonlinear_stderr ):

    output_path = Path( input_file ).parent

    rewl_substring = input_file[ input_file.find("REWL_") : input_file.find(".job") ]

    full_output = str(output_path) + "nonlinear_observables-" + rewl_substring

    header_string = ""
    for string in header_lines:
        header_string += string.rstrip()
    header_string += "\n# Intensive Nonlinear Observable Names by Column"
    for label in nonlinear_labels:
        header_string += "\n#    " + label
    header_string += "\n"

    np.savetxt( full_output + ".job_mean", nonlinear_avgs, delimiter = "  ", newline="\n", header=header_string, comments="" )
    np.savetxt( full_output + ".job_stderr", nonlinear_stderr, delimiter = "  ", newline="\n", header=header_string, comments="" )

    return None


def main():

    args = setup_args()

    labels, header_lines, avg_data, stderr_data = read_in_data( args.averages, args.stderr )

    nonlinear_labels, nonlinear_avgs, nonlinear_stderr = None, None, None

    if "Ising" in args.model:

        nonlinear_labels, nonlinear_avgs, nonlinear_stderr = calculate_Ising_observables( labels, avg_data, stderr_data )


    write_output_data( args.averages, nonlinear_labels, header_lines, nonlinear_avgs, nonlinear_stderr )

    return 0

if __name__ == "__main__":
    main()
