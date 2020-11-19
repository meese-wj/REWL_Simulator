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
class susceptibility:
    def __init__(self):
        return
    def avg( self, m, m2, size, temp ):
        return ( m2 - size * m ** 2 ) / temp
    def stderr( self, m_err, m2_err, m, size, temp ):
        return np.sqrt( m2_err ** 2 + (2 * size * m * m_err) ** 2 ) / temp

# Calculate the Binder cumulant from an
# intensive one-component order parameter
class scalar_binder_cumulant:
    def __init__(self):
        return
    def avg(self, m2, m4, size):
        return 1. - m4 / (3. * size * m2 ** 2)
    def stderr(self, m2_err, m4_err, m2, m4, size ):
        first_term = ( m4_err / (3 * size * m2 ** 2) ) ** 2
        second_term = ( 2 * m4 * m2_err / (3 * size * m2 ** 3) ) ** 2
        return np.sqrt( first_term + second_term )

# Calculate the Binder cumulant from an
# intensive two-component order parameter
# Calculate the Binder cumulant from an
# intensive one-component order parameter
class two_component_binder_cumulant:
    def __init__(self):
        return
    def avg(self, m2, m4, size):
        return 2. * ( 1. - m4 / (2. * size * m2 ** 2) )
    def stderr(self, m2_err, m4_err, m2, m4, size ):
        first_term = ( m4_err / (size * m2 ** 2) ) ** 2
        second_term = ( 2 * m4 * m2_err / (size * m2 ** 3) ) ** 2
        return np.sqrt( first_term + second_term )

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", help="Model type.", type=str)
    parser.add_argument("averages", help="Path to job-averaged canonical observable averages.", type=str)
    parser.add_argument("stderr", help="Path to job-averaged canonical observable standard error.", type=str)
    parser.add_argument("Lsize", help="Linear system size to consider.", type=str)

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
            cutoff = line
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

    num_obs = len(nonlinear_labels) - 1

    susc = susceptibility()
    binder = scalar_binder_cumulant()

    nonlinear_avgs = np.zeros((avg_data.shape[0], 1 + num_obs))
    nonlinear_stderr = np.zeros(nonlinear_avgs.shape)

    # Set the temperature
    nonlinear_avgs[:,T_idx] = avg_data[:,T_idx]
    nonlinear_stderr[:,T_idx] = stderr_data[:,T_idx]

    # Calculate the susceptibility and error
    nonlinear_avgs[:,1] = susc.avg( avg_data[:, m_idx], avg_data[:, m2_idx], Nfloat, avg_data[:, T_idx] )
    nonlinear_stderr[:,1] = susc.stderr( stderr_data[:, m_idx], stderr_data[:, m2_idx], avg_data[:, m_idx], Nfloat, avg_data[:, T_idx] )

    # Calculate the Binder cumulant and error
    nonlinear_avgs[:,2] = binder.avg( avg_data[:, m2_idx], avg_data[:, m4_idx], Nfloat )
    nonlinear_stderr[:,2] = binder.stderr( stderr_data[:, m2_idx], stderr_data[:, m4_idx], avg_data[:, m2_idx], avg_data[:, m4_idx], Nfloat )

    return nonlinear_labels, nonlinear_avgs, nonlinear_stderr

def calculate_Ashkin_Teller_observables( labels, Lsize, avg_data, stderr_data, dim = 2 ):

    Nfloat = float(Lsize) ** dim
    T_idx = labels.index("Temperature")
    sig_idx = labels.index("Sigma Mag")
    sig2_idx = labels.index("Sigma Mag2")
    sig4_idx = labels.index("Sigma Mag4")

    tau_idx = labels.index("Tau Mag")
    tau2_idx = labels.index("Tau Mag2")
    tau4_idx = labels.index("Tau Mag4")

    ord_idx = labels.index("Order Parameter")
    ord2_idx = labels.index("Order Parameter2")
    ord4_idx = labels.index("Order Parameter4")

    nem_idx = labels.index("Nematicity")
    nem2_idx = labels.index("Nematicity2")
    nem4_idx = labels.index("Nematicity4")

    nonlinear_labels = [ "1: Temperature",
                         "2: Sigma Susceptibility", "3: Sigma Binder Cumulant",
                         "4: Tau Susceptibility", "5: Tau Binder Cumulant",
                         "6: Nematicity Susceptibility", "7: Nematicity Binder Cumulant",
                         "8: Susceptibility", "9: Binder Cumulant" ]


    num_obs = len(nonlinear_labels) - 1

    susc = susceptibility()
    scalar_binder = scalar_binder_cumulant()
    spinor_binder = two_component_binder_cumulant()

    nonlinear_avgs = np.zeros((avg_data.shape[0], 1 + num_obs))
    nonlinear_stderr = np.zeros(nonlinear_avgs.shape)

    # Set the temperature
    nonlinear_avgs[:,T_idx] = avg_data[:,T_idx]
    nonlinear_stderr[:,T_idx] = stderr_data[:,T_idx]

    # Calculate the sigma susceptibility and error
    nonlinear_avgs[:,1] = susc.avg( avg_data[:, sig_idx], avg_data[:, sig2_idx], Nfloat, avg_data[:, T_idx] )
    nonlinear_stderr[:,1] = susc.stderr( stderr_data[:, sig_idx], stderr_data[:, sig2_idx], avg_data[:, sig_idx], Nfloat, avg_data[:, T_idx] )

    # Calculate the sigma Binder cumulant and error
    nonlinear_avgs[:,2] = scalar_binder.avg( avg_data[:, sig2_idx], avg_data[:, sig4_idx], Nfloat )
    nonlinear_stderr[:,2] = scalar_binder.stderr( stderr_data[:, sig2_idx], stderr_data[:, sig4_idx], avg_data[:, sig2_idx], avg_data[:, sig4_idx], Nfloat )

    # Calculate the tau susceptibility and error
    nonlinear_avgs[:,3] = susc.avg( avg_data[:, tau_idx], avg_data[:, tau2_idx], Nfloat, avg_data[:, T_idx] )
    nonlinear_stderr[:,3] = susc.stderr( stderr_data[:, tau_idx], stderr_data[:, tau2_idx], avg_data[:, tau_idx], Nfloat, avg_data[:, T_idx] )

    # Calculate the tau scalar_binder cumulant and error
    nonlinear_avgs[:,4] = scalar_binder.avg( avg_data[:, tau2_idx], avg_data[:, tau4_idx], Nfloat )
    nonlinear_stderr[:,4] = scalar_binder.stderr( stderr_data[:, tau2_idx], stderr_data[:, tau4_idx], avg_data[:, tau2_idx], avg_data[:, tau4_idx], Nfloat )

    # Calculate the nem susceptibility and error
    nonlinear_avgs[:,5] = susc.avg( avg_data[:, nem_idx], avg_data[:, nem2_idx], Nfloat, avg_data[:, T_idx] )
    nonlinear_stderr[:,5] = susc.stderr( stderr_data[:, nem_idx], stderr_data[:, nem2_idx], avg_data[:, nem_idx], Nfloat, avg_data[:, T_idx] )

    # Calculate the nem scalar_binder cumulant and error
    nonlinear_avgs[:,6] = scalar_binder.avg( avg_data[:, nem2_idx], avg_data[:, nem4_idx], Nfloat )
    nonlinear_stderr[:,6] = scalar_binder.stderr( stderr_data[:, nem2_idx], stderr_data[:, nem4_idx], avg_data[:, nem2_idx], avg_data[:, nem4_idx], Nfloat )

    # Calculate the order susceptibility and error
    nonlinear_avgs[:,7] = susc.avg( avg_data[:, ord_idx], avg_data[:, ord2_idx], Nfloat, avg_data[:, T_idx] )
    nonlinear_stderr[:,7] = susc.stderr( stderr_data[:, ord_idx], stderr_data[:, ord2_idx], avg_data[:, ord_idx], Nfloat, avg_data[:, T_idx] )

    # Calculate the order Binder cumulant and error
    nonlinear_avgs[:,8] = spinor_binder.avg( avg_data[:, ord2_idx], avg_data[:, ord4_idx], Nfloat )
    nonlinear_stderr[:,8] = spinor_binder.stderr( stderr_data[:, ord2_idx], stderr_data[:, ord4_idx], avg_data[:, ord2_idx], avg_data[:, ord4_idx], Nfloat )

    return nonlinear_labels, nonlinear_avgs, nonlinear_stderr

def write_output_data( input_file, nonlinear_labels, header_lines, nonlinear_avgs, nonlinear_stderr ):

    output_path = Path( input_file ).parent

    rewl_substring = input_file[ input_file.find("REWL_") : input_file.find(".job") ]

    full_output = str(output_path) + "/nonlinear_observables-" + rewl_substring

    header_string = ""
    for ldx in range(0, len(header_lines)):
        header_string += header_lines[ldx].rstrip()
        if ldx != len(header_lines) - 1:
            header_string += "\n"
    header_string += "\n# Intensive Nonlinear Observable Names by Column"
    for label in nonlinear_labels:
        header_string += "\n#    " + label
    header_string += "\n#"

    np.savetxt( full_output + ".job_mean", nonlinear_avgs, delimiter = "  ", newline="\n", header=header_string, comments="" )
    np.savetxt( full_output + ".job_stderr", nonlinear_stderr, delimiter = "  ", newline="\n", header=header_string, comments="" )

    return None


def main():

    args = setup_args()

    labels, header_lines, avg_data, stderr_data = read_in_data( args.averages, args.stderr )

    nonlinear_labels, nonlinear_avgs, nonlinear_stderr = None, None, None

    if "Ising" in args.model:

        nonlinear_labels, nonlinear_avgs, nonlinear_stderr = calculate_Ising_observables( labels, args.Lsize, avg_data, stderr_data )

    elif "Ashkin_Teller" in args.model:

        nonlinear_labels, nonlinear_avgs, nonlinear_stderr = calculate_Ashkin_Teller_observables( labels, args.Lsize, avg_data, stderr_data )


    write_output_data( args.averages, nonlinear_labels, header_lines, nonlinear_avgs, nonlinear_stderr )

    return 0

if __name__ == "__main__":
    main()
