'''
This file is to plot the splitting
in the pseudo-transition temperatures
for the Baxter susceptibility and
other susceptibilities.
'''

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')  # THIS IS REQUIRED FOR WSL2
import matplotlib.pyplot as plt
from parse_file_header import collect_labels
from parse_and_plot_observables import collect_observables_and_data
from multiple_couplings import *
from bisect import bisect_left

output_path = "Figures"
baxter_string = "Baxter"
susc_string = "Susceptibility"
data_file_stem = "nonlinear_observables"
observable_marker = "Intensive Nonlinear Observable"

Tmin_value = 0.25
Tmax_value = 4

capsize = 2
markersize = 5
marker_edge_color = "black"
marker_edge_width = 1


def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("coupling_symbol", help = "The constant value(s) to parse through. May be a space separated list.", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling(s). May be a space separated list.", type = str)
    parser.add_argument("Tc", default = None, help = "Value of the clean Tc for vertical lines", type = str)

    return parser.parse_args()

def gather_susc_labels( susc_str, labels ):
    '''
    Search through the data labels and find
    all that have Susceptibilty in them
    '''
    Tindex = -1
    susc_labels = []
    for ldx in range(len(labels)):
        if "Temperature" in labels[ldx]:
            Tindex = ldx
        elif susc_str in labels[ldx]:
            susc_labels.append( (ldx, labels[ldx]) )
    return Tindex, susc_labels

def find_pseudo_Tc( Tmin, Tmax, temperatures, susceptibility ):
    '''
    Find the peak temperature in the susceptibility.
    Return the index and the temperature of the maximum.
    '''
    Tmin_idx = bisect_left( temperatures, Tmin )
    Tmax_idx = bisect_left( temperatures, Tmax )
    max_index = Tmin_idx + np.argmax(susceptibility[Tmin_idx:Tmax_idx])
    return (max_index, temperatures[max_index])

def find_all_Tc( T_label_index, Tmin, Tmax, susc_labels, data_tuples ):
    '''
    Iterate through the data tuples and extract all
    the pseudo Tcs and susceptibility maxima.

    susc_labels is a list of tuples of the label index
    and the label string
    '''
    sifter_values = np.zeros((1, len(data_tuples)))
    susc_values = np.zeros((len(susc_labels), len(data_tuples)))
    susc_errs = np.zeros((len(susc_labels), len(data_tuples)))
    Tc_values = np.zeros((len(susc_labels), len(data_tuples)))

    for idx in range(len(data_tuples)):
        sifter_values[0,idx] = float( data_tuples[idx][0] )

        for ldx  in range(len(susc_labels)):
            index, Tc_val = find_pseudo_Tc( Tmin, Tmax, data_tuples[idx][1][:, T_label_index], data_tuples[idx][1][:, susc_labels[ldx][0] ] )
            susc_values[ldx, idx] = data_tuples[idx][1][index, susc_labels[ldx][0]]
            susc_errs[ldx, idx] = data_tuples[idx][2][index, susc_labels[ldx][0]]
            Tc_values[ldx, idx] = Tc_val
    print(Tc_values)
    return sifter_values, susc_values, susc_errs, Tc_values

def plot_Tc_data( plot_directory, sifter_coupling, clean_Tc, coupling_tuples, susc_labels, sifter_values, Tc_values ):
    '''
    Plot the Tc values for all susceptibilities.
    '''
    fig, ax = plt.subplots(1,1)
    for ldx in range(len(susc_labels)):
        label = susc_labels[ldx][1]
        if label == "Susceptibility":
            label = "Two-Component " + label
        ax.plot( Tc_values[ldx,:], sifter_values[0,:], label = r"%s" % label,
                 marker = "o", ms = markersize, mec = marker_edge_color,
                 mew = marker_edge_width, ls = None )

    ymin, ymax = ax.get_ylim()
    const_h_values = np.linspace(ymin, ymax, 10)
    const_Tc = float(clean_Tc) + 0. * const_h_values
    ax.plot( const_Tc, const_h_values, color = "gray", lw = 1, ls = "dashed", label = r"Clean $T_c = %s$" % clean_Tc )
    ax.set_ylim( (ymin, ymax) )

    if sifter_coupling == "h":
        ax.set_ylabel(r"Random Field Strength $%s$" % sifter_coupling)
    ax.set_xlabel(r"Susceptibility Pseudo-$T_c$")
    ax.legend()
    ax.set_title(r"Ashkin_Teller2d_RFAT_Baxter%s" % latex_couplings(coupling_tuples) )

    plotname = "Susceptibility Pseudo-Tc Splitting.png"
    plt.savefig( plot_directory + "/" + plotname )
    plt.close()

    return

def plot_susc_data( plot_directory, sifter_coupling, clean_Tc, coupling_tuples, susc_labels, sifter_values, susc_values, susc_errs ):
    '''
    Plot the max susceptibilities values for all susceptibilities.
    '''
    fig, ax = plt.subplots(1,1)
    for ldx in range(len(susc_labels)):
        label = susc_labels[ldx][1]
        if label == "Susceptibility":
            label = "Two-Component " + label
        lines = ax.plot( sifter_values[0,:], susc_values[ldx,:], label = r"%s" % label,
                         marker = "o", ms = markersize, mec = marker_edge_color,
                         mew = marker_edge_width, ls = None )
        ax.errorbar( sifter_values[0,:], susc_values[ldx,:], yerr = susc_errs[ldx,:],
                     label = None, color = "None", ecolor = lines[-1].get_color(),
                     marker = "o", ms = markersize, mec = marker_edge_color,
                     mew = marker_edge_width, mfc = lines[-1].get_color(),
                     ls = None, capsize=capsize )

    if sifter_coupling == "h":
        ax.set_xlabel(r"Random Field Strength $%s$" % sifter_coupling)
    ax.set_ylabel(r"Peak Susceptibility")
    ax.legend()
    ax.set_title(r"Ashkin_Teller2d_RFAT_Baxter%s" % latex_couplings(coupling_tuples) )

    plotname = "Peak Susceptibility Splitting.png"
    plt.savefig( plot_directory + "/" + plotname )
    plt.close()

    return

def main():

    args = setup_args()

    sifter_coupling, coupling_tuples = parse_couplings( args.coupling_symbol, args.coupling_value )

    # Exit the program if the coupling tuples are null
    if len(coupling_tuples) == 0:
        return
    if sifter_coupling != "h":
        print("\nSifter coupling is not h. Exiting.\n")
        return

    plot_directory = check_for_output(sifter_coupling, coupling_tuples, output_path)

    labels, data_tuples = collect_observables_and_data( data_file_stem, observable_marker, sifter_coupling, coupling_tuples )

    Tindex, susc_labels = gather_susc_labels( susc_string, labels )
    print(susc_labels,"\n", len(data_tuples))

    sifter_values, susc_values, susc_errs, Tc_values = find_all_Tc( Tindex, Tmin_value, Tmax_value, susc_labels, data_tuples )

    plot_Tc_data( plot_directory, sifter_coupling, args.Tc, coupling_tuples, susc_labels, sifter_values, Tc_values )
    plot_susc_data( plot_directory, sifter_coupling, args.Tc, coupling_tuples, susc_labels, sifter_values, susc_values, susc_errs )


    return

if __name__=="__main__":
    main()
