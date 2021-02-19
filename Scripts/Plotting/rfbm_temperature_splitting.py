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
fontsize = 12


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

def find_Tc_error( max_index, temperatures, susceptibility, susc_err ):
    '''
    Define the Tc error as the half width of the peak
    in the susceptibility where the width is taken
    along the line at
        susceptibility(Tc) - susceptibility_error(Tc).
    '''
    lower_susc = susceptibility[max_index] - susc_err[max_index]
    if lower_susc <= 0.:
        # If the lower susceptibility is negative
        # then return the largest reasonable
        # interval.
        return temperatures[max_index] - temperatures[0]

    lower_temp, upper_temp = None, None
    for i in range(1, min(max_index, len(temperatures) - max_index)):
        if susceptibility[max_index - i] - lower_susc < 0. and lower_temp == None:
            lower_temp = temperatures[max_index - i]
        if susceptibility[max_index + i] - lower_susc < 0. and upper_temp == None:
            upper_temp = temperatures[max_index + i]
        if lower_temp != None and upper_temp != None:
            break;
    # Return half the range of temperatures
    # as the error.
    return 0.5 * ( upper_temp - lower_temp )


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
    Tc_errors = np.zeros((len(susc_labels), len(data_tuples)))

    for idx in range(len(data_tuples)):
        sifter_values[0,idx] = float( data_tuples[idx][0] )

        for ldx  in range(len(susc_labels)):
            index, Tc_val = find_pseudo_Tc( Tmin, Tmax, data_tuples[idx][1][:, T_label_index], data_tuples[idx][1][:, susc_labels[ldx][0] ] )
            susc_values[ldx, idx] = data_tuples[idx][1][index, susc_labels[ldx][0]]
            susc_errs[ldx, idx] = data_tuples[idx][2][index, susc_labels[ldx][0]]
            Tc_values[ldx, idx] = Tc_val
            Tc_errors[ldx, idx] = find_Tc_error( index, data_tuples[idx][1][:, T_label_index],
                                                 data_tuples[idx][1][:, susc_labels[ldx][0] ], data_tuples[idx][2][:, susc_labels[ldx][0] ]  )
    print(Tc_values, "\n", Tc_errors)
    return sifter_values, susc_values, susc_errs, Tc_values, Tc_errors

def plot_Tc_data( plot_directory, sifter_coupling, clean_Tc, coupling_tuples, susc_labels, sifter_values, Tc_values, Tc_errors, ashkin_teller ):
    '''
    Plot the Tc values for all susceptibilities.
    '''
    fig, ax = plt.subplots(1,1)
    for ldx in range(len(susc_labels)):
        label = susc_labels[ldx][1]
        if ashkin_teller and label == "Susceptibility":
            label = "Two-Component " + label
        lines = ax.plot( Tc_values[ldx,:], sifter_values[0,:], label = r"%s" % label,
                         marker = "o", ms = markersize, mec = marker_edge_color,
                         mew = marker_edge_width, ls = None )
        ax.errorbar( Tc_values[ldx,:], sifter_values[0,:], xerr = Tc_errors[ldx,:],
                     label = None, color = "None", ecolor = lines[-1].get_color(),
                     marker = "o", ms = markersize, mec = marker_edge_color,
                     mew = marker_edge_width, mfc = lines[-1].get_color(),
                     ls = None, capsize=capsize )

    ymin, ymax = ax.get_ylim()
    const_h_values = np.linspace(ymin, ymax, 10)
    const_Tc = float(clean_Tc) + 0. * const_h_values
    ax.plot( const_Tc, const_h_values, color = "gray", lw = 1, ls = "dashed", label = r"Clean $T_c = %s$" % clean_Tc )
    ax.set_ylim( (ymin, ymax) )

    if sifter_coupling == "h":
        ax.set_ylabel(r"Random Field Strength $%s$" % sifter_coupling, fontsize = fontsize)
    ax.set_xlabel(r"Susceptibility Pseudo-$T_c$", fontsize = fontsize)
    ax.legend(fontsize = fontsize)
    if ashkin_teller:
        ax.set_title(r"Ashkin_Teller2d_RFAT_Baxter%s" % latex_couplings(coupling_tuples), fontsize = fontsize )
    else:
        ax.set_title(r"Ising2d_RFIM%s" % latex_couplings(coupling_tuples), fontsize = fontsize )

    plotname = "Susceptibility Pseudo-Tc Splitting.png"
    plt.savefig( plot_directory + "/" + plotname )
    plt.close()

    return

def plot_Tc_data_inverted( plot_directory, sifter_coupling, clean_Tc, coupling_tuples, susc_labels, sifter_values, Tc_values, Tc_errors, ashkin_teller ):
    '''
    Plot the Tc values for all susceptibilities.
    '''
    fig, ax = plt.subplots(1,1)
    for ldx in range(len(susc_labels)):
        label = susc_labels[ldx][1]
        if ashkin_teller and label == "Susceptibility":
            label = "Two-Component " + label
        lines = ax.plot( sifter_values[0,:], Tc_values[ldx,:], label = r"%s" % label,
                         marker = "o", ms = markersize, mec = marker_edge_color,
                         mew = marker_edge_width, ls = None )
        ax.errorbar( sifter_values[0,:], Tc_values[ldx,:], yerr = Tc_errors[ldx,:],
                     label = None, color = "None", ecolor = lines[-1].get_color(),
                     marker = "o", ms = markersize, mec = marker_edge_color,
                     mew = marker_edge_width, mfc = lines[-1].get_color(),
                     ls = None, capsize=capsize )

    xmin, xmax = ax.get_xlim()
    const_h_values = np.linspace(xmin, xmax, 10)
    const_Tc = float(clean_Tc) + 0. * const_h_values
    ax.plot( const_h_values, const_Tc, color = "gray", lw = 1, ls = "dashed", label = r"Clean $T_c = %s$" % clean_Tc )
    ax.set_xlim( (xmin, xmax) )

    if sifter_coupling == "h":
        ax.set_xlabel(r"Random Field Strength $%s$" % sifter_coupling, fontsize = fontsize)
    ax.set_ylabel(r"Susceptibility Pseudo-$T_c$", fontsize = fontsize)
    ax.legend(fontsize = fontsize)
    if ashkin_teller:
        ax.set_title(r"Ashkin_Teller2d_RFAT_Baxter%s" % latex_couplings(coupling_tuples), fontsize = fontsize )
    else:
        ax.set_title(r"Ising2d_RFIM%s" % latex_couplings(coupling_tuples), fontsize = fontsize )

    plotname = "Susceptibility Pseudo-Tc Splitting inverted.png"
    plt.savefig( plot_directory + "/" + plotname )
    plt.close()

    return



def plot_susc_data( plot_directory, sifter_coupling, clean_Tc, coupling_tuples, susc_labels, sifter_values, susc_values, susc_errs, ashkin_teller ):
    '''
    Plot the max susceptibilities values for all susceptibilities.
    '''
    fig, ax = plt.subplots(1,1)
    for ldx in range(len(susc_labels)):
        label = susc_labels[ldx][1]
        if ashkin_teller and label == "Susceptibility":
            label = "Order Parameter " + label
        lines = ax.plot( sifter_values[0,:], susc_values[ldx,:], label = r"%s" % label,
                         marker = "o", ms = markersize, mec = marker_edge_color,
                         mew = marker_edge_width, ls = None )
        ax.errorbar( sifter_values[0,:], susc_values[ldx,:], yerr = susc_errs[ldx,:],
                     label = None, color = "None", ecolor = lines[-1].get_color(),
                     marker = "o", ms = markersize, mec = marker_edge_color,
                     mew = marker_edge_width, mfc = lines[-1].get_color(),
                     ls = None, capsize=capsize )

    if sifter_coupling == "h":
        ax.set_xlabel(r"Random Field Strength $%s$" % sifter_coupling, fontsize = fontsize )
    ax.set_ylabel(r"Peak Susceptibility", fontsize = fontsize )
    ax.legend(fontsize = fontsize)
    if ashkin_teller:
        ax.set_title(r"Ashkin_Teller2d_RFAT_Baxter%s" % latex_couplings(coupling_tuples), fontsize = fontsize )
    else:
        ax.set_title(r"Ising2d_RFIM%s" % latex_couplings(coupling_tuples), fontsize = fontsize )

    plotname = "Peak Susceptibility Splitting.png"
    plt.savefig( plot_directory + "/" + plotname )
    plt.close()

    return

def plot_anomalous_susc( plot_directory, Tindex, susc_labels, sifter_coupling, clean_Tc, coupling_tuples, labels, data_tuples ):
    '''
    Plot the anomalous diagmetic susceptibility given by
    chi_a = chi_order_parameter - chi_sigma - chi_tau
    '''
    order_idx = 0
    sigma_idx = 0
    tau_idx = 0
    for i in range(len(susc_labels)):
        if "Susceptibility" == susc_labels[i][1]:
            order_idx = susc_labels[i][0]
        elif "Sigma" in susc_labels[i][1]:
            sigma_idx = susc_labels[i][0]
        elif "Tau" in susc_labels[i][1]:
            tau_idx = susc_labels[i][0]

    error_inc = 100
    fig, ax = plt.subplots(1,1)
    for i in range(len(data_tuples)):
        anom_susc_values = data_tuples[i][1][:,order_idx] - ( data_tuples[i][1][:,sigma_idx] + data_tuples[i][1][:,tau_idx] )
        anom_errs = np.sqrt( data_tuples[i][2][:,order_idx] ** 2 +  data_tuples[i][2][:,sigma_idx] ** 2 + data_tuples[i][2][:,tau_idx] ** 2 )
        lines = ax.plot( data_tuples[i][1][:,Tindex], anom_susc_values, label = r"$%s = %.2f$" % (sifter_coupling, float(data_tuples[i][0]) ) )
        ax.errorbar( data_tuples[i][1][::error_inc,Tindex], anom_susc_values[::error_inc], yerr = anom_errs[::error_inc],
                     label = None, color = "None", ecolor = lines[-1].get_color(),
                     marker = "o", ms = markersize, mec = marker_edge_color,
                     mew = marker_edge_width, mfc = lines[-1].get_color(),
                     ls = None, capsize=capsize )



    epsilon_range = .25
    ax.set_xlim( float(clean_Tc) * (1-epsilon_range), float(clean_Tc) * (1 + epsilon_range) )

    ymin, ymax = ax.get_ylim()
    const_Tc = float(clean_Tc) + 0. * np.linspace(ymin, ymax, 10)
    ax.plot( const_Tc, np.linspace(ymin, ymax, 10), color = "gray", lw = 1, ls = "dashed", label = r"Clean $T_c = %s$" % clean_Tc )
    ax.set_ylim( (ymin, ymax) )

    ax.set_xlabel("%s" % labels[Tindex], fontsize = fontsize)
    ax.set_ylabel("Anomalous Susceptibility", fontsize = fontsize)
    ax.set_title(r"Ashkin_Teller2d_RFAT_Baxter%s" % latex_couplings(coupling_tuples), fontsize = fontsize )
    ax.legend()

    plotname = "Anomalous Susceptibility.png"
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
    ashkin_teller = False
    if "Tau Susceptibility" in labels:
        ashkin_teller = True

    Tindex, susc_labels = gather_susc_labels( susc_string, labels )
    print(susc_labels,"\n", len(data_tuples))

    ''' Convert susceptibility into correlator by multiplying by temperature
    for i in range(len(data_tuples)):
        temperatures = data_tuples[i][1][:,Tindex]
        for si in range(len(susc_labels)):
            s = susc_labels[si][0]
            data_tuples[i][1][:,s] = data_tuples[i][1][:,s] * temperatures
            data_tuples[i][2][:,s] = data_tuples[i][2][:,s] * temperatures
    '''

    sifter_values, susc_values, susc_errs, Tc_values, Tc_errors = find_all_Tc( Tindex, Tmin_value, Tmax_value, susc_labels, data_tuples )

    if ashkin_teller:
        plot_anomalous_susc(  plot_directory, Tindex, susc_labels, sifter_coupling, args.Tc, coupling_tuples, labels, data_tuples )

    plot_Tc_data( plot_directory, sifter_coupling, args.Tc, coupling_tuples, susc_labels, sifter_values, Tc_values, Tc_errors, ashkin_teller )
    plot_Tc_data_inverted( plot_directory, sifter_coupling, args.Tc, coupling_tuples, susc_labels, sifter_values, Tc_values, Tc_errors, ashkin_teller )
    plot_susc_data( plot_directory, sifter_coupling, args.Tc, coupling_tuples, susc_labels, sifter_values, susc_values, susc_errs, ashkin_teller )


    return

if __name__=="__main__":
    main()
