# Plot the observables
# for many system sizes at some fixed
# coupling.

import argparse
import numpy as np
from parse_file_header import collect_labels
import matplotlib as mpl
mpl.use('Agg')  # THIS IS REQUIRED FOR WSL2
import matplotlib.pyplot as plt
import os

#data_file_stem = "self_averaged_observables"
#observable_marker = "Intensive Observable"
output_path = "Figures"

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_name", help = "Model being studied", type = str)
    parser.add_argument("data_file_stem", help = "Observable type being plotted", type = str)
    parser.add_argument("observable_marker", help = "Line in data header before observable labels", type = str)
    parser.add_argument("coupling_symbol", help = "The constant value to parse through", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling", type = str)
    parser.add_argument("--Tc", default = None, help = "Optional value of Tc for vertical lines", type = str)

    return parser.parse_args()

def check_for_output(coupling_symbol, coupling_value):

    key_string = coupling_symbol + "-" + coupling_value
    if not os.path.isdir( os.getcwd() + "/" + output_path + "_" + key_string ):
        os.mkdir( os.getcwd() + "/" + output_path + "_" + key_string )

    return os.getcwd() + "/" + output_path + "_" + key_string

def find_string_value( string_type, file_string ):

    start = file_string.find(string_type + "-") + len(string_type + "-")
    end = start + file_string[start:].find("_")

    return file_string[start:end]


def collect_observables_and_data( data_file_stem, observable_marker, coupling_symbol, coupling_value, comment = "#"):

    labels = []

    data_tuples = []

    key_string = coupling_symbol + "-" + ("%.6f" % float(coupling_value))
    print("\nKey String:", key_string)

    for fl in os.listdir( os.getcwd() ):
        if not os.path.isdir( fl ) and ( data_file_stem in fl and key_string in fl ):

            if len(labels) == 0:
                labels = collect_labels( fl, observable_marker, comment )

            Lvalue = find_string_value("L", fl)

            data = np.loadtxt(fl, delimiter = "  ", dtype = float, comments = comment)

            print(fl, Lvalue)

            data_tuples.append( (Lvalue, data) )



    print(labels)

    data_tuples.sort(key = lambda tup: int(tup[0]) )

    return labels, data_tuples

# Plot the probability density as functions of energy at a fixed temperature
def plot_probability_density( model_name, data_file_stem, coupling_string, coupling_value, labels, data_tuples, plot_directory, Tc_val = None ):

    xlabel = labels[0]
    key_string = coupling_string + " = " + "%.3f" % float(coupling_value)

    fig, ax = plt.subplots(1,1)

    for Ldx in range(0, len(data_tuples)):

        Lvalue = data_tuples[Ldx][0]
        Lfloat = float(Lvalue)

        # TODO: This will break in higher dimensions
        Nfloat = Lfloat ** 2

        extensive_energy = Nfloat * data_tuples[Ldx][1][:,0]
        extensive_logdos = Nfloat * data_tuples[Ldx][1][:,1]

        # Find the exponent of the probability density
        exponent = extensive_logdos - extensive_energy / float(Tc_val)

        # Normalize it by its maximum
        exponent -= np.max(exponent)

        # Find the partition function
        partition = np.sum( np.exp( exponent ) )

        # Find the probability density as a function
        # of the extensive energy
        density = np.exp( exponent ) / partition


        # Rescale it by N to plot the normalized
        # density along the intensive energy axis
        ax.plot( data_tuples[Ldx][1][:,0], Nfloat * density, label = r"$L = %s$" % Lvalue )

    # Set xlim
    # ax.set_xlim([-1.55,-1.25])

    ax.set_xlabel(r"Energy per Site $[E/N]$", fontsize = 12)
    ax.set_ylabel(r"Probability Density $[N\cdot E^{-1}]$", fontsize = 12)
    ax.legend(fontsize = 10)

    ax.set_title(model_name + " at $T_c = %s$: " % (Tc_val) + key_string, fontsize = 12)

    plotname = "%s" % ("probability density vs energy.png")
    plt.savefig(plot_directory + "/" + plotname)

    plt.close()


def plot_data_tuples( model_name, data_file_stem, coupling_string, coupling_value, labels, data_tuples, plot_directory, Tc_val = None ):

    xlabel = labels[0]
    key_string = coupling_string + " = " + "%.3f" % float(coupling_value)

    for lbl in range(1, len(labels)):

        print("\nPlotting %s vs %s" % (labels[lbl], xlabel))
        fig, ax = plt.subplots(1,1)

        for Ldx in range(0, len(data_tuples)):

            Lvalue = data_tuples[Ldx][0]
            ax.plot(data_tuples[Ldx][1][:,0], data_tuples[Ldx][1][:,lbl], label = r"$L = %s$" % Lvalue)

        ymin, ymax = ax.get_ylim()
        if Tc_val != None and Tc_val != "":
            if "microcanonical" not in data_file_stem:
                ax.set_xlim([0, 2 * float(Tc_val)])
                #epsilon_range = 0.025
                #ax.set_xlim([(1 - epsilon_range) * float(Tc_val), (1 + epsilon_range) * float(Tc_val)])
                ax.plot( float(Tc_val) + 0. * np.linspace(0,1,10), ymin + (ymax - ymin) * np.linspace(0,1,10), color = "gray", lw = 1, ls = "dashed", label = r"$T_c = %s$" % Tc_val )
            elif "microcanonical" in data_file_stem and lbl == 1:
                plot_probability_density( model_name, data_file_stem, coupling_string, coupling_value, labels, data_tuples, plot_directory, Tc_val )

        ax.set_ylim([ymin, ymax])

        ax.set_xlabel(xlabel, fontsize = 12)
        ax.set_ylabel(labels[lbl], fontsize = 12)
        ax.set_title(model_name + ": " + key_string, fontsize = 12)
        ax.legend(fontsize = 12)

        plotname = "%s" % (labels[lbl] + "_vs_" + xlabel + ".png")

        plt.savefig(plot_directory + "/" + plotname)

        plt.close()

def main():

    args = setup_args()

    plot_directory = check_for_output(args.coupling_symbol, args.coupling_value)

    labels, data_tuples = collect_observables_and_data( args.data_file_stem, args.observable_marker, args.coupling_symbol, args.coupling_value )

    plot_data_tuples( args.model_name, args.data_file_stem, args.coupling_symbol, args.coupling_value, labels, data_tuples, plot_directory, args.Tc )



    return None

if __name__ == "__main__":
    main()
