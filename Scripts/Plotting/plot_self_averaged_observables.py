# Plot the self-averaged observables
# for many system sizes at some fixed
# coupling.

import argparse
import numpy as np
from parse_file_header import collect_labels
import matplotlib as mpl
mpl.use('Agg')  # THIS IS REQUIRED FOR WSL2
import matplotlib.pyplot as plt
import os

data_file_stem = "self_averaged_observables"
observable_marker = "Intensive Observable"
output_path = "Figures"

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_name", help = "Model being studied", type = str)
    parser.add_argument("coupling_symbol", help = "The constant value to parse through", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling", type = str)

    return parser.parse_args()

def check_for_output():

    if not os.path.isdir( os.getcwd() + "/" + output_path ):
        os.mkdir( os.getcwd() + "/" + output_path )

def find_string_value( string_type, file_string ):

    start = file_string.find(string_type + "-") + len(string_type + "-")
    end = start + file_string[start:].find("_")

    return file_string[start:end]


def collect_observables_and_data(coupling_symbol, coupling_value, comment = "#"):

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

def plot_data_tuples( model_name, coupling_string, coupling_value, labels, data_tuples ):

    xlabel = labels[0]
    key_string = coupling_string + " = " + "%.3f" % float(coupling_value)

    for lbl in range(1, len(labels)):

        print("\nPlotting %s vs %s" % (labels[lbl], xlabel))
        fig, ax = plt.subplots(1,1)

        for Ldx in range(0, len(data_tuples)):

            Lvalue = data_tuples[Ldx][0]
            ax.plot(data_tuples[Ldx][1][:,0], data_tuples[Ldx][1][:,lbl], label = "L = %s" % Lvalue)

        ax.set_xlabel(xlabel, fontsize = 12)
        ax.set_ylabel(labels[lbl], fontsize = 12)
        ax.set_title(model_name + ": " + key_string, fontsize = 12)
        ax.legend(fontsize = 12)

        plotname = "%s" % (labels[lbl] + "_vs_" + xlabel + ".png")

        plt.savefig(output_path + "/" + plotname)

        plt.close()

def main():

    args = setup_args()

    check_for_output()

    labels, data_tuples = collect_observables_and_data( args.coupling_symbol, args.coupling_value )

    plot_data_tuples( args.model_name, args.coupling_symbol, args.coupling_value, labels, data_tuples )



    return None

if __name__ == "__main__":
    main()
