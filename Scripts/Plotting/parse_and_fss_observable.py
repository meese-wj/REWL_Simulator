# Plot the finite size scaling of one
# specific observable.

import argparse
import numpy as np
from parse_file_header import collect_labels
import matplotlib as mpl
mpl.use('Agg')  # THIS IS REQUIRED FOR WSL2
import matplotlib.pyplot as plt
import os

output_path = "Figures"
dimension = 2

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_name", help = "Model being studied", type = str)
    parser.add_argument("fss_observable", help = "Which observable to undergo finite size scaling", type = str)
    parser.add_argument("data_file_stem", help = "Observable type being plotted", type = str)
    parser.add_argument("observable_marker", help = "Line in data header before observable labels", type = str)
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


def collect_observables_and_data( data_file_stem, fss_observable, observable_marker, coupling_symbol, coupling_value, comment = "#"):

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

    if fss_observable not in labels:
        print("\nERROR. %s is not found in the observable labels. Exiting.\n", fss_observable)
        return None, None

    data_tuples.sort(key = lambda tup: int(tup[0]) )

    return labels, data_tuples

def plot_data_tuples( model_name, fss_observable, coupling_string, coupling_value, labels, data_tuples ):

    xlabel = "System Size L"
    key_string = coupling_string + " = " + "%.3f" % float(coupling_value)

    # Find the index of the FSS variable
    lbl = labels.index( fss_observable )

    # Now store the max values in a list
    Lvalues = []
    max_obs = []
    for Ldx in range(0, len(data_tuples)):
        Lvalues.append( float( data_tuples[Ldx][0] ) )
        max_obs.append( np.max( data_tuples[Ldx][1][:,lbl] ) )

    Lvalues = np.array(Lvalues)
    max_obs = np.array(max_obs)

    print("\nPlotting FSS max{ %s } vs %s" % (labels[lbl], xlabel))
    fig, ax = plt.subplots(1,1)

    # Now fit on a loglog plot
    fit_scaling = np.polyfit( np.log(Lvalues), np.log(max_obs), 1 )
    fit_predictions = np.poly1d( fit_scaling )

    label_string = ""
    if fss_observable == "Specific Heat":
        label_string = "$c_V \sim t^{-\\alpha},\; \\alpha = %.3f$" % fit_scaling[0]
        print("\nSpecific Heat alpha = %.3f\n" % fit_scaling[0])
    elif fss_observable == "Susceptibility":
        label_string = "$\chi \sim t^{-\gamma},\; \gamma = %.3f$" % fit_scaling[0]
        print("\nSusceptibility gamma = %.3f\n" % fit_scaling[0])

    # Finally Plot the results
    ax.scatter( Lvalues, max_obs, label = None )
    ax.plot( Lvalues, np.exp( fit_predictions( np.log(Lvalues) ) ), color = "red", ls = "dashed", label = r"FSS Fit: %s" % label_string )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(xlabel, fontsize = 12)
    ax.set_ylabel("Peak %s" % labels[lbl], fontsize = 12)
    ax.set_title(model_name + ": " + key_string, fontsize = 12)
    ax.legend(fontsize = 12)

    plotname = "%s" % ("Peak " + labels[lbl] + " vs " + xlabel + ".png")

    plt.savefig(output_path + "/" + plotname)

    plt.close()

def main():

    args = setup_args()

    check_for_output()

    labels, data_tuples = collect_observables_and_data( args.data_file_stem, args.fss_observable, args.observable_marker, args.coupling_symbol, args.coupling_value )

    if labels == None or data_tuples == None:
        return None

    plot_data_tuples( args.model_name, args.fss_observable, args.coupling_symbol, args.coupling_value, labels, data_tuples )



    return None

if __name__ == "__main__":
    main()
