# Plot the observables for some fixed
# couplings over a range of changed
# values (e.g. energy vs temperature
# for may different system sizes).

import argparse
import numpy as np
from parse_file_header import collect_labels
from multiple_couplings import *
import matplotlib as mpl
mpl.use('Agg')  # THIS IS REQUIRED FOR WSL2
import matplotlib.pyplot as plt
import os
from bisect import bisect_left

#data_file_stem = "self_averaged_observables"
#observable_marker = "Intensive Observable"
output_path = "Figures"
error_incrementer = 60
capsize = 2
markersize = 5
marker_edge_color = "black"
marker_edge_width = 1

def setup_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_name", help = "Model being studied", type = str)
    parser.add_argument("data_file_stem", help = "Observable type being plotted", type = str)
    parser.add_argument("observable_marker", help = "Line in data header before observable labels", type = str)
    parser.add_argument("coupling_symbol", help = "The constant value(s) to parse through. May be a space separated list.", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling(s). May be a space separated list.", type = str)
    parser.add_argument("--Tc", default = None, help = "Optional value of Tc for vertical lines", type = str)

    return parser.parse_args()

def collect_observables_and_data( data_file_stem, observable_marker, sifter_coupling, coupling_tuples, comment = "#"):

    labels = []

    data_tuples = []

    coupling_string = get_coupling_string( coupling_tuples, isfloat = True )
    print("\nCoupling String:", coupling_string, data_file_stem, coupling_tuples)

    for fl in os.listdir( os.getcwd() ):
        if not os.path.isdir( fl ) and ( data_file_stem in fl and couplings_in_file(coupling_tuples, fl) and "stderr" not in fl ):

            if len(labels) == 0:
                labels = collect_labels( fl, observable_marker, comment )
                print(labels)

            sifter = find_string_value(sifter_coupling, fl)
            print(sifter)

            data = np.loadtxt(fl, delimiter = "  ", dtype="float64", comments = comment)
            errs = np.zeros((0,0))

            if "job_mean" in fl:
                err_file = fl[ : fl.find("job_mean") ] + "job_stderr"
                errs = np.loadtxt( err_file, delimiter = "  ", dtype="float64", comments = comment )

            print(fl, sifter)

            data_tuples.append( (sifter, data, errs) )

    print(labels)

    data_tuples.sort(key = lambda tup: float(tup[0]) )

    return labels, data_tuples

# Find successive crossing temperature from Binder cumulants
# or correlation lengths
def crossing_temperatures( model_name, data_file_stem, coupling_tuples, label, label_idx, data_tuples, plot_directory, Tmin, Tmax ):

    possible_labels = [ "Binder", "Correlation Length" ]

    if len(data_tuples) < 3:
        return

    perform_analysis = False
    for pl in possible_labels:
        if pl in label:
            perform_analysis = True

    if not perform_analysis:
        return

    Lvalues = []
    Crossings = []

    for ldx in range(len(data_tuples) - 1):

        Lvalues.append( data_tuples[ldx][0] )
        # Find the closest indices of Tmin and Tmax
        Tmin_idx = bisect_left( data_tuples[ldx][1][:,0], Tmin )
        Tmax_idx = bisect_left( data_tuples[ldx][1][:,0], Tmax )
        # Calculate the distance between the data
        # at L and the next largest L
        # diff_values = np.abs( data_tuples[ldx + 1][1][Tmin_idx:Tmax_idx, label_idx] - data_tuples[ldx][1][Tmin_idx:Tmax_idx, label_idx] )
        diff_values = data_tuples[ldx + 1][1][Tmin_idx:Tmax_idx, label_idx] - data_tuples[ldx][1][Tmin_idx:Tmax_idx, label_idx]
        sign_values = np.sign(diff_values)
        sign_idx = 0
        for idx in range(1, len(sign_values)):
            if sign_values[idx] != sign_values[idx -1]:
                sign_idx = idx
                break

        Crossings.append( data_tuples[ldx][1][ Tmin_idx + idx, 0] )

    Lvalues = np.array(Lvalues)
    Crossings = np.array(Crossings)

    # Now plot the successive crossings vs L
    fig, ax = plt.subplots(1,1)
    ax.plot( Lvalues, Crossings,
             color = "None", marker = "o", mfc = "C0",
             ms = markersize, mec = marker_edge_color, mew = marker_edge_width )


    ax.set_xlabel(r"System Size $L$", fontsize = 12)
    ax.set_ylabel(r"%s Crossing Temperature" % label, fontsize = 10)

    midpoint = 0.5 * ( Tmin + Tmax )
    ax.set_ylim([(1-0.1)*Tmin, (1+0.05)*Tmax])

    plottitle = model_titles(model_name) + latex_couplings(coupling_tuples) + " with $T \in [%.3f, %.3f]$" % (Tmin, Tmax)
    ax.set_title(r"%s" % plottitle, fontsize = 12)

    plotname = label + " Crossing Temperatures.png"
    plt.savefig( plot_directory + "/" + plotname )

    plt.close()

# Plot the probability density as functions of energy at a fixed temperature
def plot_probability_density( model_name, data_file_stem, coupling_tuples, labels, data_tuples, plot_directory, Tc_val = None, plt_err = False ):

    xlabel = labels[0]
    key_string = latex_couplings(coupling_tuples)

    auxiliary_microcanonical = "Order Parameter2"

    if auxiliary_microcanonical in labels:

        fig, ax = plt.subplots(1,2, sharey=True)

        for Ldx in range(0, len(data_tuples)):

            sifter = data_tuples[Ldx][0]
            sFloat = float(sifter)

            # TODO: This will break in higher dimensions
            Nfloat = sFloat ** 2

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
            # density = np.exp( exponent ) / partition
            density = np.exp( exponent )

            # Rescale it by N to plot the normalized
            # density along the intensive energy axis
            lines1 = ax[1].plot( data_tuples[Ldx][1][:,0], density, label = r"$L = %s$" % sifter )

            order2_col = labels.index("Order Parameter2")
            #ax[0].plot( data_tuples[Ldx][1][:,order2_col] / Nfloat, Nfloat * density, label = r"$L = %s$" % Lvalue )
            lines0 = ax[0].plot( data_tuples[Ldx][1][:,order2_col] / Nfloat, density, label = r"$L = %s$" % sifter )

            if plt_err and data_tuples[Ldx][2].shape != (0,0):
                energy_err = Nfloat * data_tuples[Ldx][2][:,0]
                logdos_err = Nfloat * data_tuples[Ldx][2][:,1]
                # Ignore the error in the max density since 1/density(max_exponent) -> 0
                density_err = density * np.sqrt( logdos_err ** 2. + ( energy_err / float(Tc_val) ) ** 2. )
                ax[0].errorbar(data_tuples[Ldx][1][::error_incrementer, order2_col] / Nfloat, density[::error_incrementer],
                               xerr = data_tuples[Ldx][2][::error_incrementer, 0] / Nfloat, yerr = density_err[::error_incrementer],
                               ms = markersize, mec = marker_edge_color, mew = marker_edge_width, mfc = lines[-1].get_color(),
                               color = "None", marker = "o", ecolor = lines[-1].get_color(),
                               fmt='none', capsize=capsize, label = None)
                ax[1].errorbar(data_tuples[Ldx][1][::error_incrementer, 0], density[::error_incrementer],
                               xerr = data_tuples[Ldx][2][::error_incrementer, 0], yerr = density_err[::error_incrementer], errorevery=error_incrementer,
                               color = "None", marker = "o", ecolor = lines[-1].get_color(),
                               ms = markersize, mec = marker_edge_color, mew = marker_edge_width, mfc = lines[-1].get_color(),
                               fmt='none', capsize=capsize, label = None)

        # Set xlim
        ax[1].set_xlim([-5,-3])

        ax[1].set_xlabel(r"Energy per Site $[E/N]$", fontsize = 12)
        ax[0].set_xlabel(r"%s" % labels[order2_col], fontsize = 12)
        #ax[0].set_ylabel(r"Probability Density $[N\cdot E^{-1}]$", fontsize = 12)
        ax[0].set_ylabel(r"Scaled Probability Density $[E^{-1}]$", fontsize = 12)
        ax[0].legend(fontsize = 10)

        if "RF" in model_name:
            fig.suptitle(model_titles(model_name) + " at $T_c^{(0)} = %s$: " % (Tc_val) + key_string, fontsize = 12)
        else:
            fig.suptitle(model_titles(model_name) + " at $T_c = %s$: " % (Tc_val) + key_string, fontsize = 12)

        plotname = "%s" % ("probability density vs energy.png")
        plt.savefig(plot_directory + "/" + plotname)

    else:
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
            partition = np.sum( np.exp( exponent ), dtype="float64" )

            # Find the probability density as a function
            # of the extensive energy
            # density = np.exp( exponent ) / partition
            density = np.exp( exponent )

            # Rescale it by N to plot the normalized
            # density along the intensive energy axis
            # ax.plot( data_tuples[Ldx][1][:,0], Nfloat * density, label = r"$L = %s$" % Lvalue )
            lines = ax.plot( data_tuples[Ldx][1][:,0], density, label = r"$L = %s$" % Lvalue )

            if plt_err and data_tuples[Ldx][2].shape != (0,0):
                energy_err = Nfloat * data_tuples[Ldx][2][:,0]
                logdos_err = Nfloat * data_tuples[Ldx][2][:,1]
                # Ignore the error in the max density since 1/density(max_exponent) -> 0
                density_err = density * np.sqrt( logdos_err ** 2. + ( energy_err / float(Tc_val) ) ** 2. )
                ax.errorbar(data_tuples[Ldx][1][::error_incrementer, 0], density[::error_incrementer],
                            xerr = data_tuples[Ldx][2][::error_incrementer, 0], yerr = density_err[::error_incrementer],
                            color = "None", marker = "o", ecolor = lines[-1].get_color(),
                            ms = markersize, mec = marker_edge_color, mew = marker_edge_width, mfc = lines[-1].get_color(),
                            fmt='none', capsize=capsize, label = None)

        # Set xlim
        ax.set_xlim([-2.0,-1.25])

        ax.set_xlabel(r"Energy per Site $[E/N]$", fontsize = 12)
        #ax.set_ylabel(r"Probability Density $[N\cdot E^{-1}]$", fontsize = 12)
        ax.set_ylabel(r"Scaled Probability Density $[E^{-1}]$", fontsize = 12)
        ax.legend(fontsize = 10)

        if "RF" in model_name:
            fig.suptitle(model_titles(model_name) + " at $T_c^{(0)} = %s$: " % (Tc_val) + key_string, fontsize = 12)
        else:
            fig.suptitle(model_titles(model_name) + " at $T_c = %s$: " % (Tc_val) + key_string, fontsize = 12)

        plotname = "%s" % ("probability density vs energy.png")
        plt.savefig(plot_directory + "/" + plotname)



    plt.close()


def get_y_range( yvalues, xvalues, xmin, xmax, extension=0.05 ):

    # Find the values of x that are closest to xmin and xmax
    # xvalues must be sorted!
    xindex_min = bisect_left( xvalues, xmin )
    xindex_max = bisect_left( xvalues, xmax )

    # First prune the yvalues for any nans which break Python
    # prune_yvalues = yvalues[ np.isnan(yvalues) != True ]

    # Now find the min and max of y on this range
    ymin = np.nanmin( yvalues[xindex_min:xindex_max] )
    ymax = np.nanmax( yvalues[xindex_min:xindex_max] )

    # Set the ylimits to be some fractional extension of
    # the yrange
    yrange = ymax - ymin
    return ymin - extension * yrange, ymax + extension * yrange

def plot_data_tuples( model_name, data_file_stem, sifter_coupling, coupling_tuples, labels, data_tuples, plot_directory, Tc_val = None ):

    xlabel = labels[0]
    key_string = latex_couplings(coupling_tuples)

    for lbl in range(1, len(labels)):

        print("\nPlotting %s vs %s" % (labels[lbl], xlabel))
        fig, ax = plt.subplots(1,1)

        epsilon_range = .1
        #xmin, xmax, plt_ymin, plt_ymax = 0, 0, 0, 0
        xmin, xmax, plt_ymin, plt_ymax = None, None, None, None
        if epsilon_range != None and Tc_val != None and Tc_val != "":
            xmin, xmax = (1 - epsilon_range) * float(Tc_val), (1 + epsilon_range) * float(Tc_val)
            xmin, xmax = 1.5, 1.25 * float(Tc_val)

        for Ldx in range(0, len(data_tuples)):

            sifter = data_tuples[Ldx][0]
            lines = ax.plot(data_tuples[Ldx][1][:,0], data_tuples[Ldx][1][:,lbl], label = r"%s" % pretty_label_string(sifter_coupling, sifter))

            if data_tuples[Ldx][2].shape != (0,0):
                ax.errorbar(data_tuples[Ldx][1][::error_incrementer, 0], data_tuples[Ldx][1][::error_incrementer, lbl],
                            xerr = data_tuples[Ldx][2][::error_incrementer, 0], yerr = data_tuples[Ldx][2][::error_incrementer, lbl],
                            color = "None", marker = "o", ecolor = lines[-1].get_color(),
                            ms = markersize, mec = marker_edge_color, mew = marker_edge_width, mfc = lines[-1].get_color(),
                            ls=None, capsize=capsize, label = None)

            if epsilon_range != None and Tc_val != None and Tc_val != "":
                if "microcanonical" not in data_file_stem and "Counts" not in labels[lbl]:
                    test_ymin, test_ymax = get_y_range( data_tuples[Ldx][1][:,lbl], data_tuples[Ldx][1][:,0], xmin, xmax )
                    #print("Plot Min = ", plt_ymin, "  Test Min = ", test_ymin)
                    #print("Plot Max = ", plt_ymax, "  Test Max = ", test_ymax)
                    #if abs(test_ymin) > abs(plt_ymin):
                    if plt_ymin == None or test_ymin < plt_ymin:
                        plt_ymin = test_ymin
                    #if abs(test_ymax) > abs(plt_ymax):
                    if plt_ymax == None or test_ymax > plt_ymax:
                        plt_ymax = test_ymax

        if epsilon_range != None and Tc_val != None and Tc_val != "":
            if "microcanonical" not in data_file_stem and "Counts" not in labels[lbl]:
                ax.set_ylim([plt_ymin, plt_ymax])

        if Tc_val != None and Tc_val != "":
            if "microcanonical" not in data_file_stem:
                ax.set_xlim([0, 2 * float(Tc_val)])
                if epsilon_range != None:
                    ax.set_xlim([xmin, xmax])

                ymin, ymax = ax.get_ylim()
                if "RF" in model_name:
                    ax.plot( float(Tc_val) + 0. * np.linspace(0,1,10), ymin + (ymax - ymin) * np.linspace(0,1,10), color = "gray", lw = 1, ls = "dashed", label = r"$T_c^{(0)} = %s$" % Tc_val )
                else:
                    ax.plot( float(Tc_val) + 0. * np.linspace(0,1,10), ymin + (ymax - ymin) * np.linspace(0,1,10), color = "gray", lw = 1, ls = "dashed", label = r"$T_c = %s$" % Tc_val )
                ax.set_ylim([ymin, ymax])
            elif sifter_coupling == "L" and "microcanonical" in data_file_stem and lbl == 1:
                plot_probability_density( model_name, data_file_stem, coupling_tuples, labels, data_tuples, plot_directory, Tc_val )

        if sifter_coupling == "L" and xmin != None and xmax != None and "nonlinear" in data_file_stem:
                crossing_temperatures( model_name, data_file_stem, coupling_tuples, labels[lbl], lbl, data_tuples, plot_directory, (1-0.1)*float(Tc_val), (1+0.05)*float(Tc_val) )

        ax.set_xlabel(xlabel, fontsize = 12)
        ax.set_ylabel(labels[lbl], fontsize = 12)
        ax.set_title( model_titles(model_name) + key_string, fontsize = 12)
        ax.legend(fontsize = 12)

        plotname = "%s" % (labels[lbl] + "_vs_" + xlabel + ".png")

        plt.savefig(plot_directory + "/" + plotname)

        plt.close()

def main():

    args = setup_args()

    sifter_coupling, coupling_tuples = parse_couplings( args.coupling_symbol, args.coupling_value )

    # Exit the program if the coupling tuples are null
    if len(coupling_tuples) == 0:
        return

    plot_directory = check_for_output(sifter_coupling, coupling_tuples, output_path)

    labels, data_tuples = collect_observables_and_data( args.data_file_stem, args.observable_marker, sifter_coupling, coupling_tuples )

    plot_data_tuples( args.model_name, args.data_file_stem, sifter_coupling, coupling_tuples, labels, data_tuples, plot_directory, args.Tc )

    return None

if __name__ == "__main__":
    main()
