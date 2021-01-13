# Plot the finite size scaling of one
# specific observable.

import argparse
import numpy as np
import scipy.optimize as opt
from parse_file_header import collect_labels
import multiple_couplings
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
    parser.add_argument("coupling_symbol", help = "The constant value(s) to parse through. May be a space separated list.", type = str)
    parser.add_argument("coupling_value",  help = "Value of the coupling(s). May be a space separated list.", type = str)

    return parser.parse_args()

def collect_observables_and_data( data_file_stem, fss_observable, observable_marker, coupling_tuples, comment = "#"):

    labels = []

    data_tuples = []

    coupling_string = get_coupling_string( coupling_tuples, isfloat = True )
    print("\nCoupling String:", coupling_string)

    for fl in os.listdir( os.getcwd() ):
        if not os.path.isdir( fl ) and ( data_file_stem in fl and couplings_in_file(coupling_tuples, fl) and "stderr" not in fl ):

            if len(labels) == 0:
                labels = collect_labels( fl, observable_marker, comment )

            Lvalue = find_string_value("L", fl)

            data = np.loadtxt(fl, delimiter = "  ", dtype = float, comments = comment)
            error = np.zeros((0,0))
            if "job_mean" in fl:
                error_string = fl[:fl.find(".job_mean")] + ".job_stderr"
                error = np.loadtxt( error_string, delimiter = "  ", dtype = "float64", comments = comment )

            print(fl, Lvalue)

            data_tuples.append( (Lvalue, data, error) )



    print(labels)

    if fss_observable not in labels:
        print("\nERROR. %s is not found in the observable labels. Exiting.\n", fss_observable)
        return None, None

    data_tuples.sort(key = lambda tup: int(tup[0]) )

    return labels, data_tuples

def plot_data_tuples( model_name, fss_observable, coupling_tuples, labels, data_tuples, plot_directory ):

    xlabel = "System Size L"
    key_string = latex_couplings(coupling_tuples)

    # Find the index of the FSS variable
    Tlbl = labels.index( "Temperature" )
    lbl = labels.index( fss_observable )

    # Now store the max values in a list
    Lvalues = []
    max_obs, max_obs_err = [], []
    Tc_obs = []
    for Ldx in range(0, len(data_tuples)):
        Lvalues.append( float( data_tuples[Ldx][0] ) )
        max_obs.append( np.max( data_tuples[Ldx][1][:,lbl] ) )
        Tc_obs.append( data_tuples[Ldx][1][ np.argmax( data_tuples[Ldx][1][:,lbl] ), Tlbl ] )

        if data_tuples[Ldx][2].shape != (0,0):
            max_obs_err.append(  data_tuples[Ldx][2][ np.argmax( data_tuples[Ldx][1][:,lbl] ), lbl ] )

    Lvalues = np.array(Lvalues)
    max_obs = np.array(max_obs)
    max_obs_err = np.array(max_obs_err)
    Tc_obs = np.array(Tc_obs)

    print("\nPlotting FSS max{ %s } vs %s" % (labels[lbl], xlabel))
    fig, ax = plt.subplots(1,1)

    # Now fit on a loglog plot
    fit_scaling, covariance = None, None
    if data_tuples[Ldx][2].shape != (0,0):
        fit_scaling, covariance = opt.curve_fit( lambda x,a,b: a * x + b, np.log(Lvalues), np.log(max_obs), sigma = max_obs_err )
    else:
        fit_scaling, covariance = opt.curve_fit( lambda x,a,b: a * x + b, np.log(Lvalues), np.log(max_obs) )

    smooth_L = np.linspace(np.min(Lvalues), np.max(Lvalues), 1000)
    fit_predictions = fit_scaling[0] * np.log( smooth_L ) + fit_scaling[1]
    errors = np.sqrt( np.diag( covariance ) )

    label_string = ""
    if fss_observable == "Specific Heat":
        label_string = "$c_V \sim L^{\\alpha/\\nu},\; \\alpha/\\nu = %.3f \pm %.3f$" % (fit_scaling[0], errors[0])
        print("\nSpecific Heat alpha/nu = %.3f +/- %.3f\n" % (fit_scaling[0], errors[0]) )
    elif fss_observable == "Susceptibility":
        label_string = "$\chi \sim L^{\gamma/\\nu},\; \gamma/\\nu = %.3f \pm %.3f$" % (fit_scaling[0], errors[0])
        print("\nSusceptibility gamma/nu = %.3f +/- %.3f\n" % (fit_scaling[0], errors[0]) )

    # Finally Plot the results
    lines = ax.plot( smooth_L, np.exp( fit_predictions), color = "red", ls = "dashed", label = r"FSS Fit: %s" % label_string )
    lines = ax.scatter( Lvalues, max_obs, label = None )
    if data_tuples[Ldx][2].shape != (0,0):
        lines = ax.errorbar( Lvalues, max_obs, yerr=max_obs_err, color = "None", ecolor = "C0", ms = 5, marker = "o", mfc = "C0", mec = "black", mew = 1,  ls=None, label = None, capsize=2)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(xlabel, fontsize = 12)
    ax.set_ylabel("Peak %s" % labels[lbl], fontsize = 12)
    ax.set_title(model_name + ": " + key_string, fontsize = 12)
    ax.legend(fontsize = 12)

    plotname = "%s" % ("Peak " + labels[lbl] + " vs " + xlabel + ".png")

    plt.savefig(plot_directory + "/" + plotname)

    # Now plot Tc on loglog
    if len(data_tuples) > 4:

        Tfit_scaling, Tcovariance = opt.curve_fit( lambda x,a,b,p: a + b * x ** p, 1./Lvalues, Tc_obs )
        smooth_over_L = np.linspace(0, 1/np.min(Lvalues), 1000)
        Tfit_predictions = Tfit_scaling[0] + Tfit_scaling[1] * (smooth_over_L) ** Tfit_scaling[2]
        Terrors = np.diag(Tcovariance)

        nu, nu_err = 1/Tfit_scaling[2], Terrors[2]/Tfit_scaling[2]**2

        label_string = "$T_c(L) = T_c(\infty) + aL^{-1/\\nu}$ \n$T_c(\infty) = %.3f \pm %.3f$ \n$ \\nu = %.3f \pm %.3f$" % ( Tfit_scaling[0], Terrors[0], nu, nu_err )
        ylabel = labels[lbl] + " Pseudo $T_c$"
        if fss_observable == "Specific Heat":
            print("\nSpecific Heat Tc(inf) = %.3f +/- %.3f\n" % (Tfit_scaling[0], Terrors[0]) )
            print("\nSpecific Heat nu = %.3f +/- %.3f\n" % (nu, nu_err) )
        elif fss_observable == "Susceptibility":
            print("\nSusceptibility Tc(inf) = %.3f +/- %.3f\n" % (Tfit_scaling[0], Terrors[0]) )
            print("\nSusceptibility nu = %.3f +/- %.3f\n" % (nu, nu_err) )

        fig, ax = plt.subplots(1,1)
        lines = ax.plot( smooth_over_L ** Tfit_scaling[2], Tfit_predictions, color = "red", ls = "dashed", label = r"%s" % label_string )
        lines = ax.plot( ( 1./Lvalues ) ** Tfit_scaling[2], Tc_obs, color = "None", marker = "o", mfc = "C0", mec = "black", mew = 1, label = None )

        ax.set_xlabel(r"Scaled System Size $L^{-1/\nu}$", fontsize = 12)
        ax.set_ylabel(ylabel, fontsize = 12)
        ax.set_title(model_name + ": " + key_string, fontsize = 12)
        ax.legend( fontsize = 10 )

        plotname = labels[lbl] + " pseudo Tc scaling.png"
        plt.savefig(plot_directory + "/" + plotname)

    else:
        print("\nCannot obtain " + labels[lbl] + " pseudo Tc due to too few system sizes.\n")


    plt.close()

def main():

    args = setup_args()

    coupling_tuples = parse_couplings( args.coupling_symbol, args.coupling_value )

    # Exit the program if the coupling tuples are null
    if len(coupling_tuples) == 0:
        return

    plot_directory = check_for_output(coupling_tuples)

    labels, data_tuples = collect_observables_and_data( args.data_file_stem, args.fss_observable, args.observable_marker, coupling_tuples )

    if labels == None or data_tuples == None:
        return None

    if len(data_tuples) < 3:
        print("\nCannot perform FSS of %s due to too few system sizes. Exiting.\n" % args.fss_observable)
        return

    plot_data_tuples( args.model_name, args.fss_observable, coupling_tuples, labels, data_tuples, plot_directory )



    return None

if __name__ == "__main__":
    main()
