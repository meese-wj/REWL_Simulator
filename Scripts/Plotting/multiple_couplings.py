# Use this file to house the multiple couplings
# parsing functions.

import os

def parse_couplings(coupling_symbols, coupling_values):
    # Strip the space-separated lists and combine into
    # a list of tuples.
    # This assumes the couplings will be listed as
    #   coupling_list = "sifter coupling1 ..."
    # and the values are
    #   values_list = "value1 ..."
    coupling_list = coupling_symbols.split()
    sifter_coupling = coupling_list[0]
    coupling_list = coupling_list[1:]
    values_list = coupling_values.split()
    print(coupling_list, values_list)
    if len(coupling_list) != len(values_list):
        print("\nCoupling mismatch.\nCoupling list = ", coupling_symbols, "\nCoupling values = ", coupling_values)
        return []
    tuples = []
    for i in range(len(coupling_list)):
        tuples.append( (coupling_list[i], values_list[i]) )
    return sifter_coupling, tuples

def get_key_string(coupling_symbol, coupling_value, isfloat = False):
    if coupling_symbol != "L" and isfloat:
        return coupling_symbol + "-" + ("%.6f" % float(coupling_value))
    else:
        return coupling_symbol + "-" + coupling_value

def get_coupling_string( coupling_tuples, isfloat = False ):
    string = get_key_string( coupling_tuples[0][0], coupling_tuples[0][1], isfloat )
    for i in range(1,len(coupling_tuples)):
        string += "_" + get_key_string( coupling_tuples[i][0], coupling_tuples[i][1], isfloat )
    return string

def get_coupling_string_bash( sifter, coupling_tuples ):
    # Return the coupling string as bash
    # would see it because I'm lazy
    string = sifter + " " + coupling_tuples[0][0]
    for i in range(1, len(coupling_tuples)):
        string += " " + coupling_tuples[i][0]
    string += "-" + coupling_tuples[0][1]
    for i in range(1, len(coupling_tuples)):
        string += " " + coupling_tuples[i][1]
    return string

def latex_couplings( coupling_tuples ):
    string = ": $%s = %s$" % (coupling_tuples[0][0], coupling_tuples[0][1])
    for i in range(1,len(coupling_tuples)):
        string += ", $%s = %s$" % (coupling_tuples[i][0], coupling_tuples[i][1])
    return string

def check_for_output(sifter_coupling, coupling_tuples, output_path):
    # Create a figure directory in the data folder
    #coupling_string = get_coupling_string(coupling_tuples)
    coupling_string = get_coupling_string_bash( sifter_coupling, coupling_tuples )
    if not os.path.isdir( os.getcwd() + "/" + output_path + "_" + coupling_string ):
        os.mkdir( os.getcwd() + "/" + output_path + "_" + coupling_string )

    return os.getcwd() + "/" + output_path + "_" + coupling_string

def find_string_value( string_type, file_string ):
    # Determine if the string is in the file name
    start = file_string.find(string_type + "-") + len(string_type + "-")
    dat_ender = file_string[start:].find(".dat.")
    underscore_ender = file_string[start:].find("_")
    ender = dat_ender
    if ender > underscore_ender:
        ender = underscore_ender
    end = start + ender

    return file_string[start:end]

def couplings_in_file( coupling_tuples, file_string ):
    # Return true if the couplings are all
    # in the file string name
    present = True
    for cdx in range(len(coupling_tuples)):
        key_string = get_key_string( coupling_tuples[cdx][0], coupling_tuples[cdx][1], isfloat = True )
        present = present and (key_string in file_string)
    return present

def pretty_label_string(coupling, value_string):
    # Return a nice looking label string
    if coupling == "L":
        return "$%s = %d$" % (coupling, int(value_string))
    else:
        return "$%s = %.2f$" % (coupling, float(value_string))
