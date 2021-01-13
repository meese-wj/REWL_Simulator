# Use this file to house the multiple couplings
# parsing functions.

import os

def parse_couplings(coupling_symbols, coupling_values):
    # Strip the space-separated lists and combine into
    # a list of tuples.
    coupling_list = coupling_symbols.split()
    values_list = coupling_values.split()
    if len(coupling_list) != len(values_list):
        print("\nCoupling mismatch.\nCoupling list = ", coupling_symbols, "\nCoupling values = ", coupling_values)
        return []
    tuples = []
    for i in range(len(coupling_list)):
        tuples.append( (coupling_list[i], values_list[i]) )
    return tuples

def get_key_string(coupling_symbol, coupling_value, isfloat = False):
    if isfloat:
        return coupling_symbol + "-" + ("%.6f" % float(coupling_value))
    else:
        return coupling_symbol + "-" + coupling_value

def get_coupling_string( coupling_tuples, isfloat = False ):
    string = get_key_string( coupling_tuples[0][0], coupling_tuples[0][1], isfloat )
    for i in range(1,len(coupling_tuples)):
        string += "_" + get_key_string( coupling_tuples[i][0], coupling_tuples[i][1], isfloat )
    return string

def latex_couplings( coupling_tuples ):
    string = ": $%s = %s$" % (coupling_tuples[0][0], coupling_tuples[0][1])
    for i in range(1,len(coupling_tuples)):
        string += ", $%s = %s$" % (coupling_tuples[i][0], coupling_tuples[i][1])
    return string

def check_for_output(coupling_tuples):
    # Create a figure directory in the data folder
    coupling_string = get_coupling_string(coupling_tuples)
    if not os.path.isdir( os.getcwd() + "/" + output_path + "_" + coupling_string ):
        os.mkdir( os.getcwd() + "/" + output_path + "_" + coupling_string )

    return os.getcwd() + "/" + output_path + "_" + coupling_string

def find_string_value( string_type, file_string ):
    # Determine if the string is in the file name
    start = file_string.find(string_type + "-") + len(string_type + "-")
    end = start + file_string[start:].find("_")

    return file_string[start:end]

def couplings_in_file( coupling_tuples, file_string ):
    # Return true if the couplings are all
    # in the file string name
    coupling_string = get_coupling_string(coupling_tuples, isfloat = True)
    return ( coupling_string in file_string )


