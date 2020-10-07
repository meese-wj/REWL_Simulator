import argparse

def set_up_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", help = "Which data file to parse", type = str)
    parser.add_argument("observable_marker", help = "The line to parse for", type = str)

    return parser.parse_args()

def concatenate_strings( string_list, starter = 2, spacer = " " ):
    new_string = ""

    for piece_dx in range(starter, len(string_list)):
        new_string += string_list[piece_dx]
        if piece_dx != len(string_list) - 1:
            new_string += spacer

    return new_string

def collect_labels(data_file, observable_marker, comment = "#"):

    labels = []
    header_lines = []

    data_file = open(data_file, "r")

    data_lines = data_file.readlines()

    line = 0
    slicer = 0
    while data_lines[line][0] == comment:
        header_lines.append( data_lines[line] )
        if observable_marker in data_lines[line]:
            slicer = line
        line += 1

    # Look through the file header
    for header_dx in range(slicer + 1, len(header_lines)):
        split_line = header_lines[header_dx].split()
        if len(split_line) > 1:
            labels.append( concatenate_strings(split_line) )

    data_file.close()
    return labels

def main():

    args = set_up_args()
    collect_labels( args.data_file, args.observable_marker )

    return None

if __name__ == "__main__":
    main()
