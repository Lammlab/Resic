"""Filter_pileup_by_reference_base_set

Usage:
    filter_pileup_by_reference_base_set.py param
    filter_pileup_by_reference_base_set.py example
    filter_pileup_by_reference_base_set.py <pileup_filename> <output_filename> <reference_nucleotide>...
    filter_pileup_by_reference_base_set.py -h | --help

Options:
    -h --help   show this screen
            
"""

##############################################################################################################################
# Author:  Roni Haas
# Refactored by: Clara Frydman
# Main Goal: Receives a pileup file and a nucleotide, and filters the pileup by reference nucleotide. 
##############################################################################################################################

import sys
from docopt import docopt


class IllegalArgument(Exception):
    def __init__(self):
        self.message = "Please enter the name of pileup file, output file and at least one nucleotide for reference base"

    def __str__(self):
        return self.message

    def __repr__(self):
        return self.message


def filter_pileup_by_reference_base_set(pileup, output, listBases):
    """
    :param pileup:  the pile up we want to filter
    :param listBases: a list of ref bases we want to filter by
    :param output: the output file that will contain lines with only the bases in listBases as ref base
    """
    with open(pileup, "r") as file1:
        line = file1.readline()
        with open(output, "w") as fileout:
            while line is not None and line != "":
                lineParts = line.split()
                if lineParts[2] in listBases:
                    fileout.write(line)
                line = file1.readline()


def print_params():
    params_string = """
    Parameters:

    """
    print(params_string)


def print_example():
    example_string = """
    Example: The following execution: 
    
    """
    print(example_string)


if __name__ == "__main__":

    args = docopt(__doc__)

    if args['param']:
        print_params()
        sys.exit()

    if args['example']:
        print_example()
        sys.exit()

    input_pileup = args['<pileup_filename>']
    output_filename = args['<output_filename>']
    reference_nucleotides = args['<reference_nucleotide>']

    if not reference_nucleotides or reference_nucleotides == []:
        raise ValueError("Must be at least one nucleotide to filter by.")

    for reference_nucleotide in reference_nucleotides:
        if reference_nucleotide.lower() not in ['a', 'c', 'g', 't']:
            raise ValueError("Reference nucleotides must be one of: A, T, G, C, a, t, c, g.")

    filter_pileup_by_reference_base_set(input_pileup, output_filename, reference_nucleotides)

    """
    input : python filter_pileup_by_reference_base_set.py pileup_filename output_filename ref_bases_list
            example : python filter_pileup_by_reference_base_set.py exmple.pileup output.pileup A G C

            pileup_filename - the file we want to filter from by ref bases - to keep the lines with those ref bases
            output_filename - the output file for the result filter
            ref_bases_list - a list of ref bases we want to keep in the output file

    output : the lines from the pileup_filename file that contain one of the ref_bases_list as ref base.

    """

    # if len(sys.argv)<3 : raise IllegalArgument
    # pileup = sys.argv[1]
    # output = sys.argv[2]

    # numberOfBases = len(sys.argv) -3
    # if numberOfBases == 0: 
    #     raise IllegalArgument
    # if numberOfBases > 0:
    #     one = sys.argv[3]
    #     if numberOfBases > 1:
    #         two = sys.argv[4]
    #         if numberOfBases > 2:
    #             three = sys.argv[5]
    #             if numberOfBases > 3:
    #                 four = sys.argv[6]
    #                 if numberOfBases > 4:
    #                     five = sys.argv[7]
    #                     filter_pileup_by_reference_base_set(pileup,output,[one,two,three,four,five])
    #                 else:
    #                     filter_pileup_by_reference_base_set(pileup, output,[one, two, three, four])
    #             else:
    #                 filter_pileup_by_reference_base_set(pileup,output, [one, two, three])
    #         else:
    #             filter_pileup_by_reference_base_set(pileup,output, [one, two])
    #     else: 
    #         filter_pileup_by_reference_base_set(pileup,output, [one])
