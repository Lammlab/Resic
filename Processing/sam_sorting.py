""" sam_sorting.

Usage:
  sam_sorting.py param
  sam_sorting.py example
  sam_sorting.py <sam_in_file> <sam_out_file>
  sam_sorting.py (-h | --help)

Options:
  -h --help     Show this screen.

"""
from Utility.generators_utilities import key_sorted_gen
from docopt import docopt
import re


# Exceptions

class IllegalFile(Exception):
    def __init__(self):
        self.message = "File doesn't exists"

    def __str__(self):
        return self.message

    def __repr__(self):
        return self.message


class IllegalFileWrapper(IllegalFile):
    def __init__(self):
        self.message = "File doesn't exists"

    def __str__(self):
        return self.message

    def __repr__(self):
        return self.message


"""
Function Object which used by the interval finder algorithm, that will extract the fields of the range_start and range_end from the sorted_sam file
:param line: the current line in the file
:return: (chromosome name,start range,end range)
"""


def sam_extract_range(line):
    split_line = re.split('\t|\n', line)
    if (line[0] == "@"):
        return (str(float("-inf")), 0, 0)
    elif (split_line[2] == "*"):  # The misaligned lines will be in the end of the output file
        return ("~", 0, 0)
    else:
        return (split_line[2], int(split_line[3]), int(split_line[3]) + len(split_line[9]))


"""
:param sam_filename: the filename of the sam file
:return: output file called ${out_filename} that sorted mainly according to the chromosome number, and secondly according to the interval range
"""


def sam_sorted(sam_filename, out_filename):
    # Sorting the file
    with open(out_filename, "w") as output_file, \
            open(sam_filename, "r") as file:
        """
        file_lines = file.readlines()
        for line in file_lines:
            if line[0] == "@":
                output_file.write(line)
            else:
                break
        for line in sorted(filter(lambda x: x[0] != "@", file_lines), key = sam_extract_range):
            output_file.write(line)
        """

        for line in key_sorted_gen(sam_extract_range, file):
            output_file.write(line)


def parameters_description():
    string = "The parameters are:" + "\n" \
             + "sam_in_file: The name of the sam file that you would like to sort" + "\n" \
             + "sam_out_file: The name of the output (the sorted sam file)"
    print(string)


def example_description():
    input_string = "The input is:\n" \
                   + "NB500948:31:HJYCCBGXX:1:11101:26701:10507\t0\tchrI\t19959533\t255\t5M\t*\t0\t0\tCCGAC\tAAAAA\tXA:i:0\tMD:Z:123\tNM:i:0\n" \
                   + "NB500948:31:HJYCCBGXX:1:11101:11478:10492\t16\tchrI\t15070824\t255\t12M\t*\t0\t0\tAGCTCACGTTGA\tEEEEEEE/E/EE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                   + "NB500948:31:HJYCCBGXX:1:11101:11478:10492\t16\tchrI\t15063627\t255\t10M\t*\t0\t0\tAGCTCACGTT\tEEEAAEAEEE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n"

    output_string = "The output is:\n" \
                    + "NB500948:31:HJYCCBGXX:1:11101:11478:10492\t16\tchrI\t15063627\t255\t10M\t*\t0\t0\tAGCTCACGTT\tEEEAAEAEEE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                    + "NB500948:31:HJYCCBGXX:1:11101:11478:10492\t16\tchrI\t15070824\t255\t12M\t*\t0\t0\tAGCTCACGTTGA\tEEEEEEE/E/EE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                    + "NB500948:31:HJYCCBGXX:1:11101:26701:10507\t0\tchrI\t19959533\t255\t5M\t*\t0\t0\tCCGAC\tAAAAA\tXA:i:0\tMD:Z:123\tNM:i:0\n"

    print(input_string)
    print(output_string)


"""

"""
if __name__ == "__main__":
    arguments = docopt(__doc__)

    if arguments['param']:
        parameters_description()

    if arguments['example']:
        example_description()

    if (arguments['<sam_in_file>'] and arguments['<sam_out_file>']):
        sam_sorted(arguments['<sam_in_file>'], arguments['<sam_out_file>'])
