"""merge sam
Usage:
	merge_sams.py param
	merge_sams.py example
	merge_sams.py <in_sam> <in_sam> ... -o out_sam
	merge_sams.py -h

Options:
-h, --help                              Show this screen.
-o output_file, --out output_file       Name of output file


"""
from docopt import docopt
from Utility.generators_utilities import fileGetChunks
from Utility.Samfile_class import split_sam
from Processing.sam_sorting import sam_sorted
from Utility.temp_files import new_temp_file
import os
import logging

"""
@param: in_list, list of sam file paths

The function add newline (if necessary)
"""


def add_newline(in_list):
    is_changed = False

    for in_sam in in_list:
        with open(in_sam, "a+b") as fp:
            fp.seek(-2, 2)
            bytes = fp.read(2)
            last_chars = [chr(b) for b in bytes]

            # In this case there isn't a new line at the end of the file, and therefore the we need to add newline
            if (set(last_chars) == {'\n', '\r'} or last_chars[1] == '\n' or last_chars[1] == '\r') == False:
                is_changed = True
                with open(in_sam, "a+") as fp_write:
                    fp_write.write("\n")

    # Check wheather we need to warn the user about editing a file (only if we had to add newline in at least one file)
    if is_changed == True:
        logging.warning("Changes to the file are done")


"""
@param: ih_list, list of sam file paths
@param: out_sam, path to write merged sam to
"""


def merge_sams(in_list, out_sam):
    (headers, alligns, misalligns) = ([], [], [])
    # split all files
    for in_sam in in_list:
        head, all, misall = split_sam(in_sam)
        for list, file in zip([headers, alligns, misalligns], [head, all, misall]):
            list.append(file)

    unsorted_out = new_temp_file()

    # write all headers than all alligned than all un alligned
    with open(unsorted_out, "w") as out_fp:
        for part_list in [headers, alligns, misalligns]:
            for part in part_list:
                with open(part, "r") as fp:
                    for chunk in fileGetChunks(fp):
                        out_fp.writelines(chunk)
                os.remove(part)

    sam_sorted(unsorted_out, out_sam)
    os.remove(unsorted_out)

    return


def parameters_description():
    string = "The parameters are:" + "\n" \
             + "Infiles: a list of sam file names\n" \
             + "-o Outfile: an \'o\' flag followed by the path of the desired output filename.\n " \
             + "If non is given the output is printed to the stdout channel "
    print(string)


def example_description():
    input_sam_1 = "@HD\tVN:1.6\tbbb\n" \
                  + "@SQ\tSN:ref\tLN:45_2\n" \
                  + "seq1\t16\tchrX\t200\t255\t10M\t*\t0\t0\tAGCTCACGTT\tEEEAAEAEEE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                  + "ERROR2\t0\t*\t19959533\t255\t5M\t*\t0\t0\tCCGAC\tAAAAA\tXA:i:0\tMD:Z:123\tNM:i:0\n"

    input_sam_2 = "@HD\tVN:1.6\taaa\n" \
                  + "@SQ\tSN:ref\tLN:45_1\n" \
                  + "seq1\t16\tchrX\t100\t255\t10M\t*\t0\t0\tAGCTCACGTT\tEEEAAEAEEE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                  + "ERROR1\t0\t*\t19959533\t255\t5M\t*\t0\t0\tCCGAC\tAAAAA\tXA:i:0\tMD:Z:123\tNM:i:0\n"

    input_sam_3 = "@HD\tVN:1.6\tccc\n" \
                  + "@SQ\tSN:ref\tLN:45_3\n" \
                  + "seq1\t16\tchrY\t100\t255\t10M\t*\t0\t0\tAGCTCACGTT\tEEEAAEAEEE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                  + "ERROR3\t0\t*\t19959533\t255\t5M\t*\t0\t0\tCCGAC\tAAAAA\tXA:i:0\tMD:Z:123\tNM:i:0\n"

    output_sam = "@HD\tVN:1.6\tbbb\n" \
                 + "@SQ\tSN:ref\tLN:45_2\n" \
                 + "@HD\tVN:1.6\taaa\n" \
                 + "@SQ\tSN:ref\tLN:45_1\n" \
                 + "@HD\tVN:1.6\tccc\n" \
                 + "@SQ\tSN:ref\tLN:45_3\n" \
                 + "seq1\t16\tchrX\t100\t255\t10M\t*\t0\t0\tAGCTCACGTT\tEEEAAEAEEE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                 + "seq1\t16\tchrX\t200\t255\t10M\t*\t0\t0\tAGCTCACGTT\tEEEAAEAEEE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                 + "seq1\t16\tchrY\t100\t255\t10M\t*\t0\t0\tAGCTCACGTT\tEEEAAEAEEE\tXA:i:0\tMD:Z:53G7C19A68\tNM:i:3\n" \
                 + "ERROR2\t0\t*\t19959533\t255\t5M\t*\t0\t0\tCCGAC\tAAAAA\tXA:i:0\tMD:Z:123\tNM:i:0\n" \
                 + "ERROR1\t0\t*\t19959533\t255\t5M\t*\t0\t0\tCCGAC\tAAAAA\tXA:i:0\tMD:Z:123\tNM:i:0\n" \
                 + "ERROR3\t0\t*\t19959533\t255\t5M\t*\t0\t0\tCCGAC\tAAAAA\tXA:i:0\tMD:Z:123\tNM:i:0\n"

    print("first sam is:\n")
    print(input_sam_1)
    print("second sam is:\n")
    print(input_sam_2)
    print("third sam is:\n")
    print(input_sam_3)
    print("merged sam is:\n")
    print(output_sam)


if __name__ == "__main__":
    arguments = docopt(__doc__)
    if arguments['param']:
        parameters_description()

    if arguments['example']:
        example_description()

    add_newline(arguments['<in_sam>'])
    merge_sams(arguments['<in_sam>'], arguments['--out'])
