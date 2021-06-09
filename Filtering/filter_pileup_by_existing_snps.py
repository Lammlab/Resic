"""filter_pileup_by_existing_snp.py

The script receives pileup file, cns file(/pileup file) and an output filename.
it filter out the lines from the pileup file that overlapped with the cns file if there
was a change in the read string (one at least of atcgATCG and the other IUPAC nucleotide code in the reads string ) of the cns file.
by defult it keeps the lines that didnt overlap, if flag --rm-non-overlapping is added
it filter those lines out.


Usage:
    snprationotincs.py param
    snprationotincs.py example
    snprationotincs.py <pileup_input> <cns_file> <output_file> [<output_file_rest>] [--rm-non-overlapping]
    snprationotincs.py -h | --help

Options:
    -h --help                       Show this screen.
    -r --rm-non-overlapping         Remove the lines that from pileup_input that didnt match no line in cns_file

"""


# the cns_file arg can be a pileup file or a Bsortedcns file #

from Utility.parallel_generator import *
from Utility.Pileup_class import Pileup_line
import itertools
from Utility.Bsotredcns_class import Bsortedcns_line
from Utility.generators_utilities import *
from docopt import docopt
import os


'''def get_pos_and_id(line):
    """
    @param line : a pileup line object
    return : ref_base , gene_pos of the pileup line
    """

    return line.reference_id,line.gene_pos'''


'''def is_with_any_change(pileup_line): 
    """
    @param pileup_line: a pileup line (object)
    return : true if in the reads string there is one or more of a c g t A C G T , else false
    """
    #all the possible IUPAC nucleotide code
    nuc_changes = ["a","c","g","t","r","y","s","w","k","m","b","d","h","v","n"]
    all_nuc_changes = [letter.upper() for letter in nuc_changes] + nuc_changes

    if any(change in pileup_line.reads() for change in all_nuc_changes):
        return True
    else : 
        return False'''


def filter_pileup_by_ref_file_comparing(pileup_input, ref_file, output_file,rest_output,ref_type=Pileup_line, remove_non_matched = False):
    """
    @param pileup_input: pileup file name
    @param ref_file: pileup/Bsortedcns file name
    @param output_file: output file name with path
    @param rest_output: the name of the filterout file
    @param ref_type: "pileup" or "Bsortedcns" file
    @param remove_non_matched: if remove_non_matched = true it filter lines that didnt overlap out, else keeps them.
    action : receives pileup file, cns file(/pileup file) and an output file.
                it filter out the lines from the pileup file that overlapped with the cns file if there 
                was a change in the read string (one at least of atcgATCG in the reads string ).
                by defult it keeps the lines that didnt overlap, if remove_non_matched = true 
                it filter those lines out.
                the filter- out lines are printed to a file named output_file (without ending ) + "filter_out."+ (the ending)
    """


    #change_functor = is_with_any_change

    with open(pileup_input, "r") as in_file, open(ref_file,"r") as ref, open(output_file, "w") as out, open(rest_output,"w") as out2 :
        
        input_gen = class_generator(Pileup_line,file=in_file)
        ref_gen = class_generator(ref_type,file=ref)

        get_pos_and_id = lambda x : (x.reference_id,x.gene_pos)

        parallel_gen = parallel_generator([input_gen,ref_gen],[get_pos_and_id,get_pos_and_id])

        for [list_from_in_file, list_from_ref_file] in parallel_gen:
                # when there is a line that matched the ref - check the changes, 
                # if not changed then print the line to output
            if (list_from_in_file is not None) and (list_from_ref_file is not None): # if there is overlap
                if not(list_from_ref_file[0].is_with_any_change()): # when there is no change in the reads of the ref file (all "," or ".")
                    out.write(str(list_from_in_file[0]) + "\n")
                else : 
                    out2.write(str(list_from_in_file[0]) + "\n")
            elif (list_from_in_file is not None) and not(remove_non_matched):  # when not matched, just print (only if requested to print the non matched)
                out.write(str(list_from_in_file[0]) + "\n")
            elif (list_from_in_file is not None): # the removed files are printed to the filter_out file
                out2.write(str(list_from_in_file[0]) + "\n")


def print_params():
    pass


def print_example():
    pass


if __name__ == '__main__':
    import sys
    args = docopt(__doc__)
    if args["param"]:
        print_params()
        sys.exit()

    if args["example"]:
        print_example()
        sys.exit()

        # Extracting the parameters
    pileup_input = args['<pileup_input>']
    ref_file = args['<cns_file>']
    output_file = args['<output_file>']
    remove = args['--rm-non-overlapping']
    rest = args['<output_file_rest>']

    if ref_file.split(".")[-1] == "Bsortedcns":
        ref_type1 = Bsortedcns_line
    else:
        ref_type1 = Pileup_line
    if not(rest):
        filename, file_extension = os.path.splitext(output_file)
        rest_output = filename + ".filter_out" + file_extension
    else: rest_output = rest 

    filter_pileup_by_ref_file_comparing(pileup_input,ref_file,output_file,rest_output,ref_type=ref_type1,remove_non_matched=remove)