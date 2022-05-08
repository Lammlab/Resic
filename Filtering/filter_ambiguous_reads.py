"""filter_ambiguous_reads

	Usage:
		filter_ambiguous_reads param
		filter_ambiguous_reads example
		filter_ambiguous_reads <fastq_file> <out_file>
		filter_ambiguous_reads -h | --help

	Options:
		-h --help   Show this screen
"""

#######################################################################################################################
# Author:  Roni Haas
# Main goal: Takes a fastq file and by using biopython library performs several filtering steps for
# ambiguous reads, reflected by the nucleotides appearance in the read, quality scores, %N and more.
# The output is a new fastq file with only high quality reads
#######################################################################################################################

from docopt import docopt
from Bio import SeqIO
from math import log
import sys


def filter_by_single_nucleotide_appearance(read_sequance, upper_appearance_thresh=60, lower_appearance_thresh=10, \
                                           N_max_appearance_thresh=10):
    """
    This function filters based on the % appearance of each nucleotide in the read
    :param read_sequance: the seq of the read taken from the fastq file
    :param upper_appearance_thresh: the max % of each nucleotide out of all read's bp
    :param lower_appearance_thresh: the min % of each nucleotide out of all read's bp
    :param N_max_appearance_thresh: the max % of Ns out of all read's bp
    :return True if meets all requirements, False if not :
    """
    nucl_changes = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    upper_string = read_sequance.upper()
    # counts the number of shows for each nucleotide or N
    for nucl in list(upper_string):
        try:
            nucl_changes[nucl] = nucl_changes[nucl] + 1
        except KeyError:
            return False
    # replase the total number of shows into the % of shows out of the read len
    for key, value in nucl_changes.items():
        nucl_changes[key] = value * 100 / len(read_sequance)
    # test if the read meets the requirements
    return all(x <= upper_appearance_thresh for x in (list(nucl_changes.values())[:-1])) and \
           all(x >= lower_appearance_thresh for x in list(nucl_changes.values())[:-1]) and \
           nucl_changes['N'] < N_max_appearance_thresh


def ave_qual(quals):
    """
    This function calculates the average basecall quality of a read. calculates the average error
    probability and convert average back to Phred scale
    NOTE: Phred scores are log probabilities, so simply taking an average of those is wrong.
    :param quals: a list of integer qualities for all read's bp
    :return: the average error probability
    """
    if quals:
        return -10 * log(sum([10 ** (q / -10) for q in quals]) / len(quals), 10)
    else:
        return None


def filter_by_quality(quality_string, quality_thrsh=25):
    """
    :param quality_string: a list of qualities for all read's bp
    :param quality_thrsh: the max quality_thrsh aloud
    :return: True if avg quality meets all requirements
    """
    return ave_qual(quality_string) > quality_thrsh


def filter_by_long_stretches_repeats(seq_string, stretche_max_len=20):
    """
    :param seq_string: the read sequence
    :param stretche_max_len: the max shows in a row for a single nucleotide
    :return: True if meets thr requirement, False if no
    """

    list_of_nuc = ["A", "C", "G", "T"]
    for nuc in list_of_nuc:
        counter = 0
        # for the first match
        prev_nuc = nuc
        for nuc_from_seq in seq_string:
            if nuc_from_seq == nuc and prev_nuc == nuc:
                counter += 1
                prev_nuc = nuc_from_seq
                if counter >= stretche_max_len:
                    return False
            elif nuc_from_seq == nuc and prev_nuc != nuc:
                counter = 0
                prev_nuc = nuc_from_seq
                continue
            elif nuc_from_seq != nuc:
                prev_nuc = nuc_from_seq
                continue
    return True


def extract_from_fastq(fq, output_fq):
    """
    Takes a fastq file, examines each read using all the above functions, and writes to a
    new file the non-ambiguous reads
    :param fq: the fastq file
    :param output_fq: the output fastq file after filtering
    """
    input_iterator = SeqIO.parse(fq, "fastq")
    # gos over each record and tests if the read meets the requirements
    short_iterator = (rec for rec in input_iterator if filter_by_quality(rec.letter_annotations["phred_quality"]) \
                      and filter_by_single_nucleotide_appearance(rec.seq) and filter_by_long_stretches_repeats(rec.seq))
    # writes to a new file after the conversion to a fastq format
    SeqIO.write(short_iterator, output_fq, "fastq")


def print_params():
    string = "The parameters are:" + "\n" \
             + "fastq_file: a standard fastq file" + "\n" \
             + "out_file: the desired output fastq name after filtering according to the reads quality"
    print(string)


def print_example():
    string_for_print = "Input fastq file. First read is of low quality (in this case >10% N in the seq):" + "\n" \
                       + "@HISEQ:79:C7V66ANXX:4:1101:1069:2041 1:N:0:\n" \
                       + "CNNGCACCATCANNAGATCGGAANNNNNCGTATGCCGNNNNCTGCTTGAA\n" \
                       + "+\n" \
                       + "BBBBFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFFF\n" \
                       + "@HISEQ:79:C7V66ANXX:4:1101:1153:2070 1:N:0:\n" \
                       + "CTGGCACCATCAATAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAA\n" \
                       + "+\n" \
                       + "<BBBBFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFF\n" \
                       + "Output fastq file:\n" \
                       + "@HISEQ:79:C7V66ANXX:4:1101:1153:2070 1:N:0:\n" \
                       + "CTGGCACCATCAATAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAA\n" \
                       + "+\n" \
                       + "<BBBBFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFF\n"
    print(string_for_print)


if __name__ == "__main__":
    args = docopt(__doc__)
    if args['param']:
        print_params()
        sys.exit()

    if args['example']:
        print_example()
        sys.exit()
    fastq_file = args['<fastq_file>']
    output = args['<out_file>']
    extract_from_fastq(fastq_file, output)
