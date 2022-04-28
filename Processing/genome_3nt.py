"""
Functions to sort/combine nucleotide pairs/combine lines of multiline formats, and to revert combined formats to regular form

Usage:
    genome_3nt.py FASTA_FILENAME FASTA_OUTPUT FASTQ_FILENAME FASTQ_OUTPUT
Options:
    -FASTA_FILENAME fasta filename
"""

from Utility.generators_utilities import class_generator
from Utility.parallel_generator import parallel_generator
import subprocess
from Utility.Fastq_class import Fastq
from Utility.Samfile_class import SamLine
import re
from docopt import docopt
import itertools
from Utility.multiline_sort import multiline_sort
import os
import shutil


def pre(fasta_file, fastq_file, out_fasta_file, out_fastq_file, nt_replacement):
    if type(nt_replacement) is list:
        nt_replacement = "".join(nt_replacement)

    if fasta_file == out_fasta_file:
        fasta_file_temp = fasta_file + "_temp_copy.fasta"
        shutil.copyfile(fasta_file, fasta_file_temp)
        fasta_nt = fasta_to_3nt(nt_replacement, fasta_file_temp, output_file=out_fasta_file)
        os.remove(fasta_file_temp)
    else:
        fasta_nt = fasta_to_3nt(nt_replacement, fasta_file, output_file=out_fasta_file)

    # logging.info(f"==exec==\nmultiline_to_3nt({2},{'~'}, {nt_replacement}, {fasta_file},output_file={out_fasta_file})")
    # logging.info(f"==exec==\nmultiline_to_3nt({2},chr(127), {nt_replacement}, {fastq_file},output_file={out_fastq_file})")
    if fastq_file is not None:
        fastq_nt = multiline_to_3nt(4, chr(127), 2, nt_replacement, fastq_file, output_file=out_fastq_file)
    else:
        fastq_nt = None

    return fastq_nt, fasta_nt


def post(fastq_file, sam_file, out_fastq_file, out_sam_file, **kwargs):
    # logging.info(f"==exec==\nrevert_fastq_to_4nt({fastq_file}, kwargs['fastq_lib'],output_filename=out_fastq_file)multiline_to_3nt({2},chr(127), {fastq_file},output_file={out_fastq_file})")
    if fastq_file is not None:
        fastq_ntr = revert_fastq_to_4nt(fastq_file, kwargs['fastq_lib'], output_filename=out_fastq_file)
    else:
        fastq_ntr = None
    # logging.info(f"==exec==\nmultiline_to_3nt({2},chr(127), {fastq_file},output_file={out_fastq_file})")
    sam_ntr = revert_sam_to_4nt(sam_file, output_filename=out_sam_file, fastq_lib=kwargs['fastq_lib'])

    return fastq_ntr, sam_ntr


# nt_pairings = {'AG': 'A', 'GA': 'A', 'CTU': 'U', 'AC': 'C', 'CA': 'C', 'CG': 'G', 'GC': 'G'}
nt_pairings = {"{x}{y}".format(x=x, y=y): min(x, y) for x, y in itertools.product("AGCT", repeat=2)}


def fasta_to_3nt(nt_replacement, filename, output_file=None):
    '''
    alternative function to multiline_to_3nt for fasta files with multiple lines for the sequence, instead of a single
    line with the entire sequence
    :param nt_replacement: nucleotides to combine (ex: 'AG', 'CT')
    :param filename: filename for the input
    :param output_file: filename for the output, if desired
    :return: output filename
    '''
    if output_file is None:
        splited = filename.split(".")
        output_file = ".".join(splited[:-1]) + "_nt." + splited[-1]
    with open(filename, 'r'), open(output_file, 'w'):
        # TODO: does case matter?
        # if type(nt_replacement) is list:
        # nt_replacement = "".join(nt_replacement)
        # lower_case_nt = nt_replacement.lower()
        lower_case_replacement = nt_pairings[nt_replacement].lower()
        # print(lower_case_replacement)
        command = f"cat {filename} | sed -r '/^[ \t]*$/d' | sed '/^[NAGCTnagct]/s/[{nt_replacement}]/{nt_pairings[nt_replacement]}/g' " f" | sed '/^[NAGCTnagct]/s/{[nt_replacement.lower()]}/{lower_case_replacement}/g' > {output_file} "
        # command = f"cat {filename} | sed -r '/^[ \t]*$/d' | sed '/[AGCTagct]*/s/[{nt_replacement}]/{nt_pairings[nt_replacement]}/ig' > {output_file} "
        # ^ also works
        # print(command)
        try:
            subprocess.run(command, shell=True)
        except Exception as e:
            raise e
    return output_file


def multiline_to_3nt(number_of_lines, new_delim, sequence_line, nt_replacement, multi_filename=None, output_file=None):
    '''

    :param number_of_lines: number of lines in each entry
    :param new_delim: any character that is not within the alphabet of possible characters that are used in the format
    :param sequence_line: the number of the line containing the sequence to be modified
    :param nt_replacement: the nucleotides to be combined
    :param multi_filename: name of the file containing the data
    :param output_file: optional filename to write output to. if not specified output goes to input filename with .nt extension
    :return:
    '''

    if output_file is None:
        splited = multi_filename.split(".")
        output_file = ".".join(splited[:-1]) + "_nt." + splited[-1]
    with open(multi_filename, 'r'), open(output_file, 'w'):
        # TODO: does case matter?
        if type(nt_replacement) is list:
            nt_replacement = "".join(nt_replacement)
        # command = f"cat {multi_filename} | sed -r '/^[ \t]*$/d' |paste -d \"{new_delim}\" {' -'*int(number_of_lines)} \
        # | awk -F{new_delim} 'BEGIN {{OFS = \"{new_delim}\"}}; {{gsub(/[{nt_replacement}]/,\"{nt_pairings[nt_replacement]}\",${sequence_line})}};{{print}};'\
        # | sed 's/{1}/\\n/g' > {output_file} "
        # | awk -F{new_delim} 'BEGIN {{OFS = \"{new_delim}\"}}; {{gsub(/[tolower({nt_replacement})]/,tolower(\"{nt_pairings[nt_replacement]}\"),${sequence_line})}};{{print}}' \

        # OLD command with number variable replacement rather than format string replacement.
        command = "cat {0} | sed -r '/^[ \t]*$/d' |paste -d \"{1}\" {2}\
        | awk -F{1} 'BEGIN {{OFS = \"{1}\"}}; {{${5} = toupper(${5}); gsub(/[{3}]/,\"{4}\",${5})}};{{print}}' |" \
                  "sed 's/{1}/\\n/g' > {6} ".format(multi_filename, new_delim, ' -' * int(number_of_lines),
                                                    nt_replacement, nt_pairings[nt_replacement], sequence_line,
                                                    output_file)
        # print(command)

        # bash command for fastqs, AT/at to A/a
        # "cat smallseq.fastq | sed -r '/^[ \t]*$/d' |paste -d "~" - - - -| awk -F~ 'BEGIN {OFS="~"};{gsub(/[AT]/,"A",$2);gsub(/[at]/,"a",$2)};{print}' | sed 's/~/\n/g' >smallseq.nt.fastq"

        try:
            subprocess.run(command, shell=True)
        except Exception as e:
            raise e
    return output_file


def revert_fastq_to_4nt(fastq_3nt_filename, fastq_4nt_filename, output_filename=None):
    '''

    :param fastq_3nt_filename: negative reads from alignment of 3nt fastq library
    :param fastq_4nt_filename: original 4nt fastq library
    :param output_filename:
    :return:
    '''
    sorted_fastq_3nt = multiline_sort(4, '~', 1, fastq_3nt_filename)
    sorted_fastq_4nt = multiline_sort(4, '~', 1, fastq_4nt_filename)
    if output_filename is None:
        splited = fastq_3nt_filename.split(".")
        output_filename = ".".join(splited[:-1]) + "_ntR." + splited[-1]

    with open(sorted_fastq_3nt, 'r') as fastq_3nt, open(sorted_fastq_4nt, 'r') as fastq_4nt:
        gen_3nt = class_generator(Fastq, file=fastq_3nt, number_of_lines=4)
        gen_4nt = class_generator(Fastq, file=fastq_4nt, number_of_lines=4)
        get_fastq_id = lambda fq: fq.strings[0]
        with open(output_filename, 'w') as out:
            p_gen = parallel_generator([gen_3nt, gen_4nt], [get_fastq_id, get_fastq_id])
            for [fq_3, fq_4] in p_gen:
                if (fq_3 is not None) and (fq_4 is not None):
                    fq_3[0].strings[1] = fq_4[0].strings[1]
                    out.write(re.sub('\n\n', '\n', str(fq_3[0])))

    return output_filename


def reverse_multiline(number_of_lines, new_delim, sequence_line, multi_filename, output_file=None):
    '''

    :param number_of_lines: number of lines in each entry
    :param new_delim: any character that is not within the alphabet of possible characters that are used in the format
    :param sequence_line: the number of the line containing the sequence to be modified
    :param nt_replacement: the nucleotides to be combined
    :param multi_filename: name of the file containing the data
    :param output_file: optional filename to write output to. if not specified output goes to input filename with .nt extension
    :return:
    '''

    if output_file is None:
        splited = multi_filename.split(".")
        output_file = ".".join(splited[:-1]) + "_nt." + splited[-1]
    with open(multi_filename, 'r'), open(output_file, 'w'):
        # OLD command with number variable replacement rather than format string replacement.
        command = "cat {0} | sed -r '/^[ \t]*$/d' |paste -d \"{1}\" {2}\
        | awk -F{1} 'BEGIN {{OFS = \"{1}\"}};{{cmd=\"echo \"${3}\" | rev\"; cmd | getline reversed_line ; close(cmd) ;gsub(${3},reversed_line,${3})}};{{print}}' |" \
                  "sed 's/{1}/\\n/g' > {4} ".format(multi_filename, new_delim, ' -' * int(number_of_lines),
                                                    sequence_line, output_file)
        # print(command)

        # bash command for fastqs, AT/at to A/a
        # "cat smallseq.fastq | sed -r '/^[ \t]*$/d' |paste -d "~" - - - -| awk -F~ 'BEGIN {OFS="~"};{gsub(/[AT]/,"A",$2);gsub(/[at]/,"a",$2)};{print}' | sed 's/~/\n/g' >smallseq.nt.fastq"

        try:
            subprocess.run(command, shell=True)
        except Exception as e:
            raise e
    return output_file


def complement_multiline(number_of_lines, new_delim, sequence_line, multi_filename, output_filename=None):
    '''
    Complement the sequence line of a multilined format file. Nucleotides are switched to their appropriate base pair (A->T, C->G, etc)
    :param number_of_lines: number of lines in each entry
    :param new_delim: any character that is not within the alphabet of possible characters that are used in the format
    :param sequence_line: the number of the line containing the sequence to be modified
    :param nt_replacement: the nucleotides to be combined
    :param multi_filename: name of the file containing the data
    :param output_file: optional filename to write output to. if not specified output goes to input filename with .nt extension
    :return:
    '''
    # TODO: only works on fasta now
    if output_filename is None:
        splited = multi_filename.split(".")
        output_filename = ".".join(splited[:-1]) + "_complement." + splited[-1]
    with open(multi_filename, 'r') as fasta_file, open(output_filename, 'w') as out:
        translation = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
        # TODO: use "tr AGCTagct TCGAtcga" command to flip the nucleotides
        fasta_gen = class_generator(str, file=fasta_file)
        for line in fasta_gen:
            if line[0] != '>':
                line = line.translate(translation)
            out.write(line.strip("\n") + "\n")
    return output_filename


def reverse_complement_sequence(sequence):
    '''
    reverse complement of sequence of nucleotides
    :return: reversed sequence
    '''
    reverse_dict = {'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', 'n': 'n'}
    reversed_sequence = ''
    for nt in reversed(sequence):
        reversed_sequence += reverse_dict[nt]
    return reversed_sequence


def revert_samline_to_4nt(samline, fastq_entry):
    '''
    replace the sequence field of a SAM file line with the sequence from a fastq entry

    :param samline: samline to revert. Of type SamLine, from Utility.Samfile_class
    :param fastq_entry: fastq entry to use as reference
    :return: the updated SamLine object
    '''
    fastq_sequence = fastq_entry.get_seq().strip("\n")
    if samline.flags == 16:
        fastq_sequence = reverse_complement_sequence(fastq_sequence)
    samline.sequence = fastq_sequence
    return samline


def flip_sam_flags(samfile, flag_to_replace, new_flag, output_filename=None):
    if output_filename is None:
        splited = samfile.split(".")
        output_filename = ".".join(splited[:-1]) + "_flipped." + splited[-1]

    # create temp file if input and output filenames are the same
    if samfile == output_filename:
        temp_name = samfile + "temp"
        shutil.copyfile(samfile, temp_name)
        samfile = temp_name

    with open(samfile, 'r') as sam_file, open(output_filename, 'w') as out:
        header_gen = class_generator(str, skip_condition=(lambda line: line[0] != '@'), file=sam_file)
        for line in header_gen:
            out.write(str(line))
    with open(samfile, 'r') as sam_file, open(output_filename, 'a') as out:
        sam_gen = class_generator(SamLine, skip_condition=(lambda line: line[0] == '@'), file=sam_file)
        for line in sam_gen:
            if line.flags == flag_to_replace:
                line.flags = new_flag  # samclass turns flag to string using str(new_flag), turning(for example) 16->20
            out.write(str(line) + '\n')

    if samfile == output_filename:
        os.remove(temp_name)
    return output_filename


def sort_sam_alphabetic(samfile):
    splited = samfile.split(".")
    sorted_name = ".".join(splited[:-1]) + "_sorted." + splited[-1]
    command = "cat {0} | sed -r '/^[ \t]*$/d' | LC_ALL=C sort -k1,1n > {1}".format(samfile, sorted_name)
    try:
        subprocess.run(command, shell=True)
    except Exception as e:
        raise (e)
    return sorted_name


def revert_sam_to_4nt(aligned_3nt_sam_filename, output_filename=None, sorted=False, **kwargs):
    '''

    :param aligned_3nt_sam_filename: filename of sam format file, which should contain sam lines that correspond to
    a library of reads in fastq format to revert the sam file according to.
    :param output_filename: filename to write output to
    :param kwargs: fastq_lib: the original library of sequence reads that was used to align the sam file
    :return: None, in case of error
    '''
    if 'fastq_lib' not in kwargs:
        print("no reference library provided")
        return None
    if output_filename is None:
        splited = aligned_3nt_sam_filename.split(".")
        output_filename = ".".join(splited[:-1]) + "_ntR." + splited[-1]

    sequence_library_filename = kwargs['fastq_lib']
    if not sorted:
        # sort library of fastq reads according to sequence id
        sorted_lib_filename = multiline_sort(4, "~", 1, sequence_library_filename)

        # sort sam file according to query sequence id
        sorted_sam_filename = sort_sam_alphabetic(aligned_3nt_sam_filename)
    else:
        sorted_lib_filename = sequence_library_filename
        sorted_sam_filename = aligned_3nt_sam_filename

    # create temp file to avoid overwriting input with output if same input and output filenames
    if output_filename == aligned_3nt_sam_filename:
        temp_name = sorted_sam_filename + "temp_copy"
        shutil.copyfile(sorted_sam_filename, temp_name)
        sorted_sam_filename = temp_name

    # add the header lines
    with open(sorted_sam_filename, 'r') as sam_file, open(output_filename, 'w') as out:
        header_gen = class_generator(str, skip_condition=(lambda line: line[0] != '@'), file=sam_file)
        for line in header_gen:
            out.write(str(line))
    # add the reverted sam lines
    with open(sorted_sam_filename, 'r') as sam_file, open(sorted_lib_filename, 'r') as lib_file, open(output_filename,
                                                                                                      'a') as out:
        aligned_3nt_gen = class_generator(SamLine, skip_condition=(lambda line: line[0] == '@'), file=sam_file)
        seq_library_gen = class_generator(Fastq, file=lib_file, number_of_lines=4)

        # TODO: update functors when fastq_class/samfile_class are updated
        get_sam_qname = lambda x: x.seq_id.split(" ")[0]
        get_fastq_id = lambda x: x.strings[0][1:].split(" ")[0]
        p_gen = parallel_generator([aligned_3nt_gen, seq_library_gen], [get_sam_qname, get_fastq_id])
        for [samlines, fastq_entry] in p_gen:
            # print("samlines: ",samlines,"fastq_entry: ",fastq_entry)
            if (samlines is not None and fastq_entry is not None):
                for samline in samlines:
                    reverted_samline = revert_samline_to_4nt(samline, fastq_entry[0])
                    out.write(str(reverted_samline) + "\n")
    if output_filename == aligned_3nt_sam_filename:
        os.remove(temp_name)
    return output_filename


# ignore main, not for use
if __name__ == '__main__':
    arguments = docopt(__doc__)
    try:
        multiline_to_3nt(4, '~', 2, 'AG', arguments['FASTQ_FILENAME'], arguments['FASTQ_OUTPUT'])
        multiline_to_3nt(2, '~', 2, 'AG', arguments['FASTA_FILENAME'], arguments['FASTA_OUTPUT'])
    except Exception as e:
        print("error")
