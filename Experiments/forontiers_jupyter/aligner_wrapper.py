import Utility.threading_tree_utils
import subprocess
import re
import os
import sys
import shutil
from contextlib import nullcontext
import logging
import  sys
import doctest

class AlignerWrapper:
    """
    This class defines the general definition of an aligner, using a provided format string, and allows the execution
    of an alignment using defined aligner.

    format string should include:

    aligner : name of the aligner to be used
    {flags}
    {lib} | {lib_1} {lib_2}  : lib for single-end reads, lib_1 and lib_2 for paired end reads
    {reference_file}
    {positive_alignment}
    {negative_library}
    {log} : optional log file

    format string for building index can also be included:

    index_builder : index building program to be used
    {flags}
    {reference_file}

    index_output_format_list: list of format strings that describe which files constiture the aligner index
    {reference_file}

    example format string for alignment:

    "aligner {flags} {lib_1} {lib_2} {reference_file} --output_sam {positive_alignment} --nomatch {negative_library} 2>{log}"

    example usage:
    # format_string =" new_aligner -ref {reference_file} -in {lib} {flags} --pos-out {positive_alignment} --neg {negative_library}"
    # new_aligner = AlignerWrapper(format_string)
    # new_aligner.align(lib = 'lib.fastq', flags = '', reference= 'gen.fasta', pos='pos.sam', neg='neg.sam', paired_end=False ) #logging=False?

    #doctest example:
   >>> format_string = "bowtie {flags} {lib} --positive {positive_alignment} --un {negative_library} {reference_file} --un fastq >> bowtie 2>&1"
   >>> aligner_wrap = AlignerWrapper(format_string, "bowtie-build {flags} {reference}")
   >>> aligner_wrap.align(flags='-m 3 -n 2', reference='ref.fa', indexing_flags='-f 2', pos='pos.sam', neg='neg.fastq', lib='reads.fq', paired_end=False, build_index=True)

    """

    def __init__(self, format_string, index_format_string=None,index_output_format_list=None):
        '''

        :param format_string:
        '''
        self.format_string = format_string
        self.index_format_string = index_format_string
        self.index_output_format_list = index_output_format_list
        self.logging = False
        aligner, builder = self.parse_format_string(format_string, index_format_string)
        try:
            if shutil.which(aligner) is None:
                raise AlignerDoesntExistError
            if builder is not None:
                if shutil.which(builder) is None:
                    raise BuilderDoesNotExistError
        except AlignerWraperError as e:
            print(e.message)
            return


    def parse_format_string(self, format_string, index_format_string=None):
        '''
        The function extracts the aligner, and flags, and whether the read libraries are single or paired end

        :param format_string:
        :return:
        '''
        builder = None
        try:
            if index_format_string is not None:
                builder_string = index_format_string.split()
                builder = builder_string[0]
                if re.search('{reference}', index_format_string) is None:
                    raise FormatStringError("Error: no reference for indexing")
                if re.search('{flags}', index_format_string) is None:
                    raise FormatStringError("Error: no flags for indexing")
            split_string =  format_string.split()
            aligner = split_string[0]

            lib1 = re.search('{lib_1}', format_string)
            lib2 = re.search('{lib_2}', format_string)
            if re.search('{lib}', format_string) is not None:
                if lib2 is not None or lib1 is not None:
                    raise FormatStringError("Error: no library selected")
            elif lib2 is None or lib1 is None:
                raise FormatStringError("Error: no library selected")
            if re.search('{reference_file}', format_string) is None:
                raise FormatStringError("Error: no reference file")
            if re.search('{positive_alignment}', format_string) is None: #TODO: make pos and neg optional?
                raise FormatStringError("Error: no positive alignment")
            if re.search('{negative_library}', format_string) is None:
                raise FormatStringError("Error: no negative library")
        except FormatStringError as e:
            print(e.message)
            return
        if re.search('{log}', format_string) is not None:
            self.logging = True

        return aligner, builder


    def align(self, reference, pos, neg, flags='', indexing_flags='', lib=None, lib_1=None, lib_2=None, log='', paired_end=False,
              locking_context=nullcontext(), build_index=False):
        '''
        :param flags: The flags the aligner will be run with
        :param reference: Reference nucleotides for reads to be aligned to
        :param pos: Positive alignment file
        :param neg: Negative alignments
        :param lib: Library of reads, for single-end case
        :param lib_1: First library of reads, for paired-end case
        :param lib_2: Second library of reads, for paired-end case
        :param paired_end: Set to True if using paired-end data
        :param locking_context: Use locking_context from threading_tree_utils.py if you need synchronization
        :param build_index: Set to True to build indexes for the reference library
        :param logging: Set to True to enable logging #TODO: is param still relevant?
        :return:
        '''
        #TODO: check validity of provided files
        #build_indexes
        #if not os.path.exists(fasta_file + ".fai"): #from bowtie_wrapper. TODO: Important?
        if build_index == True:
            with locking_context:
                # check if index already exists
                if not self.is_index_built(reference):
                    if self.index_format_string is None:
                        raise FormatStringError("Error: Index building not defined")
                        return
                    command = self.index_format_string.format(reference=reference, flags=indexing_flags)
                    self.execute_command(command)

        #align
        if paired_end:
            command = self.format_string.format(lib_1=lib_1, lib_2=lib_2, flags=flags, reference=reference, positive_alignment=pos, negative_library=neg, log=log )
        else:
            command = self.format_string.format(lib=lib, flags=flags, reference_file=reference, positive_alignment=pos, negative_library=neg, log=log)

        self.execute_command(command) #locking context?

    def is_index_built(self,reference):
        if self.index_format_string is None:
            raise FormatStringError("Index building not defined, no point in looking for indexing files")
            return
        index_files=[idx_f.format(reference_file=reference) for idx_f in self.index_output_format_list]

        all_indexes_exist=all([os.path.isfile(idx_f) for idx_f in index_files])

        return all_indexes_exist


    def execute_command(self, command):
        try:
            #print(command)
            subprocess.check_output(command, shell=True)  # this will continue when the call is done
        except subprocess.CalledProcessError as e:
            print(e.output)
            raise e

"""
Error classes for the Aligner Wrapper 
"""
class AlignerWraperError(Exception):
        """Parent class for other errors"""
        pass

class AlignerDoesntExistError(AlignerWraperError):
    def __init__(self):
        self.message = "Aligner does not exist"
        super().__init__(self.message)
        pass
    pass

class BuilderDoesNotExistError(AlignerWraperError):
    def __init__(self):
        self.message = "Index building tool does not exist"
        super().__init__(self.message)
        pass
    pass
"""
class for format string input errors.
:param error_message: error message to be displayed when error is raised
"""
class FormatStringError(AlignerWraperError):
    def __init__(self, error_message):
        self.message = error_message
        super().__init__(self.message)
        pass
    pass

if __name__ == '__main__':
    #doctest.testmod()
    main()


'''

def bowtie_wrapper(fastq_file: str, fasta_file: str, flags: str, output_name, locking_context):
    """

    :param fastq_file: the fastq file we want to align to fasta
    :param fasta_file: the fasta file we want to align to
    :param flags: the flags we want to use in the bowtie command
    :param output_name: the output name without the ending of file for the sam and unaligned fastqs
    :param locking_context: allow the funck to set lock for a critical code erea
    :return: the sam file name that is the alignment result , the fastq file name that contains the lines that didn't
    align to the fasta
    """

    if not os.path.exists(fasta_file + ".fai"):
        # creating index for the fasta file
        with locking_context:  # locking the case where another thread is also indexing at the same time
            if not os.path.exists(fasta_file + ".fai"):  # check if still the indexing was not made while we waited
                command = "bowtie-build {fasta} {fasta}".format(fasta=fasta_file)
                try:
                    subprocess.check_output(command, shell=True)  # this will continue when the call is done
                except subprocess.CalledProcessError as e:
                    print(e.output)
                    raise e

    # command2 = "samtools faidx " + fasta_file

    command = "bowtie {flags_} --sam {fasta_file_} {fastq_file_} {output}.sam  --un {output}.fastq >> " \
              "{output}.bowtie 2>&1".format(flags_=flags, fasta_file_=fasta_file, fastq_file_=fastq_file,
                                            output=output_name)
    try:
        subprocess.check_output(command, shell=True)  # this will continue when the call is done
    except subprocess.CalledProcessError as e:
        print(e.output)
        raise e

    return output_name + ".sam", output_name + ".fastq"

'''