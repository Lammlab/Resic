from Experiments.forontiers_jupyter.aligner_wrapper import AlignerWrapper
import Utility.threading_tree_utils



def bowtie_wrapper(fastq_file: str, fasta_file: str, flags: str, output_name, locking_context):
    '''
    :param fastq_file: the fastq file we want to align to fasta
    :param fasta_file: the fasta file we want to align to
    :param flags: the flags we want to use in the bowtie command
    :param output_name: the output name without the ending of file for the sam and unaligned fastqs
    :param locking_context: allow the funck to set lock for a critical code erea
    :return: the sam file name that is the alignment result , the fastq file name that contains the lines that didn't
    align to the fasta
    '''
    aligner_wrapper = AlignerWrapper("bowtie {flags} --sam {reference_file} {lib} {positive_alignment}.sam --un {negative_library}.fastq >> {log}.bowtie 2>&1", "bowtie-build {reference} {reference} {flags}")
    aligner_wrapper.align(reference=fasta_file, pos=output_name,neg=output_name, flags=flags, lib=fastq_file, log=output_name, build_index=True, locking_context=locking_context)

    return output_name + ".sam", output_name + ".fastq"


if __name__ == '__main__':
    positive_fastq =  '/home/lammlab/Documents/frontiers_test/files_for_test/test_test_mRNA.fastq'
    reference_library =  "/home/lammlab/Documents/frontiers_test/files_for_test/ws220-genes_expanded.fasta"
    negative_fastqs =  '/home/lammlab/Documents/frontiers_test/files_for_test/test_test_negative.fastq'
    pos,neg = bowtie_wrapper(fastq_file=positive_fastq, fasta_file=reference_library, flags=" -m 2 -n 3",)


