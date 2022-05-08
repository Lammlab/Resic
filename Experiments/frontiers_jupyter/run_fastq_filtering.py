"""
To run Resic, edit the real_pipe function (below) and run this file.
This file is the runfile for Resic's fastq_filtering step, in which input fastq files are filtered, filtering out
reads that meet the following criteria: one or more nucleotides represent over 60% or under 10% of the read sequence,
more than 10% Ns (when a base call could not be done), average Phred quality score < 25,
or more than 20 repeats of a single nucleotide in a row
"""

import logging
import os
import shutil
import Filtering.filter_ambiguous_reads as filter_read
import Processing.include_exclude_alignment as in_ex_align
from Experiments.frontiers_jupyter.pipe_utils import all_genome_pair_combinations, genome_3nt_all_combination_spec,\
    genome_3nt_factory
from Experiments.frontiers_jupyter.directory_structure_definer import DirectoryStructure


class PipeTester():
    """
    @param@ snp_database: list of file locations of snp database(s) in vcf format
    """

    def __init__(self, root_dir, positive_fastqs, negative_fastqs, fasta, graph_dict, spec_dict, group_dict, aligner,
                 parallel_limit, disable_parallel, skip_existing_files=False, snp_database=[]):
        self.root_dir = root_dir
        self.positive_fastqs = positive_fastqs
        self.negative_fastqs = negative_fastqs
        self.snp_database = snp_database
        self.fasta = fasta
        self.graph_dict = graph_dict
        self.spec_dict = spec_dict
        self.group_dict = group_dict
        self.aligner = aligner
        self.parallel_limit = parallel_limit
        self.disable_parallel = disable_parallel
        self.skip_existing_files = skip_existing_files
        # create directory structure
        self.dirstruct = DirectoryStructure(root_dir)

    def remove_files(self):
        shutil.rmtree(self.root_dir)

    def filter_fastq(self):

        root_dir = self.root_dir
        reference_library = self.fasta
        positive_fastqs = self.positive_fastqs
        negative_fastqs = self.negative_fastqs
        spec_dict = self.spec_dict
        graph_dict = self.graph_dict
        dirstruct = self.dirstruct
        skip_existing_files = self.skip_existing_files
        aligner = self.aligner

        for fastq in positive_fastqs + negative_fastqs:
            try:
                fastq_name = fastq.split("/")[-1]
            except:
                fastq_name = fastq

            temp = os.path.join(root_dir, fastq_name + "2")
            shutil.copyfile(fastq, temp)
            filter_read.extract_from_fastq(temp, fastq)
            os.remove(temp)


def Resic_graph_config(bowtie_parrallel=1):
    aligner = in_ex_align.bowtie_wrapper

    norep_align = f" -p {bowtie_parrallel} -m 2 -l 50 -n 3 --chunkmbs 200 --strata --best"
    repetative_align = f" -p {bowtie_parrallel} -m 20 -l 50 -n 3 --chunkmbs 200 --strata --best --max reads_exceeding_the_-m_limit.fastq"
    logging.error(f"{norep_align}")

    # building dict
    # no special processing
    spec_dict = dict()
    spec_dict["norep"] = (norep_align, None, None)
    spec_dict["rep"] = (repetative_align, None, None)

    # adding non rep and rep 3nt genome
    for base_pairs_name, pre, post in genome_3nt_all_combination_spec():
        spec_dict["norep_hyper_" + base_pairs_name] = (norep_align, pre, post)

    for base_pairs_name, pre, post in genome_3nt_all_combination_spec():
        spec_dict["rep_hyper_" + base_pairs_name] = (repetative_align, pre, post)

    # generate empty dict
    graph_dict = {name: [] for name in spec_dict.keys()}
    for key, fathers in graph_dict.items():
        # non rep is starting node
        if key == "norep":
            continue

        # rep is son of norep
        if key == "rep":
            fathers.append("norep")
            continue
        # if we reached here, the name has a base_pair
        key_split = key.split('_')
        base_pair_name = '_'.join(key_split[-2:])
        # hyper norep are all sons of rep
        if "norep_hyper_" in key:
            fathers.append("rep")
            continue

        # hyper rep is son of hyper norep with same bps
        if "rep_hyper_" in key:
            fathers.append("norep_hyper_" + base_pair_name)
            continue

    group_dict = {
        "norep": {"norep"},
        "rep": {"rep"},
        "norep_hyper": {"norep_hyper_" + "_".join(name) for name in all_genome_pair_combinations()},
        "rep_hyper": {"rep_hyper_" + "_".join(name) for name in all_genome_pair_combinations()},
    }
    return aligner, spec_dict, graph_dict, group_dict


def naive_3nt_graph_config():
    '''
    naive definition of graph for testing. **not for usage in a real pipeline.**
    '''
    aligner = in_ex_align.bowtie_wrapper

    norep_align = " -p 3 -m 2 -n 3"
    repetative_align = " -p 3 -m 100 -n 3 "

    # building dict
    # 3nt pre proccessing with A-G

    pre, post = genome_3nt_factory('A', 'G')
    spec_dict = dict()
    spec_dict["norep"] = (norep_align, pre, post)
    spec_dict["rep"] = (repetative_align, pre, post)

    graph_dict = {name: [] for name in spec_dict.keys()}

    for key, fathers in graph_dict.items():
        # non rep is starting node
        if key == "norep":
            continue

        # rep is son of norep
        if key == "rep":
            fathers.append("norep")
            continue

    group_dict = {
        "norep": {"norep"},
        "rep": {"rep"},
    }
    return aligner, spec_dict, graph_dict, group_dict


def real_pipe():
    # this is the main function.
    # all you need to do is to add the paths to your files and run this file

    # this option reverse the negative and positive fastq samples. It is used for quality control
    # do NOT change it for a regular run.
    reverse = False

    # the output directory 
    root_dir = "/Data2/users/clara/resic_split/root_dir"

    # no need to add a path for a regular run.
    if reverse:
        root_dir += ""

    # the directory in which the fastq input files are located
    data_dir = "/Data2/users/clara/resic_fix/root_dir"
    # a path to a fasta file to use as a reference genome (either indexed or not indexed)
    reference_library = "/Data2/reference_genomes/ws220/ws220-allchro_with_chr_M.fasta"

    # a list of fastq sample names for which you want to detect editing sites
    positive_fastqs = ['N2_E_27_Orna_illumina_rep_1.united.fastq.collapsed',
                       'N2_E_27_Orna_illumina_rep_2.united.fastq.collapsed']

    # Optional-  a list of fastq sample names that will be used to exclude non-editing sites.
    # these can be DNA reads or mutant strains that are lacking the editing mechanism.               
    negative_fastqs = ['BB21_E_27_Orna_illumina.united.fastq.collapsed',
                       'BB21_E_29_Alla_illumina.united.fastq.collapsed']

    # Optional - a dbSNP vcf file, currently WITH NO HEADERS, to exclude SNPs.
    # remove all lines that start with #
    snp_database = []

    # Note: currently we support using either a snp database or negative files and not both at once.

    if not os.path.exists(root_dir):
        os.makedirs(root_dir)

    for fastq in positive_fastqs + negative_fastqs:
        if not os.path.exists(os.path.join(root_dir, fastq)):
            shutil.copyfile(os.path.join(data_dir, fastq), os.path.join(root_dir, fastq))

    positive_fastqs = [os.path.join(root_dir, fastq) for fastq in positive_fastqs]

    negative_fastqs = [os.path.join(root_dir, fastq) for fastq in negative_fastqs]

    if reverse:
        positive_fastqs, negative_fastqs = (negative_fastqs, positive_fastqs)

    # Resic's ful graph
    aligner, spec_dict, graph_dict, group_dict = Resic_graph_config(bowtie_parrallel=12)

    # naive is non rep and rep
    # aligner,spec_dict,graph_dict,group_dict=naive_graph_config()
    # naive 3nt is non rep and rep with AG 3nt pre and post
    # aligner,spec_dict,graph_dict,group_dict=naive_3nt_graph_config()

    test_pipe = PipeTester(root_dir, positive_fastqs, negative_fastqs, reference_library,
                           graph_dict, spec_dict, group_dict, aligner, parallel_limit=6, disable_parallel=False,
                           skip_existing_files=True, snp_database=snp_database)

    """
    You can comment function that you want to skip. All functions together creats RECIC
    Thresholds can be modified.
    """
    # comment the test_pipe.remove_files if you want to keep previous files in the output directory

    # test_pipe.remove_files()

    test_pipe.filter_fastq()



if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    real_pipe()
