"""
To run resic, edit the real_pipe function (below) and run this file.
This file is the runfile for the editing percent summary creation step, in which files
containing statistics about the different editing rates in each library will be produced.

"""

import logging
import os
import shutil
import itertools
import Processing.include_exclude_alignment as in_ex_align
from Processing.analyze_editing_percent import analyse_multiple_editing_percent_files
from Experiments.frontiers_jupyter.pipe_utils import all_genome_pair_combinations, genome_3nt_all_combination_spec
from Experiments.frontiers_jupyter.directory_structure_definer import \
DirectoryStructure, Stages, ConcensusStage, EditTypeStage

global_list_colors = ['#d12325', '#7C00BE', '#55E5D8', '#41C109', '#DEF407', '#9F3A7D', '#3C2C0D', '#089E4C', '#F69804',
                      '#2482BD', '#A0A269', '#CC1D74', '#1DA6CC', '#d6adc6', '#018439', '#3939a7', '#37462a', '#11dacd',
                      '#ab6b5f', '#1eb1cf', '#1a48a7']


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
        self.Disable_parallel = disable_parallel
        self.skip_existing_files = skip_existing_files
        # create directory structure
        self.dirstruct = DirectoryStructure(root_dir)

    def remove_files(self):
        """
        This function cleans the root output dir
        :return: None
        """
        shutil.rmtree(self.root_dir)


    def editing_percent_plot_generate_summaries(self, editing_min_thresh=0, editing_max_thresh=100, noise_thresh=1,
                                                read_thresh=2):
        reference_library = self.fasta
        positive_fastqs = self.positive_fastqs
        spec_dict = self.spec_dict
        graph_dict = self.graph_dict
        dirstruct = self.dirstruct
        skip_existing_files = self.skip_existing_files
        group_dict = self.group_dict

        nodes = list(graph_dict.keys())
        # calculate editing percent and editing site distribution summary for each lib and node
        in_pileups = [dirstruct.pathName(fastq, node, Stages.consensus, ConcensusStage.filtered) for node, fastq
                      in itertools.product(nodes, positive_fastqs)]
        out_pileups = [dirstruct.pathName(fastq, node, Stages.editing_type_count, EditTypeStage.edit_percent_pileup) for
                       node, fastq
                       in itertools.product(nodes, positive_fastqs)]
        out_summaries = [dirstruct.pathName(fastq, node, Stages.editing_type_count, EditTypeStage.file_summary) for
                         node, fastq
                         in itertools.product(nodes, positive_fastqs)]

        analyse_multiple_editing_percent_files(in_pileups, out_pileups, out_summaries, total_summary_file=None,
                                               add_headers=True, summary_only=False, min_editing=editing_min_thresh,
                                               max_editing=editing_max_thresh,
                                               max_noise=noise_thresh, min_reads=read_thresh, edit_tag='edited',
                                               parallel_limit=self.parallel_limit,
                                               Disable_parallel=self.Disable_parallel)


def Resic_graph_config(bowtie_parrallel=1):
    """
    Configure parameters for resic run's alignment graph
    :param bowtie_parrallel: number of threads to run aligner with (-p parameter, default: 1)
    :return:
    """
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


def real_pipe():
    # this is the main function.
    # all you need to do is to add the paths to your files and run this file

    # this option reverse the negative and positive fastq samples. It is used for quality control
    # do NOT change it for a regular run.
    reverse = False

    # the output directory 
    root_dir = ""

    # no need to add a path for a regular run.
    if reverse:
        root_dir += ""

    # the directory in which the fastq input files are located
    data_dir = ""
    # a path to a fasta file to use as a reference genome (either indexed or not indexed)
    reference_library = ""

    # a list of fastq sample names for which you want to detect editing sites
    positive_fastqs = []

    # Optional-  a list of fastq sample names that will be used to exclude non-editing sites.
    # these can be DNA reads or mutant strains that are lacking the editing mechanism.               
    negative_fastqs = []

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
    You can comment function that you want to skip. All functions together create RESIC
    Thresholds can be modified.
    """
    # comment the test_pipe.remove_files if you want to keep previous files in the output directory
    # test_pipe.remove_files()

    test_pipe.editing_percent_plot_generate_summaries(editing_min_thresh=30, editing_max_thresh=99, noise_thresh=3,
                                                      read_thresh=2)


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    real_pipe()
