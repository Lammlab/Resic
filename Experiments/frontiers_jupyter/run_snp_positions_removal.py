"""
To run resic, edit the real_pipe function (below) and run this file.
This is the snp_positions_removal runfile, in which positions that are know SNP locations (present in the input SNP
database) will be filtered out.
"""

import logging
import os
import shutil
import Processing.include_exclude_alignment as in_ex_align
from Filtering.filter_pileup_by_multiple_existing_snps import snp_algebra, snp_detect, snp_algebra_from_vcf_file
from Experiments.frontiers_jupyter.pipe_utils import all_genome_pair_combinations, genome_3nt_all_combination_spec
from Experiments.frontiers_jupyter.directory_structure_definer import DirectoryStructure, Stages, PileupStage
from Experiments.frontiers_jupyter.parallel_commands import parallel_commands
from Utility.Pileup_class import Pileup_line
from Utility.generators_utilities import class_generator
from Utility.parallel_generator import parallel_generator


class PipeTester():
    """
    @param@ snp_database: list of file locations of snp database(s) in vcf format
    """

    def __init__(self, root_dir, positive_fastqs, negative_fastqs, fasta, graph_dict, spec_dict, group_dict, aligner,
                 parallel_limit, Disable_parallel, skip_existing_files=False, snp_database=[]):
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
        self.Disable_parallel = Disable_parallel
        self.skip_existing_files = skip_existing_files
        # create directory structure
        self.dirstruct = DirectoryStructure(root_dir)

    def remove_files(self):
        shutil.rmtree(self.root_dir)


    def snp_removal_test(self):

        reference_library = self.fasta
        positive_fastqs = self.positive_fastqs
        negative_fastqs = self.negative_fastqs
        snp_database = self.snp_database
        spec_dict = self.spec_dict
        graph_dict = self.graph_dict
        dirstruct = self.dirstruct
        skip_existing_files = self.skip_existing_files

        # all negative pileups from all nodes
        negative_pileups = [dirstruct.pathName(neg, node, Stages.pileup_generation, PileupStage.sorted_pileup) \
                            for neg in negative_fastqs for node in graph_dict.keys()]

        # combine all negatives into single filtered negative, to save runtime in snp_removal
        # parameters for parallel generator

        # TODO: I added the if:
        neg_obj_list = [open(p) for p in negative_pileups if os.path.isfile(p)]
        # neg_obj_list = [open(p) for p in negative_pileups]
        neg_gen_list = []
        get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)

        not_snp_database = True
        if snp_database != []:
            not_snp_database = False

        if not_snp_database:
            # filter out bad lines from negatives
            filtered_negatives = self.root_dir + "/" + "ALL_NEGATIVES_UNFILTERED.pileup"
            for file_obj in neg_obj_list:
                neg_gen_list.append(class_generator(Pileup_line, file=file_obj))
            with open(filtered_negatives, 'w') as out:
                parallel_gen = parallel_generator(neg_gen_list, [get_pos_and_id for i in range(len(neg_gen_list))])
                for parallel_line_list in parallel_gen:
                    snp = False
                    for pileup_list in parallel_line_list:
                        if pileup_list is not None:
                            if pileup_list[0].is_with_any_change():
                                snp = True

                                break
                    if not snp:
                        line_to_write = next(item for item in parallel_line_list if item is not None)
                        out.write(str(line_to_write[0]) + "\n")
            for obj in neg_obj_list:
                obj.close()

        commands = []
        for fastq in positive_fastqs:
            for node in graph_dict.keys():
                pileup_name = dirstruct.pathName(fastq, node, Stages.read_threshold, need_suffix=True)
                filtered_pileup_name = dirstruct.pathName(fastq, node, Stages.snp_removal, need_suffix=True)

                if skip_existing_files and os.path.isfile(filtered_pileup_name):
                    logging.info(
                        f"SKIP pileup snp removal filtering for file {pileup_name} since {filtered_pileup_name} already exists")
                    continue

                # TODO: I added this:
                if not os.path.isfile(pileup_name):
                    logging.info(
                        f"SKIP pileup snp removal filtering for file {pileup_name} since {pileup_name} doesn't exist / was not created")
                    continue

                # snp algebre interface is (pos,neg_list,snp_detector_func,out_put_name,is_input_sorted,not_snp_database)
                if not_snp_database:
                    command = [snp_algebra, pileup_name, [filtered_negatives], snp_detect, filtered_pileup_name, True,
                               False]
                else:
                    command = [snp_algebra_from_vcf_file, pileup_name, snp_database[0], filtered_pileup_name, True]
                commands.append(command)
        parallel_commands(commands, parallel_limit=self.parallel_limit, Disable_parallel=self.Disable_parallel)

        # snp algebre interface is (pos,neg_list,snp_detector_func,out_put_name,is_input_sorted)


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
                           graph_dict, spec_dict, group_dict, aligner, parallel_limit=6, Disable_parallel=False,
                           skip_existing_files=True, snp_database=snp_database)

    """
    You can comment function that you want to skip. All functions together creats RECIC
    Thresholds can be modified.
    """
    # comment the test_pipe.remove_files if you want to keep previous files in the output directory
    # test_pipe.remove_files()

    test_pipe.snp_removal_test()


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    real_pipe()
