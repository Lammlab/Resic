"""
To run resic, edit the real_pipe function (below) and run this file.
This is the runfile for the pileup_creation step, in which sam files are sorted and converted into corresponding
pileup files.
"""

import logging
import os
import shutil
from Processing.merge_sams import merge_sams
import subprocess
import Processing.include_exclude_alignment as in_ex_align
from Experiments.frontiers_jupyter.pipe_utils import all_genome_pair_combinations, genome_3nt_all_combination_spec
from Experiments.frontiers_jupyter.directory_structure_definer import DirectoryStructure, Stages, AlignStage, \
PileupStage
from Experiments.frontiers_jupyter.parallel_commands import parallel_commands


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


    def pileup_creation_test(self):

        reference_library = self.fasta
        positive_fastqs = self.positive_fastqs
        negative_fastqs = self.negative_fastqs
        spec_dict = self.spec_dict
        graph_dict = self.graph_dict
        dirstruct = self.dirstruct
        skip_existing_files = self.skip_existing_files

        def sam_to_pileup(sam_name, fasta_name, bam_name, sorted_bam_name, pileup_name, sorted_pileup_name):
            # TODO: changed this:
            # command=(   f"oldsamtools view -bS {sam_name} > {bam_name} && " +
            #             f"oldsamtools sort {bam_name} {sorted_bam_name} && " +
            #             f"oldsamtools mpileup -f {fasta_name} {sorted_bam_name}.bam > {pileup_name} &&" + #| tail -n +3 
            #             f"sed -r '/^[\t]*$/d' <{pileup_name} | sort -k1,1 -k2,2n -o {sorted_pileup_name} "
            #             )
            command = (f"samtools view -bS {sam_name} > {bam_name} && " +
                       f"samtools sort -o {sorted_bam_name} {bam_name} && " +
                       f"samtools mpileup -f {fasta_name} {sorted_bam_name} > {pileup_name} &&" +  # | tail -n +3
                       f"sed -r '/^[\t]*$/d' <{pileup_name} | sort -k1,1 -k2,2n -o {sorted_pileup_name} "
                       )
            try:
                outout = subprocess.check_output(command, shell=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"pileup conversions of {sam_name} failed")
                raise e

            return

        def sams_to_pileup(fasta_name, sam_name, bam_name, sorted_bam_name, pileup_name, sorted_pileup_name,
                           antisense_sam, combined_sam_name):
            merge_sams([sam_name, antisense_sam], combined_sam_name)
            sam_to_pileup(combined_sam_name, fasta_name, bam_name, sorted_bam_name, pileup_name, sorted_pileup_name)
            return

        commands = []
        for fastq in positive_fastqs + negative_fastqs:
            for node in graph_dict.keys():
                sam_name = dirstruct.pathName(fastq, node, Stages.graph_aligner, AlignStage.post_sam, True)
                bam_name = dirstruct.pathName(fastq, node, Stages.pileup_generation, PileupStage.bam, True)
                sorted_bam_name = dirstruct.pathName(fastq, node, Stages.pileup_generation, PileupStage.sorted_bam,
                                                     True)
                pileup_name = dirstruct.pathName(fastq, node, Stages.pileup_generation, PileupStage.pileup, True)
                sorted_pileup_name = dirstruct.pathName(fastq, node, Stages.pileup_generation,
                                                        PileupStage.sorted_pileup, True)
                antisense_sam = dirstruct.pathName(fastq, node, Stages.graph_aligner, AlignStage.antisense_post_sam,
                                                   True)
                combined_sam_name = dirstruct.pathName(fastq, node, Stages.pileup_generation, PileupStage.combined_sam,
                                                       True)
                fasta = reference_library

                if skip_existing_files and os.path.isfile(sorted_pileup_name):
                    logging.info(f"SKIP pileup creation for file {sam_name} since {sorted_pileup_name} already exists")
                    continue
                # TODO: I added this if 10.10
                # if node == "norep" or node == "rep":
                if os.path.isfile(sam_name) and (node == "norep" or node == "rep"):
                    command = [sam_to_pileup, sam_name, fasta, bam_name, sorted_bam_name, pileup_name,
                               sorted_pileup_name]
                else:
                    if (node != "norep" and node != "rep") and os.path.isfile(antisense_sam) and os.path.isfile(
                            sam_name):
                        command = [sams_to_pileup, fasta, sam_name, bam_name, sorted_bam_name, pileup_name,
                                   sorted_pileup_name, antisense_sam, combined_sam_name]
                    else:
                        logging.info(f"SKIP pileup creation for file {sam_name}, file was not created")
                        continue

                commands.append(command)

        parallel_commands(commands, parallel_limit=self.parallel_limit, Disable_parallel=self.Disable_parallel)


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
    root_dir = "/Data2/users/clara/resic_split/root_dir"

    # no need to add a path for a regular run.
    if reverse:
        root_dir += ""

    # the directory in which the fastq input files are located
    data_dir = "/Data2/users/berta/Resic_test_output"
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
                           graph_dict, spec_dict, group_dict, aligner, parallel_limit=6, Disable_parallel=False,
                           skip_existing_files=True, snp_database=snp_database)

    """
    You can comment function that you want to skip. All functions together creats RECIC
    Thresholds can be modified.
    """
    # comment the test_pipe.remove_files if you want to keep previous files in the output directory
    # test_pipe.remove_files()

    test_pipe.pileup_creation_test()


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    real_pipe()
