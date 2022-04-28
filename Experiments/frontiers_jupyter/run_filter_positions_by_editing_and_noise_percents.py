"""
To run resic all you need to do is to edit the real_pipe function (below) and run this file.
"""

import logging
import os
import shutil
import Processing.include_exclude_alignment as in_ex_align
from Processing.analyze_editing_percent import filter_pileup_by_categories
from Experiments.frontiers_jupyter.pipe_utils import all_genome_pair_combinations, genome_3nt_all_combination_spec
from Experiments.frontiers_jupyter.directory_structure_definer import DirectoryStructure, Stages
from Experiments.frontiers_jupyter.parallel_commands import parallel_commands


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

    # def filter_fastq(self):

    #     root_dir = self.root_dir
    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     negative_fastqs = self.negative_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files
    #     aligner=self.aligner

    #     for fastq in positive_fastqs + negative_fastqs:
    #         try:
    #             fastq_name = fastq.split("/")[-1]
    #         except:
    #             fastq_name = fastq

    #         temp = os.path.join(root_dir, fastq_name + "2")
    #         shutil.copyfile(fastq, temp)
    #         filter_read.extract_from_fastq(temp, fastq)
    #         os.remove(temp)

    # def include_exclude_alignment_test(self):
    #     # calling the function that run's the analysis (alignments)
    #     root_dir = self.root_dir
    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     negative_fastqs = self.negative_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files
    #     aligner=self.aligner

    #     for fastq in positive_fastqs + negative_fastqs:
    #         in_ex_align.do_include_exclude_alignment(fastq, reference_library
    #             ,spec_dict, graph_dict, dirstruct, aligner, skip_existing_files)

    # def pileup_creation_test(self):

    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     negative_fastqs = self.negative_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files

    #     def sam_to_pileup(sam_name,fasta_name,bam_name,sorted_bam_name,pileup_name,sorted_pileup_name):
    #         # TODO: changed this:
    #         # command=(   f"oldsamtools view -bS {sam_name} > {bam_name} && " +
    #         #             f"oldsamtools sort {bam_name} {sorted_bam_name} && " +
    #         #             f"oldsamtools mpileup -f {fasta_name} {sorted_bam_name}.bam > {pileup_name} &&" + #| tail -n +3 
    #         #             f"sed -r '/^[\t]*$/d' <{pileup_name} | sort -k1,1 -k2,2n -o {sorted_pileup_name} "
    #         #             )
    #         command=(   f"samtools view -bS {sam_name} > {bam_name} && " +
    #         f"samtools sort -o {sorted_bam_name} {bam_name} && " +
    #         f"samtools mpileup -f {fasta_name} {sorted_bam_name} > {pileup_name} &&" + #| tail -n +3 
    #         f"sed -r '/^[\t]*$/d' <{pileup_name} | sort -k1,1 -k2,2n -o {sorted_pileup_name} "
    #         )
    #         try:
    #             outout=subprocess.check_output(command,shell=True)
    #         except subprocess.CalledProcessError as e:
    #             logging.error(f"pileup conversions of {sam_name} failed")
    #             raise e

    #         return
    #     def sams_to_pileup(fasta_name,sam_name,bam_name,sorted_bam_name,pileup_name,sorted_pileup_name,antisense_sam, combined_sam_name):
    #         merge_sams([sam_name, antisense_sam], combined_sam_name)
    #         sam_to_pileup(combined_sam_name,fasta_name,bam_name,sorted_bam_name,pileup_name,sorted_pileup_name)
    #         return

    #     commands=[]
    #     for fastq in positive_fastqs+negative_fastqs:
    #         for node in graph_dict.keys():
    #             sam_name=dirstruct.pathName(fastq,node,Stages.graph_aligner,AlignStage.post_sam,True)
    #             bam_name=dirstruct.pathName(fastq,node,Stages.pileup_generation,PileupStage.bam,True)
    #             sorted_bam_name=dirstruct.pathName(fastq,node,Stages.pileup_generation,PileupStage.sorted_bam,True)
    #             pileup_name=dirstruct.pathName(fastq,node,Stages.pileup_generation,PileupStage.pileup,True)
    #             sorted_pileup_name=dirstruct.pathName(fastq,node,Stages.pileup_generation,PileupStage.sorted_pileup,True)
    #             antisense_sam=dirstruct.pathName(fastq,node,Stages.graph_aligner,AlignStage.antisense_post_sam,True)
    #             combined_sam_name=dirstruct.pathName(fastq,node,Stages.pileup_generation,PileupStage.combined_sam,True)
    #             fasta= reference_library

    #             if skip_existing_files and os.path.isfile(sorted_pileup_name):
    #                 logging.info(f"SKIP pileup creation for file {sam_name} since {sorted_pileup_name} already exists")
    #                 continue
    #             # TODO: I added this if 10.10
    #             #if node == "norep" or node == "rep":
    #             if os.path.isfile(sam_name) and (node == "norep" or node == "rep"):
    #                 command=[sam_to_pileup,sam_name,fasta,bam_name,sorted_bam_name,pileup_name,sorted_pileup_name]
    #             else:
    #                 if (node != "norep" and node != "rep") and os.path.isfile(antisense_sam) and os.path.isfile(sam_name):
    #                     command=[sams_to_pileup,fasta,sam_name,bam_name,sorted_bam_name,pileup_name,sorted_pileup_name,antisense_sam,combined_sam_name]
    #                 else:
    #                     logging.info(f"SKIP pileup creation for file {sam_name}, file was not created")
    #                     continue

    #             commands.append(command)

    #     parallel_commands(commands,parallel_limit=self.parallel_limit,disable_parallel=self.disable_parallel)

    # def filter_nochange_test(self):
    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files

    #     def filter_no_change(pileup,filtered_pileup):
    #         filter_pileup_by_categories(pileup,filtered_pileup,None,True,None,None)
    #         return

    #     commands=[]
    #     for fastq in positive_fastqs:
    #         for node in graph_dict.keys():
    #             pileup_name=dirstruct.pathName(fastq,node,Stages.pileup_generation,PileupStage.sorted_pileup,True)
    #             filtered_pileup_name=dirstruct.pathName(fastq,node,Stages.no_change,need_suffix=True)

    #             if skip_existing_files and os.path.isfile(filtered_pileup_name):
    #                 logging.info(f"SKIP pileup no change filtering for file {pileup_name} since {filtered_pileup_name} already exists")
    #                 continue

    #             # Todo I added this:
    #             if not os.path.isfile(pileup_name):
    #                 logging.info(f"SKIP pileup no change filtering for file {pileup_name} since {pileup_name} was not created / does not exist")
    #                 continue

    #             command=[filter_no_change,pileup_name,filtered_pileup_name]
    #             commands.append(command)

    #     parallel_commands(commands,parallel_limit=self.parallel_limit,disable_parallel=self.disable_parallel)

    # def filter_readthreshold_test(self, threshold=1):

    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files

    #     def filter_by_threshold(pileup,filtered_pileup,thresh):
    #         filter_pileup_by_categories(pileup,filtered_pileup,thresh,None,None,None)
    #         return

    #     commands=[]
    #     for fastq in positive_fastqs:
    #         for node in graph_dict.keys():
    #             pileup_name=            dirstruct.pathName(fastq,node,Stages.no_change,need_suffix=True)
    #             filtered_pileup_name=   dirstruct.pathName(fastq,node,Stages.read_threshold,need_suffix=True)

    #             if skip_existing_files and os.path.isfile(filtered_pileup_name):
    #                 logging.info(f"SKIP pileup read threshold filtering for file {pileup_name} since {filtered_pileup_name} already exists")
    #                 continue
    #             # Todo: I added this:
    #             if not os.path.isfile(pileup_name):
    #                 logging.info(f"SKIP pileup read threshold filtering for file {pileup_name} since {pileup_name} does not exist / was not created")
    #                 continue

    #             command=[filter_by_threshold,pileup_name,filtered_pileup_name,threshold]
    #             commands.append(command)

    #     parallel_commands(commands,parallel_limit=self.parallel_limit,disable_parallel=self.disable_parallel)

    # def snp_removal_test(self):

    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     negative_fastqs = self.negative_fastqs
    #     snp_database = self.snp_database
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files

    #      # all negative pileups from all nodes
    #     negative_pileups=[dirstruct.pathName(neg,node,Stages.pileup_generation,PileupStage.sorted_pileup) \
    #                              for neg in negative_fastqs for node in graph_dict.keys()]

    #     #combine all negatives into single filtered negative, to save runtime in snp_removal
    #     #parameters for parallel generator

    #     # TODO: I added the if:
    #     neg_obj_list = [open(p) for p in negative_pileups if os.path.isfile(p)]
    #     #neg_obj_list = [open(p) for p in negative_pileups]
    #     neg_gen_list = []
    #     get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)

    #     not_snp_database = True
    #     if snp_database != []:
    #         not_snp_database = False

    #     if not_snp_database:
    #         #filter out bad lines from negatives
    #         filtered_negatives = self.root_dir +"/" +"ALL_NEGATIVES_UNFILTERED.pileup"
    #         for file_obj in neg_obj_list:
    #             neg_gen_list.append(class_generator(Pileup_line, file = file_obj))
    #         with open(filtered_negatives, 'w') as out:
    #             parallel_gen = parallel_generator(neg_gen_list, [get_pos_and_id for i in range(len(neg_gen_list))])
    #             for parallel_line_list in parallel_gen: 
    #                 snp = False
    #                 for pileup_list in parallel_line_list:
    #                     if pileup_list is not None:
    #                        if pileup_list[0].is_with_any_change():
    #                             snp = True

    #                             break
    #                 if not snp:
    #                     line_to_write = next(item for item in parallel_line_list if item is not None)
    #                     out.write(str(line_to_write[0]) + "\n")
    #         for obj in neg_obj_list:
    #             obj.close()

    #     commands=[]
    #     for fastq in positive_fastqs:
    #         for node in graph_dict.keys():
    #             pileup_name=            dirstruct.pathName(fastq,node,Stages.read_threshold,need_suffix=True)
    #             filtered_pileup_name=   dirstruct.pathName(fastq,node,Stages.snp_removal,need_suffix=True)

    #             if skip_existing_files and os.path.isfile(filtered_pileup_name):
    #                 logging.info(f"SKIP pileup snp removal filtering for file {pileup_name} since {filtered_pileup_name} already exists")
    #                 continue

    #             # TODO: I added this:
    #             if not os.path.isfile(pileup_name):
    #                 logging.info(f"SKIP pileup snp removal filtering for file {pileup_name} since {pileup_name} doesn't exist / was not created")
    #                 continue

    #             # snp algebre interface is (pos,neg_list,snp_detector_func,out_put_name,is_input_sorted,not_snp_database)
    #             if not_snp_database:
    #                 command=[snp_algebra,pileup_name,[filtered_negatives],snp_detect,filtered_pileup_name,True, False]
    #             else:
    #                 command=[snp_algebra_from_vcf_file,pileup_name,snp_database[0],filtered_pileup_name,True]
    #             commands.append(command)
    #     parallel_commands(commands,parallel_limit=self.parallel_limit,disable_parallel=self.disable_parallel)

    #             # snp algebre interface is (pos,neg_list,snp_detector_func,out_put_name,is_input_sorted)

    def editing_percent_test(self, editing_min_threshold=30, editing_max_threshold=99, noise_threshold=3,
                             editing_read_thresh=2):

        reference_library = self.fasta
        positive_fastqs = self.positive_fastqs
        spec_dict = self.spec_dict
        graph_dict = self.graph_dict
        dirstruct = self.dirstruct
        skip_existing_files = self.skip_existing_files

        def filter_editperc(pileup, filtered_pileup, edit_min_thresh, edit_max_thresh, noise_thresh):
            filter_pileup_by_categories(pileup, filtered_pileup, None, None, edit_min_thresh, edit_max_thresh,
                                        noise_thresh, editing_read_thresh)
            return

        commands = []
        for fastq in positive_fastqs:
            for node in graph_dict.keys():
                pileup_name = dirstruct.pathName(fastq, node, Stages.snp_removal, need_suffix=True)
                filtered_pileup_name = dirstruct.pathName(fastq, node, Stages.editing_percent, need_suffix=True)

                if skip_existing_files and os.path.isfile(filtered_pileup_name):
                    logging.info(
                        f"SKIP pileup editing percent filtering for file {pileup_name} since {filtered_pileup_name} already exists")
                    continue
                # Todo: I added this if:
                if not os.path.isfile(pileup_name):
                    logging.info(
                        f"SKIP pileup editing percent filtering for file {pileup_name} since {pileup_name} does not exist / was not created")
                    continue

                command = [filter_editperc, pileup_name, filtered_pileup_name, editing_min_threshold,
                           editing_max_threshold, noise_threshold]
                commands.append(command)

        parallel_commands(commands, parallel_limit=self.parallel_limit, disable_parallel=self.disable_parallel)

    # def hyper_non_relevant_editing_site_test(self):
    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files

    #     all_pileup_files_from_last_stage = []
    #     for fastq in positive_fastqs:
    #         for node in graph_dict.keys():
    #             #all_pileup_files_from_last_stage.append(dirstruct.pathName(fastq,node,Stages.editing_percent,need_suffix=True))
    #             # Todo: I added this if:
    #             if os.path.isfile(dirstruct.pathName(fastq,node,Stages.editing_percent,need_suffix=True)):
    #                 all_pileup_files_from_last_stage.append(dirstruct.pathName(fastq,node,Stages.editing_percent,need_suffix=True))
    #             else:
    #                 logging.info(f"SKIP pileup hyper non-relevant filtering for file {dirstruct.pathName(fastq,node,Stages.editing_percent,need_suffix=True)} since this file does not exist / was not created")

    #     filter_hyper_non_relevant_editing_sites(all_pileup_files_from_last_stage)

    # def unique_site_test(self):
    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files

    #     commands=[]
    #     for fastq in positive_fastqs:
    #         pileup_files_from_all_nodes = []
    #         for node in graph_dict.keys():
    #             pileup_files_from_all_nodes.append(dirstruct.pathName(fastq,node,Stages.editing_percent,need_suffix=True))
    #         write_unique_sites(pileup_files_from_all_nodes, sorted_input=True)

    # def concensus_test(self,consensus_threshold=0.5):

    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files

    #     commands=[]
    #     for node in graph_dict.keys():
    #         pileups=            [dirstruct.pathName(fastq,node,Stages.editing_percent_unique) for fastq in positive_fastqs ]
    #         filtered_pileups=   [dirstruct.pathName(fastq,node,Stages.concensus,ConcensusStage.filtered) for fastq in positive_fastqs ]
    #         concensus_pileup=   dirstruct.pathName(None,node,Stages.concensus,ConcensusStage.concensus)

    #         if skip_existing_files and \
    #                 all([os.path.isfile(file) for file in filtered_pileups+[concensus_pileup]]):

    #             logging.info(
    #                 f"SKIP pileup concensus filtering for file {node} since {filtered_pileups + [concensus_pileup]} already exists")
    #             continue

    #         # format for filter by concensus is (in_pileups,thresh,out_pileups,concensus_pileup,is_sorted)
    #         command=[filter_by_consensus,pileups,consensus_threshold,filtered_pileups,concensus_pileup,True]
    #         commands.append(command)

    #     parallel_commands(commands,parallel_limit=self.parallel_limit,disable_parallel=self.disable_parallel)

    # def site_loss_plot_test(self):

    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files
    #     group_dict=self.group_dict

    #     dict_colors = calculate_dictinary_colors_for_site_loss()

    #     commands=[]
    #     node_names=list(graph_dict.keys())
    #     for fastq in positive_fastqs:
    #         command=[site_loss_by_group_plot,fastq,node_names,group_dict,dirstruct,dict_colors]
    #         commands.append(command)

    #     parallel_commands(commands,parallel_limit=self.parallel_limit,disable_parallel=self.disable_parallel)

    # def editing_percent_plot_generate_summaries(self, editing_min_thresh=0, editing_max_thresh=100, noise_thresh=1, read_thresh=2):

    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files
    #     group_dict=self.group_dict

    #     nodes=list(graph_dict.keys())
    #     # calculate editing percent and editing site distribution summary for each lib and node
    #     in_pileups=[dirstruct.pathName(fastq,node,Stages.concensus,ConcensusStage.filtered) for node,fastq
    #                 in itertools.product(nodes,positive_fastqs) ]
    #     out_pileups=[dirstruct.pathName(fastq,node,Stages.editing_type_count,EditTypeStage.edit_percent_pileup) for node,fastq
    #                 in itertools.product(nodes, positive_fastqs)]
    #     out_summaries=[dirstruct.pathName(fastq,node,Stages.editing_type_count,EditTypeStage.file_summary) for node,fastq
    #                 in itertools.product(nodes, positive_fastqs)]

    #     analyse_multiple_editing_percent_files(in_pileups, out_pileups, out_summaries, total_summary_file=None,
    #                      add_headers=True, summary_only=False, min_editing=editing_min_thresh, max_editing=editing_max_thresh,
    #                      max_noise=noise_thresh, min_reads=read_thresh, edit_tag='edited', parallel_limit=self.parallel_limit,
    #                      disable_parallel=self.disable_parallel)

    # def editing_percent_plot_generate_plot(self):
    #     reference_library = self.fasta
    #     positive_fastqs = self.positive_fastqs
    #     spec_dict = self.spec_dict
    #     graph_dict = self.graph_dict
    #     dirstruct = self.dirstruct
    #     skip_existing_files=self.skip_existing_files
    #     group_dict = self.group_dict

    #     dict_colors = calculate_dictinary_colors_for_editing_percent(dirstruct, positive_fastqs, group_dict)

    #     # Generate group summary files and plots
    #     commands = []
    #     # input : the pileups after the intire thing - then do editing percent anaysis
    #     for fastq in positive_fastqs:
    #         command = [editing_type_count_by_group_plot,fastq,group_dict,dirstruct,dict_colors]
    #         commands.append(command)

    #     parallel_commands(commands,parallel_limit=self.parallel_limit,disable_parallel=self.disable_parallel)


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


# def naive_3nt_graph_config():
#     '''
#     naive definition of graph for testing. **not for usage in a real pipeline.**
#     '''
#     aligner=in_ex_align.bowtie_wrapper

#     norep_align = " -p 3 -m 2 -n 3"
#     repetative_align = " -p 3 -m 100 -n 3 "

#     # building dict
#     # 3nt pre proccessing with A-G

#     pre,post=genome_3nt_factory('A','G')
#     spec_dict = dict()
#     spec_dict["norep"] = (norep_align, pre, post)
#     spec_dict["rep"] = (repetative_align, pre, post)

#     graph_dict = {name: [] for name in spec_dict.keys()}

#     for key, fathers in graph_dict.items():
#         # non rep is starting node
#         if key == "norep":
#             continue

#         # rep is son of norep
#         if key == "rep":
#             fathers.append("norep")
#             continue

#     group_dict = {
#         "norep": {"norep"},
#         "rep": {"rep"},
#     }
#     return aligner,spec_dict,graph_dict,group_dict


# def calculate_dictinary_colors_for_site_loss():
#     dict_colors = {}
#     Pileup_stages_to_plot = ['pileup_generation_sorted_pileup', 'no_change', 'read_threshold', 'snp_removal',
#                              'editing_percent', 'concensus_filtered']
#     for index, column in enumerate(Pileup_stages_to_plot):
#         dict_colors[column] = global_list_colors[index]
#     return dict_colors


# def calculate_dictinary_colors_for_editing_percent(dirstruct, positive_fastqs, group_dict):
#     dict_colors = {}
#     # get editing percent pileup and summary file names
#     editing_percent_pileups = [dirstruct.pathName(positive_fastqs[0], node, Stages.editing_type_count,
#                                                   EditTypeStage.edit_percent_pileup) for node in group_dict['norep']]
#     summary_files = [dirstruct.pathName(positive_fastqs[0], node, Stages.editing_type_count, EditTypeStage.file_summary)
#                      for node in group_dict['norep']]
#     # calculate aggregate distribution
#     aggregate_counts, count_summary, pileup_length = editing_site_count_per_type(editing_percent_pileups,
#                                                                                  summary_files)
#     for index, column in enumerate(aggregate_counts.columns):
#         dict_colors[column] = global_list_colors[index]
#     return dict_colors


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
                           graph_dict, spec_dict, group_dict, aligner, parallel_limit=6, disable_parallel=False,
                           skip_existing_files=True, snp_database=snp_database)

    """
    You can comment function that you want to skip. All functions together creats RECIC
    Thresholds can be modified.
    """
    # comment the test_pipe.remove_files if you want to keep previous files in the output directory
    # test_pipe.remove_files()
    # test_pipe.filter_fastq()

    # test_pipe.include_exclude_alignment_test()

    # test_pipe.pileup_creation_test()

    # test_pipe.filter_nochange_test()

    # test_pipe.filter_readthreshold_test(threshold=2)

    # test_pipe.snp_removal_test()

    test_pipe.editing_percent_test(editing_min_threshold=30, editing_max_threshold=99, noise_threshold=3,
                                   editing_read_thresh=2)

    # test_pipe.hyper_non_relevant_editing_site_test()

    # test_pipe.unique_site_test()

    # test_pipe.concensus_test(consensus_threshold=0.5)

    # test_pipe.site_loss_plot_test()

    # test_pipe.editing_percent_plot_generate_summaries(editing_min_thresh=30, editing_max_thresh=99, noise_thresh=3, read_thresh=2)

    # test_pipe.editing_percent_plot_generate_plot()


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    real_pipe()
