from typing import List,Dict
from collections import namedtuple
import subprocess
from Utility.parallel_generator import parallel_generator
from Utility.generators_utilities import class_generator
from Utility.Fastq_class import Fastq
import multiprocessing
import itertools
import os
import re
import logging
import Utility
from contextlib import nullcontext
from Experiments.forontiers_jupyter.directory_structure_definer import DirectoryStructure ,Stages,AlignStage
from Experiments.forontiers_jupyter.aligner_wrapper import AlignerWrapper
from Processing.genome_3nt import complement_multiline, flip_sam_flags, reverse_multiline

Fathers_list = List[str]
Specifications = namedtuple('Specifications', ['id', 'flags', 'pre_function', 'post_function'])
Dict_specifications = Dict[Specifications, Fathers_list]

# default aligner
bowtie_wrapper = AlignerWrapper("bowtie {flags} --sam {reference_file} {lib} {positive_alignment} --un {negative_library} >> {log}.bowtie 2>&1", "bowtie-build {reference} {reference} {flags}",
                                index_output_format_list=["{reference_file}"+suff for suff in [".1.ebwt",".2.ebwt",".3.ebwt",".4.ebwt",".rev.1.ebwt",".rev.2.ebwt"]])


def do_include_exclude_alignment(fastq_file: str, fasta_file: str, dict_data, dict_graph,
                                 dir_struct:DirectoryStructure, aligner:AlignerWrapper, skip_existing_files=False):
    # gaol:Combines dict_data and dict_graph information into one dict
    # from dict_data(spec_dict from the main tool): key and value are given now as the key with fields under
    # "Specifications" 
    # from dict_graph: the value will be the value in the new dict
    spec_dict = dict()
    for id1 in dict_graph.keys():
        spec_dict[Specifications(id1, dict_data[id1][0], dict_data[id1][1], dict_data[id1][2])] = dict_graph[id1]

    include_exclude_alignment(fastq_file, fasta_file, spec_dict, dir_struct=dir_struct,aligner=aligner, skip_existing_files=skip_existing_files)

def include_exclude_alignment(fastq_file: str, fasta_file: str, dict_spec: Dict_specifications,
                                 dir_struct:DirectoryStructure, aligner:AlignerWrapper, skip_existing_files=False):
    """

    :param fastq_file: the path to the fastq file we want to work on
    :param fasta_file: the path to the fasta file we want to work on
    :param dict_spec: a dict which its values are named tuples that each tuple contains an id string, a string of flags
    , a pre processing function and a post processing function . and its values are the id of the stage dependencies (fathers)
    :param action_func: this is the main action function we want our file to go through (after the pre-processing func
    and before the post processing func
    :param dir_struct : a DirectorySturcture objects that gives the correct paths to give to files
    :return: the name of the last fastq file created in the graph (the non - aligned sequances of the last graph step)
    """

    # sending a thread for each vertex in "tree"
    threads = []
    thread_tree = Utility.threading_tree_utils.ThreadingTree(dict_spec)
    #thread_tree: the key is the information for each editing type (ID), the valuse is the thread for each vertex
    for spec, _ in dict_spec.items():
        t = multiprocessing.Process(target=do_pre_align_post, args=(fastq_file, fasta_file, spec.id, thread_tree, aligner,
                                                             dir_struct, skip_existing_files))
        t.start()
        threads.append(t)

    for thread in threads:
        thread.join()

    return thread_tree.get_filenames_of_leafs()


def do_pre_align_post(fastq, fasta_file, id_step, tree_threads, aligner:AlignerWrapper, dirstruct:DirectoryStructure
                      ,skip_existing_files=False):
    """
    :param fastq: the original fastq file name
    :param fasta_file: the fasta file name
    :param id_step: the id of the curr stage to run
    :param tree_threads: the object of ThreadingTree that controls the specification dict and the events
    :param aligner: the aligner wrapper object to use to perform alignment
    :param dirstruct: The Directory sturcture object that generates pathnames for this RESIK run
    :param skip_existing_files: Wether to use existing files instead of recomputing them.
        Used for manual failure recovery

    :return: the name of the negative fastq file (the lines that didn't aligned)
    """

    with tree_threads.threading_context(id_step):  # this is the inner event handler

        # get misaligned fastqs from fathers
        fastqs_list = tree_threads.get_previous_steps_file_names(id_step)
        
        # if all father fastq are None, means that no sequences are available for alignment
        if len(fastqs_list) >=1 and all((prev_fastq == None for prev_fastq in fastqs_list)) :
            logging.info(f"file {fastq} stage {id_step} couldn't be preformed - all the sequences were aligned in the steps before it")
            tree_threads.set_filename_of_finished_alignment(id_step, None)
            return None, None

        # removing those that are none
        fastqs_list=[fas for fas in fastqs_list if fas != None]

        logging.info(f" performing pre-align-post for file {fastq} stage {id_step}  based on {str(fastqs_list)}")
 
        _, spec = tree_threads.get_fathers_and_spec(id_step)

        # if more than 1 father, have to merge them
        if len(fastqs_list) > 1:
            #TODO we need to be able to skip this if recover is true too
            output=dirstruct.pathName(fastq,id_step,Stages.graph_aligner,AlignStage.in_fastq)
            union_fastqs(fastqs_list, output)
            in_fastq = output
        # if only one input, take it
        elif len(fastqs_list) == 1:
            in_fastq = fastqs_list[0]
        # if no input from graph we are a root node, in which case
        # the original fastq is out input
        else:
            in_fastq = fastq

        #set variable to true if on a hyper stage (e.g.: rep_hyper_A_C)
        hyper = not (id_step == "rep" or id_step == "norep")
            
        # define file names
        pre_fasta = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.pre_fasta) #reference file after converting to 3nt
        pre_fastq = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.pre_fastq) #read library after converting to 3nt
        aligned_sam = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.aligned_sam) #aligned reads from alignment of pre_fastq to pre_fasta
        misaligned_fastq = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.misaligned_fastq) #misaligned reads from alignment of antisense_fastq to antisense_fasta
        post_sam = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.post_sam) #aligned sam after reverting to 4nt
        post_fastq = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.post_fastq) #misaligned fastq reads after reverting to 4nt
        if hyper == True:
            antisense_fastq = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.antisense_fastq)#antisense fastq, misaligned reads from alignment of pre_fastq to pre_fasta
            antisense_fasta = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.antisense_fasta)#antisense fasta, with nts complemented(A->T, C->G, etc) then converted to 3nt
            antisense_sam = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.antisense_sam)#antisense sam, alignment of misaligned reads to antisense fasta, with 0 flags flipped to 16
            antisense_post_sam = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.antisense_post_sam)#antisense sam after post function
            sense_filtered_sam = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.sense_filtered_sam) #aligned sam after reverting to 4nt
            antisense_flipped_sam = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.antisense_flipped_sam)#antisense sam after flipping mapped reads to antisense reads (samline flags from 0 to 16)
            antisense_filtered_sam = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.antisense_filtered_sam)#antisense sam after flipping mapped reads to antisense reads (samline flags from 0 to 16)
            sense_misaligned_fastq = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.sense_misaligned_fastq)#misaligned reads from sense alignment
        #decide if any of the stages need skipping based on skip_exisitng_files

        # if no pre function, give the naive pre function (identity pre) and no point in skipping
        skip_pre = False
        if spec.pre_function is None:
            spec = spec._replace(pre_function = Naive_pre)
        else:
            # setting skip pre true only if skip exisiting file and aligned output exist
            if os.path.isfile(pre_fasta) and os.path.isfile(pre_fastq):
                if hyper: 
                    if os.path.isfile(antisense_fasta):
                        skip_pre = skip_existing_files
                else:
                    skip_pre = skip_existing_files


        # if no post function, give the naive pre function (identity pre) and no point in skipping
        # if no post need to have the alinger stage output the post sam
        skip_post=False
        if spec.post_function is None:
            spec = spec._replace(post_function = Naive_post)
            aligned_sam = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.post_sam)
            misaligned_fastq = dirstruct.pathName(fastq, id_step, Stages.graph_aligner, AlignStage.post_fastq)
        else:
            # setting skip post true only if skip exisiting file and post output exist
            if os.path.isfile(post_sam) and os.path.isfile(post_fastq):
                if hyper: 
                    if os.path.isfile(antisense_sam) and os.path.isfile(antisense_post_sam):
                        skip_post = skip_existing_files
                else:
                    skip_post = skip_existing_files
                
        #TODO: skipping innacurate, since align and post stage currently share antisense fastq
        
        #setting skip allign true only if skip exisiting file and aligned output exist
        skip_align = False
        if os.path.isfile(aligned_sam) and os.path.isfile(misaligned_fastq):
            if hyper:
                if os.path.isfile(antisense_sam) and os.path.isfile(antisense_fastq):
                    skip_align = skip_existing_files
            else:
                skip_align = skip_existing_files
                    
        ## computing pre
        if skip_pre :
            logging.info(f"skipping pre calculation of {id_step} on {fasta_file} and {in_fastq} ")
        else:
            logging.info(f"computing pre calculation of {id_step} on {fasta_file} and {in_fastq} ")
            # interface expected is f(in_fasta,in_fastq,out_fasta,out_fastq) where arguemtns are file names
            pre_fastq, pre_fasta = spec.pre_function(fasta_file, in_fastq,pre_fasta,pre_fastq)
            if hyper == True:
                complement_multiline(2,'~',2,fasta_file,output_filename=antisense_fasta)#TODO: complement before or after 3nt?                 
                spec.pre_function(antisense_fasta, None, antisense_fasta, None)

        # computing alignment
        if skip_align :
            logging.info(f"skipping alignment  of {id_step} on {pre_fasta} and {pre_fastq} ")
            if hyper == True:
                logging.info(f"skipping alignment  of {id_step} on {antisense_fasta} and {misaligned_fastq} ")
        else:
            # get a lock that will be acquired when building index, since we cant have siblings build index together.
            locking_context = tree_threads.locking_context()
            secondary_locking_context = tree_threads.secondary_locking_context()
            logging.info(f"computing alignment of {id_step} on {pre_fasta} and {pre_fastq} ")

            if hyper == True:
                aligner.align(lib=pre_fastq, reference=pre_fasta,flags=spec.flags,
                        pos=aligned_sam, neg=sense_misaligned_fastq,locking_context=locking_context,build_index=True, log=aligned_sam)
                logging.info(f"computing antisense alignment of {id_step} on {antisense_fasta} and {misaligned_fastq} ")
                reverse_multiline(4,chr(127),2,sense_misaligned_fastq, antisense_fastq)
                aligner.align(lib=antisense_fastq, reference=antisense_fasta,flags=spec.flags,
                        pos=antisense_sam, neg=misaligned_fastq,locking_context=secondary_locking_context,build_index=True, log=antisense_sam)
                        
            else:
                aligner.align(lib=pre_fastq, reference=pre_fasta,flags=spec.flags,
                        pos=aligned_sam, neg=misaligned_fastq,locking_context=locking_context,build_index=True, log=aligned_sam)


        #post calculation
        if skip_post :
            logging.info(f"skipping post calculation of {id_step} on {aligned_sam} and {misaligned_fastq} ")

        else:
            logging.info(f"computing post calculation of {id_step} on {aligned_sam} and {misaligned_fastq} ")
            # interface expected is f(in_fastq,in_sam,out_fastq,out_sam, **kwargs) where arguemtns are file names
            
            if hyper == True:
            
                command =   f"samtools view -h -F 4 {aligned_sam} | samtools view -h -F 16 > {sense_filtered_sam} &&" +\
                f"samtools view -h -F 4 {antisense_sam} | samtools view -h -F 16 > {antisense_filtered_sam}"
                try:
                    subprocess.check_output(command, shell=True)
                except Exception as e:
                    raise e         
                
                spec.post_function(misaligned_fastq, sense_filtered_sam, post_fastq, post_sam, fastq_lib=fastq)
                flip_sam_flags(antisense_filtered_sam, 0, 16, output_filename=antisense_flipped_sam)
                spec.post_function(None, antisense_flipped_sam,None,antisense_post_sam, fastq_lib=fastq)

            else:
                post_fastq, post_sam = spec.post_function(misaligned_fastq, aligned_sam,post_fastq,post_sam, fastq_lib=fastq)
        
        
        logging.info(f"setting output of step {id_step} as {post_sam} and {post_fastq}")
        # letting tree threads now what the negative fastq of this step is
        tree_threads.set_filename_of_finished_alignment(id_step, post_fastq)

        return

def Naive_pre(in_fasta,in_fastq,pre_fasta,pre_fastq,**kwargs):
    """
    Pre processing function for identity pre processing.
    It ignores the out_paths given and returns the original fasta and fastq given
    :param in_fasta: input fasta
    :param in_fastq: input fastq
    :param pre_fasta: desired name for preprocessed fasta
    :param pre_fastq: desired name for preprocessed fasta
    :return: names of preproccessed fasta and fastq as a 2-tuple
    """
    return (in_fastq,in_fasta)#fastq,fasta expected in rest of code

def Naive_post(in_fastq,in_sam,post_fastq,post_sam,**kwargs):
    """
    Pre processing function for identity pre processing.
    It ignores the out_paths given and returns the original fasta and fastq given
    :param in_fastq: input fastq, the reads that failed to align
    :param in_sam: input sam the reads that did align
    :param post_fastq: desired name for prostprocessed fastq
    :param post_sam: desired name for preprocessed sam
    :return: names of post_processed fastq and sam as a 2-tuple
    """
    return (in_fastq,in_sam)

def get_id(fastq: Fastq):
    return fastq.get_id()


def union_fastqs(list_fastqs, output):
    """
    take a list of fastqs, and output the returns a name of file that contains the union (set)
    :param list_fastqs: the name of fastq we want to union into one
    :return: the name of the new intersected fastq file
    """
    # preparing input :

    list_opened_files = [open(file1, "r") for file1 in list_fastqs if file1 is not None]

    list_gens = [class_generator(Fastq, file=file2, number_of_lines=4) for file2 in list_opened_files]
    list_func = [item for item in itertools.repeat(get_id, len(list_gens))]

    with open(output, "w") as output:
        # from each equivalent class (same fastq line) we need only one item.
        for list_of_equivalence_classes in parallel_generator(list_gens, list_func):
            for item in list_of_equivalence_classes:
                if item is not None:
                    output.write(str(item[0]))
                    break

    for file in list_opened_files:
        file.close()
