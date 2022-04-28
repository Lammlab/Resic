"""
This basically deifnes how the output directory of resic will look like.

"""

import os 
import logging
import traceback
from enum import Enum
from Experiments.forontiers_jupyter.pipe_utils import all_genome_pair_combinations

class Stages(Enum):
    graph_aligner=1
    pileup_generation=2
    no_change=3
    read_threshold=4
    snp_removal=5
    snp_removal_from_db=6
    editing_percent=7
    editing_percent_unique=8
    concensus=9
    site_loss=10
    editing_type_count=11

class AlignStage(Enum):
    in_fastq=0
    in_fasta=1
    pre_fastq=2
    pre_fasta=3
    aligned_sam=4
    misaligned_fastq=5
    post_sam=6
    post_fastq=7
    antisense_fastq=8
    antisense_fasta=9
    antisense_sam=10
    antisense_filtered_sam=11
    antisense_flipped_sam=12
    antisense_post_sam=13
    sense_filtered_sam=14
    sense_misaligned_fastq=15



class PileupStage(Enum):
    bam=0
    sorted_bam=1
    pileup=2
    sorted_pileup=3
    combined_sam=4 

class ConcensusStage(Enum):
    filtered=0
    concensus=1


class SiteLossStage(Enum):
    file_summary=0
    group_summary=1
    plot=2

class EditTypeStage(Enum):
    edit_percent_pileup=0
    file_summary=1
    group_count_summary=2
    group_distribution_summary=3
    plot=4

class DirectoryStructure():
    def __init__(self,root_dir,graph_nodes=None):
        self.root_dir=root_dir
        self.graph_nodes=graph_nodes
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
        return

    def pathName(self,lib_name,node_name,stage:Stages,aux=None,need_suffix=True):

        dirname=self.dirName(lib_name,node_name,stage,aux)
        filename=self.fileName(lib_name,node_name,stage,aux,need_suffix)

        path=os.path.join(dirname,filename)

        #if directories dont exist, create them
        try:
            os.makedirs(dirname)
        except os.error :
            pass

        return path

    def fileName(self,lib_name,node_name,stage:Stages,aux=None,need_suffix=True):

        if lib_name != None: #lib specific files
            file_dir,file_name_with_ext=os.path.split(lib_name)

            file_name = os.path.splitext(file_name_with_ext)[0]
        else:
            file_name="global"


        if aux is None:
            name=f"{file_name}_{node_name}_{stage.name}"
        else:
            name=f"{file_name}_{node_name}_{stage.name}_{aux.name}"


        # TODO compute file type
        if need_suffix:
            if stage == Stages.graph_aligner:
                if 'fasta' in aux.name:
                    name += '.fasta'
                elif 'fastq' in aux.name:
                    name += '.fastq'
                elif 'sam' in aux.name:
                    name += '.sam'
            # figures
            elif stage == Stages.site_loss or stage == Stages.editing_type_count:
                if "plot" in aux.name:
                    name += ".tiff"
                else:
                    name += '.csv'
            elif stage == Stages.pileup_generation and 'bam' in aux.name:
                name += '.bam'
            elif stage == Stages.pileup_generation and 'sam' in aux.name:
                name += '.sam'
            else:
                name += '.pileup'

        return name
    
    def dirName(self,lib_name,node_name,stage:Stages,aux=None):
        # figures
        if stage in [Stages.site_loss,Stages.editing_type_count]:
            dir = os.path.join(self.root_dir,'figures',stage.name)

        # filtering products
        elif stage in [Stages.pileup_generation, Stages.no_change, 
            Stages.read_threshold, Stages.snp_removal, Stages.snp_removal_from_db,Stages.editing_percent
            , Stages.concensus] or \
             ( stage == Stages.graph_aligner and (aux == AlignStage.post_sam or aux == AlignStage.antisense_post_sam)) :
            dir = os.path.join(self.root_dir, 'raw_data', node_name, stage.name)

        # alignment by products
        elif stage == Stages.graph_aligner:
            dir = os.path.join(self.root_dir, 'raw_data', 'intermediate_alignment_files')
        
        elif  stage in [Stages.editing_percent_unique]:
            dir = os.path.join(self.root_dir,'raw_data', node_name, "editing_percent")

        else:
            raise(AttributeError(" The following arguments are not valid names \n " +
                " ".join([str(i) for i in [lib_name,node_name,stage,aux]]) +
                    '\n'
            ))

        return dir

    def print_struct(self):
        return None


        x="-root directory \n" \
               "\t -raw_data \n" \
               "\t\t -norep \n" \
               "\t\t\t -graph_aligner \n" \
               "\t\t\t\t {fastqInput1_norep_graph_aligner_in_fastq.fastq} {fastqInput2_norep_graph_aligner_aligned_sam.sam}\n" \
               "\t\t\t -pileup_generation\n" \
               "\t\t\t\t {fastqInput1_norep_pileup_generation.bam} {fastqInput1_norep_pileup_generation.pileup}" \
               "\t\t\t -no_change\n" \
               "\t\t\t -read_threshold\n" \
               "\t\t\t -snp_removal\n" \
               "\t\t\t -editing_percent\n" \
               "\t\t\t -concensus\n" \
               "\t\t\t -site_loss\n" \
               "\t\t\t -editing_type_count\n" \
               "\t\t -norep_hyper_AC \n" \
               "\t\t\t -graph_aligner \n" \
               "\t\t\t\t {fastqInput1_norep_graph_aligner_in_fastq.fastq} {fastqInput2_norep_graph_aligner_aligned_sam.sam}\n" \
               "\t\t\t -pileup_generation\n" \
               "\t\t\t\t {fastqInput1_norep_hyper_AC_pileup_generation.bam} {fastqInput1_norep_hyper_AC_pileup_generation.pileup}" \
               "\t\t\t -no_change\n" \
               "\t\t\t -read_threshold\n" \
               "\t\t\t -snp_removal\n" \
               "\t\t\t -editing_percent\n" \
               "\t\t\t -concensus\n" \
               "\t\t\t -site_loss\n" \
               "\t\t\t -editing_type_count\n" \
               "\t\t -rep_hyper_AG \n" \
               "\t\t\t -graph_aligner \n" \
               "\t\t\t\t {fastqInput1_rep_hyper_AG_graph_aligner_in_fastq.fastq} {fastqInput2_rep_hyper_AG_graph_aligner_aligned_sam.sam}\n" \
               "\t\t\t -pileup_generation\n" \
               "\t\t\t\t {fastqInput1_rep_hyper_AG_pileup_generation.bam} {fastqInput1_rep_hyper_AG_pileup_generation.pileup}" \
               "\t\t\t -no_change\n" \
               "\t\t\t -read_threshold\n" \
               "\t\t\t -snp_removal\n" \
               "\t\t\t -editing_percent\n" \
               "\t\t\t -concensus\n" \
               "\t\t\t -site_loss\n" \
               "\t\t\t -editing_type_count\n" \
               "\t\t -rep_transcriptome_GT \n" \
               "\t\t\t -graph_aligner \n" \
               "\t\t\t\t {fastqInput1_rep_transcriptome_GT_graph_aligner_in_fastq.fastq} {fastqInput2_rep_transcriptome_GT_graph_aligner_aligned_sam.sam}\n" \
               "\t\t\t -pileup_generation\n" \
               "\t\t\t\t {fastqInput1_rep_transcriptome_GT_pileup_generation.bam} {fastqInput1_rep_transcriptome_GT_pileup_generation.pileup}" \
               "\t\t\t -no_change\n" \
               "\t\t\t -read_threshold\n" \
               "\t\t\t -snp_removal\n" \
               "\t\t\t -editing_percent\n" \
               "\t\t\t -concensus\n" \
               "\t\t\t -site_loss\n" \
               "\t\t\t -editing_type_count\n" \
               "\t -figures"



