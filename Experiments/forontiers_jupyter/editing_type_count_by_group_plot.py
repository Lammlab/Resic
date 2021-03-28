"""editing_type_count_by_group
    Generating a editing site distribution plot for RESIC
"""
from Experiments.forontiers_jupyter.bar_utils import stacked_bar
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import logging
from typing import Dict,List
from os import listdir
from os.path import isfile, join, dirname
from Utility.generators_utilities import group_generator
from Utility.count_utils import get_line_count
import pandas as pd
from Experiments.forontiers_jupyter.directory_structure_definer import DirectoryStructure,Stages,ConcensusStage,EditTypeStage


def editing_type_count_by_group_plot(lib_name,group_dict:Dict,dirstruct:DirectoryStructure,dict_colors):
    """
    Generates a plot showing the editing percent type distribution of a certain group of nodes
    :param lib_name: name of original library file
    :param group_dict: dict of form {group_name:list_of_nodes_in_group}
    :param dirstruct: DirectoryStructure object for this RESIK run
    :param dict_colors: Dictionary of the nucleotides combination to his color bin
    :return: Populates The directory structure with:
                * A summary file showing the editing type distribution per node for each  group
                * A summary file showing the aggregate editing type distribution for all nodes for each group
                * A plot showing the aggregate distribution of editing sites across all groups
    """

    group_counts=dict()

    # get aggregate counts per group
    for group_name,group_nodes in group_dict.items():
        # get editing percent pileup and summary file names
        editing_percent_pileups=[dirstruct.pathName(lib_name,node,Stages.editing_type_count,EditTypeStage.edit_percent_pileup)
                             for node in group_nodes ]
        summary_files=[dirstruct.pathName(lib_name,node,Stages.editing_type_count,EditTypeStage.file_summary)
                             for node in group_nodes ]

        # calculatte aggregate distribution
        aggregate_counts,count_summary,pileup_length=editing_site_count_per_type(editing_percent_pileups,summary_files)
        # save it for plot
        group_counts[group_name]=aggregate_counts

        #output aggregate counts to file
        aggregate_summary_file=dirstruct.pathName(lib_name,group_name,Stages.editing_type_count,EditTypeStage.group_distribution_summary)
        count_summary.to_csv(aggregate_summary_file)
        #output counts per file to file
        group_summary_file=dirstruct.pathName(lib_name,group_name,Stages.editing_type_count,EditTypeStage.group_count_summary)
        count_summary.to_csv(group_summary_file)

    # generating the plot
    try:
        plt.figure()
        group_names=[name for name in group_dict.keys()]
        data=pd.concat(aggregate_counts for aggregate_counts in group_counts.values())

        data.index=group_names
        data=data.transpose()

        plt_res, axes = stacked_bar(data,  show_values=True, value_format="{:.3f}",
                    y_label="Percent of sites",size_plot=[18,20],use_dataframe=True,throw_zeros=True,dict_colors=dict_colors)

        #puts the ledgends outside of the plot
        plt_res.subplots_adjust(right=0.62)
        plt_res.legend(loc='center left',bbox_to_anchor=(1, 0.5),handles=axes[::-1])

        output_path = dirstruct.pathName(lib_name,None,Stages.editing_type_count,EditTypeStage.plot)
        plt_res.savefig(output_path)
        plt_res.show()
    except:
        logging.exception("edit plot failed")



def editing_site_count_per_type(pileup_files, summary_files)->pd.DataFrame:
    """
    Generates a dataframe whos rows are the pileup files and columns are the editing type and cell i,j
        contains the number of sites with editing type j in pileup i
    :param pileup_files: list of path_names for pileup files
    :param summary_files: list of summary file for each pileup in pileup_files, according to analyze editing percent
    :return: (total_dist,summary,length_of_pileups) where
        total_dist: sum of counts per editing change type across all pileups
        summary: dataframe as described above
        length_of_pileups: vector of number of lines in each pileup file
    """

    # get vector of file length
    length_of_pileups=[get_line_count(pileup) for pileup in pileup_files]

    # get editing type percentages per file
    summary_dfs=[load_summary_file_to_df(summary) for summary in summary_files]

    #make length of pileups a vector
    length_of_pileups=np.array(length_of_pileups).reshape([-1,1])
    #concatenate summaries into one big summary
    summary=pd.concat(summary_dfs)

    # multiply each row in summary by its file length to get absolute site nubmers
    absolute_lengths=length_of_pileups*summary.values

    summary.iloc[:,:]=np.array(absolute_lengths,dtype=np.int32)

    #sum up groups
    total_dist=pd.DataFrame(summary.sum(axis=0)).transpose()
    total_length=absolute_lengths.sum()

    return total_dist,summary,length_of_pileups

def load_summary_file_to_df(summary_file):
    #Currently assume only one file summarised in summary file
    #TODO fix this

    #get first two lines, throw first one with originalk pileup name
    with open(summary_file,"r") as fp:
        lines=fp.readlines()
        if len(lines)==0:
            logging.warning(f"{summary_file} has zero lines")
            return None

        _,edit_percents=lines[:2]



    #split by tabs
    edit_percents.strip('\n')
    edit_percents=re.split('\t|:',edit_percents)

    # splits to key,val pairs
    keys,vals=zip(*group_generator(edit_percents,n=2))
    vals=np.array([float(v) for v in vals])
    #put into a df
    file_name=os.path.basename(summary_file)
    df=pd.DataFrame(vals,index=keys,columns=[file_name])
    df=df.transpose()
    return df

if __name__=="__main__":
    pass
