"""site_loss_by_group_plot
The script receives an output file path for the plot and a list of directories for each group a directory that will
contain a directories of all the stages necessary for the stacked plot. the folders of the original pileups before all
the steps name should be received in an original_folder_name parameter. all the other steps folders should have the same
name for all of the groups folders.
Usage:
    site_loss_by_group_plot.py <output> <original_folder_name> <dir_path>...

"""


'''
example structure:

the folders as dir_path list : group1 , group2 , group3 , group4 

in group1 (the same in all other group2-4 folders) :
original , step1 , step2 , step3 ..
(each of these confines the group results in pileup files after each filtering step)


'''
from  Experiments.forontiers_jupyter.bar_utils import stacked_bar
from os import listdir
from os.path import isfile, join, dirname
from Utility.generators_utilities import getLineFromChunk
import pandas as pd
import itertools
import logging
from Utility.count_utils import  column_delta_df,get_line_count
from Experiments.forontiers_jupyter.directory_structure_definer import DirectoryStructure,Stages \
    ,AlignStage,PileupStage,ConcensusStage,SiteLossStage


list_order_of_folders = ["consensus_filterout","editing_sites_filterout","snp_removed","reads_threshold_filterout","no_change_filterout"]


def make_data_frames(dict_group_with_steps, dict_files_with_steps, output,list_order_of_folders):
    """
    creates data frames of the values that created the plots
    :param dict_group_with_steps: a dict of group as key and a dict of (step - key , lines num - value) as value
    :param dict_files_with_steps: a dict of filename as key and a dict of (step - key , lines num - value) as value
    :param output: the output folder to print the table to
    """

    data1 = dict()
    data2 = dict()

    filename1 = join(output, "files_size_by_step_table.csv")
    filename2 = join(output, "groups_size_by_step_table.csv")

    data1["steps"] = list_order_of_folders
    data2["steps"] = list_order_of_folders

    for filename, dict_steps in dict_files_with_steps.items():
        data1[filename] = [dict_steps[step] for step in list_order_of_folders]

    for group, dict_steps in dict_group_with_steps.items():
        data2[group] = [dict_steps[step] for step in list_order_of_folders]
    open(filename1,"w").close()
    open(filename2,"w").close()
    pd.DataFrame(data1).to_csv(filename1, sep=',')
    pd.DataFrame(data2).to_csv(filename2, sep=',')

def line_count_dict_to_csv(lc_dict,columns_sort_by=None)->pd.DataFrame:
    """

    :param lc_dict: dict of form {(key,stage,aux):line_count}
    :param columns_sort_by: list of [(stage,aux),..] type that defines how we should order the columns
    :return: pd Dataframe df such that df["key","stage_aux"] = line_count
    """

    def stage_aux_name(stage,aux):
        if aux == None:
            return f"{stage.name}"
        else:
            return f"{stage.name}_{aux.name}"


    # make lists of unique row names and col names
    rows=list({key for key,_,_ in lc_dict.keys()})
    columns=list({(stage,aux) for _,stage,aux in lc_dict.keys()})

    # if given column order to sort by, sort it
    if columns_sort_by != None:
        def sort_by_index(i):
            # sorts by index of columns sort by
            return columns_sort_by.index(i)

        columns=sorted(columns,key=sort_by_index)

    #change columns into strings
    columns=[stage_aux_name(stage,aux) for stage,aux in columns]

    # make dataframe with col and row names
    df=pd.DataFrame(index=rows,columns=columns)

    # put data in df
    for (key,stage,aux),value in lc_dict.items():
        df.loc[key,stage_aux_name(stage,aux)]=value



    return df

def site_loss_by_group_plot(lib_name,node_names,group_dict,dirstruct:DirectoryStructure,dict_colors):
    """
    generates plot of site loss across the different stages of analysis of a library
    :param lib_name: the original library name
    :param node_names: list of nodes in graph aligner
    :param group_dict: Dict of (group_name,set_of_nodes) That shows how to aggregate the nodes
    :param dirstruct: Directory structure of the project
    :param dict_colors: Dictionary of color bar for each combinations of nucleotides
    :return: puts size summary files per node and per group as well as plot
    """

    # list of Stage,auxStage pairs to plot
    Pileup_stages_to_plot=[
        (Stages.pileup_generation,PileupStage.sorted_pileup),
        (Stages.no_change, None),
        (Stages.read_threshold, None),
        (Stages.snp_removal, None),
        (Stages.editing_percent, None),
        (Stages.concensus, ConcensusStage.filtered)
    ]

    # mapping between node,stage,aux_stage and the file that contains it
    file_names_dict={}
    for node,(stage,aux) in itertools.product(node_names,Pileup_stages_to_plot):
        file_names_dict[node,stage,aux]=dirstruct.pathName(lib_name,node,stage,aux)

    # mapping betweeen node,stage,aux and line count
    line_count_dict={}
    for key,file_name in file_names_dict.items():
        line_count_dict[key]=get_line_count(file_name)

    # mapping between group_name,stage,aux to total line count of the group
    group_line_count_dict={}
    for group_name,group_nodes_list in group_dict.items():
        for stage,aux in Pileup_stages_to_plot:
            line_count_sum=sum(line_count_dict[node,stage,aux] for node in group_nodes_list)
            group_line_count_dict[group_name,stage,aux]=line_count_sum


    line_count_df=line_count_dict_to_csv(line_count_dict,columns_sort_by=Pileup_stages_to_plot)
    group_line_count_df=line_count_dict_to_csv(group_line_count_dict,columns_sort_by=Pileup_stages_to_plot)

    lc_filename=dirstruct.pathName(lib_name,None,Stages.site_loss,SiteLossStage.file_summary)
    lc_group_filename=dirstruct.pathName(lib_name,None,Stages.site_loss,SiteLossStage.group_summary)
    loss_plot_filename=dirstruct.pathName(lib_name,None,Stages.site_loss,SiteLossStage.plot)

    line_count_df.to_csv(lc_filename)
    group_line_count_df.to_csv(lc_group_filename)

    try:
        row_delta_histogram(loss_plot_filename,group_line_count_df,dict_colors)
    except:
        logging.exception("site loss plot failed")


def row_delta_histogram(file_name,df:pd.DataFrame,dict_colors):
    # reorders colums from smallest to biggest (asusming the oposite ordering)
    df=df[df.columns[::-1]]

    data=column_delta_df(df)
    # move all column names in according to the delta calculation (0 column meaningful after the delta and will be removed)
    data.columns=[data.columns[0],data.columns[0],data.columns[1],data.columns[2],data.columns[3],data.columns[4]]

    #Takes only the first delta (pileup gen minus no change)
    initial_counts=data.iloc[:,-1]
    #remove the 6st column that is meaningful after the delta
    other_counts=data.iloc[:,1:]

    #maximum sum or rows (TODO:take the max of the delta. we need to total reads?)
    max_y= other_counts.sum(axis=1).max()
    # leave some white space at top of plot
    max_y=int(1.1*max_y)

    show_values=True
    y_label="Number_of_sites"

    #transpose data
    other_counts=other_counts.transpose()

    plt_res, axes = stacked_bar(other_counts,
                                y_label=y_label,size_plot=[20,5],show_values=show_values,
                                use_dataframe=True,throw_zeros=True,dict_colors=dict_colors) #
    item_last_data = axes[-1]

    # adding the number of total lines of the original pileups
    for index ,bar in enumerate(item_last_data):
        w, h = bar.get_width(), bar.get_height()
        plt_res.text(bar.get_x() + w / 2, bar.get_y() + h*1.1,'%d' % int (initial_counts[index]), ha="center",va="center")

    #puts the ledgends outside of the plot
    plt_res.subplots_adjust(right=0.62)
    plt_res.legend(loc='center left',bbox_to_anchor=(1, 0.5),handles=axes[::-1])

    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    plt_res.savefig(file_name)

    plt_res.show()



if __name__=="__main__":

    from docopt import docopt

    args = docopt(__doc__)

    # Extracting the parameters
    list1 = args['<dir_path>']
    output = args['<output>']
    original_folder_name = args['<original_folder_name>']

    site_loss_by_group_plot(list1, output, original_folder_name)

