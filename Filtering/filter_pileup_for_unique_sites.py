

########################################################################################################################
# Author:  Roni Haas
# Main goal: Takes several pileup files and for each one prints to a new file only unique sites
# (that is, only sites that are not in any of the other pileup files)
########################################################################################################################

import shutil
import Utility.generators_utilities as gen_util
from Utility.Pileup_class import Pileup_line
from Utility.parallel_generator import parallel_generator
from Utility.multiline_sort import multiline_sort_pileup

def get_candidate_nucl(pileup_line):
    clean_str=pileup_line.reads_string
    clean_str = [char for char in clean_str if not char in ('^', '$')]
    clean_str = ''.join(clean_str)
    sense_string = clean_str.replace('^','')
    sense_string = clean_str.replace('$','')
    sense_string = clean_str.replace(',','.')
    # TODO before use assumed letters are only upper case. I fix this assumption here
    sense_string_upper = sense_string.upper()
    # Find the candidate nucleotide reads
    nucl_A, nucl_C, nucl_T, nucl_G = 0, 0, 0, 0
    nucl_changes = {"A": 0, "C": 0, "G": 0, "T": 0, "a": 0, "c": 0, "g": 0, "t": 0}
    for nucl in list(sense_string):
        if (nucl=='A' or nucl=='C' or nucl=='G' or nucl=='T' or nucl=='a' or nucl=='c' or nucl=='g' or nucl=='t'):
            nucl_changes[nucl] = nucl_changes[nucl]+1
        else:
            continue
        # get the maximal nucleous change, key is the value of the dict
    (candidate_nucl, candidate_nucl_reads) = max(nucl_changes.items(), key=lambda x: x[1])
    return(candidate_nucl)

def write_unique_sites_doing_nothing(positive_pileup_list, sorted_input=True):

    """
    :param pileup_filename: positive pileup to be filtered
    :param sorted_input: binary flag, set to False if input pileups are not sorted
    :return
    """

    # sorting pileups if necessary and preparing filenames
    sorted_positive_pileups = []
    if not sorted_input:
        for pileup in positive_pileup_list:
            multiline_sort_pileup(1, '~~~', 1, pileup, 2, pileup + "_sorted.pileup")
            sorted_positive_pileups.append(pileup + "_sorted.pileup")
    else:
        sorted_positive_pileups = positive_pileup_list

    #create output files for each pileup
    output_positive_pileups = []
    for pile in sorted_positive_pileups:
        splited = pile.split(".")
        shutil.copyfile(pile, ".".join(splited[:-1]) + "_unique." + splited[-1])
        output_positive_pileups.append(".".join(splited[:-1]) + "_unique." + splited[-1])

    return output_positive_pileups

def write_unique_sites(positive_pileup_list, sorted_input=True):
    """
    :param pileup_filename: positive pileup to be filtered
    :param sorted_input: binary flag, set to False if input pileups are not sorted
    :return
    """

    # sorting pileups if necessary and preparing filenames
    sorted_positive_pileups = []
    if not sorted_input:
        for pileup in positive_pileup_list:
            multiline_sort_pileup(1, '~~~', 1, pileup, 2, pileup + "_sorted.pileup")
            sorted_positive_pileups.append(pileup + "_sorted.pileup")
    else:
        sorted_positive_pileups = positive_pileup_list

    #create output files for each pileup
    output_positive_pileups = []
    for pile in sorted_positive_pileups:
        splited = pile.split(".")
        output_positive_pileups.append(".".join(splited[:-1]) + "_unique." + splited[-1])

    # parameters for parallel generator
    pos_obj_list =  [open(pile) for pile in sorted_positive_pileups]
    out_obj_list = [open(out,'w') for out in output_positive_pileups]

    pos_gen_list = []
    get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
    # logic of the function
    try:
        for file_obj in pos_obj_list:
            pos_gen_list.append(gen_util.class_generator(Pileup_line, file=file_obj))

        parallel_gen = parallel_generator([*pos_gen_list],
                                              [get_pos_and_id for i in range(len([*pos_gen_list]))])

        for listlist in parallel_gen:  # listlist: list of lists, each entry in listlist is a list of items from one generator
            for pileup in listlist:
                if pileup is not None:
                    if pileup[0].reference_id == "NC_000001.10" and pileup[0].gene_pos == 564559:
                        print(pos_obj_list, listlist)
            true_index_list = [index for index, value in enumerate(listlist) if value]
            if len(true_index_list) == 1:
                #print(true_index_list)
                out_obj_list[true_index_list[0]].write(str(listlist[true_index_list[0]][0]) + "\n")
            # for conditions that in hyper there are 1 site that were mapped in 2 different times (A-G (A_G file)and 
            # A-g(C_T file)). True for non stranded sittuations. In these cases we leave only once the site randomly
            if len(true_index_list) == 2:
                if get_candidate_nucl(listlist[true_index_list[0]][0]).upper() == get_candidate_nucl(listlist[true_index_list[1]][0]).upper():
                    out_obj_list[true_index_list[1]].write(str(listlist[true_index_list[1]][0]) + "\n")
            if len(true_index_list) == 3:
                if get_candidate_nucl(listlist[true_index_list[0]][0]).upper() == get_candidate_nucl(listlist[true_index_list[1]][0]).upper() and \
                get_candidate_nucl(listlist[true_index_list[2]][0]).upper() == get_candidate_nucl(listlist[true_index_list[0]][0]).upper():
                    out_obj_list[true_index_list[1]].write(str(listlist[true_index_list[1]][0]) + "\n")

        
    except Exception as e:
        print(e)

    for out_file in out_obj_list:
            out_file.close()

    for file in pos_obj_list:
        file.close()

    return output_positive_pileups

if __name__ == '__main__':

	pass

