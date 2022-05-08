#######################################################################################################################
# Author:  Roni Haas
# Main Goal: Due to the Resic algorithm, it is possible for a non-hyper editing site to be identified as a hyper
# editing site.
# Namely, for every two distinct nucleotides X and Y under the 3nt genome scheme, nucleotide mismatches other than hyper
# X to Y or Y to X may be recorded, due to conditions created by the 3nt scheme. 
# Therefore, this script is needed to filter the hyper editing files to include only X -> Y or Y -> X changes.
######################################################################################################################

import os
import shutil


def get_candidate_nucl(pileup_line):
    clean_str = pileup_line.split()[4]
    clean_str = [char for char in clean_str if not char in ('^', '$')]
    clean_str = ''.join(clean_str)
    sense_string = clean_str.replace('^', '')
    sense_string = clean_str.replace('$', '')
    sense_string = clean_str.replace(',', '.')

    # Find the candidate nucleotide reads
    nucl_changes = {"A": 0, "C": 0, "G": 0, "T": 0, "a": 0, "c": 0, "g": 0, "t": 0}
    for nucl in list(sense_string):
        if nucl in nucl_changes.keys():
            nucl_changes[nucl] = nucl_changes[nucl] + 1
        else:
            continue
        # get the maximal nucleous change, key is the value of the dict
    (candidate_nucl, candidate_nucl_reads) = max(nucl_changes.items(), key=lambda x: x[1])
    return (candidate_nucl)


def filter_for_specific_node_XtoY_editing_sites(pileup_file, node_name):
    from_nuc = node_name.split("_")[2]
    to_nuc = node_name.split("_")[3]
    nuc_pair_list = {"A": "T", "C": "G", "T": "A", "G": "C"}

    original_pileup_file = os.path.join(pileup_file + "_before_" + from_nuc + "_" + to_nuc + "_filtering")
    shutil.copyfile(pileup_file, original_pileup_file)

    with open(original_pileup_file, 'r') as pileup:
        with open(pileup_file, 'w') as new_pileup:
            for line in pileup:
                reference_nucl = line.split()[2]
                reference_nucl = reference_nucl.upper()
                if reference_nucl == from_nuc and get_candidate_nucl(line) == to_nuc or \
                        reference_nucl == to_nuc and get_candidate_nucl(line) == from_nuc or \
                        reference_nucl == nuc_pair_list[to_nuc] and get_candidate_nucl(line) == nuc_pair_list[
                    from_nuc].lower() or \
                        reference_nucl == nuc_pair_list[from_nuc] and get_candidate_nucl(line) == nuc_pair_list[
                    to_nuc].lower():
                    new_pileup.write(str(line))

    return pileup_file


def filter_hyper_non_relevant_editing_sites(pileup_filename_list):
    for file in pileup_filename_list:
        is_hyper = False
        node_name = file.split("/")[-3]
        try:
            is_hyper = (node_name.split("_")[1] == "hyper")
        except:
            is_hyper == False
        if is_hyper:
            filter_for_specific_node_XtoY_editing_sites(file, node_name)
