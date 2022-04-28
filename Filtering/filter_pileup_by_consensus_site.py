
"""
Filter a list of pileups to contain only sites that where seen in over k% of the pileups

Usage:
	filter_pileup_by_consensus_site --percentage=K [-u|--unsorted] [--consensus=CONSENSUS --dir=DIR PILEUPS...
	filter_pileup_by_consensus_site params
	filter_pileup_by_consensus_site example
	filter_pileup_by_consensus_site (-h|--help)

Options:
	--percentage=K	the percentage of pileup files a line must be included in to be within the consensus. 0<k<1
	--unsorted	adding this flag will sort the input files before filtering them
	--consensus=CONSENSUS	option name of file to write consensus file to
	--dir	directory to write the filtered pileups to
	PILEUPS	the filenames of the pileup files to be filtered

"""


##############################################################################################################################
# Author:  Roni Haas
# Main Goal: This script filters a list of pileup files so that they only contain sites covered in more than a given
# percent of the pileup files.
##############################################################################################################################


from Utility.generators_utilities import class_generator
from Utility.parallel_generator import parallel_generator
from Utility.Pileup_class import Pileup_line
from Processing.pileup_sorting import pileup_sort
from docopt import docopt
import tempfile
import os
import shutil


def get_candidate_nucl(pileup_line):
	clean_str=pileup_line.split()[4]
	clean_str = [char for char in clean_str if not char in ('^', '$')]
	clean_str = ''.join(clean_str)
	sense_string = clean_str.replace('^','')
	sense_string = clean_str.replace('$','')
	sense_string = clean_str.replace(',','.')
	# TODO before use assumed letters are only upper case. I fix this assumption here
	#sense_string_upper = sense_string.upper()
	# Find the candidate nucleotide reads
	nucl_changes = {"A": 0, "C": 0, "G": 0, "T": 0, "a": 0, "c": 0, "g": 0, "t": 0}
	for nucl in list(sense_string):
		if nucl in nucl_changes.keys():
			nucl_changes[nucl] = nucl_changes[nucl]+1
		else:
			continue
		# get the maximal nucleous change, key is the value of the dict
	(candidate_nucl, candidate_nucl_reads) = max(nucl_changes.items(), key=lambda x: x[1])
	return(candidate_nucl)


def filter_for_specific_node_XtoY_editing_sites(pileup_file, node_name):
	from_nuc = node_name.split("_")[2]
	to_nuc = node_name.split("_")[3]
	nuc_pair_list = {"A":"T", "C":"G", "T":"A", "G":"C"}

	original_pileup_file = os.path.join(pileup_file + "_before_" + from_nuc + "_" + to_nuc + "_filtering")
	shutil.copyfile(pileup_file, original_pileup_file)

	with open (original_pileup_file,'r') as pileup:
		with open (pileup_file,'w') as new_pileup:
			for line in pileup:
				reference_nucl = line.split()[2]
				reference_nucl = reference_nucl.upper()
				if reference_nucl == from_nuc and get_candidate_nucl(line) == to_nuc or \
				 reference_nucl == to_nuc and get_candidate_nucl(line) == from_nuc or \
				 reference_nucl == nuc_pair_list[to_nuc] and get_candidate_nucl(line) == nuc_pair_list[from_nuc].lower() or \
				 reference_nucl == nuc_pair_list[from_nuc] and get_candidate_nucl(line) == nuc_pair_list[to_nuc].lower():
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

def filter_by_consensus(pileup_filename_list, k,filtered_output_list,concensus_file,sorted_input=True):
	#filter_hyper_non_relevant_editing_sites(pileup_filename_list)
	create_consensus_file(pileup_filename_list, k, concensus_file, sorted_input)
	filter_pileups_by_consensus(pileup_filename_list, concensus_file, sorted_input, filtered_output_list)



def create_consensus_file(pileup_filename_list, k, consensus_filename, sorted_input=True):
	'''
	:param pileup_filename_list: list of pileup files to create a consensus list of pileup lines out of
	:param k: the percentage of pileup files a line must be included in to be within the consensus. 0<k<1
	:param consensus_filename: name of output file
	:param sorted_input: binary flag, set to False if input pileups are not sorted
	:return: consensus file
	'''
	sorted_pileups = []
	if not sorted_input:
		for pileup in pileup_filename_list:
			pileup_sort(pileup, pileup+"_sorted.pileup")
			sorted_pileups.append(pileup+"_sorted.pileup")
	else:
		sorted_pileups = pileup_filename_list

	pileup_files = [open(pileup, 'r') for pileup in sorted_pileups]

	try:
		generators = [class_generator(Pileup_line, file=file) for file in pileup_files]
		get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
		parallel_gen = parallel_generator(generators, [get_pos_and_id for i in range(len(generators))])
		number_of_pileups = len(pileup_files)
		with open(consensus_filename, 'w') as out:
			for pileuplist_list in parallel_gen:
				hits = sum(x is not None for x in pileuplist_list) #number of non-None lists, the number of pileups which have the entry
				if (hits/number_of_pileups) >= k:
					pileup_list = next(item for item in pileuplist_list if item is not None) #first non-None pileup, all
															# lines in the list from parallel_gen should have the same line
					out.write(str(pileup_list[0]) + "\n")
	finally:
		for file in pileup_files:
			file.close()
		return consensus_filename


def filter_pileups_by_consensus(pileup_list,consensus_file, sorted_input=True, output_list=None):
	'''

	:param pileup_list: list of pileup files to be filtered
	:param consensus_file: consensus pileup, pileups from pileup list will be filtered so that only lines that are in
	the consensus pileup remain in each pileup file
	:param sorted_input: binary flag, set to False if input pileups are not sorted
	:param output_dir: directory to write output to
	:return: list of filtered pileup files
	'''
	# if not sorted, sort
	if not sorted_input:
		sorted_pileups = []
		for pileup in pileup_list:
			file_path,extension=os.path.splitext(pileup)
			out_name=f"{file_path}_sorted.{extension}"
			pileup_sort(pileup, out_name)
			sorted_pileups.append(out_name)
	else:
		sorted_pileups=pileup_list

	#if no output list, make conensus files in the same place as the input file
	if output_list is None:
		output_list=[]
		for pileup in sorted_pileups:
			file_path,extension=os.path.splitext(pileup)
			out_name=f"{file_path}_kfiltered.{extension}"
			output_list.append(out_name)


	for (pileup,pileup_out) in zip(sorted_pileups, output_list):
		with open(pileup, 'r') as pileup_object, open(pileup_out, 'w') as out, open(consensus_file,'r') as consensus:
			pileup_gen = class_generator(Pileup_line, file=pileup_object)
			con_gen = class_generator(Pileup_line, file=consensus)
			get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
			parallel_gen = parallel_generator([pileup_gen, con_gen], [get_pos_and_id, get_pos_and_id])
			for [pile, con] in parallel_gen:
				if pile is not None and con is not None:
					out.write(str(pile[0])+'\n')
	return output_list


if __name__ == '__main__':
	arguments = docopt(__doc__)

	if arguments['-u'] or arguments['--unsorted']:
		sorted = False
	else:
		sorted = True


	pileups = arguments['PILEUPS']
	kpercent = arguments['--percentage']
	output_dir = arguments['--dir']
	consensus = arguments['--consensus']
	if consensus is None:
		consensus = tempfile.TemporaryFile()
	consensus_file = create_consensus_file(pileups, kpercent, consensus, sorted)
	filter_pileups_by_consensus(pileups, consensus_file, sorted, output_dir)

