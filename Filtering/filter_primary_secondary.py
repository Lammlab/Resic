"""Filter Primary Secondary

this script returns 2 filenames that are with ending ".secondary.fastq" and ".primary.fastq"
the secondary one contains the lines with seqs that their sizes are between 19 and 23
the primary one contains the lines with seqs that their sizes are between 25 and 30

Usage:
	filter_primary_secondary.py param
	filter_primary_secondary.py example
	filter_primary_secondary.py <fastq_filename>
	filter_primary_secondary.py -h | --help

Options:
	-h --help	show this screen
			
"""


from Filtering.fastq_seq_size_filtering import *
import os
from docopt import docopt

def filterToPrimarySecondaryLibs(fastqFilename):
	"""
	:param fastqFilename: the fastq we want to create primary and secondary libs from
	:return: the primary and secondary filenames that we filtered by primary siRNA 25-30nt ,secondary siRNA 19-23nt
	"""

	open(fastqFilename,'r').close()


	fastqName = fastqFilename.split(".")
	fastqName.pop(len(fastqName)-1)
	mainName = ".".join(fastqName)
	primary = mainName + ".primary.fastq"
	secondary = mainName + ".secondary.fastq"
	open(primary,"w").close()
	open(secondary,"w").close()
	non="non_temp.txt"
	open(non,"w").close()

	filterFastqFileBySize(fastqFilename, 25, 30,primary, non)	#primaryRange = (25,30)	
	filterFastqFileBySize(fastqFilename, 19, 23,secondary, non)	#secondaryRange = (19,23)
	
	os.remove(non)

	return primary, secondary




if __name__ == "__main__" :
	'''
	input : this script receives a fastq file that we would like to filter into 2 files by the size of the seq 
	output : this script returns 2 filenames that are with ending ".secondary.fastq" and ".primary.fastq"
	         the secondary one contains the lines with seqs that their sizes are between 19 and 23
	         the primary one contains the lines with seqs that their sizes are between 25 and 30
	'''
	import sys
	arguments = docopt(__doc__)
	fastq_filename = arguments['<fastq_filename>']
	primary , secondary = filterToPrimarySecondaryLibs(fastq_filename)

