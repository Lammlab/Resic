"""pileup_sorting.

Usage:
  pileup_sorting.py param
  pileup_sorting.py example
  pileup_sorting.py <pileup_in_file> <pileup_out_file>
  pileup_sorting.py (-h | --help)

Options:
  -h --help     Show this screen.

"""
from docopt import docopt
from Utility.Pileup_class import Pileup_line
from Utility.generators_utilities import key_sorted_gen, class_generator



#key for sorting according to id and position of each line
#reciving pileup line
def get_key_tuple_from_pileup(line: Pileup_line):
	return line.ref(),line.pos()[0]

def pileup_sort(in_pileup, out_pileup):
	with open(in_pileup , "r") as file_pileup:
		gen_pileup = class_generator(Pileup_line, file = file_pileup)
		gen_sorted_pileup=key_sorted_gen(get_key_tuple_from_pileup,file=file_pileup,gen=gen_pileup)
		open(out_pileup, "w").close()

		for line in gen_sorted_pileup :
			with open(out_pileup,"a+") as file1:
				file1.write(line)
				

def parameters_description():
	pass

def example_description():
	pass

if __name__ == "__main__":
	arguments = docopt(__doc__)

	if arguments['param']:
		parameters_description()

	if arguments['example']:
		example_description()

	if (arguments['<pileup_in_file>'] and arguments['<pileup_out_file>']):
		pileup_sort(arguments['<pileup_in_file>'],arguments['<pileup_out_file>'])



