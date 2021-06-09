#Exceptions: Fastq class
class MyException(Exception):
	def __init__(self):
		self.message = ''
	def __repr__(self):
		return self.message

class InValidFastQFile(MyException):
	def __init__(self):
		self.message=''
class NumOfLinesNotDivisibleBy4(InValidFastQFile):
	def set_message(self, num_of_lines):
		self.message = 'FastQ file has ' + str(
			num_of_lines) + ' which is not divisible by 4 (the number of lines of each sequence)'
		return self
class InValidSequence(InValidFastQFile):
	def set_message(self, approx_line):
		self.message = 'Invalid Sequence of FastQ File: Approximate line = ' + str(approx_line)
		return self
class FirstLineNotStartWithAt(InValidSequence):
	def set_message(self , line):
		self.message += 'Line ' + str(line) + ' does not start with \'@\''
		return self
class ThirdLineNotStartWithPlus(InValidSequence):
	def set_message(self,line):
		self.message += 'Line ' + str(line) + ' does not start with \'+\''
		return self
class SecondAndForthLineNotSameSize(InValidSequence):
	def set_message(self, line):
		self.message += 'Lines ' + str(line+1) + ',' + str(line+3) + ' are not of the same length'
		return self

class WrongUsage_FastQClass(MyException):
	def __init__(self):
		self.message = 'Wrong Usage of FastQ Class:'
class InValidCuttingIndices(WrongUsage_FastQClass):
	def __init__(self):
		self.message = 'Invalid Filtering Indices'

class Fastq:
	"""
	first string - identifier
	second string - sequence
	third string - identifier (may be different than first line)
	forth line - confidence (of each nucleotide in the sequence)
	"""

	def __init__(self, strings):
		"""
		:param strings: The 4 lines (as a list of strings) which composes each fastq sequence
		"""
		if strings[0][0] != "@":
			raise FirstLineNotStartWithAt
		if strings[2][0] != "+":
			raise ThirdLineNotStartWithPlus
		if len(strings[1]) != len(strings[3]):
			raise SecondAndForthLineNotSameSize
		strings = [item.strip("\n") for item in strings]
		self.strings = strings
	def __repr__(self):
		return "\n".join(self.strings) + "\n"
	#TODO: removed newlines between lines, remember to fix
	def cut_seq(self, start_index, end_index):
		"""
		:return: FastQ object, but the sequence and the quality score are cut according to the indices
		"""
		if not isinstance(start_index,int) or not isinstance(end_index,int) or end_index <= start_index:
			raise InValidCuttingIndices
		temp = Fastq([self.strings[0] ,
					 self.strings[1][start_index:end_index] ,
					 self.strings[2] ,
					 self.strings[3][start_index:end_index]])
		return temp

	def get_seq(self):
		return self.strings[1]
	def get_id(self):
		return self.strings[0].split("@")[1]


	def get_id(self):
		return self.strings[0]

	def __getitem__(self, slice):
		"""
		:param slice: A slice [start:end] of 2 integers specifying a range in the sequence
		:return: a FastQ sequence cut in this manner: 1st,3rd lines remain the same; 2nd,4th lines are cut according to [start:end]
		"""
		try:
			return self.cut_seq(slice.start,slice.stop)
		except InValidCuttingIndices as e:
			raise e
	def __len__(self):
		"""
		:return: Length of nucleutide sequence
		"""
		return len(self.strings[1])