#TODO refactor

class IllegalArgument(Exception):
    def __init__(self):
        self.message = "you should enter name of pileup file, output file and at least one Nucleotide for ref base"

    def __str__(self):
        return self.message

    def __repr__(self):
        return self.message


def filter_pileup_by_reference_base_set(pileup,output,listBases):
    """
    :param pileup:  the pile up we want to filter
    :param listBases: a list of ref bases we want to filter by
    :param output: the output file that will contain lines with only the bases in listBases as ref base
    """
    with open(pileup,"r") as file1:
        line = file1.readline()
        with open(output,"w") as fileout:
            while line is not None and line != "":
                lineParts = line.split()
                if lineParts[2] in listBases:
                    fileout.write(line)
                line=file1.readline()



if __name__ == "__main__" :

    """
    input : python filter_pileup_by_reference_base_set.py pileup_filename output_filename ref_bases_list
            example : python filter_pileup_by_reference_base_set.py exmple.pileup output.pileup A G C

            pileup_filename - the file we want to filter from by ref bases - to keep the lines with those ref bases
            output_filename - the output file for the result filter
            ref_bases_list - a list of ref bases we want to keep in the output file

    output : the lines from the pileup_filename file that contain one of the ref_bases_list as ref base.

    """
    import sys
    if len(sys.argv)<3 : raise IllegalArgument
    pileup = sys.argv[1]
    output = sys.argv[2]
    
    numberOfBases = len(sys.argv) -3
    if numberOfBases == 0: 
        raise IllegalArgument
    if numberOfBases > 0:
        one = sys.argv[3]
        if numberOfBases > 1:
            two = sys.argv[4]
            if numberOfBases > 2:
                three = sys.argv[5]
                if numberOfBases > 3:
                    four = sys.argv[6]
                    if numberOfBases > 4:
                        five = sys.argv[7]
                        filter_pileup_by_reference_base_set(pileup,output,[one,two,three,four,five])
                    else:
                        filter_pileup_by_reference_base_set(pileup, output,[one, two, three, four])
                else:
                    filter_pileup_by_reference_base_set(pileup,output, [one, two, three])
            else:
                filter_pileup_by_reference_base_set(pileup,output, [one, two])
        else: 
            filter_pileup_by_reference_base_set(pileup,output, [one])




