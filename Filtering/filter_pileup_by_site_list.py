"""Filter pileup by site list

    Usage:
        Filter_pileup_by_site_list param
        Filter_pileup_by_site_list example
        Filter_pileup_by_site_list <pileup_file> <out_file> <site_list>
        Filter_pileup_by_site_list -h | --help

    Options:
        -h --help   Show this screen
"""



import sys
from docopt import docopt
import re



def filter_pileup_by_site_list(pileup,output,listSites):
    """
    :param pileup:  the pile up we want to filter
    :param listSites: a list of sites - (chr , pos) - we want to filter by
    :param output: the output file that will contain lines with only the sites in listSites as ref base
    """
    with open(pileup,"r") as file1:
        line = file1.readline()
        with open(output,"w") as fileout:
            while line is not None and line != "":
                lineParts = line.split()
                if (lineParts[0],lineParts[1]) in listSites:
                    fileout.write(line)
                line=file1.readline()



def param_description():
    """
    Description for the script
    """
    print("The parameters are\n" +
          "pileup_file: the pileup file you want to filter\n" +
          "out_file: the name of the output file\n" +
          "site_list: the name of the file where each line is (char, pos) that we want to filter by\n")


def example_description():
    """
    Usage example for the script
    """
    print("Pileup input:\n" +
          "seq1\t272\tA\t24\t,.$.....,,.,.G...G,..^+.\t<<<<<<<<<<<<<<<<<<<<<\t\n" +
          "seq1\t273\tC\t23\t,.....,,.,A,..,.,AAA\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t274\tT\t24\t,.$.....,,.,.,,,,.,..^+.\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t275\tT\t23\t,.....,,AAA,GG..,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t276\tT\t24\t,.$.....CCC,.,,,,.,..^+.\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t277\tT\t23\t,.....,,***.,.,..,.,..A\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t278\tT\t23\t,.....,,.,.AA.,.,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t279\tT\t23\t,...AA,,A,A,....,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t280\tG\t24\t,.$..A.AA,.AAA,,A.,..^+.\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t281\tT\t23\t,.....,,.,.,,,AA,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t282\tT\t24\t,.$.....,,.,..,,,.,..^+.\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t283\tT\t23\t,.....,AAA.,CCCC,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t284\tT\t23\t,.....,,.,.,,,,.,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t285\tG\t23\tAAAAAACCCCC,,,,.,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t286\tT\t24\t,.$.....+3ABC,,.,A,,,,AA.....\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t287\tT\t23\t,.....,-2AA,.,.,..,.,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t288\tT\t24\t,.$AAAAAAAAA.AAAAAAA^+.\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t289\tG\t23\tTTTTTTTTTT..TTTTTTTTT\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq1\t290\tC\t23\tGG$GGGGGGGGGGGGGGGGG\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
          "seq2\t291\tC\t23\t...$GGGGGG....A......\t<<<<<<<<<<<<<<<<<<<<\t" +
          "Sites list file:\n" +
          "seq1, 273\n seq1, 274\n seq1, 282\n seq1, 283\n seq1, 284\n seq1, 288\n seq1, 290\n seq2, 291\n" +
          "Calling the script:\n" +
          "python filter_pileup_by_site_list pileup_file out_file site_list\n" +
          "And the output:\n" +
          "seq1\t273\tC\t23\t,.....,,.,A,..,.,AAA\t<<<<<<<<<<<<<<<<<<<<\t" +"\n"+
		  "seq1\t274\tT\t24\t,.$.....,,.,.,,,,.,..^+.\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
		  "seq1\t282\tT\t24\t,.$.....,,.,..,,,.,..^+.\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
		  "seq1\t283\tT\t23\t,.....,AAA.,CCCC,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
		  "seq1\t284\tT\t23\t,.....,,.,.,,,,.,..A\t<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
		  "seq1\t288\tT\t24\t,.$AAAAAAAAA.AAAAAAA^+.\t<<<<<<<<<<<<<<<<<<<<<\t" + "\n" +
		  "seq1\t290\tC\t23\tGG$GGGGGGGGGGGGGGGGG\t<<<<<<<<<<<<<<<<<<<<\t"+ "\n" +
		  "seq2\t291\tC\t23\t...$GGGGGG....A......\t<<<<<<<<<<<<<<<<<<<<\t")



"""
:param pileup:  the pile up we want to filter
:param output: the output file that will contain lines with only the sites in listSites as ref base
:param listSites: file where each line is (char, pos) that we want to filter by
:output: output file will contain only the lines in the given (chr, pos) positions
"""

if __name__ == "__main__" :
    args = docopt(__doc__)
    if args["param"]:
        param_description()
        sys.exit()

    if args["example"]:
        example_description()
        sys.exit()

    pileup = args['<pileup_file>']
    output = args['<out_file>']
    listSites = args['<site_list>']
    listPos = []
    with open(listSites,"r") as file1:
        for item in file1:
            iList=re.split( ",|\n" , item)
            if len(iList)==1:
                iList=item.split("\t")
            listPos.append((iList[0],iList[1]))
    filter_pileup_by_site_list(pileup,output,listPos)

