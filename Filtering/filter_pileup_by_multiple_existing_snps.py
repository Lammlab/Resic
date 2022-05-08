import Utility.generators_utilities as gen_util
from Utility.Pileup_class import Pileup_line
from Utility.Vcf_class import VcfClass
from Utility.parallel_generator import parallel_generator
from Processing.pileup_sorting import pileup_sort
from Utility.multiline_sort import multiline_sort_pileup


def snp_detect(pileup_line: Pileup_line):
    return pileup_line.is_with_any_change()


def snp_algebra_from_vcf_file(pileup_filename, vcf_file, output_file=None, sorted_input=True):
    """
    :param pileup_filename: positive pileup to be filtered
    :param vcf_file: vcf file that is a database containing known SNPs
    :param snp_detector: binary functor that receives a pileup line and returns True if there is an snp in the line
    :param output_file: file to write output to
    :param sorted_input: binary flag, set to False if input positive pileups are not sorted
    :param keep_pos_in_neg: binary flag, set to False if we want positions in the negative pileup to be
    left out of the output file
    :return
    """

    # sorting pileups if necessary and preparing filenames
    if not sorted_input:
        multiline_sort_pileup(1, '~~~', 1, pileup_filename, 2, pileup_filename + "_sorted.pileup")
        sorted_positive_pileup = pileup_filename + "_sorted.pileup"
    else:
        sorted_positive_pileup = pileup_filename
    # with open(sorted_positive_pileup, "r") as sorted_positive_pileup:
    #   print("sorted_positive_pileup",sorted_positive_pileup.readlines())
    if output_file is None:
        splited = sorted_positive_pileup.split(".")
        output_file = ".".join(splited[:-1]) + "_no-snp." + splited[-1]
    get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)
    # logic of the function
    try:
        with open(vcf_file, 'r') as vcf_obj, open(sorted_positive_pileup, 'r') as pileup_obj, \
                open(output_file, 'w') as out:

            pos_gen = gen_util.class_generator(Pileup_line, file=pileup_obj)
            vcf_gen = gen_util.class_generator(VcfClass, file=vcf_obj)
            parallel_gen = parallel_generator([pos_gen, vcf_gen],
                                              [get_pos_and_id for i in range(2)])
            for listlist in parallel_gen:  # listlist: list of lists, each entry in listlist is a list of items from one generator
                if (listlist[0] is not None and listlist[
                    1] is None):  # listlist[0]: list of pileup_lines(should be 1 line) from the input pileup
                    out.write(str(listlist[0][0]) + "\n")
    except Exception as e:
        print(e)
    return output_file


def snp_algebra(pileup_filename, negative_pileup_list, snp_detector=snp_detect, output_file=None,
                sorted_input=True, keep_pos_in_neg=True):
    """

    :param pileup_filename: positive pileup to be filtered
    :param negative_pileup_list: list of negative pileups to filter positive pileup by
    :param snp_detector: binary functor that receives a pileup line and returns True if there is an snp in the line
    :param output_file: file to write output to
    :param sorted_input: binary flag, set to False if input pileups are not sorted
    :param keep_pos_in_neg: binary flag, set to False if we want positions in the negative pileup to be
    left out of the output file
    :return
    """

    # sorting pileups if necessary and preparing filenames
    sorted_negative_pileups = []
    if not sorted_input:
        pileup_sort(pileup_filename, pileup_filename + "_sorted.pileup")
        sorted_positive_pileup = pileup_filename + "_sorted.pileup"
        for pileup in negative_pileup_list:
            pileup_sort(pileup, pileup + "_sorted.pileup")
            sorted_negative_pileups.append(pileup + "_sorted.pileup")
    else:
        sorted_positive_pileup = pileup_filename
        sorted_negative_pileups = negative_pileup_list
    if output_file is None:
        splited = sorted_positive_pileup.split(".")
        output_file = ".".join(splited[:-1]) + "_no-snp." + splited[-1]

    # parameters for parallel generator
    neg_obj_list = [open(pile) for pile in sorted_negative_pileups]
    neg_gen_list = []
    get_pos_and_id = lambda x: (x.reference_id, x.gene_pos)

    # logic of the function
    try:
        for file_obj in neg_obj_list:
            neg_gen_list.append(gen_util.class_generator(Pileup_line, file=file_obj))
        with open(sorted_positive_pileup, 'r') as pileup_obj, open(output_file, 'w') as out:
            pos_gen = gen_util.class_generator(Pileup_line, file=pileup_obj)
            parallel_gen = parallel_generator([pos_gen, *neg_gen_list],
                                              [get_pos_and_id for i in range(len([pos_gen, *neg_gen_list]))])
            for listlist in parallel_gen:  # listlist: list of lists, each entry in listlist is a list of items from one generator
                if listlist[0] is not None:  # listlist[0]: list of pileup_lines(should be 1 line) from the input pileup
                    no_snp = True
                    position_is_in_neg = False
                    for pileup_list in listlist[1:]:
                        if pileup_list is not None:
                            position_is_in_neg = True
                            if snp_detector(pileup_list[0]):
                                no_snp = False
                    if keep_pos_in_neg is True:
                        if (no_snp is True) and (position_is_in_neg is True):
                            out.write(str(listlist[0][0]) + "\n")
                    else:  # get here if keep_pos_in_neg flag is False
                        if no_snp is True:
                            out.write(str(listlist[0][0]) + "\n")

    except Exception as e:
        print(e)

    finally:
        for fd in neg_obj_list:
            fd.close()
        return output_file
