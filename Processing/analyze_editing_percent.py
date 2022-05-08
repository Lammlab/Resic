"""analyze_editing_percent
Usage:
    analyze_editing_percent.py param
    analyze_editing_percent.py example
    analyze_editing_percent.py --reads_min=<thresh> --edit_min=<percent> --noise_max=<percent> [ --debug --tag=<tag_string> ] -s  <summary_output_filename> <pileup_files>...
    analyze_editing_percent.py --reads_min=<thresh> --edit_min=<percent> --noise_max=<percent> [ --debug --tag=<tag_string> --add_csv_headers --out_dir=<out> ]   <summary_output_filename> <pileup_files>...
    analyze_editing_percent.py -h | --help

Options:
    -h --help                       Show this screen.
    -r --reads_min=<thresh>         Minimum ammount of reads requires to flag a site as an editing site
    -e --edit_min=<percent>         Minimum percent of reads required from the most frequent bp change to flag as an editing site
    -n --noise_max=<percent>        Maximum percent of other bp changes to flag as an editing site
    -o --out_dir=<out>              Directory for putting output pileup files relevant if changing the pileup files default value is dir of the summary file       
    -s --summarise_only             Do not output new pileups. only produce the summary file Cannot be combined with --tag or --add_csv_headers             
    -t --tag=<tag_string>           A tag to add to sites that where considered editing sites
    -v --add_csv_headers            Out put the pileups with a csv style header
    -d --debug                      use for debug options - will not run the commands in parallel, but one by one
"""

from Utility.generators_utilities import class_generator
from Utility.Pileup_class import Pileup_line
import itertools
import os
from docopt import docopt
import sys
from Experiments.frontiers_jupyter.parallel_commands import parallel_commands


def is_pileup_line_edited(pileup_line: Pileup_line, read_thresh_hold=1, editing_min_percent_threshold=0.,
                          editing_max_percent_threshold=100, noise_percent_threshold=100, const_tag='',
                          editing_read_thresh=0):
    """

    adds editing percentage and noise percentage to pileup line, tags as editing site if passed thresholds
    :param: pileup_line, pileup line
    :param: read_thresh_hold, minimum number of reads required for a valid editing site
    :param: editing_min_percent_threshold, the minimum value of valid editing_percent
    :param: editing_max_percent_threshold, the maximum value of valid editing_percent
    :param: noise_percent_threshold, the maximum value of valid noise_percent
    :param: const_tag, string representing constant tag that will be added for lines that are considered as editing site
    :param add_percents: whether to add percents or to leave lines unchanged
    :return: (editing_type,pileup_line)  where editing_type is a tuple of the form (from_base,to_base)
    and pileupline is the oirginal line with the new fields in
    """
    clean_str = pileup_line.clean_base_string
    sense_string = clean_str.replace(',', '.')
    # TODO before use assumed letters are only upper case. I fix this assumption here
    sense_string = sense_string.upper()

    # Finds the candidate nucleotide reads
    nucl_A, nucl_C, nucl_T, nucl_G = 0, 0, 0, 0
    nucl_changes = {"A": 0, "C": 0, "G": 0, "T": 0}

    for nucl in list(sense_string):
        if (nucl == 'A' or nucl == 'C' or nucl == 'G' or nucl == 'T'):
            nucl_changes[nucl] = nucl_changes[nucl] + 1
        else:
            continue
    # get the maximal nucleous change, key is the value of the dict
    (candidate_nucl, candidate_nucl_reads) = max(nucl_changes.items(), key=lambda x: x[1])

    if candidate_nucl_reads == 0:
        editing_percent = 0.0
        noise_percent = 0.0
    else:
        editing_percent = round(float(100 * (candidate_nucl_reads / pileup_line.base_count)), 3)
        # Find the noise reads
        # all ACTGs that are not candidate_nucl
        noise_reads = 0

        for nucl in list(pileup_line.reads_string.upper()):
            if (nucl != candidate_nucl) and (nucl in ['A', 'C', 'T', 'G']):
                noise_reads += 1
        noise_percent = round(float(100 * (noise_reads / pileup_line.base_count)), 3)
    if is_reads_threshold_match(pileup_line, read_thresh_hold) and (editing_percent >= editing_min_percent_threshold
                                                                    and editing_percent <= editing_max_percent_threshold) \
            and (noise_percent <= noise_percent_threshold) and (editing_percent is not 0.0) \
            and (candidate_nucl_reads >= editing_read_thresh):
        is_editing_site = True
    else:
        is_editing_site = False
        '''
        if (pileup_line._gene_pos == "10517276"):
            print("pileupline",pileup_line, "read_thresh_hold",read_thresh_hold, "editing_percent",editing_percent, 
                "noise_percent", noise_percent, "candidate_nucl_reads", candidate_nucl_reads)
            print("editing_min_percent_threshold",editing_min_percent_threshold, "editing_max_percent_threshold",
                    editing_max_percent_threshold, "noise_percent_threshold",noise_percent_threshold, "editing_read_thresh", editing_read_thresh)
        '''
    # change pileup line
    pileup_line.tags["editing_percent"] = str(editing_percent)
    pileup_line.tags["noise_percent"] = str(noise_percent)
    if (const_tag != '') and is_editing_site:
        pileup_line.tags["const_tag"] = const_tag

    # calculate editing type
    if is_editing_site:
        # TODO patch change to make sure we dont get lowercase reference or nucl change
        editing_type = (pileup_line.reference.upper(), candidate_nucl.upper())
        # editing_type = (pileup_line.reference, candidate_nucl)
    else:
        editing_type = 'unchanged'

    return editing_type, pileup_line


def is_reads_threshold_match(pileup_line, treshold):
    return pileup_line.base_count >= treshold


def filter_pileup_by_categories(pileup_filename, output, reads_threshold=None, any_change=None,
                                editing_min_percent=None,
                                editing_max_percent=None, noise_percent=None, editing_read_thresh=0):
    """
    this function filtered a pileup file by reads threshold *or* by any change (editing percent > 0) *or* by editing
    percent threshold and noise_percent threshold
    :param pileup_filename: pileup file path
    :param output: the output path for the filtered pileup
    :param reads_threshold: if we want to filter by reads threshold - put a number in this var
    :param any_change: if we want to filter by any change (more then 0 editing percent) put true in this var
    :param editing_min/max_percent: if we want to filter by editing percent threashold, put threshold here
    :param noise_percent: if we want to filter by noise_percent threashold, put threshold here
    """

    with open(pileup_filename) as pileup, open(output, "w") as out:
        for line in class_generator(Pileup_line, file=pileup):

            if reads_threshold is not None and is_reads_threshold_match(line, reads_threshold):
                out.write(str(line) + "\n")

            if any_change:
                editing_type, pileup_line = is_pileup_line_edited(line)
                if editing_type != 'unchanged':
                    out.write(str(line) + "\n")

            if editing_min_percent is not None:
                if all([editing_max_percent, noise_percent]):
                    editing_type, pileup_line = is_pileup_line_edited(line,
                                                                      editing_min_percent_threshold=editing_min_percent,
                                                                      editing_max_percent_threshold=editing_max_percent,
                                                                      noise_percent_threshold=noise_percent,
                                                                      editing_read_thresh=editing_read_thresh)
                elif editing_max_percent is not None:
                    editing_type, pileup_line = is_pileup_line_edited(line,
                                                                      editing_min_percent_threshold=editing_min_percent,
                                                                      editing_max_percent_threshold=editing_max_percent,
                                                                      editing_read_thresh=editing_read_thresh)
                elif noise_percent_threshold is not None:
                    editing_type, pileup_line = is_pileup_line_edited(line,
                                                                      editing_min_percent_threshold=editing_min_percent,
                                                                      noise_percent_threshold=noise_percent,
                                                                      editing_read_thresh=editing_read_thresh)
                else:
                    editing_type, pileup_line = is_pileup_line_edited(line,
                                                                      editing_min_percent_threshold=editing_min_percent,
                                                                      editing_read_thresh=editing_read_thresh)

                if editing_type != 'unchanged':
                    out.write(str(line) + "\n")


def get_header_line(const_tag=''):
    if const_tag == '':
        optional_fields_list = ['editing_min_percent', 'noise_percent']
    else:
        optional_fields_list = ['editing_min_percent', 'noise_percent', 'const_tag']

    columns_list = ['Reference Id', 'Gene Position', 'Reference', 'Base Count', 'Reads String', 'Quality String']
    columns_list += optional_fields_list
    columns_list.append('Misc Tags')

    header_line = '\t'.join(columns_list)
    return header_line


def analyse_editing_percent(pileup_file, out_file, summary_file=None, add_headers=False, summary_only=False,
                            min_editing=0,
                            max_noise=100, min_reads=1, edit_tag=''):
    """
    analyses pileup file editing sites and summarises it
    @param pileup_file: input pileup file name
    @param out_file: name of file to write new lines
    @param summary_file: File to put the summary_string
    @param add_headers: Boolean, whether to add csv like headers to new pileups
    @param summary_only: Boolean, whether to only generate summary file and not generate edited pileup files
    @param min_editing: minimal editing threshold, percent
    @param max_noise:  maximal noise threshold, percent
    @param min_reads: minimal reads per site to be considered an editing site, int
    @param edit_tag: tag to add to sites that are classified as editing sites
    @param kwargs:
    @return: dict with histogram of site editing types, side effect creates out_file a pile up file with tags
    indicating noise, editing percent and wether a given site was classified as an editing site
    """

    # with open in read and out write

    with open(pileup_file, 'r') as in_fp, \
            open(out_file, 'w') as out_fp:

        # "r","y","s","w","k","m","b","d","h","v","n"
        # intitalize zero counts in summary dict
        # summary_dict = {(nuc1, nuc2): 0 for nuc1, nuc2 in itertools.product('acgtACGT', repeat=2) if nuc1.upper() != nuc2.upper()}
        # TODO why do we need lower case references?

        summary_dict = {(nuc1, nuc2): 0 for nuc1, nuc2 in itertools.product('ACGTN', repeat=2) if
                        nuc1.upper() != nuc2.upper()}
        # add RYSWKMBDHV to string if needed

        summary_dict['unchanged'] = 0

        if add_headers:
            out_fp.writelines([get_header_line(edit_tag)])
            out_fp.write("\n")

            list_headers = ["editing_min_percent", "noise_percent"]
            if edit_tag:
                list_headers.append("const_tag")

        pileup_gen = class_generator(Pileup_line, file=in_fp)
        line_num = -1

        for line_num, pileup_line in enumerate(pileup_gen):
            edit_type, new_pileup_line = is_pileup_line_edited(pileup_line,
                                                               read_thresh_hold=min_reads,
                                                               editing_min_percent_threshold=float(min_editing),
                                                               noise_percent_threshold=float(max_noise),
                                                               const_tag=edit_tag)
            summary_dict[edit_type] += 1

            # add empty tag if we make it into a csv so the columns will align
            if add_headers and edit_type == 'unchanged' and edit_tag != '':
                new_pileup_line.tags["const_tag"] = ''

            # if not summarise only, print the pileup line
            if not summary_only:
                # out_fp.write("\n")
                if not add_headers:
                    out_fp.writelines([str(new_pileup_line)])
                    out_fp.write("\n")
                else:
                    out_fp.writelines([new_pileup_line.line_to_csv_with_short_tags(list_headers)])
                    out_fp.write("\n")

        # dont skip file, create empty summary
        # if line_num == -1 :
        #    open(summary_file, "w").close()
        #    return #when the file is empty - skip the file!
        ##out_fp.write("\n")

        # make summary counts into percentage
        total_line_number = line_num + 1

        summary_dict_sub = {key: float(val) / total_line_number if total_line_number != 0 else 0 for key, val in
                            summary_dict.items()}

        # printing the summary format to an individual file
        with open(summary_file, "w") as file1:

            summer_str = individual_summary_format(pileup_file, summary_dict_sub)

            file1.write(summer_str)

        return summer_str


def individual_summary_format(filename, summary_dict):
    """

    @param filename: the pileup file name, the dict refers to
    @param summary_dict: the summary dict of the pileup file
    @return: a string that represent the summary dict :  pileup_filename\nchange:percentage\tchange:percentage\t...
    """

    summer_str = ""
    for change, percentage in summary_dict.items():
        if summer_str == "":
            summer_str = str(change) + ":" + str(percentage)
        else:
            summer_str = summer_str + "\t" + str(change) + ":" + str(percentage)
    return filename + "\n" + summer_str


def individual_summary_format_to_dict(lines):
    """

    @param lines: the lines of the summary file
    @return: a dict that represent the lines (with the pileup name inside)
    """
    summary_dict = dict()
    summary_dict["name of file"] = lines.pop(0).split('\n')[0]

    for line in lines[0].split("\t"):
        change_and_percentage = line.split(":")
        summary_dict[change_and_percentage[0]] = change_and_percentage[1]
    return summary_dict


def summary_format(summary_dict):
    """
    formats summary dicts into desired format
    :param summary_dict: dictionary with percents per
    :return: str of dict in desired format
    """
    return str(summary_dict.items())


def analyse_multiple_editing_percent_files(pileup_files, out_pileups, summary_files, total_summary_file=None,
                                           add_headers=False,
                                           summary_only=False, min_editing=0.0, max_editing=100.0, max_noise=100.0,
                                           min_reads=1,
                                           edit_tag=None, parallel_limit=2, Disable_parallel=True):
    """
    analyses pileup file and summarises it
    :param pileup_files: list of pileup_files
    :param out_pileups: list of path names for pileups with editing percent
    :param summary_file: list of path names for summary of all pileups
    :param total_summary_file: list of path name for summary of all pileups
    :param add_headers: Boolean, whether to add csv like headers to new pileups
    :param summary_only: Boolean, whether to only generate summary file and not generate edited pileup files
    :param min_editing: minimal editing threshold, percent (with the dot : 0.01 = 1%)
    :param max_noise: maximal noise threshold, percent
    :param min_reads: minimal reads per site to be considered an editing site, int
    :param edit_tag: tag to add to sites that are classified as editing sites
    :param parallel_limit: cap on number of parallel processes to run (passes to parallel_commands)
    :param Disable_parallel: Bool wether to do processes one by one or run in parallel (passes to parallel_commands)
    :return:
    """

    # collecting commands for parallel
    commands = []
    for in_file, out_file, summary_file in zip(pileup_files, out_pileups, summary_files):
        # False,False map to add_headers and summary_only
        command = [analyse_editing_percent, in_file, out_file, summary_file,
                   add_headers, summary_only, min_editing, max_noise, min_reads, edit_tag]
        commands.append(command)

    # run parallel editing percent + summary printing for the pileup files
    parallel_commands(commands, parallel_limit=parallel_limit, Disable_parallel=Disable_parallel)

    # if total summary was requested,
    # get the summary strings to write to the summary output file
    if total_summary_file != None:
        lines_list = get_full_summary_str(summary_files)

        with open(total_summary_file, "w") as sumfile:
            sumfile.write("\n".join(lines_list))

    return


def get_full_summary_str(out_files):
    """

    @param out_files: the list of the names of the pileup output files (with the editing percent)
    @return: string contains the lines of the required summary file
    """
    # creating a dict of change and percent from each file summary
    keys_set = dict()
    summary_list = []
    for file_summer in out_files:
        with open(file_summer, "r") as file1:
            lines = file1.readlines()
        if len(lines) == 0:
            continue
        summary_dict_file = individual_summary_format_to_dict(lines)
        for key in summary_dict_file.keys():
            keys_set[key] = 1
        summary_list.append(summary_dict_file)
    lines_list = []

    # preparing header

    keys_set.pop("name of file", None)
    header_line = "name of file"
    keys_list = list(keys_set.keys())
    for header in keys_list:
        header_line = header_line + "\t" + str(header)
    lines_list.append(header_line)
    for dict_summer in summary_list:
        line = dict_summer["name of file"]

        for change in keys_list:
            if change in dict_summer.keys():
                line = line + "\t" + dict_summer[change]  # when the base change is in that summeray file
            else:
                line = line + "\t "  # when the change is not there - put empty stub
        lines_list.append(line)
    return lines_list


def print_params():
    string = "The parameters are:" + "\n" \
             + "-editing_min_percent_threshold: the minimum editing percent that should be for the line that is considered as editing site.\n" \
             + "-editing_max_percent_threshold the maximum editing percent that should be for the line that is considered as editing site.\n" \
             + "-noise_percent_threshold: the maximum noise percent (\"5\" = 5%) that should be for the line that is considered as editing site.\n" \
             + "-out_filename: the name of the statistic output file (file that contain for each file the percent of lines that are considered as editing site).\n" \
             + "-pileup_files: list of pileup filenames.\n" \
             + "-out: Option to specify the name of the output directory (the default is the current directory).\n" \
             + "-const_tag: Option to add a constant tag string to each line in the file that is considered editing site.\n" \
             + "-noadd: Flag indicating that the pileup files already have 2 extra columns of editing_percent and noise_percent, and therefore I will not add these 2 extra columns.\n" \
             + "-addheaders: Choose this flag if you would like to add headers to each pileup file (used in order to parse the pileup to csv format easily).\n"
    print(string)


def print_example():
    input_pileup_str = "seq1\t272\tA\t24\t,.$....G,,.,.G...G,..^+.\t<<<<<<<<<<<<<<<<<<<<<\n" \
                       + "seq1\t273\tC\t23\t,.....,,.,A,..,.,AAA\t<<<<<<<<<<<<<<<<<<<<\n" \
                       + "seq1\t275\tT\t23\t,.....,,AAA,GG..,..A\t<<<<<<<<<<<<<<<<<<<<\n" \
                       + "seq1\t275\tT\t23\t,.....,,,,,,....,..,\t<<<<<<<<<<<<<<<<<<<<\n" \
                       + "seq1\t288\tT\t24\t,.$AAAAAAAAA.AAAAAAA\t<<<<<<<<<<<<<<<<<<<<<\n" \
                       + "seq1\t285\tG\t23\tAAAAAACCCCC,,,,.,..A\t<<<<<<<<<<<<<<<<<<<<\n"

    print("The pileup file is:\n" + input_pileup_str + "\n")

    output_pileup_str = "seq1\t272\tA\t24\t,.$....G,,.,.G...G,..^+.\t<<<<<<<<<<<<<<<<<<<<<\tediting_percent:Z:12.5\tnoise_percent:Z:0.0\n" \
                        + "seq1\t273\tC\t23\t,.....,,.,A,..,.,AAA\t<<<<<<<<<<<<<<<<<<<<\tediting_percent:Z:20.0\tnoise_percent:Z:0.0\n" \
                        + "seq1\t275\tT\t23\t,.....,,AAA,GG..,..A\t<<<<<<<<<<<<<<<<<<<<\tediting_percent:Z:20.0\tnoise_percent:Z:10.0\n" \
                        + "seq1\t275\tT\t23\t,.....,,,,,,....,..,\t<<<<<<<<<<<<<<<<<<<<\tediting_percent:Z:0.0\tnoise_percent:Z:0.0\n" \
                        + "seq1\t288\tT\t24\t,.$AAAAAAAAA.AAAAAAA\t<<<<<<<<<<<<<<<<<<<<<\tediting_percent:Z:80.0\tnoise_percent:Z:0.0\n" \
                        + "seq1\t285\tG\t23\tAAAAAACCCCC,,,,.,..A\t<<<<<<<<<<<<<<<<<<<<\tediting_percent:Z:35.0\tnoise_percent:Z:25.0\n"

    print("The output pileup file (with editing percent and noise percent) is:\n" + output_pileup_str)

    statistic_output_str = "\tName of File\tPercent of Valid Lines\n" \
                           + "0\tpileup_filename\t50.0\n"

    print(
        "The statistic output with editing_min_percent_threshold=30 and noise_percent_threshold=50 is:\n" + statistic_output_str)


if __name__ == '__main__':
    args = docopt(__doc__)
    if args["param"]:
        print_params()
        sys.exit()

    if args["example"]:
        print_example()
        sys.exit()

        # Extracting the parameters
    pileups_names = args['<pileup_files>']

    count_threshold = int(args['--reads_min'])
    editing_min_percent_threshold = float(args['--edit_min'])
    noise_percent_threshold = float(args['--noise_max'])

    if count_threshold < 1:
        raise ValueError("count minimum cannot be nonpositive")
    # Check that the range of the percent is correct
    if editing_min_percent_threshold < 0 and editing_min_percent_threshold > 100 or noise_percent_threshold < 0 or noise_percent_threshold > 100 or editing_max_percent_threshold < 0 or editing_max_percent_threshold > 100:
        raise ValueError("In editing or noise, Percent is not valid")

    summary_filename = args['<summary_output_filename>']

    out_dir = args['--out_dir']
    if out_dir is None:
        out_dir = os.path.dirname(summary_filename)
        out_pileups = [pile + "_edited" for pile in pileups_names]
    else:
        summary_filename = os.path.join(out_dir, summary_filename)
        pileup_base_names = [os.path.basename(pile) for pile in pileups_names]
        out_pileups = [os.path.join(out_dir, pile) for pile in pileup_base_names]

    # if out dir doesnt exist make it
    os.makedirs(out_dir, exist_ok=True)

    # give all summary pileups the same sumamry
    summary_filename = [summary_filename] * (len(pileups_names))

    const_tag = args['--tag']
    summarise_only = args['--summarise_only']
    need_to_add_headers = args['--add_csv_headers']

    debug_op = args['--debug']
    '''if debug_op is None:
        debug_op = False
    else:
        debug_op = True'''

    analyse_multiple_editing_percent_files(pileups_names, out_pileups,
                                           summary_files=summary_filename,
                                           add_headers=need_to_add_headers,
                                           summary_only=summarise_only,
                                           min_editing=editing_min_percent_threshold,
                                           max_editing=editing_max_percent_threshold,
                                           max_noise=noise_percent_threshold,
                                           min_reads=count_threshold,
                                           edit_tag=const_tag,
                                           )
