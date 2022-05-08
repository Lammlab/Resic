import subprocess
import os
import logging


def multiline_sort(number_of_lines, new_delim, field_number, multi_filename, secondary_field=None, output_file=None,
                   overwrite=False):
    """
    the script sorts a file of any multi-line format (SAM, FASTQ, FASTA, etc...), according to a chosen line in the format.
    ** function is case_sensitive

    :param number_of_lines: number of lines in each entry of the format
    :param new_delim: any character that is not within the alphabet of possible characters that are used in the format
    :param field_number the number of the field/line by which to sort
    :param multi_filename: name of the file containing the data
    :param secondary_field: optional secondary field to sort by
    :param output_file: optional filename to write output to. if not specified output goes to input filename with .sorted extension
    :param overwrite: Boolean flag, if True overwrites outputfile if file already exists. If false and output_feil exists, dont change it
    :return: output in the provided filename with .nt extension
    """

    # make default path if one wasnt supplied
    if output_file is None:
        splited = multi_filename.split(".")
        output_file = ".".join(splited[:-1]) + "_sorted." + splited[-1]

    if (not overwrite) and os.path.isfile(output_file):
        logging.warning(f"Not sorting {multi_filename} since {output_file} already exists")
        return output_file

    if secondary_field is None:
        secondary_field = field_number
    with open(multi_filename, 'r'), open(output_file, 'w'):
        command = "cat {0} | sed -r '/^[ \t]*$/d' | paste -d \"{1}\" {2} | LC_ALL=C sort  -t '~' -k{3},{3}n -k{4},{4}n | sed 's/{1}/\\n/g' " \
                  "> {5}".format(multi_filename, new_delim, ' -' * int(number_of_lines), field_number, secondary_field,
                                 output_file)

        # TODO: sort -t '~' should be sort -t {new_delim} ???
        # TODO: tests for sorting by secondary field
        # print(command, '\n')
        try:
            subprocess.run(command, shell=True).returncode
        except Exception as e:
            raise e
    return output_file


# 'combine multi-line into one line seperated by chosen delimiter'
# 'sort according to selected field'
# 'split back into lines'

### Tests are in genome_3nt_test. TODO: move sort tests to seperate file

def multiline_sort_pileup(number_of_lines, new_delim, field_number, multi_filename, secondary_field=None,
                          output_file=None, overwrite=False):
    # the script sorts a file of any multi-line format (SAM, FASTQ, FASTA, etc...), according to a chosen line in the format.
    # ** function is case_sensitive
    # :param number_of_lines: number of lines in each entry of the format
    # :param new_delim: any character that is not within the alphabet of possible characters that are used in the format
    # :param field_number the number of the field/line by which to sort
    # :param multi_filename: name of the file containing the data
    # :param secondary_field: optional secondary field to sort by
    # :param output_file: optional filename to write output to. if not specified output goes to input filename with .sorted extension
    # :param overwrite: Boolean flag, if True overwrites outputfile if file already exists. If false and output_feil exists, dont change it
    # :return: output in the provided filename with .nt extension

    # make default path if one wasnt supplied
    if output_file is None:
        splited = multi_filename.split(".")
        output_file = ".".join(splited[:-1]) + "_sorted." + splited[-1]

    if (not overwrite) and os.path.isfile(output_file):
        logging.warning(f"Not sorting {multi_filename} since {output_file} already exists")
        return output_file

    if secondary_field is None:
        secondary_field = field_number
    with open(multi_filename, 'r'), open(output_file, 'w'):
        command = "cat {0} |sed -r '/^[ \t]*$/d'| LC_ALL=C sort -k{3},{3}V -k{4},{4}n" \
                  "> {5}".format(multi_filename, new_delim, ' -' * int(number_of_lines), field_number, secondary_field,
                                 output_file)

        # TODO: sort -t '~' should be sort -t {new_delim} ???
        # TODO: tests for sorting by secondary field
        # print(command, '\n')
        try:
            subprocess.run(command, shell=True).returncode
        except Exception as e:
            raise e
    return output_file


# 'combine multi-line into one line seperated by chosen delimiter'
# 'sort according to selected field'
# 'split back into lines'

def multiline_to_oneline(number_of_lines, new_delim, multi_filename, output_file=None):
    '''
    convert multiline format to single line format with chosen delimiter

    :param number_of_lines: number of lines in each entry of the format
    :param new_delim: any character that is not within the alphabet of possible characters that are used in the format
    :param multi_filename: name of the file containing the data
    :return: output_file, if none specified, multi_filename with .oneline extension
    '''
    if output_file is None:
        splited = multi_filename.split(".")
        output_file = ".".join(splited[:-1]) + "_oneline." + splited[-1]
    with open(multi_filename, 'r'), open(output_file, 'w'):
        command = "cat {0} | sed -r '/^[ \t]*$/d' | paste -d \"{1}\" {2} > {3}".format(multi_filename, new_delim,
                                                                                       ' -' * int(number_of_lines),
                                                                                       output_file)
        try:
            subprocess.run(command, shell=True)
        except Exception as e:
            raise e
    return output_file
