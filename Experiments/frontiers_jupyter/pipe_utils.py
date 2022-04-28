import os
from tkinter import filedialog
from tkinter import *
import itertools
import subprocess
from functools import partial
import Processing.genome_3nt as genome_3nt


# TODO note that orchestration is done in functions with clear explanation
def all_genome_pair_combinations():
    yield from itertools.combinations("ACGT", 2)


def genome_3nt_factory(from_nt, to_nt):
    """
    generates the pre and post 3nt genome processing functions for a certain NT pairing
    returns 2 function objects
    """
    pre = partial(genome_3nt.pre, nt_replacement=[from_nt, to_nt])
    post = partial(genome_3nt.post, nt_replacement=[from_nt, to_nt])
    return pre, post


# todo - this is a stub function - need to remove it and replace it with the calling of the real one
def transcriptom_func(one, two):
    print("placeholder: transcriptome")
    return two, one


def genome_3nt_all_combination_spec():
    """
    returns a list of 3-ples of the form ('X_Y',pre,post) where pre and post are the
    3nt genome preprocessing and postprocessing functions of X to Y genome mapping
    for all combinations of 2 different nucleotides X,Y
    """
    three_nt_spec = []
    for from_nt, to_nt in all_genome_pair_combinations():
        name = "%s_%s" % (from_nt, to_nt)
        pre, post = genome_3nt_factory(from_nt, to_nt)
        three_nt_spec.append((name, pre, post))

    return three_nt_spec


# file selection screen function
def files_selector():
    """
    returns the filenames list the user picked in the popup window
    """
    root = Tk()
    filenames = filedialog.askopenfilenames(initialdir=os.getcwd(), title="Select files",
                                            filetypes=(("all files", "*.*"), ("fastq files", "*.fastq"),
                                                       ("pileup files", "*.pileup"), ("fasta files", "*.fasta")))
    filenames_list = root.tk.splitlist(filenames)
    root.destroy()
    return list(filenames_list)


def file_selector():
    """
    returns the filename the user picked in the popup window
    """
    root = Tk()
    filename = filedialog.askopenfilename(initialdir=os.getcwd(), title="Select files",
                                          filetypes=(("all files", "*.*"), ("fastq files", "*.fastq"),
                                                     ("pileup files", "*.pileup"), ("fasta files", "*.fasta")))
    root.destroy()
    return filename


def folder_selector():
    """
    returns the folder the user picked in the popup window
    """
    root = Tk()
    folder_selected = filedialog.askdirectory()
    root.destroy()
    return folder_selected


def print_structure(startpath):
    """
    prints the directory structure from startpath
    """
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * level
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))


def call_process(command):
    res = subprocess.call(command, shell=True)
    if res:
        print("error")
