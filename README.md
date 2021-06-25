# RESIC A tool for comprehensive adenosine to inosine RNA Editing Site Identification and Classification

This is the code repository of the [RESIC](https://www.biorxiv.org/content/10.1101/2021.04.11.439401v1) Paper.

RESIC is  an efficient pipeline that combines several approaches for the detection and classification of RNA editing sites. The pipeline can be used for all organisms and can use any number of RNA-sequencing datasets as input. RESIC provides: 

1. The detection of editing sites in both repetitive and non-repetitive genomic regions
2. The identification of hyper-edited regions
3. Optional exclusion of polymorphism sites to increase reliability, based on DNA, and ADAR-mutant RNA sequencing datasets, or SNP databases.

RESIC implements a graph aligner module that enables modular composition of new alignment schemes. 



### Installing 

resic requires:
* [Bowtie version 0.12.7](https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download)
* [Samtools 1.7](https://github.com/samtools/samtools/releases/tag/1.7)

```bash
# download the repo
git clone https://github.com/Lammlab/Resic
cd Resic

# install python dependencies
pip install requirements.txt

# install samtools and bowtie

# add modules to PYTHONPATH
echo 'export PYTHONPATH="$PYTHONPATH:${HOME}/Resic"' > ~/.bashrc

```

### Running

Edit the function `real_pipe` in the following file and then run it:

```bash
python Experiments/forontiers_jupyter/resic_run_file.py
```

