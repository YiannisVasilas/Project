#!/usr/bin/evn python
"""A Python script to run Velveth and Velvetg on sequences in fastq
format. Kmer length is set to 19.
python velvet.py <filename> <insertlength>
Given:
-- Fastq file from sequences. 
-- The length of the inserted sequence.
Return:
-- Contigs in fasta format. File is named "contigs.fa".
"""

from sys import argv
import subprocess
import tempfile

__author__ = "Ioannis"
__date__ = "15th November 2018"
__version__ = "v1"
__email__ = ""

def runVelveth(workdir, file1):
    cmd = "velveth {} 19 -fastq -short {}".format(workdir, file1)
    res = subprocess.check_call(cmd, shell=True)
    return res

def runVelvetg(workdir, ins_length):
    cmd = "velvetg {} -ins_length {} -exp_cov auto -min_contig_lgth 150".format(workdir, ins_length)
    res = subprocess.check_call(cmd, shell=True)
    return res

if __name__ == "__main__":
    file1 = argv[1]
    ins_length = 407
    workdir = tempfile.mkdtemp()
    runVelveth(workdir, file1)
    runVelvetg(workdir, ins_length)
    output = workdir + "/contigs.fa"
    cmd = "mv {} .".format(output)
    subprocess.check_call(cmd, shell=True)

