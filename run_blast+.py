#!/usr/bin/evn python
"""A Python script to run a blast on mylist.
mylist is a set of tissuenames which canbe found at certain locations.
mylist = ["SRR122239", "SRR122243", "SRR122244", "SRR122245", "SRR122251",
          "SRR122252", "SRR122253", "SRR122254"]
run with:
python run_blast+.py [p|x|n]

Given:
p   blastp
x   blastx
n   blastn

Return:
output files are named

blast?_output_{0}.txt

with ? being the option and {0} being the string in the list

WARNING: Files are created and truncated without a message. Therefore check
whether the files exist before running this script!
"""

from sys import argv
import subprocess


def runblastp(num_threads, db, outfmt, query, out):
    cmd = "blastp -num_threads {} -db {} -outfmt {} -query {} -out {}".format(num_threads, db, outfmt, query, out)
    res = subprocess.check_call(cmd, shell=True)
    print res

def runblastx(num_threads, db, outfmt, query, out):
    cmd = "blastx -num_threads {} -db {} -outfmt {} -query {} -out {}".format(num_threads, db, outfmt, query, out)
    res = subprocess.check_call(cmd, shell=True)
    print res

def runblastn(num_threads, db, outfmt, query, out):
    cmd = "blastn -num_threads {} -db {} -outfmt {} -query {} -out {}".format(num_threads, db, outfmt, query, out)
    res = subprocess.check_call(cmd, shell=True)
    print res

if __name__ == "__main__":
    mylist = ["SRR122239", "SRR122243",
              "SRR122244", "SRR122245",
              "SRR122251", "SRR122252",
              "SRR122253", "SRR122254"]
    num_threads = 12
    outfmt = 7
    
    if argv[1] == 'p':
        for item in mylist:
            print "running blastp for {}".format(item)
            db = "/home/leest008/progs_nobackup/cro_scaffolds.min_200bp.proteindb"
            query = "/local/data/BIF30806_2016_2/project/groups/Fire_Breathing_Rubber_Duckies/project_mapping/\
            {0}_trimmed_MINLENGTH20_unmapped/transcripts.fa.transdecoder_dir/longest_orfs.pep".format(item)
            out = "blastp_output_{0}.txt".format(item)
            runblastp(num_threads, db, outfmt, query, out)
           
    elif argv[1] == 'x':
        for item in mylist:
            print "running blastx for {}".format(item)
            db = "/home/leest008/progs_nobackup/unisprot.db"
            query = "/local/data/BIF30806_2016_2/project/groups/Fire_Breathing_Rubber_Duckies/project_mapping/{0}_trimmed_MINLENGTH20_unmapped/transcripts.fa".format(item)
            out = "blastx_output_{0}.txt".format(item)
            runblastx(num_threads, db, outfmt, query, out)
           
    elif argv[1] == 'n':
        for item in mylist:
            print "running blastn for {}".format(item)
            db = "/home/leest008/progs_nobackup/cro_scaffolds.min_200bp"
            query = "/local/data/BIF30806_2016_2/project/groups/Fire_Breathing_Rubber_Duckies/project_mapping/\
            {0}_trimmed_MINLENGTH20_unmapped/transcripts.fa".format(item)
            out = "blastn_output_{0}.txt".format(item)
            runblastn(num_threads, db, outfmt, query, out)
