#! /usr/bin/env python3

"""
Written by: Marieke van den Hurk
"""

from sys import argv
import subprocess
import os.path

def sam_to_sortedbam(infile):
    """Runs samtools to convert SAM file to sorted BAM file
    
    infile --- string; name of the input SAM file
    """
    
    outfile = infile[:-4]
    
    command = "samtools view -bS {} | samtools sort - {}".format(infile, outfile)
    subprocess.check_call(command, shell = True)

if __name__ == "__main__":
    samfilename = argv[1]
    sam_to_sortedbam(samfilename)

