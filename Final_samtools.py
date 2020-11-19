#!/usr/bin/env python
"""
Marijn Schipper, Ioannis Vasilas, Mayke Schelhaas and Anna Pijnacker

Script for samtools sorting and from sam to bam 
"""
from sys import argv
import subprocess
import os.path


def remove_file(file_name):
    """Removes a specified 
    
    Keyword arguments
        file_name: string -- Name of the file to be removed
    Returns
        remove samfile
    """
    cmd = "rm {}".format(file_name)
    run = subprocess.check_output(cmd, shell =True)


def run_samtools(samfile, sam_output):
    """Runs samtools v1.3.1 to convert SAM to BAM
    
    Keyword arguments
        samfile -- output of hisat2, in samformat
    Returns:
        new_output -- name of outputfile in bam extension
    """
    sam_name = samfile.split(".")[0]
    new_output = "{}.bam".format(sam_name)
    if not os.path.exists(sam_output):
        cmd = "mkdir {}".format(sam_output)
        subprocess.check_call(cmd, shell = True)
    if os.path.exists(new_output): 
        return new_output
    cmd = "samtools view  -S -b {}  | samtools sort  -o ./{}/{}  ".format(\
     samfile, sam_output ,new_output)
    run = subprocess.check_output(cmd, shell = True)
    remove_file(samfile)
    return new_output
    
    

    
    
