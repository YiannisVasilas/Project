#!/usr/bin/env python3

""" Author: The_Parsers (Vasilas Ioannis )
1/12/2017
Usage: Running samtools
  
"""
from sys import argv
import subprocess
import os.path


def run_samtools(directory = ""): Must but the locacion of directory
    cmd = "for sam in %s*.sam ; do samtools view -bS $sam | samtools sort \
- $sam.bam ; done"%(directory)
    print cmd
    output_check = subprocess.check_output(cmd, shell=True)
    call_check = subprocess.check_call(cmd, shell=True)
    return call_check #must be 0


 
