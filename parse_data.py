#!/usr/bin/env python3

"""
Author: Yan Wang
Parse the sample data file, extract data filenames and settings
"""

from sys import argv
import subprocess
import os.path

def parse_data_set(samples):
    '''parse the file of samples need to be analysed
    return the file name of reference genome and dict of groups {sample_name:settings} 
    
    samples: file, file of samples need to be analysed
    '''
    group = {}
    for line in samples:
        if not line.strip():
            continue
        else:
            if line.startswith('reference:'):
                ref_fasta = line.split()[1]
                ref_gff3 = line.split()[2]
            elif not line.startswith('#'):
                parts = line.split()
                group[parts[0]] = parts[1:]
    return ref_fasta,ref_gff3,group
    
if __name__ == '__main__':
    with open(argv[1]) as data_file:
        ref_fa, ref_gff, groups = parse_data_set(data_file)
        print(ref_fa,ref_gff)
        print(groups)
        
