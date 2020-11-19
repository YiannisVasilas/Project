#!usr/bin/env python3
"""
Author: Joeri, Marieke

"""
# import statements
from sys import argv
import os
import subprocess

# functions
def unzipper(in_fn):
    """
    runs gunzip tool, unzips fastq.gz file  
      
    Keyword arguments:
        in_fn: string, name of file to be unzipped
    Returns:
        string, filename of unzipped file
        
    """
    print("unzipping file...")
    out_fn = in_fn[:-3]
    if not os.path.exists(in_fn):
        raise ValueError("{} does not exist".format(in_fn))
    if os.path.exists(out_fn):
        print('... tool not run, {} already exsists.'.format(out_fn))
        return out_fn
    cmd = 'gunzip {}'.format(in_fn)
    subprocess.check_call(cmd, shell = True)
    return out_fn

def fastq_splitter(input_fn):
    out_fn_1eft = input_fn[:-6] + '_left.fastq'
    out_fn_right = input_fn[:-6] + '_right.fastq'
    with open(input_fn, 'r') as file_object:
        
        
        
    

if __name__ == "__main__":
    # get filename
    filename = argv[1]
    
    # unzip file
    unzipped_filename = unzipper(filename)
    print(unzipped_filename)
    
    # parse and split fastq into two files
    fastq_splitter(unzipped_filename)
    
    
