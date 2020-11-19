#!usr/bin/env python3
"""
Authors: Joeri, Marieke

script that splits fastq.gz file with concatenated reads to 2 files
input: fastq.gz file
output two fastq.gz files for left and right reads

Functions:
    unzipper: unzips .gz file
    fastq_splitter: splits concatenated fastq file
    zipper: zips fastq file to fastq.gz
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
    
    if os.path.exists(out_fn):
        print('unzipping tool not run, {} already exists.'.format(out_fn))
        return out_fn
    cmd = 'gunzip {}'.format(in_fn)
    subprocess.check_call(cmd, shell = True)
    return out_fn

def fastq_splitter(input_fn, read_length = 100):
    """ splits fastq file with concatenated reads to left and right file
    
    Keyword arguments:
        input_fn: string, name of input fastq file
        read_length: integer, length of one part of the paired read
    Returns:
        tuple with names of both output files.
    """
    out_fn_left = input_fn[:-6] + '_left.fastq'
    out_fn_right = input_fn[:-6] + '_right.fastq'
    
    with open(input_fn, 'r') as input_file,\
         open(out_fn_left, 'w') as left_file,\
         open(out_fn_right, 'w') as right_file:
             for line in input_file:
                 if line.startswith('@') or line.startswith('+'):
                     left_file.write(line)
                     right_file.write(line)
                 else:
                     left_file.write(line[:read_length]+'\n')
                     right_file.write(line[-read_length-1:])
    return (out_fn_left, out_fn_right)
        
def zipper(in_fn):
    """
    runs gzip tool, zips fastq.gz file  
      
    Keyword arguments:
        in_fn: string, name of file to be zipped
    Returns:
        string, filename of zipped file
    """
    print("zipping file...")
    out_fn = in_fn + '.gz'
    if not os.path.exists(in_fn):
        raise ValueError("{} does not exist".format(in_fn))
    if os.path.exists(out_fn):
        print('zipping tool not run, {} already exsists.'.format(out_fn))
        return out_fn
    cmd = 'gzip {}'.format(in_fn)
    subprocess.check_call(cmd, shell = True)
    return out_fn
    

if __name__ == "__main__":
    # get filename
    filename = argv[1]
    
    # unzip file
    #unzipped_filename = unzipper(filename)
    #print(unzipped_filename)
    
    # parse and split fastq into two files
    left_out, right_out = fastq_splitter(argv[1])
    
    # zip splitted files
    zipped_left_out = zipper(left_out)
    zipped_right_out = zipper(right_out)
    
    print(zipped_left_out)
    print(zipped_right_out)
