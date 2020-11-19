#!/usr/bin/env python3

""" 
Ioannis Vasias 
Anna Pijnacker
Marijn Schipper
Fastqc 26/11/2019
FASTQ Quality 
"""

# imports
from sys import argv
import subprocess
import os.path

def loop_samples(path_to_folder):
    '''Return a list of all fastq files in the path specified
    
    Input:
        path_to_folder: The absolute path to the folder with all fastq reads
    Output:
        files: list -- list of filenames of all fastq files
    '''
    cmd = "ls {}".format(path_to_folder)
    files = subprocess.check_output(cmd, shell = True).decode("utf-8")\
        .split("\n")
    return files


def namer(input_file):
    """Determines outfilename from inputname
    input_file -- fastq file with treatment, index Seq, replicate and hpi
    Returns:
    new_name -- input_filename minus extension
    """
    outputname = str(input_file).split(".")[0].rpartition("/")[-1]
    return outputname
    
    
def fast_bias_removal(in_path, input_file, outpath_bias, outputname, first_to_keep = 14):
    """Run fastq_quality_trimmer on fastq file
    input_file -- string, filename of input FASTQ file
    threshold -- int, Quality threshold - nucleotides with lower 
        quality will be trimmed (from the end of the sequence), default = 20.
    lenth -- int, last base to keep, deafulet = 20
    Return:
    output -- a fastq file with the trimmed reads
    """
    if not os.path.exists("./{}".format(outpath_bias)):
        command1 = "mkdir {} ".format(outpath_bias)
        run = subprocess.check_output(command1,shell = True)
    new_output = "./{}/fastq_trimmed_{}.fq".format(outpath_bias, outputname)
    if os.path.exists(new_output): 
        return new_output
    command2 = 'zcat {}/{} | fastx_trimmer -f {} -o {} | gzip'\
        .format(in_path, input_file, first_to_keep, new_output)
    print(command2)
    run =subprocess.check_output(command2,shell = True)
    return new_output
    
def fast_trimming(in_file, outpath_trim, outputname, threshold =20, length = 40 ):
    """Run fastq_quality_trimmer on fastq file
    input_file -- string, filename of input FASTQ file
    threshold -- int, Quality threshold - nucleotides with lower 
        quality will be trimmed (from the end of the sequence), default = 20.
    lenth -- int, last base to keep, deafulet = 20
    Return:
    output -- a fastq file with the trimmed reads
    """
    if not os.path.exists("./{}".format(outpath_trim)):
        command1 = "mkdir {} ".format(outpath_trim)
        run = subprocess.check_output(command1,shell = True)
    new_output = "./{}/fastq_trimmed_{}.fq".format(outpath_trim, outputname)
    if os.path.exists(new_output): 
        return new_output
    command2 = 'fastq_quality_trimmer -i {} -Q 33 -t {} -o {} -l {}'\
        .format(in_file, threshold, new_output, length)
    print(command2)
    run =subprocess.check_output(command2,shell = True)
    return new_output

def faststats (in_file, outpath_stats, outputname):
    """ 
    Tool for running the fast quality stasts from Fastx-Toolkit 
        version 0.0.14  
        input: Fastq files
        output : txt files with fastq read statistics
    """
    if not os.path.exists("./{}".format(outpath_stats)):
        command = "mkdir {}".format(outpath_stats)
        run = subprocess.check_output(command,shell = True)
    new_output = "./{}/faststats_{}.txt".format(outpath_stats, outputname)
    if os.path.exists(new_output):
        return new_output
    command = "fastx_quality_stats -i {} -o {} ".format(in_file, new_output)
    run = subprocess.check_output(command, shell = True)
    return new_output
    

def fastquality (fast_stat, outputname, outpath_stats):
    
    """Tool for runnung the fast qualitu chart from Fastx_Toolkit
        version 0.0.1.4
        input = txt file from fast quality stats
        output = graph
    """
    
    if not os.path.exists("./{}".format(outpath_stats)):
        command = "mkdir {} ".format(outpath_stats)
        run = subprocess.check_output(command,shell = True)
    new_output = "./{}/qualityplot_{}.png".format(outpath_stats, outputname) 
    if os.path.exists(new_output):
        return new_output
    command = "fastq_quality_boxplot_graph.sh -i {} -o {}".format\
                (fast_stat,new_output)
    run = subprocess.check_output(command,shell = True)
    return new_output
        
def Fastqrunner(input_path, outpath_trim, outpath_stats, outpath_bias,\
 tr_setting):
    """ runs different fastx toolkit tools
    input_path = gzipped fastq files directory (str)
    outpath_trim = outpath name for trimmed files (str)
    outpath_stats = outpath name for fastq stats output (str)
    outpath_trim = outpath name for start bias trimmed files (str)
    tr_setting = True or False (bool). True toggles trimming on. False toggles
                    trimming off.
    returns:
        outpath of trimmed files (str)
    
    """
    pair_count = 0
    pair1 = ''
    if not os.path.exists(outpath_trim):
        for i in loop_samples(input_path):
            pair_count += 1
            names = namer(i)
            if names ==  '':
                continue
            if not os.path.exists('{}/fastq_trimmed_{}.fq'.format(outpath_trim, names)):
                bias = fast_bias_removal(input_path, i, outpath_bias, names)
                trim = fast_trimming(bias, outpath_trim, names)
            else: 
                trim = '{}/fastq_trimmed_{}.fq'.format(outpath_trim, names)
            if pair_count == 1:
                pair1 = trim
            elif pair_count == 2:
                if not os.path.exists('{}/fastq_trimmed_{}.fq_pairs_R2.fastq'.format(outpath_trim, names)):
                    print('running trimfix')
                    command = "python3 trimfix.py {} {}".format(pair1, trim)
                    run = subprocess.check_output(command,shell = True)
                pair_count = 0     
        stats = faststats(trim, outpath_stats, names)
        graph = fastquality(stats, names, outpath_stats)
    outpath_tr_files = './{}/'.format(outpath_trim)
    return outpath_tr_files

