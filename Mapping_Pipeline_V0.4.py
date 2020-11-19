#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

#Imports
from __future__ import division, print_function
from sys import stderr
from sys import exit as exit_program
import subprocess
import argparse
import os
import time
import re

#Classes

#Functions
def get_command_line_arguments():
    """Returns the given arguments form the shell, and shows a helpfile if necicairy
    """
    parser = argparse.ArgumentParser(description="Pipeline")
    
    parser.add_argument("-o", "--output", help="Name and path of the output bam file for kallisto", type = str, default = "out.bam")    
    
    parser.add_argument("-fq", "--fastqc", help="run fastqc, give the name and the path of the file", type = str)
    
    parser.add_argument("-p", "--preprocess", help="Use sickle (and cutadapt if -c is given) give the name and path of the input file", type = str)
    parser.add_argument("-c", "--cutadapt", help="the sequence that needs to be cut, if not given cutadapt will not be run", type = str, default = "")  
    
    parser.add_argument("-pe", "--paired", help="is the command paired end? than give the second filename", type = str, default = "")

    parser.add_argument("-k", "--kallisto", help="Run kallisto, give the name of the first file, requires -f", type = str)
    parser.add_argument("-f", "--fasta", help="name and path of the fasta file containing the transcriptome/genome", type = str)
    parser.add_argument("-l", "--length", help="required if kallisto is single end, give mean length of the fastq reads", type = float, default = 35.0)
    parser.add_argument("-s", "--std", help="required if kallisto is single end, give standard deviation of the fastq reads", type = float, default = 1.0) 

    parser.add_argument("-ort", "--ortho", help="runs orthofinder, specify working directory", type = str)
    parser.add_argument("-tree", "--tree", help="runs treemaker, specify working directory", type = str)    
    
#    parser.add_argument("-i", "--input", help="unused", type = str)

#    parser.add_argument("-n", "--number", help="unused", type = int)
    
    arguments = parser.parse_args()
    return arguments

def run_commandline(command, out_file = ""):
    """Runs xxx from the commandline and makes xxx
    
    Keyword arguments:
        in_file -- string, the name of the input file
        out_file -- string, the name of the output file, default = out.txt
    
    Returns:
        out_file -- string, the name of the output file
        
    WARNING: if the output file already exists the program will not run, if the output is specified 
    """
    if os.path.exists(out_file) == False:
        #try:
        output = subprocess.check_output(command, shell = True)
        #except subprocess.CalledProcessError as ERROR:
            #stderr.write(ERROR.stderror)
            #exit_program()
    return output

def make_single_cutadapt_command(in_file, out_file, sequence):
    """Makes the cutadapt command for single end sequences to run in the commandline
    Keyword arguments:
        in_file: str, filename in fastq format that needs to be cut
        out_file: str, the desired name of the output file
        sequence: str, the sequence to be cut
    Returns:
        command -- string, the functional command
    """
    return "~/.local/bin/cutadapt -a {} -e 0.1 -O 5 -m 15 -o {} {}".format(sequence, out_file, in_file)
    
def make_double_cutadapt_command(in_file1, in_file2, out_file1, out_file2, sequence):
    """Makes the cutadapt command for paired end sequences to run in the commandline
    Keyword arguments:
        in_file: str, filename1 in fastq format that needs to be cut
        in_file: str, filename2 in fastq format that needs to be cut
        out_file: str, the desired name of the first output file
        out_file: str, the desired name of the second output file
        sequence: str, the sequence to be cut
    Returns:
        command -- string, the functional command
    """
    return "~/.local/bin/cutadapt -a {} -A {} -e 0.1 -O 5 -m 15 -o {} -p {} {}".format(sequence, sequence, out_file1, out_file2, in_file1, in_file2)
    
def make_single_sickle_command(in_file, out_file, paired = "se", quality = 15, length = 20, scoring = "sanger"):
    """Makes the sickle command for single end sequences to run in the commandline
    Keyword arguments:
        in_file: str, the name of the input fastq file
        out_file: str, the desired name of the output file
        quality: int, the minimal quality score to be cut, default = 15 
        length: int, the minimum length to be kept,  default = 20
        scoring: int, the fastq score the file uses, default  = "sanger"
    Returns:
        command -- string, the functional command
    """
    return "sickle se -q {2} -l {3} -f {0} -t {4} -o {1}".format(in_file, out_file, quality, length, scoring)

def make_double_sickle_command(in_file1, in_file2, out_file1, out_file2, paired = "se", quality = 15, length = 20, scoring = "sanger"):
    """Makes the sickle command for double end sequences to run in the commandline
    Keyword arguments:
        in_file1: str, the name of the first input fastq file
        in_file2: str, the name of the  second input fastq file
        out_file1: str, the desired name of the first output file
        out_file2: str, the desired name of the second output file
        quality: int, the minimal quality score to be cut, default = 15 
        length: int, the minimum length to be kept,  default = 20
        scoring: int, the fastq score the file uses, default  = "sanger"
    Returns:
        command -- string, the functional command
    """
    return "sickle pe -s -q {4} -l {5} -f {0} -r {1} -t {6} -o {2} -p {3}".format(in_file1, in_file2, out_file1,out_file2, quality, length, scoring)

def make_fastqc_command(in_file):
    """Makes the fastqc command to run in the commandline
    Keyword arguments:
        in_file, str, the desired fastq filename
    Returns:
        command -- string, the functional command
    """
    return "fastqc {}".format(in_file)

def make_kallisto_index_command(in_file, out_file):
    """Makes the kallisto index command to run in the commandline
    Keyword arguments:
        in_file: str, the name of the input fasta file
        out_file: str, the desired name of the output file
    Returns:
        command -- string, the functional command
    """
    return "/local/prog/kallisto/kallisto index -i {} {}".format(out_file, in_file)

def make_single_kallisto_command(index_file, in_file, out_file, length = 35, std = 1):
    """Makes the kallisto command for single end reads to run in the commandline
    Keyword arguments:
        index_file, str, the filename of the kallisto index
        in_file: str, the name of the input fastq file
        out_file: str, the desired name of the output SAM file
        length: int, the average length of the sequnces, default = 35
        std: int, the standard deviation of length of the sequnces, default = 1
    Returns:
        command -- string, the functional command
    """
    return "/local/prog/kallisto/kallisto quant -i {0} -o ./ --single -l {3} -s {4} --pseudobam {1} > {2}".format(index_file, in_file, out_file, length, std)
    
def make_double_kallisto_command(index_file, in_file1, in_file2, out_file):
    """Makes the kallisto command for double end reads to run in the commandline
    Keyword arguments:
        index_file, str, the filename of the kallisto index
        in_file1: str, the name of the first input fastq file
        in_file2: str, the name of the  second input fastq file
        out_file: str, the desired name of the output SAM file
    Returns:
        command -- string, the functional command
    """
    return "/local/prog/kallisto/kallisto quant -i {} -o ./ --pseudobam {} {} > {}".format(index_file,in_file1, in_file2, out_file)
    
def make_samtools_sort_command(out_file, in_file):
    """Makes the samtools sort command to run in the commandline, to convert sam into bam files and sort them
    Keyword arguments:
        in_file: str, the name of the input SAM file
        out_file: str, the desired name of the output BAM file
    Returns:
        command -- string, the functional command
    """
    return "/local/prog/samtools/samtools sort -l 1 -o {} -O bam {}".format(out_file, in_file)

def make_samtools_index_command(in_file):
    """Makes the samtools index command to run in the commandline
    Keyword arguments:
        in_file: str, the name of the input BAM file
    Returns:
        command -- string, the functional command
    """
    return "/local/prog/samtools/samtools index {}".format(in_file)

def make_orthofinder_normal_command(directory):
    """Makes the orthofinder command to run in the commandline
    Keyword arguments:
        directory: str, the name and path of the directory where the fasta files for orthofinder are located, will double as working directory for orthofinder
    Returns:
        command -- string, the functional command
    """
    return "/local/data/course/project/groups/DNA_Not_Found/progs_nobackup/OrthoFinder-1.1.2/orthofinder -f {} ".format(directory)

def make_orthofinder_trees_command(directory):
    """Makes the trees from MSA from orthofinder command to run in the commandline
    Keyword arguments:
        directory: str, the name and path of the directory where the othofinder output is located, will double as working directory
    Returns:
        command -- string, the functional command
    """
    return "/local/data/course/project/groups/DNA_Not_Found/progs_nobackup/OrthoFinder-1.1.2/trees_from_MSA {}".format(directory)
"""
def timelog(start_time, process):
    time_passed = time.time() - start_time
    with open("timelog.txt", "a") as timefile:
        timefile.write("{}\t{}\n".format(process,time_passed))
    return
"""
def write_log(stamp, log):
    """writes a logfile of the output and time it took to run certain prograns
    keyword arguments:
        stamp: str, the stamp you would like to attach to the output, should be the name of the command run or START
        log: str, the standard output of the program
    """
    with open("programlog.txt", "a") as log_file:
        log_file.write("\n{}\n{}\t{}\n{}\n".format(("#"*50),stamp,time.asctime(),("#"*50)))
        log_file.write("{}".format(log))
    return

def preprocesser(adapter, in_file1, in_file2 ,out_file):
    """preprocesses the fastq files with sicle and optionally cutadapt, writes logs of what was done, removes intermediary files
    keyword arguments:
        adapter: str, the adapter sequnce to be cut with cutadapt, make it an empty string ("") if running cutadapt is not required
        in_file1: str, the name of the (first) fastq file to be preprocessed
        in_file2: str, the name of the second fastq file if the sequence is paired ended, make it an empty string ("") if the sequence is single ended
        out_file: str, the desired name of the output fastq file
    WARNING: will remove all files in the working directory that start with "temp"
    """
    #start_time = time.time()
    #timelog(start_time,"process_start")    
    
    if os.path.exists(in_file1):
        if adapter != "":
            if in_file2 == "":
                cmd = make_single_cutadapt_command(in_file1, "temp1.fq", adapter)
                in_file1 = "temp1.fq"
            else:
                cmd = make_double_cutadapt_command(in_file1, in_file2, "temp1.fq", "temp2.fq", adapter)
                in_file1 = "temp1.fq"
                in_file2 = "temp2.fq"
            output = run_commandline(cmd)
            write_log("CUTADAPT",output)
                
            #timelog(start_time,"cutadapt")
        if in_file2 == "":
            cmd = make_single_sickle_command(in_file1, "preprocess.fastq")
        else:
            cmd = make_double_sickle_command(in_file1, in_file2, "preprocess1.fastq", "preprocess2.fastq")
        output = run_commandline(cmd)
        write_log("SICKLE",output)
        if adapter != "":        
            cmd = "rm temp*.fq"
            run_commandline(cmd)
            #timelog(start_time,"sickle")
    else: 
        stderr.write(arguments.input+"file not found\n")
        exit_program()
    return

def make_bams(paired,fasta,in_file, out_file, length = 35, std = 1, index_file = "temp.index", sam_file = "temp.sam"):
    """mapps the reads in a fastq file to the genome in a fasta file, makes a log of the steps taken
    keyword arguments:
        fasta: str, name of the fasta file containing the genome (or transcriptome)
        in_file: str, the name of the (first) fastq file with the reads to be mapped
        paired: str, the name of the second fastq file if the sequence is paired ended, make it an empty string ("") if the sequence is single ended
        out_file: str, the desired name of the output bam file
        length: int, the average length of the sequnces, default = 35, for kallisto if the reads are single ended
        std: int, the standard deviation of length of the sequnces, default = 1, for kallisto if the reads are single ended
        index_file: str, the name for the intermidate index file, default = "temp.index"
        sam_file: str, the name for the intermidate sam file, default = "temp.sam"
    WARNING: will remove all files in the working directory that start with "temp"
    """
    cmd = make_kallisto_index_command(fasta, index_file)
    output = run_commandline(cmd)
    write_log("KALLISTO INDEX",output)
    if paired == "":
        cmd = make_single_kallisto_command(index_file, in_file, sam_file, length, std)
    else:
        cmd = make_double_kallisto_command(index_file, in_file, paired, sam_file)
    output = run_commandline(cmd)
    write_log("KALLISTO RUN",output)
    
    cmd = make_samtools_sort_command(out_file, sam_file)
    output = run_commandline(cmd)
    write_log("SAMTOOLS SORT",output)
    
    cmd = make_samtools_index_command(out_file)
    output = run_commandline(cmd)
    write_log("SAMTOOLS INDEX",output)
    
    cmd = "rm temp*"
    run_commandline(cmd)
    
#Main Loop
if __name__ == "__main__":

    #preprocesser()
    arguments = get_command_line_arguments()

    write_log("START", "#"*50)
    
    if arguments.fastqc != None:
        cmd = make_fastqc_command(arguments.fastqc)
        output = run_commandline(cmd)
        write_log("FASTQC", output)
        
    if arguments.preprocess != None:
        preprocesser(arguments.cutadapt, arguments.preprocess, arguments.paired, arguments.output)
    
    if arguments.kallisto != None:
        make_bams(arguments.paired,arguments.fasta,arguments.kallisto, arguments.output, arguments.length, arguments.std)

    if arguments.ortho != None:
        cmd = make_orthofinder_normal_command(arguments.ortho)
        output = run_commandline(cmd)
        write_log("KALLISTO RUN",output)
    
    if arguments.tree != None:
        cmd = make_orthofinder_trees_command(arguments.tre)
        output = run_commandline(cmd)
        write_log("KALLISTO RUN",output)
    
"""
~/.local/bin/cutadapt -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGATC -A GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGATC -e 0.1 -O 5 -m 15 -o forw_cut.fastq -p rev_cut.fastq /local/data/course/project/groups/DNA_Not_Found/SRP006011/SRR125625_1.fastq /local/data/course/project/groups/DNA_Not_Found/SRP006011/SRR125625_2.fastq

~/.local/bin/cutadapt -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGATC -A GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGATC -e 0.1 -O 5 -m 15 -o forw_cut.fastq -p rev_cut.fastq /local/data/course/project/groups/DNA_Not_Found/SRP006011/SRR125625_1.fastq /local/data/course/project/groups/DNA_Not_Found/SRP006011/SRR125625_2.fastq
"""
    
    
    
    

"""Changelog

V0.1:
First Version of the program!
made an argparser for optional arguments
made function run_commandline and make_command (currently more template functions)

V0.2
Wrote time logger
Program makes cutadapt and sickle commands and SHOULD be able to run them
made the loop that will eventually be used to run the commands for multiple files
moved the changelog to the end of the program

V0.3
rewrote everything completely
now runs all preprocesses and may run kallisto
keeps logs on what was done

V0.4
now runs orthofinder
bugfix in kallisto making indexes
bugfix in another kallisto thing
now with docstings!
"""
    
    
    
    
    
    
    
    
    
