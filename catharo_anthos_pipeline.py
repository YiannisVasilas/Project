#!/usr/bin/env python
"""

"""

from __future__ import print_function, division
from sys import argv
import sys
if sys.version_info[0]==2:
    input = raw_input
import argparse
import subprocess as sp

class ArgumentParserError(Exception):
    """Creates a new kind of error.
    """
    pass

class ThrowingArgumentParser(argparse.ArgumentParser):
     """Redefines the argparse error, to be able to catch it with try/except.
     """
     def error(self, message):
        raise ArgumentParserError(message)

def argument_parser():
    """Tries to execute the script with command line arguments.
    """
    input_file = ""
    try:
        # Creates an instance of argparse.
        argparser = ThrowingArgumentParser(prog=argv[0])
        # Adds the argument.
        argparser.add_argument("rnaseq_control", \
        help="The rnaseq_control filename, including the extension.", type=str)
        argparser.add_argument("rnaseq_meja", \
        help="The rnaseq_meja filename, including the extension.", type=str)
        argparser.add_argument("ref_genome", \
        help="The ref_genome filename, including the extension.", type=str)
        argparser.add_argument("prefix", \
        help="A prefix used for some intermediate files.", type=str)
        argparser.add_argument("proc", \
        help="The maximum number of processors to use.", type=int)
        # Parses the arguments given in the shell.
        args = argparser.parse_args()
        # Retrieves the filenames.
        rnaseq_control = args.rnaseq_control
        rnaseq_meja = args.rnaseq_meja
        ref_genome = args.ref_genome
        prefix = args.prefix
        proc = args.proc
        if not all([rnaseq_control,rnaseq_meja,ref_genome,prefix] == str):
            raise ArgumentParserError
        if not proc==int:
            raise ArgumentParserError
    except ArgumentParserError:
        print("Please make sure you provide all arguments properly.")
    return rnaseq_control, rnaseq_meja, ref_genome, prefix, str(proc)

def command_line(command):
    """Executes the command and checks if it produces an error.

    Input:
    - command: string (command to run)
    Output:
    - Any text created by the command, including errors.
    - Stops running if an error occurs.
    """
    try:
        check_output = sp.check_output(command, stderr=sp.STDOUT, shell=True)
    except sp.CalledProcessError as e:
        print("The tool created an error:"+str(e))
        sys.exit()
    return check_output

def run_commands(rnaseq_control, rnaseq_meja, ref_genome, prefix, proc):
    """Executes all commands for the Catharo anthos pipeline.

    Input:
    - rnaseq_control: string (RNAseq FASTQ file)
    - rnaseq_meja: string (RNAseq FASTQ file)
    - ref_genome: string (genome FASTQ file)
    - prefix: string (prefix for intermediate files)
    - proc: string (max. number of processors to use)
    Output:
    - Creates seperate files for the forward and reverse reads (_1 and _2).
    - Creates .sam, .bam, .gtf,
    - Writes command output to screen 
    """
    hisat_control = rnaseq_control+"_hisat"
    hisat_meja = rnaseq_meja+"_hisat"
    prefix_control = rnaseq_control.split(".")[0]
    prefix_meja = rnaseq_meja.split(".")[0]
    bam_control = prefix_control+".bam"
    bam_meja = prefix_meja+".bam"
    gtf_control = prefix_control+"_counted.gtf"
    gtf_meja = prefix_meja+"_counted.gtf"
    command = "fq_splitter_catharantos.pl {}".format(rnaseq_control)
    print(command_line(command))
    command = "fq_splitter_catharantos.pl {}".format(rnaseq_meja)
    print(command_line(command))
    command = "hisat2-build {} {}".format(ref_genome, prefix)
    print(command_line(command))
    command = "hisat2 --time -x {} -1 {} -2 {} -S {}.sam".format(prefix, \
        rnaseq_control+"_1",rnaseq_control+"_2",hisat_control)
    print(command_line(command))
    command = "hisat2 --time -x {} -1 {} -2 {} -S {}.sam".format(prefix, \
        rnaseq_meja+"_1",rnaseq_meja+"_2",hisat_meja)
    print(command_line(command))
    command = "samtools view -S -b {} > {}.temp".format(hisat_control, \
                                                        hisat_control)
    print(command_line(command))
    command = "samtools sort -@ {} {}.temp {}".format(proc, hisat_control, \
                                                      prefix_control)
    print(command_line(command))
    command = "samtools view -S -b {} > {}.temp".format(hisat_meja, \
                                                        hisat_meja)
    print(command_line(command))
    command = "samtools sort -@ {} {}.temp {}".format(proc, hisat_meja, \
                                                      prefix_meja)
    print(command_line(command))
    command = "rm *.temp"
    print(command_line(command))
    command = "stringtie -p {} -G {} -o {} {}".format(proc, ref_gff, \
                                                gtf_control, bam_control)
    print(command_line(command))
    command = "stringtie -p {} -G {} -o {} {}".format(proc, ref_gff, \
                                                      gtf_meja, bam_meja)
    print(command_line(command))
    command = "gffcompare –r {} –G –o {} {}".format(ref_gff, prefix_control, \
                                                      gtf_control)
    print(command_line(command))
    command = "gffcompare –r {} –G –o {} {}".format(ref_gff, prefix_meja, \
                                                      gtf_meja)
    print(command_line(command))
    command = "stringtie -p {} {} -G {} –B -e -o s {}".format(proc, \
                                            bam_control, ref_gff, gtf_control)
    print(command_line(command))
    command = "stringtie -p {} {} -G {} –B -e -o s {}".format(proc, bam_meja, \
                                                      ref_gff, gtf_meja)
    print(command_line(command))
    command = "python prepDE_edit.py {} {}".format(gtf_control, gtf_meja)
    print(command_line(command))

if __name__=="__main__":
    # Default files to use with the pipeline.
    rnaseq_control = "SRR646604.fastq"
    rnaseq_meja = "SRR646572.fastq"
    ref_genome = "cro_scaffolds.min_200bp.fasta"
    ref_gff = "cro_std_maker_anno.final.gff3"
    prefix = "cr_index"
    proc = str(12)
    # Checks if other files are passed as arguments.
    try:
        rnaseq_control, rnaseq_meja, ref_genome, prefix, proc = argument_parser()
    except Exception:
        print("Faulty or no additional arguments were entered. Using the \
default files as input.")
    # Runs the commands and scripts up until DESeq.
    run_commands(rnaseq_control, rnaseq_meja, ref_genome, prefix, proc)
    print("The pipeline has finished. Please continue with the DESeq script.")
    
