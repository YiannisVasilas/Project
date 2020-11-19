#!/bin/usr/env python
"""
Created on Mon Dec  5 09:15:07 2016
author: team members of PyOmics project group
Pipeline script for doing a Differential Expression analysis
"""

from __future__ import division
from sys import argv
import argparse
import subprocess
import os.path
import re



def get_arguments():
    """ Gets the arguments from the command line

    Keyword arguments:
        no parameters for this function
    """
    parser = argparse.ArgumentParser(description="Pipeline script for DE \
    analysis")
    parser.add_argument("-raw", "--raw_folder", help="folder containing \
    the raw fastq files", type=str, default="raw_samples/")
    parser.add_argument("-SEr", "--se_reads", nargs="*", help=".fastq file \
    containing single end RNA seq reads", type=str, required=False)
    parser.add_argument("-PEr", "--pe_reads", nargs="*", help=".fastq file \
    containing paired end RNA seq reads", type=str, required=False)
    parser.add_argument("-trim", "--trim_folder", help="folder containing \
    the trimmed fastq files", type=str, default="trimmed_samples/")
    parser.add_argument("-clip", "--illumina_clip", help="file containing \
    the sequences that need to be trimmed off the target sequence", 
    type=str, required=False)
    parser.add_argument("-index", "--index_filename", help="file containing \
    indices for the reference genome", type=str, required=False)
    parser.add_argument("-splice", "--splicesites", help="file containing \
    known splice sites", type=str, required=False)
    arguments = parser.parse_args()
    return arguments


def run_trimmomatic(fastq_filename, raw_folder, trim_folder, clipper_file, 
                    single = True, seed_mm=2, palin_th=30, simple_th=10, 
                    lead_val=3, trail_val=3, win_size=4, req_qual=30):
    """ Returns the trimmed reads of single or paired end data to a file
    
    Keyword arguments:
        fastq_filename -- file with  RNA-seq reads in fastq format (.fastq) 
        trimm_outfile -- string, filename for the outputfile
        clipper_file -- string, filename with adaptersequences that need to 
        be trimmed of the raw fastq reads
        seed_mm -- integer, seedMismatches: specifies maximum mismatch count
        palin_th -- integer, palindromeClipThreshold: specifies how accurate
        the match between two adapter ligated reads must be.
        simple_th -- integer, simpleClipThreshold: specifies how accurate the 
        match between any adapter must be against a read
        LEADING -- integer, specifies minimum quality required to keep a base.
        TRAILING -- integer, specifies minimum quality required to keep a base.
        win_size -- integer, windowSize:specifies the number of bases to 
        average across
        req_qual -- integer, requiredQuality: specifies the average quality 
        required
        single -- boolean, 
    Returns:
        trimm_outfile -- .fastq, name of trimmed RNA-seq reads in .fastq format
    """
    if single == True:
        # make for loop to go through every dataset
        
        trimm_outfile = "%st_%s"%(trim_folder, fastq_filename)
        # checks whether the file already exists, if not it runs the tool
        if os.path.exists(trimm_outfile):
            print "The file already exists"
            return [trimm_outfile]
        else:
            cmdSE = "TrimmomaticSE %s%s %s ILLUMINACLIP:%s%s:%d:%d:%d LEADING:%d \
TRAILING:%d SLIDINGWINDOW:%d:%d" %(raw_folder, fastq_filename, trimm_outfile, 
                                   raw_folder, clipper_file, seed_mm, palin_th, 
                                   simple_th, lead_val, trail_val, 
                                   win_size, req_qual)

            output_check = subprocess.check_output(cmdSE, shell=True)
            call_check = subprocess.check_call(cmdSE, shell=True)
#            return call_check #must be 0
            print cmdSE
            return [trimm_outfile]
    
    elif single == False:
        # make for loop to go through every dataset
        
        trimm_outfile0 = "t_%s"%(fastq_filename[0])
        trimm_outfile1 = "t_%s"%(fastq_filename[1])
        # checks whether the file already exists, if not it runs the tool
        if os.path.exists("%spaired_%s"%(trim_folder, trimm_outfile0)) \
and os.path.exists("%spaired_%s"%(trim_folder, trimm_outfile1)):
            print "The files already exists"
            return [["%spaired_%s"%(trim_folder, trimm_outfile0), \
"%spaired_%s"%(trim_folder, trimm_outfile1)]]
        else:
            cmdPE = "TrimmomaticPE %s%s %s%s %spaired_%s %sunpaired_%s \
%spaired_%s %sunpaired_%s ILLUMINACLIP:%s%s:%d:%d:%d LEADING:%d \
TRAILING:%d SLIDINGWINDOW:%d:%d" %(raw_folder, fastq_filename[0], 
                                   raw_folder, fastq_filename[1], 
                                   trim_folder, trimm_outfile0, 
                                   trim_folder, trimm_outfile0, 
                                   trim_folder, trimm_outfile1, 
                                   trim_folder, trimm_outfile1, 
                                   raw_folder, clipper_file, seed_mm, palin_th, 
                                   simple_th, lead_val, trail_val, 
                                   win_size, req_qual)

            output_check = subprocess.check_output(cmdPE, shell=True)
            call_check = subprocess.check_call(cmdPE, shell=True)
#            return call_check #must be 0
            print cmdPE
            return [["%spaired_%s"%(trim_folder, trimm_outfile0), \
"%spaired_%s"%(trim_folder, trimm_outfile1)]]


#==============================================================================
#         cmdPE = "TrimmomaticPE %s %s ILLUMINACLIP:%s:%d:%d:%d LEADING:%d \
#         TRAILING:%d SLIDINGWINDOW:%d:%d" %(fastq_filename, trimm_outfile, 
#                                            clipper_file, seed_mm, palin_th, 
#                                            simple_th, lead_val, trail_val, 
#                                            win_size, req_qual)
#                            
#         output_check = subprocess.check_output(cmdPE, shell=True)
#         call_check = subprocess.check_call(cmdPE, shell=True)
#         return call_check #must be 0
#==============================================================================

def run_hisat2(index_filename, splicesites, trimmed_input, single = True):
    """
    Returns RNA-seq reads mapped to the reference genome to a file
    
    Keyword arguments:
        index_filename: string, a file with the indices for the reference genome
        splicesites: string, a file of list of know splice sites
        trimmed_input: Files with the unpaired reads
    Command structure:
    /local/prog/hisat2/hisay2-align-s
    --wrapper basic-0
    --dta
    -x genome/cro_genome.dna
    --known-splicesite-infile genome/ssFile.table
    -S stringtie_sams/samplex.sam
    -U trimmed_SSR1820149.fastq /
    -1 trimmed_SRR1820326_1.fq -2 trimmed_samples/SRR1820326_2.fq 
    
    """
    if single == True:
        single_runcount = 0
        hisat_outfiles = []
        for i in range(len(trimmed_input)):
            single_runcount += 1
            trimmed_singles = trimmed_input[i]
            
            hisat_outfile = "stringtie_sams/sample%d.sam"%(single_runcount)
            # checking if file already exists
            if os.path.exists(hisat_outfile):
                print "The file already exists"
                return [hisat_outfile]
            else:
                cmdSE = "/local/prog/hisat2/hisat2-align-s --wrapper basic-0 \
--dta -x %s --known-splicesite-infile %s -S %s -U %s"%(index_filename, 
                                                       splicesites, hisat_outfile, 
                                                       trimmed_singles)
                print cmdSE
                hisat_outfiles.append(hisat_outfile)
                output_check = subprocess.check_output(cmdSE, shell=True)
                call_check = subprocess.check_call(cmdSE, shell=True)
#                return call_check #must be 0
        return hisat_outfiles
        
        
    elif single == False:
        pairs_runcount = 0
        hisat_outfiles = []
        for i in range(len(trimmed_input)):
            pairs_runcount += 1
            trimmed_pairs = trimmed_input[i]
            
            hisat_outfile = "stringtie_sams/samplep%d.sam"%(pairs_runcount)
            # checking if file already exists
            if os.path.exists(hisat_outfile):
                print "The file already exists"
                return [hisat_outfile]
            else:
                cmdSE = "/local/prog/hisat2/hisat2-align-s --wrapper basic-0 \
--dta -x %s --known-splicesite-infile %s -S %s -1 %s -2 %s"%(index_filename,
                                                             splicesites, 
                                                             hisat_outfile,
                                                             trimmed_pairs[0], 
                                                             trimmed_pairs[1])
                print cmdSE
                hisat_outfiles.append(hisat_outfile)
                output_check = subprocess.check_output(cmdSE, shell=True)
                call_check = subprocess.check_call(cmdSE, shell=True)
#                return call_check #must be 0
        return hisat_outfiles
                    
#==============================================================================
#         cmdPE = "hisat2 -x %s --known-splicesites-infile %s -U %s \
#       --dta-cfufflinks -S %s"%(index_filename, splicesites, trimmed_input, 
#                                 hisat_outfile)
#                                 
#        output_check = subprocess.check_output(cmdSE, shell=True)
#        call_check = subprocess.check_call(cmdSE, shell=True)
#        return call_check #must be 0
#==============================================================================

def run_samtools(directory = "stringtie_sams/"):
    cmd = "for sam in %s*.sam ; do samtools view -bS $sam | samtools sort \
- $sam.bam ; done"%(directory)
    print cmd
    output_check = subprocess.check_output(cmd, shell=True)
    call_check = subprocess.check_call(cmd, shell=True)
    return call_check #must be 0

def run_stringtie(directory = "stringtie_sams/"):
    cmd = "for bam in %s*.bam ; do /local/prog/stringtie/stringtie $bam -G genome/annots_with_introns.gff3 \
-o $bam.gtf ; done"%(directory)
    print cmd
    output_check = subprocess.check_output(cmd, shell=True)
    call_check = subprocess.check_call(cmd, shell=True)
    return call_check #must be 0
    
def run_merge_stringtie(directory = "stringtie_sams/"):
    cmd = "/local/prog/stringtie/stringtie --merge -G genome/annots_with_introns.gff3 -o merged.gtf %s*.gtf"\
%(directory)
    print cmd
    output_check = subprocess.check_output(cmd, shell=True)
    call_check = subprocess.check_call(cmd, shell=True)
    return call_check #must be 0
    
def run_ballgown(directory = "stringtie_sams/"):
    cmd = "for bam in %s*.bam ; do /local/prog/stringtie/stringtie -e -B -G merged.gtf \
-o ballgown/bg_${bam%%.*}/bg_${bam%%.*}.gtf $bam ; done"%(directory)
    print cmd
    output_check = subprocess.check_output(cmd, shell=True)
    call_check = subprocess.check_call(cmd, shell=True)
    return call_check #must be 0
    
if __name__ == "__main__":
    #Get input file names from command line
    arguments = get_arguments() 
#    print arguments

    # running the trimmomatic tool
    # only runs when se_reads contain something
    if arguments.se_reads:
        trimmed_singles = []
        for i in range(len(arguments.se_reads)):
            #adds the trimmed files to the list of trimmed_singles
            trimmed_singles += run_trimmomatic(arguments.se_reads[i],
                                               arguments.raw_folder, 
                                               arguments.trim_folder,
                                               arguments.illumina_clip,
                                               single=True)
#        print trimmed_singles
    
    if arguments.pe_reads:
        trimmed_pairs = []
        for i in range(0,len(arguments.pe_reads),2):
            #adds the trimmed files to the list of trimmed+pairs
            trimmed_pairs += run_trimmomatic([arguments.pe_reads[i], 
                                              arguments.pe_reads[i+1]], 
                                              arguments.raw_folder, 
                                              arguments.trim_folder,
                                              arguments.illumina_clip, 
                                              single=False)
#        print trimmed_pairs
    
    try:
        trimmed_singles
        singles = True
    except NameError:
        print "No trimmed singles found"
        singles = False
    if singles == True:
        single_sams = []
        single_sams = run_hisat2(arguments.index_filename, arguments.splicesites, 
                   trimmed_singles, single = True)
#        print single_sams
    
    try:
        trimmed_pairs
        pairs = True
    except NameError:
        print "No trimmed pairs found"
        pairs = False
    if pairs == True:
        paired_sams = []
        paired_sams = run_hisat2(arguments.index_filename, arguments.splicesites, 
                   trimmed_pairs, single = False)
#        print paired_sams
    
    run_samtools()
    run_stringtie()
    run_merge_stringtie()
    run_ballgown()
    
#==============================================================================
#         
#     if arguments.pe_reads==True:
#         for i in range(len(arguments.pe_reads)):
#             run_trimmomatic(arguments.pe_reads[i], arguments.illumina_clip,
#                             single=False)
#                             
#     if arguments.se_reads==True:
#         for i in range(len(arguments.se_reads)):
#             run_trimmomatic(arguments.se_reads[i], arguments.illumina_clip, 
#                             single=True)
#                             
#                     
#         
#     # run trimmomatic tool from command line
#     run_trimmomatic(arguments.se_reads[0], arguments.illumina_clip)
#==============================================================================
