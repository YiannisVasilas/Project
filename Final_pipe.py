#!/usr/bin/env python
"""
Marijn Schipper, Ioannis Vasilas, Jasper van der Rijst

Script to run all the scripts
argv[1] = path to folder with fastq files 
"""


from sys import argv
import subprocess
import os
import Fastq_module_2 as fq
import hisat as h2
import samtools as st
import htseq 
import trimfix as tfix
import gff_parser as gp


def pipe(input_path, out_trim, out_stats, out_bias, ref_genome_path, \
 path_to_gff_file, tr_setting, version_name):
    """Returns read counts from rna-seq fastq files, mapped to ref genome
     input: input_path - path to gzipped fastq files(str)
            out_trim - folder name for trim output (str)(def = out_trim)
            out_bias - folder name for bias output (str)(def = out_bias)
            out_stats - folder name for out stats (str)(def = out_stats)
            ref_genome - path to reference genome (str)
            path_to_gff_file - path to gene models .gff file (str)
            tr_set - determines whether to trim (T) or not (F) (Boolean)
        output: tab delimited txt files (htseq std output, see htseq)
    """
    #Bool switch for trimming
    if tr_setting:
        files = fq.Fastqrunner(input_path, out_trim, out_stats, out_bias,\
     tr_setting)
        #delete intermediate data
        if os.path.exists(out_bias):
            command = "rm -r {}".format(out_bias)
            run = subprocess.check_output(command,shell = True)
        files = "./{}/".format(out_trim)
    else:
        files = './fastq_links/'
    #sort different filenames into correlating group lists
    replicate_list =  h2.group_replicates(tr_setting, treatments = \
    ["Cf", "Mock1", "Mock2", "Sc"], timepoints = ["12", "24", "48"])
    #build index with Hisat 2 and ref genome
    index = h2.build_hisat2_index(ref_genome_path, version_name)
    #map and count reads for each trimmed fastq pair in input path
    for (pair1, pair2) in replicate_list:
            sam_file = h2.run_hisat2(files, pair1, pair2, index, tr_setting)
            bam_file = st.run_samtools(sam_file, out_sam)
            path_to_bam = "./{}/{}".format(out_sam, bam_file)
            htseq.run_htseq_count(path_to_bam, path_to_gff_file, out_htseq)
            #remove redundant bam files
            st.remove_file(path_to_bam)
    gp.parse_htseq_count_output(out_gcount, path_to_gff_file, "./{}/".format(out_htseq))
    return

if __name__=="__main__":
    #path to gzipped fastqfiles
    input_path = argv[1]
    tr_setting = argv[2]
    #set necessary variables
    ref_genome_path = argv[3]
    path_to_gff_file = argv[4]
    version_name = argv[5]
    if tr_setting:
        tr_add = '_trimmed'
    else:
        tr_add = '_untrimmed'
    
    #set output paths
    out_bias = 'bias'
    out_trim = 'trim'
    out_stats = 'stats_{}'.format(tr_add)
    out_sam = 'out_sam_{}{}'.format(version_name, tr_add)
    out_gcount = '{}{}'.format(version_name, tr_add)
    out_htseq = 'out_htseq_{}{}'.format(version_name, tr_add)
    
    #run pipe
    pipe(input_path, out_trim, out_stats, out_bias, ref_genome_path,\
     path_to_gff_file, tr_setting, version_name)
