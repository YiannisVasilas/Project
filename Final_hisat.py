#!/usr/bin/env python
"""
Ioannis Vasilas, Jasper van der Rijst, Marijn Schipper 

Script for Hisat2 
indexing input: ref_genome output : indexed genome
  
mapping input indexed genome trimmed reads  
"""


from sys import argv
import subprocess
import os.path
import samtools as st
import htseq

def group_replicates(tr_setting, treatments, timepoints):
    '''Yield the a list with the replicates per treatment per timepoint
    
    treatments: list -- list with the names of the treatments
    timepoints: list -- list with the timepoints
    '''
    replicate_list = []
    for treatment in treatments:
        for time in timepoints:
            if tr_setting:
                paired_files1 = ("fastq_trimmed_{}_1_{}_1.fq_pairs_R1.fastq"\
                    .format(treatment,time),\
                    "fastq_trimmed_{}_1_{}_2.fq_pairs_R2.fastq".format(treatment,time))
                replicate_list.append(paired_files1)
                paired_files2 = ("fastq_trimmed_{}_2_{}_1.fq_pairs_R1.fastq"\
                    .format(treatment,time),\
                    "fastq_trimmed_{}_2_{}_2.fq_pairs_R2.fastq".format(treatment,time))
                replicate_list.append(paired_files2)
            else:
                paired_files1 = ("{}_1_{}_1.fastq.gz"\
                .format(treatment,time),\
                "{}_1_{}_2.fastq.gz".format(treatment,time))
                replicate_list.append(paired_files1)
                paired_files2 = ("{}_2_{}_1.fastq.gz"\
                .format(treatment,time),\
                "{}_2_{}_2.fastq.gz".format(treatment,time))
                replicate_list.append(paired_files2)
    return replicate_list

def build_hisat2_index(ref_file, ref_index):
    '''Build the hisat2 index of the reference genome, return the index name
    
    ref_file: string -- path to the reference genome in fasta format
    ref_index: string -- name of the index
    '''
    if not os.path.exists(ref_index+".1.ht2"):
        cmd = 'hisat2-build %s %s' %(ref_file, ref_index)
        subprocess.check_call(cmd, shell=True)
    return ref_index

def run_hisat2(path_to_reads, pair1, pair2, ref_index, tr_setting):
    """Run hisat2 program on fastq files
    
    path_to_reads: string -- relative path to the trimmed fastq files
    pair1: list -- list of paired filenames of trimmed RNA fastq files
    pair2: list -- list of paired filenames of trimmed RNA fastq files
    ref_index: string -- name of the hisat2 index of the reference genome
    """
    m1 = "{}{}".format(path_to_reads, pair1)
    m2 = "{}{}".format(path_to_reads, pair2)
    nameparts = pair1.split("_")
    if '.gz' not in m2:
        output_name = "{}_{}_{}.sam".format(nameparts[2], nameparts[4],\
         nameparts[3])
    else:
        output_name = "{}_{}_{}.sam".format(nameparts[0], nameparts[2],\
         nameparts[1]).split('.')[0]       
    if not os.path.exists(output_name):
        print('running hisat')
        cmd = 'hisat2 -x {} -1 {} -2 {} -S {}'\
          .format(ref_index, m1, m2, output_name)
        subprocess.check_call(cmd, shell=True)
    else:
        print('skipping hisat ' + output_name) 
    return output_name
    
if __name__ == "__main__":
    treatments = ["Cf", "Mock1", "Mock2", "Sc"]
    timepoints = ["12", "24", "48"]
    for pair in group_replicates(treatments, timepoints):
        run_hisat2('.', pair[0], pair[1], "index")
    
