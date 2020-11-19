#!/usr/bin/env python
"""
Marijn Schipper

Script to run all the scripts
argv[1] = path to folder with fastq files 
"""


from sys import argv
import subprocess
import os

def symbolic_link_runner():
    """runs python script to make soft link folder in command line
    
    returns:
        folder containing soft links to path set in symbolic_link_fastq.py
    """
    if not os.path.exists('./fastq_links/'):
        cmd = "python3 symbolic_link_fastq.py {}".format(index_seq_csvfile)
        subprocess.check_call(cmd, shell = True)
    return
    
def run_pipe_diff_settings(versions):
    """runs scripts in linux command line, piping fastqfiles -> genecount
    
        versions = [[genome.gff, gene models.gff, 'genome version name'],..]
        returns:
            output in output directories specified in pipe.py
"""
    for version in versions:
        for boolean in (True, False):
            if not os.path.exists('./{}/'.format(version[2])):
                cmd = "mkdir {}".format(version[2])
                subprocess.check_call(cmd, shell = True)
            cmd = 'python3 pipe.py ./fastq_links/ {} {} {} {} '.format(\
            boolean, version[0], version[1], version[2])
            subprocess.check_call(cmd, shell = True)
            cmd = 'mv out* {}'.format(version[2])
            subprocess.check_call(cmd, shell = True)
    return
            
            

if __name__=="__main__":
    #path to gzipped fastqfiles
    index_seq_csvfile = argv[1]
    
    #make symbolic links to fastq files
    symbolic_link_runner()
    
    #create versions
    genome_path_v2 = argv[2]
    gene_models_path_v2 = argv[3]
    genome_path_v4 = argv[4]
    gene_models_path_v4 = argv[5]
    versions = [(genome_path_v2, gene_models_path_v2, '2.4_genome'),\
     (genome_path_v4, gene_models_path_v4, '4.0_genome')]
    
    #run pipeline
    run_pipe_diff_settings(versions)
    
    
    
