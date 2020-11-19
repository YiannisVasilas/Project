##!/usr/bin/env python
"""
Created on Fri Dec  2 11:14:29 2016

"""

from __future__ import division
from sys import argv
import subprocess
import logging
import os, sys



def create_index(name, genome_fasta): #removed name
    """ Creates index files for a genome fasta file using Hisat2-build 
    Checks if the index files for said fasta genome already exists, if not, builds them 
    
    Keyword arguments: 
        name -- general name that was given to the job in the original command line 
                All files will be saved under variations of this name 
        genome_fasta -- Filename (and possibly location) of the fasta fil eof the genome
        that is to be indexed 
    Returns: 
        Eight output files as made by Hisat2-build, deposited to the map 'refgenomeindices'
        in the PacWay directory. The name of these output files will be the name of the 
        genome file with '_index'behind it plus program specific adittions from Hisat2 
    """
#    name = genome_fasta.split('/')[-1]
    path_folder = '/local/data/course/project/groups/PacWay/refgenomeindices/{}_index'.format(name)
   
    cmd = "hisat2-build {} {}".format(genome_fasta, path_folder)
  
    if os.path.exists("{}.1.ht2".format(path_folder)):
        logging.info('{} is already present'.format(path_folder))
        print 'index file is present'
       
    else: 
        subprocess.check_call(cmd, shell=True)
        logging.info('running: {}'.format(cmd))
    
    return path_folder #added 
    
    
    
def assemble(name, index, RNAseq):
    """
    Maps RNA-seq reads to (the index of) a fasta genome using hisat2
    
    Keyword arguments:
        name -- general name that was given to the job in the original command line 
                All files will be saved under variations of this name 
        index -- Location and filename of the index files as made by function 'create index'
        RNAseq -- ??          
    """
    cmd = 'hisat2 -x {} -U RNA_data/{} -S refgenomemapping/{}.sam'.format(index, RNAseq, name)
    
    if os.path.exists("refgenomemapping/{}.sam".format(name)):
        logging.info('{}.sam is already present'.format(name))
        print 'sam file is present'
    
    else: 
        subprocess.check_call(cmd, shell=True)
        logging.info('running: {}'.format(cmd))
    
    return 'refgenomemapping/{}.sam'.format(name)
    
    
    
def samtobam(name, SAM): 
    """
    Makes a BAM (binary SAM) file out of a SAM file using samtools 
    
    Keyword arguments:
        name -- general name that was given to the job in the original command line 
                All files will be saved under variations of this name 
        SAM -- Name and location of the SAM file
    Returns: 
        A bam file in the refgenomemapping directory
    """
    cmd = 'samtools view -S -b {} > refgenomemapping/{}.bam'.format(SAM, name)
    cmd2 = 'samtools sort refgenomemapping/{}.bam refgenomemapping/{}.prefix'.format(name, name)    
    
    if os.path.exists("refgenomemapping/{}.bam".format(name)):
        logging.info('{}.bam is already present'.format(name))
        print 'bam file is present'
        if os.path.exists("refgenomemapping/{}.prefix.bam".format(name)):
            logging.info('{}.prefix.bam is already present'.format(name))
            print 'sorted bam file is present'
        else:
            subprocess.check_call(cmd2, shell=True)
            logging.info('running: {}'.format(cmd2))
    
    else: 
        subprocess.check_call(cmd, shell=True)
        logging.info('running: {}'.format(cmd))
        subprocess.check_call(cmd2, shell=True)
        logging.info('running: {}'.format(cmd2))
    
    return 'refgenomemapping/{}.bam'.format(name), 'refgenomemapping/{}.prefix.bam'.format(name)
  
  

def stringtie(name, sortedbam):
    """
    Runs stringtie on a sorted bam file, creating a gtf file with annotation
    from a gff3 file. 

    Keyword arguments:
        name -- general name given to the job, all files are safed under this
                name
        sortedbam -- a filename of the *.prefix.bam file to be used by Stringtie 
    """
    
    cmd = 'stringtie {} -G Catharanthus_roseus/cro_std_maker_anno.final.gff3 > refgenomemapping/{}.gtf'.format(sortedbam, name)
    if os.path.exists("refgenomemapping/{}.gtf".format(name)):
        logging.info('{}.gtf is already present'.format(name))
        print '{}.gtf file is present'.format(name)
    
    else: 
        subprocess.check_call(cmd, shell=True)
        logging.info('running: {}'.format(cmd))
    
    return 'refgenomemapping/{}.gtf'.format(name)


    
if __name__=='__main__':
 
    logging.basicConfig(filename='Loes/guided_assembly.log', level=logging.INFO)
    
#    folder = os.listdir('/local/data/course/project/groups/PacWay/RNA_data')
#    RNAfiles = []
#    for file in folder: 
#        RNAfiles += [file]
#    print RNAfiles
    
    ## all non-paired files     
#    RNAfiles =['SRR122246.fastq', 'SRR122244.fastq', 'SRR122258.fastq', 'SRR122240.fastq', 'SRR122253.fastq', 'SRR122259.fastq', 'SRR122251.fastq', 'SRR122250.fastq', 'SRR122249.fastq', 'SRR122257.fastq', 'SRR122248.fastq', 'SRR122255.fastq', 'SRR122261.fastq', 'SRR122239.fastq']
  
#    for RNA in RNAfiles: 
    RNA_seq = argv[1] 
    number = RNA_seq.split('/')[-1].split('.')[0].replace('SRR1222','')
    name = 'RNA{}'.format(number)
    index = 'refgenomeindices/cro_scaffolds.min_1000bp_index'
#       genome_fasta =    
    
#       index = create_index(name, genome_fasta)
    SAM = assemble(name, index, RNA_seq)
    BAM, sortedbam = samtobam(name, SAM)
    GTF = stringtie(name, sortedbam)
    
