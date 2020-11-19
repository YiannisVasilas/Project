#!/usr/bin/env python3

"""
Author: Yan Wang
trimming off low quality part of raw reads
"""

from sys import argv
import subprocess
import os.path



def run_trimmomatic(inp_f1, inp_f2, phred = 33,\
                    adapter = '/local/prog/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10',\
                    slidwin = '4:10', minlen = 20):
    '''return the name of result files
    
    inp_f1, inp_f2: string, the name of two files from paired end reads, 
    threads: number, the threadshold score
    phred: number, the scoring system used in fastq file
    adapter: string, the file of adapter
    slidwin: string, the parameter for sliding window trimming
    minlen: number, removes reads that fall below the specified minimal length
    '''
    if not (os.path.exists(inp_f1) and os.path.exists(inp_f1)):
        raise ValueError('input file for trimming does not exist')
    else:
        outp_f_1u = inp_f1.replace('.fastq','_unpair.fastq')
        outp_f_1p = inp_f1.replace('.fastq','_pair.fastq')
        outp_f_2u = inp_f2.replace('.fastq','_unpair.fastq')
        outp_f_2p = inp_f2.replace('.fastq','_pair.fastq')
        #path_outp = '/local/data/course/project/groups/the_parsers/' + oup_f_1u
        if not os.path.exists(outp_f_1p):
            cmd = 'java -jar /local/prog/trimmomatic/trimmomatic-0.33.jar PE -thread 30 -phred{} \
            {} {} {} {} {} {} ILLUMINACLIP:{} LEADING:3 TRAILING:3 SLIDINGWINDOW:{} MINLEN:{}'\
                  .format(phred, inp_f1, inp_f2, outp_f_1p, outp_f_1u, outp_f_2p, outp_f_2u\
                          ,adapter, slidwin, minlen)
            subprocess.check_call(cmd, shell =  True)
        return outp_f_1p, outp_f_1u, outp_f_2p, outp_f_2u
    
if __name__ == '__main__':
    raw_reads_left = argv[1]
    raw_reads_right = argv[2]
    run_trimmomatic(raw_reads_left, raw_reads_right)
