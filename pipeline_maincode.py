#!/usr/bin/env python3

"""
Author: Yan Wang
course project
main code for pipeline
"""

from sys import argv
import subprocess
import os.path

import parse_data as pd
import fastq_splitter as fs
import fastqc_gunzip as fg
import trimmomatic as tr
import hisat2 as hi
import sambam as sb
import cuffdiff as cd 

if __name__ == '__main__':
    with open(argv[1]) as datasets:
        ref_fa, ref_gff, groups = pd.parse_data_set(datasets)
    
    work_dir = str(subprocess.check_output('pwd')).strip("b'").strip("'\\n")
    
    samples = list(groups.keys())
    expr = []
    ctrl = []
    for i in range(len(samples)):
        
        unzip_raw_file = fs.unzipper(samples[i])
        read_len = int(groups[samples[i]][1])
        raw_left, raw_right = fs.fastq_splitter(unzip_raw_file, read_len) 
        
        adapter = groups[samples[i]][3]
        phred = int(groups[samples[i]][4])
        slidwin = groups[samples[i]][5]
        minlen = int(groups[samples[i]][6])
        left_pair, left_unpair, right_pair, right_unpair = tr.run_trimmomatic(raw_left, raw_right, phred, adapter, slidwin, minlen)
        
        if groups[samples[i]][2] == 'yes':
            #fastqc consumes huge time, run it if necessary
            out_dir = work_dir
            fg.run_fastqc(raw_left, raw_right, out_dir)
            fg.run_fastqc(left_pair,right_pair, out_dir)
    
        # making index
        base = work_dir + '/C_roseus_index' 
        hi.hisat2_index(ref_fa, base) 

        # run hisat
        index = 'C_roseus_index'
        hisat_out = samples[i].replace('.fastq.gz', '_mapped.sam')
        hi.hisat2(index, left_pair, right_pair, hisat_out)
        if not os.path.exists(hisat_out):
            raise ValueError('sam file not exist')
            
        # change sam to bam file
        sb.sam_to_sortedbam(hisat_out)
        bam_file = hisat_out.replace('.sam','.bam')
        if not os.path.exists(bam_file):
            raise ValueError('bam file not exist')
        
        # group the files together to run cuffdiff
        
        if groups[samples[i]][0] == 'expr':
            expr.append(bam_file)
        else:
            ctrl.append(bam_file)
        
    # use the bam files from two groups to run cuffdiff    
    inp_expr = ','.join(expr)
    inp_ctrl = ','.join(ctrl)
    sample_file = [inp_ctrl, inp_expr]
    print(sample_file)
    out_dir = 'cuffdiff_out'
    cd.run_cuffdiff(ref_gff, ref_fa, sample_file, out_dir)
    
    
    
    
    
    
    
    
    
    
    
    
