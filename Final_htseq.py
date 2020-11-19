#!usr/bin/env python3
"""
@authors: Jasper, Marijn Schipper
"""

import sys
import subprocess
import os

def run_htseq_count(bam_file, gff_file, out_htseq):
    '''runs htseq -count program in linux command line.
    bam_file = input sorted bam_file of RNAseq alignment (str)
    gff_file = input refenece gene-models gff file (str)
    out_htseq = htseq output location (str)
    returns:
        name of outputted htseq -count output (str)
    '''
    if not os.path.exists(out_htseq):
        cmd = "mkdir {}".format(out_htseq)
        subprocess.check_call(cmd, shell = True)

    output_file_name = bam_file.rpartition("/")[-1].split(".")[0]+".txt"
    if not os.path.exists(output_file_name):
        cmd = "htseq-count -f bam -i Parent {} {} > ./{}/{}"\
            .format(bam_file, gff_file, out_htseq, output_file_name)
        subprocess.check_call(cmd, shell = True)
    return output_file_name


