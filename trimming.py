#!/usr/bin/env python

"""
Author: Adzkia Salima

Script to run trimmomatic

Usage:
    python trimming.py accs_num folder_input adapter output-dir

    accs_num: list of accession number of files to be trimmed
    folder_input: folder of input files
    adapter
    



"""

from sys import argv
import subprocess
import os.path



def run_trimomatic(data_dir, output_dir, accs_num, adapter):
    """Runs trimmomatic program

    Keyword arguments:
        data_dir: input data directory
        output_dir: output data directory
        accs_num" accession number of bad reads

    """

    cmd = 'java -jar /local/prog/Trimmomatic-0.36/trimmomatic-0.36.jar\
     PE %s/%s_STARoutput/Unmapped.out.mate1 %s/%s_STARoutput/Unmapped.out.mate2 \
     %s/%s_1t.fastq %s/%s_1unpaired.fastq\
     %s/%s_2t.fastq %s/%s_2unpaired.fastq\
     ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'\
     %(data_dir, accs_num, data_dir, accs_num, output_dir, accs_num, output_dir, accs_num, output_dir, accs_num, output_dir, accs_num, adapter)

    result = subprocess.check_output(cmd, shell = True)
    subprocess.check_output('rm *unpaired*', shell = True)
   
    return result
           

   
if __name__ == "__main__":

    accs_nums = argv[1]
    data_dir = argv[2]
    adapter = argv[3]
    output_dir = argv[4]
    for accs_num in accs_nums:
      run_trimomatic(data_dir, output_dir, accs_num, adapter1)
      
##    data_dir = "/local/data/BIF30806_2016_2/project/groups/Unixcorn*power/data/denovo_data"
##    #accs_nums = ["SRR122236", "SRR122237", "SRR122238", "SRR1820326", "SRR342017", "SRR342019", "SRR342022", "SRR342023", "SRR646572", "SRR1144633","SRR648705"] #original run
##    #accs_nums = ["SRR646596","SRR646604","SRR924147","SRR646598", "SRR648707"] #missing file1
##    accs_nums1 = ["SRR342017", "SRR342019", "SRR342022", "SRR342023", "SRR648705", "SRR648707"]#adapter fail, run using adapater truseq2
##    #accs_nums = ["SRR924147", "SRR924148"] #data nyusul1
##    #accs_nums1 = ["SRR122236", "SRR122237", "SRR122238","SRR342017", "SRR342019", "SRR342022", "SRR342023"]#original
##    adapter1 = "TruSeq2-PE.fa"
##    #accs_nums2 = ["SRR646572", "SRR646596","SRR646604","SRR1144633","SRR648705","SRR648707","SRR924147", "SRR924148","SRR1820326"]#original
##    #accs_nums2 = ["SRR648705","SRR648707"] #wytze request 12/12/2016
##    adapter2 = "TruSeq3-PE.fa"
##    output_dir = "/local/data/BIF30806_2016_2/project/groups/Unixcorn*power/trimomatic_test/experiment2"
##    for accs_num in accs_nums1:
##        run_trimomatic(data_dir, output_dir, accs_num, adapter1)
##    #for accs_num in accs_nums2:
##    #    run_trimomatic(data_dir, output_dir, accs_num, adapter2)
