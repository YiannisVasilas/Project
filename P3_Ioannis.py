#!/usr/bin/env python3
"""
Author: Ioannis
Author : Anastasia
Student Number:

This python script has to parse a GenBank file and outputs a FASTA file
and an ordered table with detail statistics
"""



from sys import argv
import subprocess 
import os.path 
from operator import itemgetter

def record_finder(lines):
    """Return list of records separated by each read
    lines: open file or list of lines
    """
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("@"):
            try:
                yield curr
            except:
                pass
            curr = []
            curr.append(line.strip())
        else:
            curr.append(line.strip())
    if curr:
        yield curr
    
        
def parse_fastq (rec_lines):
    data = []
    data.append(rec_lines[0][1:])
    data.append(rec_lines[1])
    data.append(rec_lines[3])
    return data
        
def trans_ascii (data):
    quality = []
    asl = data[2]
    for char in asl:
        change= ord (char)-64
        quality.append(change)
    return quality
    
def cal_length (datalist):
    sortedlist = sorted(datalist, key = lambda x : len(x[1]))
    maxle = len(sortedlist[-1][1])
    minle = len(sortedlist[0][1])
    average = sum([len(x[2]) for x in sortedlist])/len(sortedlist)
    return maxle, minle, average
    
def aver_score(datalist):
    scores_per_position = []
    count = 0
    for tupl in datalist:
        count +=1
        sum_of_position = 0
        for element in tupl[3]:
            sum_of_position += element
           
        aver_pos = sum_of_position/ len(datalist)
        
        scores_per_position += [aver_pos]
    print(count)
    return scores_per_position

def trimming(input_file, threshold =30):
        """Run fastq_quality_trimmer on fastq file
    input_file: string, filename of input FASTQ file
    threshold: int, Quality threshold - nucleotides with lower 
    quality will be trimmed (from the end of the sequence).
    """
    output_file = "trimmed{}.fq".format(threshold)
    command = 'fast_quality_trimmer -Q64 -t {} -i {} -o{}'\
    .format(thershold,input_file,output_file)
    e =subprocess.check.output(command,shell = True)
    return output_file
    
    
def report():
    pass
    
if __name__ == "__main__":
    data = []
    for rec in record_finder(open (argv[1])):
        #print(rec)
        parse_fasta = parse_fastq(rec)
        qual = trans_ascii (parse_fasta) 
        data.append((parse_fasta[0],parse_fasta[1], parse_fasta[2], qual))
    spec_lengths = cal_length(data)
    avscore = aver_score(data)
    print (avscore)
