#!/usr/bin/env python3
"""
Author: Ioannis Vasilas

This python script has to parse a GenBank file and outputs a FASTA file
and an ordered table with detail statistics
"""

#imports

from sys import argv
import subprocess 
import os.path 


#functions

def record_finder(lines):
    """Return list of records separated by each read
    lines: open file or list of lines
    
    input: lines of an open file
    output: list of a block of one read
    """
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("@"): #separate each block of reads by @
            try:
                yield curr
            except:
                pass
            curr = []
            curr.append(line.strip())
        else:
            curr.append(line.strip())
    if curr:
        yield curr #Sandra et al. 2019
    
        
        
        
        
def parse_fastq (rec_lines):
    """it takes from the block of the read the label,
    sequence and scores
    
    input: lines of a block of a read
    output: tuple of the 3 elements
    """
    data = []
    data.append(rec_lines[0][1:])
    data.append(rec_lines[1])
    data.append(rec_lines[3])
    return data
        
        
        
def trans_ascii (data):
    """translates ascii code to score
    
    input: a tuple which its 3rd element is the ascii code
    output: phred score by position
    """
    quality = []
    asl = data[2]
    for char in asl:
        change= ord (char)-64
        quality.append(change)
        
    return quality
    
    
    
    
def cal_length (datalist):
    """calculation of the min length, max length and average length
    
    input: a data list of tuples which contains sequences as a second element
    """
    sortedlist = sorted(datalist, key = lambda x : len(x[1]))
    maxle = len(sortedlist[-1][1])
    minle = len(sortedlist[0][1])
    average = sum([len(x[2]) for x in sortedlist])/len(sortedlist)
    return minle, maxle,  average
    
    
    
    
def aver_score(datalist):
    """calculation of the average phred score per position
    
    input: a data list of tuples which contains phred scores as a fourth element
    output: a list of average phred score per position
    """
    scores_per_position = []
    
    for tupl in datalist:
        count = 0
        sum_of_position = 0
        for element in tupl[3]:
            sum_of_position += element
            count +=1
        aver_pos = sum_of_position/ count
        scores_per_position += [aver_pos]
        
    return scores_per_position




def trimming(input_file, threshold =30):
    """Run fastq_quality_trimmer on fastq file
    input_file: string, filename of input FASTQ file
    threshold: int, Quality threshold - nucleotides with lower 
    quality will be trimmed (from the end of the sequence).
    
    output: a fastq file with the trimmed reads
    """
    output_file = "trimmed{}.fq".format(threshold)
    if os.path.exists(output_file): 
        return
    command = 'fastq_quality_trimmer -Q64 -t {} -i {} -o{}'\
        .format(threshold,input_file,output_file)
    e =subprocess.check_output(command,shell = True)
    return output_file
    
    
#main
    
if __name__ == "__main__":
    #parsing, calculating lengths, calculating phred scores, average scores and creating a
    #of tuples for the non-trimmed data
    data = []
    for rec in record_finder(open (argv[1])):
        parse_fasta = parse_fastq(rec)#parsing
        qual = trans_ascii (parse_fasta) #calculating phred scores
        data.append((parse_fasta[0],parse_fasta[1], parse_fasta[2], qual))
    spec_lengths = cal_length(data)#calculating lengths
    avscore = aver_score(data)#average scores
    
    #running tool for trimming the reads
    trimmed = trimming(argv[1])
    
    
    
    
    #parsing, calculating lengths, calculating phred scores, average scores and creating a
    #of tuples for the trimmed data
    data2 = []
    for rec in record_finder(open (argv[2])):
        parse_fasta = parse_fastq(rec)#parsing
        qual = trans_ascii (parse_fasta) #calculating phred scores
        data2.append((parse_fasta[0],parse_fasta[1], parse_fasta[2], qual))
    spec_lengths2 = cal_length(data2) #calculating lengths
    avscore2 = aver_score(data2)#average scores
    
    
    
    
    #report of the results
    report = open('results.txt', 'w')
    min1,max1,avg1 = spec_lengths
    min2, max2, avg2 = spec_lengths2
    #writing the lengths for original sequences
    report.write('ORIGINAL: min={}, max={}, avg={} \n'.format\
        (min1,max1,avg1))
    
    #writing the lengths for trimmed sequences
    report.write('TRIMMED: min={}, max={}, avg={}\n'.format\
        (min2, max2, avg2))
        
    for scores in range(len(avscore)): # writing the phred scores per sequence in the txt
        report.write(str(scores)+'\t'+str('%5.2f'%(avscore[scores]))\
            +'\t'+str('%5.2f'%(avscore2[scores]))+'\t'+\
            str('%5.2f'%(abs(avscore[scores]-avscore2[scores])))+'\n')
    report.close()
