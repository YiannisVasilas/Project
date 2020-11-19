#!/usr/bin/env python3
#
#  Author : 	Ioannis va
#  P3 exercise: old exam


from sys import argv
import os.path
import subprocess

def parse_input(input_fn):
    """
    parse fastq file
    returns fastq_dict: dict with identity, sequence, quality score of fastq file
    """
    
    fastq_dict = {}
    with open(input_fn, "r") as input_file:
        identity = ""
        sequence_coming = False
        seq = []
        quality_coming = False
        quality_score = []
        for line in input_file:
            line = line.strip()
            
            if line.startswith("@"):
                identity = line
                fastq_dict[identity] = {}
                sequence_coming = True
                quality_coming = False
                seq = []
            elif quality_coming:
                quality_score.append(line) 
                quality_score = "".join(quality_score)
                fastq_dict[identity]["Quality"] = quality_score
            elif line.startswith("+"):
                sequence_coming = False
                quality_coming = True
                fastq_dict[identity]["Sequence"] = seq
                quality_score = []
            elif sequence_coming:
                seq.append(line)
                seq = ''.join(seq)
         
    return fastq_dict
    
 
def fastq_trimmer(input_fn,output):
    """
    run fastq_quality_trimmer
    """
    if not os.path.exists(output):
        cmd = 'fastq_quality_trimmer -i %s -t %s -Q %s -o %s '\
            %(input_fn, 30, 64, output)
        subprocess.check_call(cmd, shell=True)
    return 0 
        
def sequence_length(fastq_dict):
    """
    min max avg sequence length 
    
    returns min, max, avg length of the sequences of fastq file
    """
    
    seq_length = []
    for key in fastq_dict.keys():
        seq_length.append(len(fastq_dict[key]["Sequence"]))
        max_length = max(seq_length)
        min_length = min(seq_length)
        avg_length = (sum(seq_length)/len(seq_length))
        output_template = "min: {}\tmax: {}\tavg: {}\n"
        output_lengths = output_template.format(min_length,max_length,avg_length)
    
    return output_lengths
    
def avg_quality(fastq_dict):
   
    """
    avg quality in each position for all entries
    
    avg_qual: dict, position: average of qualities scores in same position
    """
    qualities = {}
    for key in fastq_dict.keys():
        qual_score = fastq_dict[key]["Quality"]
        count = 0 
        for item in qual_score:
            count +=1
            if count in qualities:
                qualities[count].append((ord(item)-64))

            else:    
                qualities[count] = [(ord(item)-64)]

    avg_qual = {}
    for key,lis in qualities.items():
        avg_qual[key] = (sum(lis)/float(len(lis)))
    
    return avg_qual

if __name__=="__main__":
    out_file = 'finaltrimmed.fq'
    #input tomatosapmle.fq
    input_sample= argv[1]
    parse_fastq = parse_input(argv[1])
    
    #seq lengths max min average
    seq_length = sequence_length(parse_fastq)
    
    
    #fastq trimmer of tomatosample.fq
    fastq_trimmer = fastq_trimmer(input_sample,out_file)
    
    #average quality of tomatosample.fq
    avg_qual = avg_quality(parse_fastq)
    
    #trimmed file trimmed.fq
    parse_trimmed = parse_input(out_file)
    
    #seq lengths max min average
    seq_length_trimmed  = sequence_length(parse_trimmed)
    
    
    #average quality of trimmed.fq
    avg_trimmed = avg_quality(parse_trimmed)
    
    #difference in average quality trimmed-tomatosample
    print("ORIGINAL:",seq_length)
    print("TRIMMED:",seq_length_trimmed)
    for i in avg_qual:
        output_template = "{}\t{}\t{}\t{}\n"
        output_file = output_template.format(i, avg_qual[i], avg_trimmed[i], avg_trimmed[i] - avg_qual[i])
        print(output_file)
