# -*- coding: utf-8 -*-
"""
C
"""


from __future__ import division
from sys import argv
import subprocess
import os, sys
import re
import logging


def parse_gtf(filename):
    """
    Extracts scaffold number, transcript ID, reference ID and FPKM value from a gtf file and
    puts them in a matrix in this exact order. Removes all double instances of data points by 
    filtering out non-unique gene ID's. 
    
    Keyword arguments: 
        filename -- Filename and location of the gtf file to be parsed 
    Returns: 
        A matrix with scaffold number, trandscript ID, reference ID and FPKM value per 
        unique data point in that order 
    """
    
    data = open(filename, 'r')
    GTF = data.read()
    data.close()
    lines = GTF.split('\n')
#    print lines 
    usefull = []
    endlist = []
    transcript_ids = []
    
    for line in lines: 
        linelist = line.split('\t')
        scaffold = linelist[0]
        ID = linelist[-1].split(' ')
        ID += [scaffold]
        usefull += [ID]
    for elem in usefull:
        FPKMvalue = None
        scaffold = elem[-1]
        reference_id = None
        for i in range(len(elem)): 
            if elem[i] == 'transcript_id': 
                transcript_id = elem[i+1]
#                print transcript_id
            elif elem[i] == 'reference_id': 
                reference_id = elem[i+1]
#                print reference_id
            elif elem[i] == 'FPKM': 
                FPKMvalue = elem[i+1]
#                print FPKMvalue
                
            if FPKMvalue != None and transcript_ids.count(transcript_id) == 0: 
                if reference_id != None: 
                    line = [scaffold, transcript_id, reference_id, FPKMvalue]
                else: 
                    line = [scaffold, transcript_id, 'unavailable', FPKMvalue]
                endlist += [line]
                transcript_ids += [transcript_id]
#    print transcript_ids
#    print endlist
    return endlist             
                
              
               

def write_output(name, data): 
    """
    Writes output from a 4x* matrix into a text file. 
    
    Keyword arguments: 
        name -- name of the job & output file 
        data -- A 4x* matrix, as made by parse_gtf 
    """
    
    if os.path.exists('denovomapping/{}.txt'.format(name)):
        logging.info('{}.txt is already present'.format(name))
        print '{}.txt is already present'.format(name)
    else: 
        output_file = open('denovomapping/{}.txt'.format(name), 'w')
        header = 'scaffold \t transcript_id \t reference_id \t FPKM \n'
        output_file.write(header) 
        
        for thing in data: 
            line = '{} \t {} \t {} \t {} \n'.format(thing[0], thing[1], thing[2], thing[3])                      
            output_file.write(line)
            
        print 'data was written to {}.txt'.format(name)

    
    
    
if __name__=='__main__':

     folder = os.listdir('/local/data/course/project/groups/PacWay/denovomapping') 
     gtffiles = []    
     for file in folder: 
         if file.endswith('.gtf') == True: 
             gtffiles += [file]
#     print gtffiles
 
     for GTF in gtffiles:    
         gtf_file = GTF
         name = str(GTF).split('/')[-1].split(".")[0]
    
         data = parse_gtf('denovomapping/{}'.format(GTF))
         write_output(name, data)
