#!/usr/bin/env python3


""" Ioannis Vasias
16/11/2019 """
#imports 

from sys import argv
 
import myvariant


#functions


def parsevariant (fil):
    """import a file and split by line
    
    output a list"""
    results = []
    for line in fil:
        
       data= line.strip().split()
       results.append(data)
    return results
    
def retrive_variants(list_var):
    """ Takes the variants form variant.info
    input = a list of varinats
    outtupt = txt file with the variant"""
   
    mv = myvariant.MyVariantInfo
    record = open("variants.txt", 'w')
    for entry in list_var:
        var = mv.getvariant('{}:g.{}{}>{}',fields).format(entry[0],entry[1],\
            entry[2],entry[3]
        record.write(var)
    record.close()
        mv = myvariant.MyVariantInfo
    record = open("variants.txt", 'w')
    
    
def parse_records (rec_file):
    """ breaks the file into records
        input = file from variant.info
        """
            
    curr = []
    for line in rec_file:
        if not line.strip():
            continue
        if line.startswith("}"):
            if curr:
                yield curr
                curr = []
        else:
            curr.append(line.strip())
    if curr:
        yield curr


def parse_high (record_list):
    """ check the record for high putative 
        input = list of lines
        output = id and high variant""",
        
    for line in records_list:
        if 'HIGH' in line
            return record_list([1][8:-2]
            
def retrive_frequences(list_var):
        """ Takes frequnces form variants
    input = a list of varinats
    outtupt = txt file with the variant"""
    
    
    
        mv = myvariant.MyVariantInfo
    record = open("frequnces.txt", 'w')
    for entry in list_var:
        var = mv.getvariant('{}:g.{}{}>{}').fields = "exac.alleles,exac.af')\
            .format(entry[0],entry[1],entry[2],entry[3])
        record.write(var)
    record.close()
    

def MAF (frequency_list):
    """ parse the  
    input: frequences
    out: tuble with the frequnces <0.001"""
    
    
     freq = []
    for line in frequency_list:
        if line startswith('af'):
            freq +=line[0]
            if freq < 0.001:
            res = (frequency_list[1][8:-2] ,freq)
    return (freq)    
    
        
    
    
    

if __name__ == "__main__":
    
    parse = parsevariant (open(argv[1],'r'))
          
    ret= retrive_variants(parse)
    highlist = []
    for record in parse_records(open(argv[1],'r'):
        p=parse_high(record)
        if p != Nione:
            highlst.append(p)
    frequence = retrive_frequences (parse)
