#!/usr/bin/env python3
"""
Author: Ioannis Vasials
Parsing, acc. numbers, organisms, seqs
"""

#imports
from sys import argv
from operator import itemgetter

#functions

def record_finder (lines):
    """ take a list from first record in file
    """
    record = []
    for line in lines:
        if not line.strip():
            continue
        elif line.startswith('//'):
            yield record
            record = []
        else:
            line = line.strip()
            record.append(line)
    if record:
        yield record
        
def parse_data(reclines):
    """ takes lines and parses the accesion number, organism and the 
    species from GenBank file
    
    Output: tuple of the three mentioned
    """
    index = 0
    rec_id = None
    rec_species = None
    for lin in reclines: #parse accesion number and organism
    
        if lin.startswith('ACCESSION'):
            line = lin.strip().split('   ')
            rec_id = line[1]
        elif lin.startswith('ORGANISM') or lin.startswith('  ORGANISM'):
            line = lin.strip().split('  ')
            rec_species = line[1]
        elif lin == 'ORIGIN':
            index = reclines.index('ORIGIN')  
			
    sequences = []
    for seq in range(index+1, len(reclines)):# parse the sequence per line
        if reclines[seq] == '//':
            break
        else:
            piece = reclines[seq].strip().split('  ')[1]
            seqs = piece.split(' ')
            sequences += seqs
    rec_seq = ''.join(sequences) #join substring to one sequence
    return rec_id, rec_species, rec_seq
    
     
def count_GC (sequence):
    """calculates the GC content of each sequence and reports a tuple of name and GCcontent

    input: a tuple of lists of names and sequences"""
    sequence = sequence.upper()
    gccount = 0
    for base in sequence:
        if base == 'G' :
            gccount +=1
        elif base == 'C':
            gccount +=1
    gccont = (gccount/len(sequence))*100
    gccont = format(gccont,'.2f')
    return gccont
    
#writing fasta file
def write_fasta(filename, tupl):
    fast = open(filename, 'w')
    for item in tupl:
        rec_id, rec_species,gccont, rec_seq = item
        fast.write('>'+rec_id+' '+rec_species+'\n')
        fast.write(rec_seq+'\n')
    fast.close() 
    
    
def write_txt(filename, tupl):	
	#writing txt file 
    txtf = open('output.txt', 'w')
    for item in tupl:
        rec_id, rec_species,gccont, rec_seq = item
        if len(rec_id) > 12:
            rec_id = rec_id.split()
            for ind in range (len(rec_id)):
                
                txtf.write('{:12s}\t{:25s}\t{}\t{}\n'.format(rec_id[ind],rec_species,gccont,len(rec_seq)))
        else:
            txtf.write('{:12s}\t{:25s}\t{}\t{}\n'.format(rec_id,rec_species,gccont,len(rec_seq)))	
    txtf.close()   


#main
if __name__=='__main__':
    
    result = []
    with open(argv[1]) as fo:
        for record in record_finder(fo):
            rec_id, rec_species, rec_seq = parse_data(record)
            gccont = count_GC (rec_seq)
            result.append((rec_id, rec_species,float(gccont), rec_seq))
    result = sorted(result, key = itemgetter(2), reverse = True)
    
    #writing fasta file
    write_fasta('fasta_output.fasta', result)
	
	#writing txt file
	write_txt('output.txt', result)
    
    
    
    
