#!/usr/bin/env python3
#
#  Author : Sotiria Milia
#  P2 exercise: Parsing

from sys import argv

def parse_genbank(input_fn):
    """
    Returns file_gen dictionary {accession number:{organism:, DNA sequence:}} 
    from genbank file

    input file : genbank format file
    """
    
    accession_dict = {}
   
    with open(input_fn,"r") as input_file:
        accession_number = "" 
        sequence_coming =False
        seq = []
        for line in input_file:
            line = line.strip() 
            if line.startswith("//"):
                sequence_coming = False
                seq = "".join(seq)
                
                accession_dict[accession_number]["Sequence"] = seq
                seq = []
            elif sequence_coming:
                line_elements = line.split()
                seq_parts = line_elements[1:]
                for seq_part in seq_parts:
                    seq.append(seq_part)
            elif line.startswith("ACCESSION"):
                line_elements = line.split()
                accession_number = line_elements[1]
                accession_dict[accession_number]={}
            elif line.startswith("ORGANISM"):
                line_elements = line.split()
                organism = line_elements[1:] 
                organism = ' '.join(organism) 
                accession_dict[accession_number]["organism"]=organism
            elif line.startswith("ORIGIN"):
                sequence_coming = True
           
    return accession_dict
            

def out_files(file_dict):
    """
    Writes two files in fasta and tab-delimited format
    Fasta format file: >accession number, organism, sequence
    Tab-delimited: accession number, organism, gc content, length of DNA sequence  

    input: dict from parse_genbank function
    
    """
        
    with open("output_fasta.txt", "w") as output_fasta, open("output_tab.txt","w") as output_tab:
        accession_list = []
        for accession_number in file_dict.keys():
            organism = file_dict[accession_number]["organism"]
            sequence = file_dict[accession_number]["Sequence"].upper()
            gc_content = str(round((float(sequence.count("G") + sequence.count("C"))/(len(sequence)))*100,2))
            sequence_length = len(sequence)
            accession_list.append((accession_number, organism, sequence, gc_content, sequence_length))
        sorted_accession_list = sorted(accession_list, key=lambda tup:tup[3], reverse = True)
        for entry in sorted_accession_list:
            fasta_header_template = ">{} {}\n"
            fasta_header = fasta_header_template.format(entry[0],entry[1])
            output_fasta.write(fasta_header)
            output_fasta.write(entry[2])
            output_fasta.write("\n")
            tab_template = "{}\t{}\t{}\t{}\n"
            tab_line = tab_template.format(entry[0],entry[1],entry[3],entry[4])
            output_tab.write(tab_line)

    return 0
    

if __name__=="__main__":
    #input genbank file
    input_file = argv[1]
    #parse genbank file
    parse_dict = parse_genbank(input_file)
    #print(parse_dict)
    #fasta and tab-delimited output
    return_code = out_files(parse_dict)
    
