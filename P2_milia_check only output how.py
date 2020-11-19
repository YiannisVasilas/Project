#!/usr/bin/env python3
#
#  Author : Sotiria Milia
#  P2 exercise: Parsing

from sys import argv
import re ##for regular expression


def parse_genbank(input_fn):
    """
    Returns file_gen dictionary {accession number:{organism:, DNA sequence:}} 
    from genbank file

    input file : genbank format file
    """
    
    file_gen = {}
    fo = open(input_file,"r").read()



   ### input_file = open(input_file,"r")
    with open(input_fn,"r") as input_file:
        accession_number = "" ##an dwsei la8os katalavainw oti einai e3aitias autou
        sequence_coming =False
        seq = []
        for line in input_file:
            line = line.strip() ###removes first spaces
            if sequence_coming:
                line_elements = line.split()
                seq_parts = line_elem[1:]
                for seq_part in seq_parts:
                    seq.append(seq_part)
            elif line.startswith("ACCESSION"):
                line_elements = line.split()
                accession_number = line_elements[1]
                accession_dict[accession_number]={}
            elif line.startswith("ORGANISM"):
                line_elements = line.split()
                organism = line_elements[1:] ### lista poy prepei na alla3w
                organism = ' '.join(organism) ###kanei th lista enwmenh me kena
                accession_dict[accession_number]["organism"]=organism
            elif line.startswith("ORIGIN"):
                sequence_coming = True
            elif line.startswith("//"):
                sequence_coming = False
                seq = "".join(seq)
                accession_dict[accession_number]["Sequence"] = seq
                seq = [] ###prepei na einai adeia h lista 3ana alliws 8a kanei append sth palia
        
    ####input_file.close()       
        return accession_dict
            
            
            
    
    entries = fo.split("LOCUS")[1:]
    
    for line in entries:
        
        accession_num = line.split("ACCESSION")[1].partition("\n")[0].strip()
        file_gen[accession_num]  = {}
        organism = line.split("ORGANISM")[1].partition("\n")[0].strip()
        file_gen[accession_num].update({organism:""})
        ###main issue was to find a way to parse the many lines that the sequence was, \
        ###I could only get one line and not multiple but finally i found the way i am using\
        ###probably there are betters ways to parse the file
        sequence = line.split("ORIGIN")[1].partition("//")[0].strip()
        sequence = re.sub('[^atgc]','',sequence) ###google search-exact, keeps only atgc(lower) letters
        file_gen[accession_num][organism]+= sequence ###i don't know why i coulnd't use update or append

        ###i checked the type(file_gen) and it says <class 'dict'> \
        ###and i don't understand why it says class and how i made that
        
    return file_gen    
    

def out_files(file_gen):
    """
    Writes two files in fasta and tab-delimited format
    Fasta format file: >accession number, organism, sequence
    Tab-delimited: accession number, organism, gc content, length of DNA sequence  

    input: dict from parse_genbank function
    
    """
    ###((RANDOM ORDER)-needs from high to low gc content which I couldn't manage to make it in any possible way
    out_1 = open("output_fasta.txt", "w")
    out_2 = open("output_tab.txt","w")
    for k in file_gen.keys():
        for key,value in file_gen[k].items():
            gc_content = str(round((float(value.count("g") + value.count("c"))/(len(value)))*100,2))
            ### Below I add the gc content to the dict in order to find a way to sort it according to that value, but \
            ###1)i couldn't find a way to sort by the second value and 2)then also to be able to write it in the files sorted
            file_gen[k][key]= value, gc_content 
            
            out_1.write(">" + k + " " + str(key) + "\n" + file_gen[k][key][0] + "\n")
            out_2.write(k + "\t"+ str(key) +"\t" +str(gc_content) +"\t"+ str(len(value)) +"\n")
            
            


if __name__=="__main__":
    #input genbank file
    input_file = argv[1]
    #parse genbank file
    parse_dict = parse_genbank(input_file)
    #fasta and tab-delimited output
    #out_files(parse_file)
    
