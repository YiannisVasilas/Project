#!/usr/bin/env python3

""" 


GFFparser 04/12/2019
GFF parser
"""

# imports
from sys import argv
import subprocess
import os.path
from symbolic_link_fastq import list_files_in_directory
from copy import deepcopy


def gff_parser(gfffile):
    """Parses transcriptname out of GFF-file
    
    Keyarguments:
        gfffile -- name of file in GFF-format
    Returns:
        mrna_dict -- dictionary with transcript as key and Parent as value
        gene_dict -- dictionary with gene as key and counts as value set to 0,
            will be added in another function
    """
    mrna_dict = {}
    gene_dict = {}
    for line in gfffile:
        if not line.startswith("#"):
            lineparts = line.strip().split("\t")
            features = lineparts[8].split(";")
            feature_dict = {}
            for feature in features:
                x = feature.split("=")
                feature_dict[x[0]]=x[1]
            if lineparts[2] == "gene":
                gene_dict[feature_dict["Name"]] = 0
            elif lineparts[2] == "mRNA":
                mrna_dict[feature_dict["ID"]] = feature_dict["Parent"][5:]
    return mrna_dict, gene_dict
           

def gene_count(gene_count_output, mrna_dict, gene_dict, htseqfile):
    """Adds transcript counts to gene in dictionary
    
    Keyarguments:
        mrna_dict -- dictionary with transcript as key and Parent as value
        gene_dict -- dictionary with gene as key and counts as value set to 0,
            will be added in another function
        htseqfile -- name of file in tsv-format
    Returns:
        gene_dict -- dictionary {gene: count}
    """
    for line in htseqfile:
        if not line.startswith("_"):
            mrna, count = line.strip().split("\t")
            gene = mrna_dict[mrna]
            gene_dict[gene] += int(count)
    return gene_dict

def write_output(gene_count_output, gene_dict, name):
    """Writes output to tabdelimited file
    
    Keyarguments:
        gene_dict -- dictionary {gene: count}
        name -- name of the file with transcript counts in tsv-format
    Returns:
        nothing, creates tsv-file
    """
    if not os.path.exists(gene_count_output):
        command = "mkdir {}".format(gene_count_output)
        run = subprocess.check_output(command,shell = True)
    new_name = "{}/gene_count_{}".format(gene_count_output,\
     name.rpartition("/")[-1])
    f = open(new_name, "w+") 
    genes = list(gene_dict.keys())
    sorted_genes = sorted(genes)
    for gene in sorted_genes:
        f.write(gene+"\t"+str(gene_dict[gene])+"\n")
    f.close()
    
def parse_htseq_count_output(gene_count_output, path_gff, path_to_htseq):
    """Parses htseq output and stores gene count as dict {gene(str):count(int)}
    
    Keyarguments:
        gene_count_output -- filename and location for output(str)
        path_gff -- path to gene model gff (str)
        path_to_htseq -- path to htseq output file(str)
        
    Returns:
        writes gene name to gene count file
    """
    with open(path_gff) as fo:
        mrna_dict, gene_dict = gff_parser(fo)
    for file_name in list_files_in_directory(path_to_htseq):
        temp_gene_dict = deepcopy(gene_dict)
        path = path_to_htseq+file_name
        with open(path) as htseq:
            gene_dict_updated = gene_count(gene_count_output, mrna_dict,\
             temp_gene_dict, htseq)
        write_output(gene_count_output, gene_dict_updated, path)
        gene_dict_updated = {}

if __name__ == "__main__":
    parse_htseq_count_output(argv[1], argv[2])
    
