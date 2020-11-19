#!/bin/usr/env python

def parse_fasta(fasta_fn):
    """
    Parses fasta file returning a dict with  labels as keys and seqs as values.

    :param fasta_fn: Filename of the FASTA file to parse
    :return: A dictionary containing labels as keys and sequences as values.
    """
    with open(fasta_fn) as fasta_file:
        fasta_list = fasta_file.read().splitlines()
        parsed_seqs = {}
        for line in fasta_list:
            if line.startswith(">"):
                label = line[1:]
                parsed_seqs[label] = ""
            else:
                parsed_seqs[label] += line
    return parsed_seqs


with open("DE_geneslist.csv") as genelist_f:
    genes = genelist_f.readlines()
    genes_list = []
    for row in genes:
        gene = row.split(',')[1]
        gene = gene[1:5] + "T" + gene[5:-2]
        print gene
        genes_list.append(gene)

genes_list = genes_list[1:]

ref_prots = parse_fasta("cro_std_maker_anno.final.pep.fasta")


with open("sig_prots.fasta", "w") as sig_prots:
    for gene in genes_list:
        sig_prots.write(">" + gene + "\n")
        sig_prots.write(ref_prots[gene] + "\n")

