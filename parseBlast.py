#!/bin/usr/env python

Cro2Gene = {}

with open(argv[1]) as blast_file:
    blastout = blast_file.readlines()
    for line in blastout:
        if line.startswith("Query= CRO"):
            cro = line[7:-1]
            tophit = True
        if line.startswith("> gnl") and tophit == True:

            genename = line[2:]
            tophit = False
            Cro2Gene[cro] = genename

outstr = ""
for key in Cro2Gene.keys():
    outstr += "{}\t{}".format(key, Cro2Gene[key])

with open("sigGeneNames.txt", "w") as outfile:
    outfile.write(outstr)
