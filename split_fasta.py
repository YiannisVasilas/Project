#!/usr/bin/env python

"""




"""

from sys import argv
import subprocess
import os.path



def split_fasta(inp):
	"""
	Void function which takes a folder containing merged paired end reads is taken as input, writes two files of originalname_1 and originalname_2 to the input folder.
	"""
	splinp = open(inp).read().split('\n')
	print('processing %s' %(inp))
	
	# N.B. assumption is made that read length L is meets 99 < L < 1000 (the [:-3] term)
	# Creates a list of lists, where each entry is of the format [Header (incl. adj. length), complete read sequence, second header/identifier (incl. adj. length), full quality sequence]
	fw =  [[splinp[x*4][:-3] + str(int(splinp[x*4][-3:])/2), splinp[x*4+1],splinp[x*4+2][:-3] + str(int(splinp[x*4][-3:].replace('=',''))/2),splinp[x*4+3]] for x in range(len(splinp)/4)]
	outf = open('%s_1.fastq' %(inp.replace('.fastq','')),'w')
	outf2 =  open('%s_2.fastq' %(inp.replace('.fastq','')),'w')
	for x in fw:
		n = int(x[0][-3:]) 
		outf.write('%s\n%s\n%s\n%s\n' %(x[0],x[1][:n],x[2],x[3][:n]))  #In addition to the headers/identifiers it prints first half of both nucleotide and quality sequence (0 to n)
		outf2.write('%s\n%s\n%s\n%s\n' %(x[0],x[1][n:],x[2],x[3][n:])) # prints second half (n to end)
