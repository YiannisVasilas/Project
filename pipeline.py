#!/usr/bin/env python

from sys import argv
import os
import subprocess

def run_sickle(inp_dir):
	"""
	Finds all files ending with _1.fastq in the input folder and tries to find its _2.fastq counterpart.
	When a suitable pair is found they are fed to sickle and the output it put in a folder called sickle_trimmed/
	"""
	fastq_files = [filename for filename in os.listdir(inp_dir) if filename.endswith('.fastq')]
	if not os.path.exists('sickle_trimmed'):
		os.makedirs('sickle_trimmed')
	
	for fname in fastq_files:
		if fname.endswith('_1.fastq'):
			in1 = fname
			in2 = fname[:-8] + '_2.fastq' #remove last characters rather than using replace to prevent replacing mid-string _1.fastq
			if in2 in fastq_files: #We have found a pair of files of the right format, proceed to feed it to sickle
				print('Trimming the pair  (%s - %s)' %(in1,in2))
				out1 = 'sickle_trimmed/%s' %('trimmed_' + in1)
				out2 = 'sickle_trimmed/%s' %('trimmed_' + in2)
				outsing = 'sickle_trimmed/%s' %(in1[:-8]+'_singles.fastq')
				#cmd = "sickle pe -f %s -r %s -t sanger -o %s -p %s -s %s" %(inp_dir + in1, inp_dir + in2,out1,out2,outsing)
				
				#CHANGE BACK TO RUN ON SERVER!!!!!===============================================================================
				cmd = "/local/prog/sickle/sickle pe -f %s -r %s -t sanger -o %s -p %s -s %s" %(inp_dir + in1, inp_dir + in2,out1,out2,outsing)
				subprocess.check_call(cmd, shell=True)
				
def run_samtools_split(inpdir):
	vfiles = [fname for fname in inpdir if fname.endswith('.bam')]
	
	if os.path.exists('unmapped/'):
		os.makedirs('unmapped')
	if os.path.exists('mapped/'):
		os.makedirs('mapped')
	
	for fil in vfiles:
		print('printing mapped and unmapped for %s' %(fil))
		cmd = 'samtools view -f 4 %s > unmapped/unmapped_%s' %(fil,file)
		subprocess.check_call(cmd, shell=True)
		cmd = 'samtools view -F 4 %s > mapped/mapped_%s' %(fil,file)
		subprocess.check_call(cmd, shell=True)
		

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

def run_fastqc(inp_dir,outp_dir):
	if not os.path.exists(outp_dir):
		os.makedirs(outp_dir)
	
	fastq_files = [filename for filename in os.listdir(inp_dir) if filename.endswith('_1.fastq') or filename.endswith('_2.fastq')]
	
	for fname in fastq_files:
		cmd = 'fastqc %s -o %s' %(inp_dir + fname, outp_dir)
		subprocess.check_call(cmd, shell=True)

def build_hisat2_index(ref_file, ref_name):
	cmd = '/local/prog/hisat2/hisat2-build %s %s' %(ref_file, ref_name)
	subprocess.check_call(cmd, shell=True)
	
def run_hisat2(inp_dir):
	if not os.path.exists('hisat2_output'):
		os.makedirs('hisat2_output')
	fastq_files = [filename for filename in os.listdir(inp_dir) if filename.endswith('_1.fastq')]
	
	for fname in fastq_files:
		cmd = '/local/prog/hisat2/hisat2 pipe_ref -1 %s -2 %s > hisat2_output/%s' %(inp_dir + fname, inp_dir + fname[:-7] + '2.fastq', fname + '.sam')
		subprocess.check_call(cmd, shell=True)

def build_kallisto_ref(ref_folder):
	cmd = "/local/prog/kallisto/kallisto index -i ref_transcripts %s" %(ref_folder)
	subprocess.check_call(cmd, shell=True)

def run_kallisto(inp_dir):
	"""
	Finds all files ending with _1.fastq and tries to find its _2.fastq counterpart.
	When a suitable pair is found they are fed to sickle and the output it put in a folder called sickle_trimmed/
	"""
	fastq_files = [filename for filename in os.listdir(inp_dir) if filename.endswith('.fastq')]
	
	ref = 'ref_transcripts'
	if not os.path.exists('Kallisto_output'):
		os.makedirs('Kallisto_output')
	outp_dir = 'Kallisto_output'
	for fname in fastq_files:
		if fname.endswith('_1.fastq'):
			in1 = fname
			in2 = fname[:-8] + '_2.fastq' #remove last characters rather than using replace to prevent replacing mid-string _1.fastq
			if in2 in fastq_files: #We have found a pair of files of the right format, proceed to feed it to sickle
				if not os.path.exists(outp_dir + fname[:-8]):
					os.makedirs(outp_dir + fname[:-8])
				
				print('processing %s and %s' %(in1,in2))
				cmd = "/local/prog/kallisto/kallisto quant -i %s -b 30 -o %s %s %s" %(ref, outp_dir + fname[:-8] , inp_dir + in1, inp_dir + in2)
				subprocess.check_call(cmd, shell=True)
	print('Done running kallisto.')

def main():
	try:
		data_folder = [argv[1] + '/',argv[1]][argv[1].endswith('/')]
		ref_file = argv[2]
		if not os.path.exists(ref_file):
			print('Reference genome file does not exist')
			exit()
	except:
		print('Please provide the folder containing the reads as an argument, followed by the location of the reference genome file.')
		exit()
	
	ref_name = ''
	
	#Split data
	if os.path.isdir(data_folder):
		fastq_files = [file_name for file_name in os.listdir(data_folder) if file_name.endswith('.fastq') ]
		for fname in fastq_files:
			print('Splitting file %s' %(fname))
			split_fasta(data_folder + fname)
	
	#Run FastQC
	print('Running fastQC for split and original files.')
	run_fastqc(data_folder,'untrimmed_quality')
	
	#Sickle Trim
	print('Starting Sickle to trim split files.')
	run_sickle(data_folder)
	
	#Run FastQC on the trimmed files
	print('Running fastQC for trimmed files.')
	run_fastqc('sickle_trimmed/','trimmed_quality')
	
	#Build Hisat2 reference index
	#print('Building Hisat2 Index')
	#build_hisat2_index(ref_file, 'pipe_ref')
	
	#Run Hisat2
	#run_hisat2(data_folder  + 'split_fastq/sickle_trimmed/')
	
	#Split mapped/unmapped reads
	#run_samtools_split('hisat2_output')
	
	#Build kallisto reference
	print('Building Kallisto index')
	build_kallisto_ref(ref_file)
	
	#Run kallisto (with bootstrapping)
	print('Running Kallisto')
	run_kallisto(data_folder + 'split_fastq/sickle_trimmed/')
	
	
	print('RNA-seq pipeline finished.')
	
	print('Attempting to run Sleuth, this may have to be done manually if an error occurs (Sleu.R).')

if __name__	== "__main__":
	main()
