import sys,getopt,os,subprocess
from multiprocessing import Pool
from functools import partial
from collections import OrderedDict

'''
Make sure the following programs are available:

Cutadapt:
pip install --user --upgrade cutadapt
export PATH=$PATH:~/.local/bin/

Cufflinks
export PATH=$PATH:/local/prog/cufflinks/
'''

# Define options for -h, -i & -o
def main(argv):
  input_file = ''
  output_file = ''
  
  try:
    opts, args = getopt.getopt(sys.argv[1:], "h:i:g:a:o:", ["help", "input_file=", "genome=", "annotation="])
  except getopt.GetoptError:
    print '\nUsage: \n$ python geneious_pipline -i <input_file, input_file, ...> \n'
    sys.exit(2)
    
  for opt, arg in opts:
    if opt == '-h':
      print '\nUsage: \n$ python geneious_pipeline -i <input_file,input_file,...> -o <output_file> \n'
      sys.exit()
    elif opt == "-i":
      inputfile = arg.split(",")
    elif opt == "-g":
      if arg.endswith(".fasta"):
          try:
             os.rename(arg, arg.split(".fasta")[0] + ".fa")
          except:
             pass
          genome = arg.split(".fasta")[0] + ".fa"
      else:
          genome = arg
    elif opt == "-a":
      annotation = arg

  return inputfile, genome, annotation

def run_cutadapt(RNAseq_file):
  """Run cutadapt program on fastq files

  RNAseq_file: string, filename of input RNA fastq file
  outfile: string, filename of output filtered RNA sequence file
  """
  outfile = "paper_filter_" + str(RNAseq_file)  
  if not os.path.exists(outfile):
      cmd = 'cutadapt %s -m %s -q %s -o %s '\
          %(RNAseq_file, 30, 20, outfile)
      subprocess.check_call(cmd, shell=True)
  return outfile

def make_index(genome,annotation):
  """Run gffread and hisat2 program on gff3 file to build an index

  genome: string, filename of reference genome
  annotation: string, filename of annotation file
  build_outfile: string, filename of build created
  gtf_outfile: string, filename of gtf_annotation
  """
  gtf_outfile = "paper_" + annotation.split(".gff3")[0] + ".gtf"
  if not os.path.exists(gtf_outfile):
      cmd = 'gffread %s -T -o %s '\
          %(annotation,gtf_outfile)
      subprocess.check_call(cmd, shell=True)
  build_outfile = genome.split(".fa")[0]
  if not os.path.exists(build_outfile + ".1.ebwt"):
      cmd = 'bowtie-build %s %s'\
          %(genome, build_outfile)
      subprocess.check_call(cmd, shell=True)
  return build_outfile, gtf_outfile

def run_tophat(RNAseq_file, index, annotation):
  """Run tophat program on fastq files

  RNAseq_files: string, filename of input RNA fastq files
  index: string, filename of index created in make_index
  annotation: string, filename of annotation created in make_index
  """
  outfolder = RNAseq_file.split(".fastq")[0] + "_tophat"
  if not os.path.exists(outfolder):
      cmd = '/local/prog/tophat/tophat -o %s -i %s -I %s --segment-length %s --GTF %s %s %s'\
          %(outfolder, 5, 19000, 15, annotation, index, RNAseq_file)
      subprocess.check_call(cmd, shell=True)
  return outfolder

def run_cufflinks(bamfolder, txt_file):
  """Run cufflinks program on bam files

  bamfolder: list, foldername containing bamfile
  txt_file: string, filename where to append all gtf output locations
  """

  txt_file = open(txt_file, 'a')
  outfile = "transcripts.gtf"
  if not os.path.exists(bamfolder + "/" + outfile):
      cmd = 'cufflinks accepted_hits.bam'
      subprocess.check_call(cmd, shell=True, cwd=bamfolder)
  txt_file.write(bamfolder + "/" + outfile + "\n")
  txt_file.close()

def merge_cufflinks(gtf_file, annotation):
  """Run cuffmerge program on gtf files

  gtf_file: string, filename of txt file containing all locations of gft_files
  annotation: string, filename of annotation created in make_index
  """

  outfile = "merged_asm"
  if not os.path.exists(outfile):
      cmd = 'cuffmerge -g %s %s'\
          %(annotation, gtf_file)
      subprocess.check_call(cmd, shell=True)

def retrieve_bam(outfolder):
  """Copy bam files to current folder for use later

  outfolder: string, foldername containing bam file
  """

  outfile = outfolder + ".bam"
  if not os.path.exists(outfile):
      cmd = 'cp -a ./%s/*hits.bam %s'\
          %(outfolder, outfile)
      subprocess.check_call(cmd, shell=True)
  return outfile

def retrieve_gtf():
  """Copy gtf file to current folder for use later
  """

  out_gtf = "merged.gtf"
  if not os.path.exists(out_gtf):
      cmd = 'cp -a ./merged_asm/merged.gtf %s'\
          %(out_gtf)
      subprocess.check_call(cmd, shell=True)
  return out_gtf

def run_cuffdiff(merged_gtf, bam_files):
  """Run cuffdiff program on bam files

  merged_gtf: string, filename merged gtf from all samples
  bam_files: list, filenames of all bam files
  """

  outfile = "genes.fpkm_tracking"
  input_bam = ' '.join(bam_files)
  if not os.path.exists(outfile):
      cmd = 'cuffdiff %s %s'\
          %(merged_gtf, input_bam)
      subprocess.check_call(cmd, shell=True)
  return outfile

def create_endTable(inputfile, fpkm_final):
  """Create abundance table from the calculated FPKM

  inputfile: list, filenames of input RNA fastq file
  fpkm_final: string, filename of fpkm file created in run_cuffdiff
  abundance_table.tsv: tab-delimeted file containing all FPKM values
  """

  data_file = open('paper_abundance_table.tsv', 'w')
  column_names = '\t'.join(inputfile)
  data_file.write('Sample\t' + column_names + "\n")
  endTable = {}
  for line in open(fpkm_final):
      line = line.strip()
      if "CRO_" in line:
          gene_id = line.split("\t")[4]
          for x in range(len(inputfile)):
             column = 9 + x * 4
             RPKM = line.split("\t")[column]
             if gene_id in endTable:
                endTable[gene_id].append(RPKM)
             else:
                endTable[gene_id] = [RPKM]
  Sorted_endTable = OrderedDict(sorted(endTable.items(), key=lambda t: t[0]))
  for gene in Sorted_endTable:
      data_file.write(gene + "\t" + '\t'.join(Sorted_endTable[gene])+ "\n")
  data_file.close()

if __name__ == "__main__":
  pool = Pool()
  inputfile,genome,annotation = main(sys.argv[1:])
  filtered_RNAseq = pool.map(run_cutadapt, inputfile)
  pool.close()
  pool.join()

  index,gtf_annotation = make_index(genome,annotation)

  pool = Pool()
  partial_run_tophat = partial(run_tophat, index=index, annotation=gtf_annotation)
  bam_outfolders = pool.map(partial_run_tophat, filtered_RNAseq)
  pool.close()
  pool.join()
 
  cufflink_out = 'cufflink_gtf.txt'
  open(cufflink_out, 'w')

  pool = Pool()
  partial_run_cufflinks = partial(run_cufflinks, txt_file=cufflink_out)
  pool.map(partial_run_cufflinks, bam_outfolders)
  pool.close()
  pool.join()

  merge_cufflinks(cufflink_out, gtf_annotation)


  pool = Pool()
  bam_files = pool.map(retrieve_bam, bam_outfolders)
  pool.close()
  pool.join()

  gtf_file = retrieve_gtf()

  fpkm_file = run_cuffdiff(gtf_file, bam_files)

  create_endTable(inputfile, fpkm_file)
