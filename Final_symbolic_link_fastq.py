#!usr/bin/env python3
"""
@authors: Jasper, Anna and Marijn Schipper
"""

import sys
import subprocess
import os

def folder_builder():
    """Makes folder named fastq_links
    """
    if not os.path.exists('./fastq_links/'):
        cmd = "mkdir fastq_links"
        subprocess.check_call(cmd, shell = True)
    return

def list_files_in_directory(path_to_folder):
    '''Return a list of all fastq files in the path specified
    
    Input:
        path_to_folder: string -- the absolute path to a folder
    Output:
        files: list -- list of filenames of all fastq files
    '''
    cmd = "ls {}".format(path_to_folder)
    files = subprocess.check_output(cmd, shell = True).decode("utf-8")\
        .split("\n")
    files = [name for name in files if not name.startswith("JE")]
    files = [name for name in files if name]
    return files

def parse_tsv(infile):
    '''Return a nested dictionary with the sample distribution
    
    Input:
        infile: open file or list of lines of the tsv overview
    Output:
        dict_folder: nested dictionary {folder : {index_seq: (sample, time)}}
    '''
    dict_folder = {}
    identifiers = ["Mock1", "Mock2", "Sc", "Cf"]
    headers = infile.readline()
    for line in infile:
        line_parts = line.strip().split("\t")
        folder = line_parts[0][0]+line_parts[0][1]+"0"+line_parts[0][2]
        seq_index = line_parts[1]
        sample = line_parts[3]
        repl = line_parts[4]
        time = line_parts[5]
        if folder not in dict_folder:
            dict_folder[folder] = {}
        if sample in identifiers:
            dict_folder[folder][seq_index] = (sample, repl, time)
    return dict_folder

def select_file(folder_dict):
    '''Yield the path and filename of the files of interest
    
    Input:
        folder_dict: dict -- nested dict {folder : {index_seq: (sample, time)}}
    Output:
        tuple -- (path, filename)
    '''
    path = "/local/data/course/project/falcarindiol_data/transcriptome/"
    subfolders = list(folder_structure.keys())
    for folder in subfolders:
        files = list_files_in_directory(path+folder+"/")
        for fastq_file in files:
            path_to_file = path+folder+"/"+fastq_file
            file_parts = fastq_file.split(".")[0].split("_")
            if (file_parts[5] in folder_structure[folder]) and \
                (file_parts[7] == 'pf'):
                yield (path_to_file, folder, fastq_file)

def make_symbollic_link(filename, folder, path_to_file, folder_dict):
    '''Make a soft link of an file
    
    Input:
        filename: string -- name of the file to make the link
        path_to_file: string -- absolute path to the file
    Output:
        A soft link under a new name in the folder fastq_links
    The name of the links consists of treatment_identifier_hpi_pair.fastq.gz
    '''
    file_parts = filename.split(".")[0].split("_")
    sample, repl, hpi = folder_structure[folder][file_parts[5]]
    link_name = "{}_{}_{}_{}.fastq.gz  "\
        .format(sample,repl, hpi, file_parts[6])
    if not os.path.exists(link_name):
        cmd = "ln -s {} ./fastq_links/{}".format(path_to_file, link_name)
        subprocess.check_call(cmd, shell = True)

if __name__ == "__main__":
    #step 1: read command line arguments
    summary_files = sys.argv[1]
    
    #step 1.5: make output directory
    folder_builder()
    
    #step 2: parse the tsv file to identify the folder structure
    with open(summary_files) as fo:
        folder_structure = parse_tsv(fo)
    
    #step3: make soft links of the selected files    
    for (path, folder, filename) in select_file(folder_structure):
        make_symbollic_link(filename, folder, path, folder_structure)

