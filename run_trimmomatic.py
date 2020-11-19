#!/usr/bin/evn python
""" Author: The_Parsers (Vasilas Ioannis )
1/12/2017
Usage:
    python trimming.py runs the trimmomatic
"""


from sys import argv
import subprocess

def run_trimmomatic(filename, output):
    cmd = "java -jar /local/prog/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 {} {} MINLEN:20".format(filename, output)
    res = subprocess.check_call(cmd, shell=True)
    return res

if __name__ == "__main__":
    list_of_files = argv[1:]
    for filename in list_of_files:
        elements = filename.split(".")
        output = elements[0] + "_MINLENGTH20." + elements[1]
        run_trimmomatic(filename, output)
        
