#!/usr/bin/env python3

"""

This python script has to parse a GenBank file and outputs a FASTA file 
and an ordered table with detail statistics
"""

from sys import argv

def record_finder_yield(lines):
    """Return list of records separated by //
    lines: open file or list of lines"""
    curr = []
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("//"):
            if curr:
                yield curr
                curr = []
        else:
            curr.append(line.strip())
    if curr:
        yield curr

def parse_data(rec_lines):
    """Return tuple of (accession ID, organism, sequence) """
    data = []
    for lists in rec_lines:
        seqstart = False
        seq = ""
        for line in lists:
            if line.startswith("ACCESSION"):
                accession = line.split()[1]
            if line.startswith("ORGANISM"):
                organism = line.partition(" ")[2].strip()
            if line.startswith("ORIGIN"):
                seqstart = True
            if seqstart == True:
                seq += "".join(line.split()[1:])
            else:
                pass
        data.append((accession, organism, seq))
    return data
    
def gc_content(seq):
        n = len(seq)
        G = 0
        C = 0

        for base in seq:
                if base in 'Gg':
                        G += 1
                elif base in 'Cc':
                        C += 1

        return 100 * (float(G + C) / n)

def sort_seq(data):
    """
    make a tuple of (accession ID, organism, sequence, GCcontant)
    return the tuple sorted by GC contant"""
    record_gc = []
    for i in range(len(data)):
        lists = list(data[i])
        lists.append(gc_content(data[i][2]))
        record_gc.append(lists)
    sorted_data = sorted(record_gc, key=lambda x: x[3], reverse=True)
    return sorted_data
    
def output_fasta(data):
    """sorted tuple with GC contant"""
    file = open("seq.fasta","w")
    for lists in data:
            output = ">" + lists[0] + " " + lists[1] + "\n" + lists[2] + "\n"
            file.write(output)
    file.close()
    
def output_report(data):
    """a sorted tuple with GC contant"""
    file = open("report.txt","w")
    for lists in data:
        out = [lists[0],lists[1],str('%.2f' % lists[3]),\
        str(len(lists[2]))]
        output = "\t".join(out) + "\n"
        file.write(output)
    file.close

if __name__ == "__main__":
    
    # parse input data
    inp_fn = argv[1]
    rec_lines = record_finder_yield(open(inp_fn))
    data = parse_data(rec_lines)
    sorted_data = sort_seq(data)
    output_fasta(sorted_data)
    output_report(sorted_data)
