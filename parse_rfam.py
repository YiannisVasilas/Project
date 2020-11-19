#!/usr/bin/env python3

"""
Author: Vasilas Ioannis 
Script to calculate the number of base pairs in the
RNA secondary structures stored in the Rfam database
"""

from sys import argv

def record_finder(lines):
    """Yield lines from a single record of file (Rfam format)

    lines: list of lines or open file object
    """
    curr = []
    for line in lines:
        if not line.stricd p):
            continue
        if line.startswith('//'):
            yield curr
            curr = []
        else:
            curr.append(line.strip())
    if curr:
        yield curr

def parse_record(rec_lines):
    """Return tuple of (ID, description, structure) for a single record

    rec_lines: list of strings, lines belonging to a single Rfam record
    """
    for line in rec_lines:
        if line.startswith('#=GF ID'):
            rec_id = line.strip().split(' ',2)[-1].strip()
        elif line.startswith('#=GF DE'):
            rec_descr = line.strip().split(' ',2)[-1].strip()
        elif 'SS_cons' in line:
            rec_struct = line.strip().split(' ',2)[-1].strip()
    return (rec_id, rec_descr, rec_struct)
    
def calc_num_basepairs(struct_str):
    """Return the number of basepairs (int) in secondary structure

    struct_str: string, secondary structure in Stockholm format
    
    Pairs are indicated by {}, <>, (), and []
    Unpaired bases are indicated by _-:, and ~
    """
    num_open = 0 
    num_close = 0
    for sym in struct_str:
        if sym in '<{[(':
            num_open += 1
        elif sym in '>}])':
            num_close += 1
        else:
            pass
    assert num_open == num_close
    return num_open


if __name__ == "__main__":

    result = [] #data structure to store the results
    # open the file
    with open(argv[1]) as fo:
        # iterate over all records
        for record in record_finder(fo):
            # get the relevant parts
            label, desc, struct = parse_record(record)
            # perform some calculation
            num_bps = calc_num_basepairs(struct)
            result.append((num_bps, label, desc))
    # print the results
    result.sort()
    for (num_bps, label, descr) in result:
        print('{0}\t{1}\t{2}\t'.format(num_bps, label, descr))

