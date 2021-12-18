from Bio import TogoWS
import argparse
import sys
import os

def summary(options):
    with open(options.input) as file_input:
        for line in file_input:
            line = line.rstrip()
            ele = line.split("\t")
            handle = TogoWS.entry("nucleotide", ele[2], field="organism")
            if ele[3] == "*": 
                continue
            data = handle.read().strip()  # ignore trailing \n
            handle.close()
            print("{read_id}\t{hit}\t{species}".format(
                     read_id = ele[0], hit=ele[2], species=data
                    ))

if __name__ == '__main__':
    ## description - Text to display before the argument help (default: none)   
    parser=argparse.ArgumentParser(description='mbmeth') 
    parser.add_argument("-i", '--input', help="Input list")
    
    options = parser.parse_args(args=None if sys.argv[1:] else ['--help']) 
    summary(options) 
