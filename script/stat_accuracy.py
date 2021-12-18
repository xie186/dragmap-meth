from Bio import TogoWS
import argparse
import sys
import os

def summary(options):
    with open(options.input) as file_input:
        for line in file_input:
            line = line.rstrip()
            ele = line.split("\t")
            if ele[2] == "*": 
                continue
            data = "FAILED"
            attempts = 0
            while attempts < 3:
                try:
                    handle = TogoWS.entry("nucleotide", ele[2], field=options.fields)
                    data = handle.read().strip()  # ignore trailing \n
                    handle.close()
                    break
                except:
                    attempts += 1
                    pass
            print("{read_id}\t{hit}\t{species}".format(
            read_id = ele[0], hit=ele[2], species=data))
                #raise("Makeblastdb failed!")

if __name__ == '__main__':
    ## description - Text to display before the argument help (default: none)   
    parser=argparse.ArgumentParser(description='mbmeth') 
    parser.add_argument("-i", '--input', help="Input list")
    parser.add_argument("-f", '--fields', help="Fields. default: origanism", default="organism",
                           choices=("organism", "taxonomy")) 
    options = parser.parse_args(args=None if sys.argv[1:] else ['--help']) 
    summary(options) 
