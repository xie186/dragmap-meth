from Bio import TogoWS
import argparse
import sys
import os

def summary(options):
    num_reads = 0
    num_correct = 0
    with open(options.input) as file_input:
        for line in file_input:
            line = line.rstrip()
            ele = line.split("\t")
            if "FAILED" in line:
                continue
            if "BAD" in ele[0]:
                num_reads += 1
                if "Bacteria" in line:
                    num_correct += 1
                #raise("Makeblastdb failed!")
            else:
                num_reads += 1
                if options.species in line:
                    num_correct += 1
                #raise("Makeblastdb failed!")

    percentage = 100* num_correct/num_reads
    print("Percente: {perc}\n".format(perc = percentage))

if __name__ == '__main__':
    ## description - Text to display before the argument help (default: none)
    parser=argparse.ArgumentParser(description='mbmeth')
    parser.add_argument("-i", '--input', help="Input list")
    parser.add_argument("-s", '--species', help="species")
    options = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    summary(options)



