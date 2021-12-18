from Bio import TogoWS
import argparse
import sys
import os
from collections import defaultdict

def stat_read(dict_read):
    for key in dict_read:
        dict_num = dict_read[key]
        for key_sec in dict_num:
            #print("{key}\t{key_sec}\t{num}".format(key=key, key_sec=key_sec, num = dict_num[key_sec]))
            pass

def summary(options):
    num_reads = 0
    num_correct = 0
    num_read_ecoli = 0
    num_correct_ecoli = 0
    dict_read = defaultdict(dict)
    with open(options.input) as file_input:
        for line in file_input:
            line = line.rstrip()
            ele = line.split("\t")
            ele[2].replace("organisms; ", "")
            info = ele[2].split(";")[options.taxonomy]            
            if "FAILED" in line:
                continue
            print(info)
            if not dict_read[ele[0]]:
                dict_read[ele[0]] = {}
            if info not in  dict_read[ele[0]]:
                dict_read[ele[0]][info] = 0
            if "BAD" in ele[0]:
                num_read_ecoli += 1
                if options.contamination in line:
                    num_correct_ecoli += 1
                #raise("Makeblastdb failed!")
                dict_read[ele[0]][info] +=1
            else:
                num_reads += 1
                if options.species in line:
                    num_correct += 1
                dict_read[ele[0]][info] +=1
                #raise("Makeblastdb failed!")
    percentage = 100* num_correct/num_reads
    perc2  = 100*num_correct_ecoli/num_read_ecoli
    perc3  = 100*(num_correct + num_correct_ecoli) / (num_reads + num_read_ecoli)
    print("Percente (reference): {perc}\n".format(perc = percentage))
    print("Percente (E coli): {perc}\n".format(perc=perc2))
    print("Percente (all): {perc}".format(perc=perc3))
    if options.collapse:
        stat_read(dict_read)

if __name__ == '__main__':
    ## description - Text to display before the argument help (default: none)
    parser=argparse.ArgumentParser(description='mbmeth')
    parser.add_argument("-i", '--input', help="Input list")
    parser.add_argument("-c", '--contamination', help="contamination", default="Escherichia coli")
    parser.add_argument("-s", '--species', help="species", default="Zea mays")
    parser.add_argument('--collapse', help="Collapse", type=bool, default=False)
    parser.add_argument('--taxonomy', help="Which taxonomy to take", type=int, default=0)
    options = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    summary(options)



