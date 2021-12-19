#!/usr/bin/env python
#-*-- coding: utf-8 -*-

import argparse
import sys
import os
import re
import gzip
import toolshed
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
import subprocess
from operator import itemgetter
from subprocess import Popen, PIPE, STDOUT
from toolshed import nopen, reader, is_newer_b
from itertools import groupby, repeat, chain, islice
import shutil

try:
    from itertools import izip
    import string
    maketrans = string.maketrans
except ImportError: # python3
    izip = zip
    maketrans = str.maketrans

SUB_HT = "buildhashtable"
SUB_DRAGMAP = "dragmap"
SUB_SUMMARY = "summary"
ID_SPLIT = "__sequence__"
TEMP_FILE = []

__version__ = "0.1.1"

def comp(s, _comp=maketrans('ATCG', 'TAGC')):
    return s.translate(_comp)

## Adapted from bwameth.py
class Bam(object):
    __slots__ = 'read flag chrom pos mapq cigar chrom_mate pos_mate tlen \
            seq qual other converted_read'.split()
    def __init__(self, args):
        for a, v in zip(self.__slots__[:11], args):
            setattr(self, a, v)
        self.other = args[11:]
        self.flag = int(self.flag)
        self.pos = int(self.pos)
        self.tlen = int(float(self.tlen))
        #setattr(self, "converted_read",  self.read)
        self.converted_read = self.read
        self.read = self.read.split(ID_SPLIT)[0]

    def __repr__(self):
        return "Bam({chr}:{start}:{read}".format(chr=self.chrom,
                                                 start=self.pos,
                                                 read=self.read)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__[:11]) \
                         + "\t" + "\t".join(self.other)

    def is_first_read(self):
        return bool(self.flag & 0x40)

    def is_second_read(self):
        return bool(self.flag & 0x80)

    def is_plus_read(self):
        return not (self.flag & 0x10)

    def is_minus_read(self):
        return bool(self.flag & 0x10)

    def is_mapped(self):
        return not (self.flag & 0x4)

    def cigs(self):
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    def cig_len(self):
        return sum(c[0] for c in self.cigs() if c[1] in
                   ("M", "D", "N", "EQ", "X", "P"))

    def left_shift(self):
        left = 0
        for n, cig in self.cigs():
            if cig == "M": break
            if cig == "H":
                left += n
        return left

    def right_shift(self):
        right = 0
        for n, cig in reversed(list(self.cigs())):
            if cig == "M": break
            if cig == "H":
                right += n
        return -right or None

    @property
    def original_seq(self):
        try:
            return self.converted_read.split(ID_SPLIT)[1]
        except:
            sys.stderr.write(repr(self.other) + "\n")
            sys.stderr.write(self.read + "\n")
            raise

    @property
    def ga_ct(self):
        return self.converted_read.split(ID_SPLIT)[2]

    def longest_match(self, patt=re.compile("\d+M")):
        return max(int(x[:-1]) for x in patt.findall(self.cigar))



def byte2str(byte_str):
    """
    To convert byte to str. b"string" -> "string"
    When the input is gzip file, line in each line 
    is byte. 
    """
    if isinstance(byte_str, bytes):
        new_str = byte_str.decode("utf-8")
    return new_str

## source: https://github.com/lh3/readfq/blob/master/readfq.py
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                #l = byte2str(l)
                #logging.warning(l)
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            #l = byte2str(l)
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                #l = byte2str(l)
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def checkX(cmd):
    for p in os.environ['PATH'].split(":"):
        if os.access(os.path.join(p, cmd), os.X_OK):
            break
    else:
        raise Exception("executable for '%s' not found" % cmd)

def wrap(text, width=100): # much faster than textwrap
    try: xrange
    except NameError: xrange = range
    for s in xrange(0, len(text), width):
        yield text[s:s+width]

def run(cmd):
    list(nopen("|%s" % cmd.lstrip("|")))


######################
def fasta_iter(fasta_name):
    fh = nopen(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        yield header, "".join(s.strip() for s in next(faiter)).upper()

def convert_fasta(ref_fasta, just_name=False):
    out_fa = ref_fasta + ".dragmap.c2t"
    if just_name:
        return out_fa
    msg = "Build hash table: c2t in %s to %s" % (ref_fasta, out_fa)
    if is_newer_b(ref_fasta, out_fa):
        logging.warning("already converted: %s\n" % msg)
        return out_fa
    logging.warning("Start to convert the reference.")
    try:
        fh = open(out_fa, "w")
        for header, seq in fasta_iter(ref_fasta):
            header = header.split(" ")[0]
            ########### Reverse ######################
            fh.write(">r%s\n" % header)

            for line in wrap(seq.replace("G", "A")):
                fh.write(line + '\n')

            ########### Forward ######################
            fh.write(">f%s\n" % header)
            for line in wrap(seq.replace("C", "T")):
                fh.write(line + '\n')
        fh.close()
    except:
        try:
            fh.close()
        except UnboundLocalError:
            pass
        os.unlink(out_fa)
        raise
    logging.warning("Done")
    return out_fa

def convert_bed(bed, just_name=False):
    out_bed = bed + ".c2t.bed"
    fh = open(out_bed, "w")
    with open(out_bed, "w") as ob, open(bed, "r") as fin:
        for row in fin:
            row_new = "r"+row
            ob.write(row_new)
    return out_bed

def dragmap_hashtable(options):
    logging.warning("Building dragmap hash table: %s\n" % options.ref_c2t)
    cmd = "dragen-os --build-hash-table true --ht-reference {ref_c2t} --output-directory {outdir}".format(ref_c2t=options.ref_c2t, outdir=options.output_dir)
    if options.mask_bed: 
        bed_c2t = convert_bed(options.mask_bed)
        cmd = cmd + " --ht-mask-bed=" + bed_c2t
    if options.prefix: 
        cmd = cmd + " --output-file-prefix=" + options.prefix
    if options.output_dir:
        try:
            os.makedirs(options.output_dir, exist_ok = True)
            logging.info("Directory '%s' created successfully" % options.output_dir)
        except OSError as error:
            logging.info("Directory '%s' can not be created" % options.output_dir)
    try:
        run(cmd)
        logging.warning("Building dragmap hash table succeeded!")
    except:
        logging.warning("Building dragmap hash table failed!")

# Adapted from: https://stackoverflow.com/a/21529243/3327344
def hook_compressed_text(filename, mode, encoding='utf8'):
    """
    #lines are byte strings and not text string if we use gzip.open by default. 
    """
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode + 't', encoding=encoding)
    #elif ext == '.bz2':
    #    import bz2
    #    return bz2.open(filename, mode + 't', encoding=encoding)
    else:
        return open(filename, mode, encoding=encoding)


def convert_read(fastq_file, conversion):
    fq = hook_compressed_text(fastq_file, 'r')
    file_prefix = os.path.splitext(fastq_file)[0]
    output_file = file_prefix + "." + conversion
    TEMP_FILE.append(output_file)
    read_i = 0 if conversion == "CT" else 1
    #char_a, char_b = ['CT', 'GA'][read_i]
    with open(output_file, "w") as outputf:
        for name, seq, qual in readfq(fq):
            conversion_seq = ""
            if conversion == "CT":
                conversion_seq = seq.replace("C", "T")
            else:
                conversion_seq = seq.replace("G", "A")
            seq_id = name + ID_SPLIT + seq + ID_SPLIT + conversion
            outputf.write("@"+seq_id+"\n"+conversion_seq+"\n+\n"+qual+"\n")
    return output_file


def rm_temp_file():
    """
    Remove intermediate files
    """
    logging.warning("Remove intermediate files: ")
    logging.warning(TEMP_FILE)
    for file_name in TEMP_FILE:
        os.remove(file_name)

## Reference: https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

## Adapterd from bwa-meth
def handle_header(line, cmd, out=sys.stdout):
    toks = line.rstrip().split("\t")
    if toks[0].startswith("@SQ"):
        sq, sn, ln = toks  # @SQ    SN:fchr11    LN:122082543
        # we have f and r, only print out f
        chrom = sn.split(":")[1]
        if chrom.startswith('r'): return
        chrom = chrom[1:]
        toks = ["%s\tSN:%s\t%s" % (sq, chrom, ln)]
    if toks[0].startswith("@PG"):
        out.write("\t".join(toks) + "\n")
        toks = ["@PG\tID:dragmap-meth\tPN:dragmap-meth\tCL:%s" % (" ".join(x.replace("\t", "\\t") for x in sys.argv))]
    out.write("\t".join(toks) + "\n")

## Adapted from bwameth
def handle_reads(alns, set_as_failed=None):

    for aln in alns:
        orig_seq = aln.original_seq
        #assert len(aln.seq) == len(aln.qual), aln.read
        # don't need this any more.
        aln.other = [x for x in aln.other if not x.startswith('YS:Z')]

        # first letter of chrom is 'f' or 'r'
        direction = aln.chrom[0]
        aln.chrom = aln.chrom.lstrip('fr')

        if not aln.is_mapped():
            aln.seq = orig_seq
            continue

        assert direction in 'fr', (direction, aln)
        aln.other.append('YD:Z:' + direction)

        #flag alignments to this strand"
        #    " as not passing QC (0x200). Targetted BS-Seq libraries are often"
        #    " to a single strand, so we can flag them as QC failures. Note"
        #    " f == OT, r == OB. Likely, this will be 'f' as we will expect"
        #    " reads to align to the original-bottom (OB) strand and will flag"
        #    " as failed those aligning to the forward, or original top (OT).
        #if set_as_failed == direction:
        #    aln.flag |= 0x200

        # here we have a heuristic that if the longest match is not 44% of the
        # sequence length, we mark it as failed QC and un-pair it. At the end
        # of the loop we set all members of this pair to be unmapped
        if aln.longest_match() < (len(orig_seq) * 0.44):
            aln.flag |= 0x200  # fail qc
            aln.flag &= (~0x2) # un-pair
            aln.mapq = min(int(aln.mapq), 1)

        mate_direction = aln.chrom_mate[0]
        if mate_direction not in "*=":
            aln.chrom_mate = aln.chrom_mate[1:]

        # adjust the original seq to the cigar
        l, r = aln.left_shift(), aln.right_shift()
        if aln.is_plus_read():
            aln.seq = orig_seq[l:r]
        else:
            aln.seq = comp(orig_seq[::-1][l:r])

    if any(aln.flag & 0x200 for aln in alns):
        for aln in alns:
            aln.flag |= 0x200
            aln.flag &= (~0x2)
    return alns

def dragmap(options):
    #### Convert reads 
    if options.read2 is not None:
        logging.warning("Converting read 1")
    else:
        logging.warning("Converting reads.")
    read1_CT = convert_read(options.read1, "CT")
    logging.warning("Done")
    read2_GA = ""
    if options.read2 is not None:
        logging.warning("Converting read 2...")
        read2_GA = convert_read(options.read2, "GA")
        logging.warning("Done!")
    #starts the pipeline with the program to convert fastqs
    cmd = ""
    
    if not options.read2:  
        cmd = "dragen-os -r {database} -1 {read1} ".format(
                 read1 = read1_CT, database = options.hashtable_dir)
    else:
        cmd = "dragen-os -r {database} -1 {read1} -2 {read2}".format(
                 read1 = read1_CT, database = options.hashtable_dir, read2 = read2_GA)
    cmd = cmd + " --num-threads " + str(options.threads)
    if options.rgid:
        cmd = cmd + " --RGID " + options.rgid + " --RGSM " + options.rgsm
    logging.warning("Start to run: " + cmd)
    sam_iter = execute(cmd)
    for line in sam_iter:
        if not line[0] == "@": break
        if options.header == "T":
            handle_header(line, cmd)
    else:
        sys.stderr.flush()
        raise Exception("bad or empty fastqs")
    sam_iter2 = (x.rstrip().split("\t") for x in chain([line], sam_iter))
    for read_name, pair_list in groupby(sam_iter2, itemgetter(0)):
        pair_list = [Bam(toks) for toks in pair_list]
        set_as_failed = None  
        for aln in handle_reads(pair_list, set_as_failed):
            sys.stdout.write(str(aln) + '\n')
    if not options.keep_temp:
        rm_temp_file()
    logging.warning("Finished!")



if __name__ == '__main__':
    ## description - Text to display before the argument help (default: none)   
    parser=argparse.ArgumentParser(description='dragmap-meth: bisulfite reads mapper using dragmap-meth') 
    subparsers = parser.add_subparsers(help='sub-command help', dest = 'command')

    ## Parameters for building hash table
    parser_buildhashtable = subparsers.add_parser(SUB_HT, help='Build hash table of a reference genome for alignent')
        #flag alignments to this strand"
        #    " as not passing QC (0x200). Targetted BS-Seq libraries are often"
    parser_buildhashtable.add_argument('-r', '--reference', metavar='reference', \
                      help='Reference in fasta format', required=True)
    parser_buildhashtable.add_argument('-o', '--output-directory', dest='output_dir', \
                      help='Output directory ', required=True)
    parser_buildhashtable.add_argument('-p', '--prefix', metavar='prefix', \
                      help='Prefix of Hash table files ', required=False)
    parser_buildhashtable.add_argument('-m', '--mask-bed', metavar='mask_bed', \
                      help = 'Bed file for masking the reference genome', required=False)

    ## Parameters for aligning reads using dragmap 
    parse_DRAGMAP = subparsers.add_parser(SUB_DRAGMAP, help='Align bisulfite reads using dragmap')
    parse_DRAGMAP.add_argument('-ht', '--hashtable', dest='hashtable_dir', \
                      help='Direcotory of hashtable ', required=True)

    parse_DRAGMAP.add_argument('-r1', '--read1', metavar='read1', \
                      help='Read1 in fastq format. Provide only -r1 if you have single end reads.', required=True)
    parse_DRAGMAP.add_argument('-r2', '--read2', metavar='read1', \
                      help='Read2 in fastq format', required=False)
    
    parse_DRAGMAP.add_argument('--RGID', dest = 'rgid', \
                      help='Read Group ID', required=False)
    parse_DRAGMAP.add_argument('--RGSM', dest = 'rgsm', \
                      help='Read Group Sample', required=False)
 
    #-num_threads 
    parse_DRAGMAP.add_argument('-t', '--threads', type=int, \
                      help='Number of threads. Default: 1', default=1)
    parse_DRAGMAP.add_argument('--header', type=str, \
                      help='Whether to output header in SAM output file. F means Faslse. T means True. Default: T', default="T", choices=("F", "T"))

    parse_DRAGMAP.add_argument('--keep_temp', dest='keep_temp', action='store_true', 
                           help="Keep temporary files. [Default: OFF]")
    parse_DRAGMAP.set_defaults(keep_temp=False)
    options = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    dict_cmd = vars(options)
    logging.warning("Subcommand is: " + dict_cmd['command'])
    if len([x for x in (options.rgid, options.rgsm) if x is not None]) == 1:
        parser.error('--RGID and --RGSM must be given together')

    ## Build hash table
    if dict_cmd['command'] == SUB_HT:
        ## return c2t converted reference genome 
        ref_c2t = convert_fasta(options.reference)
        options.ref_c2t = ref_c2t
        print(options) 
        dragmap_hashtable(options)
    if dict_cmd['command'] == SUB_DRAGMAP:
        #logging.warning(options.keep_temp)
        dragmap(options)
