# dragmap_meth (dragmap_meth.py)

Alignment of BS-Seq reads using [dragmap](https://github.com/Illumina/DRAGMAP). 

## Intro

This works for single-end reads and for **paired-end reads from the
directional protocol** (most common).

Uses the method employed by methylcoder and Bismark of *in silico*
conversion of all C's to T's in both reference and reads.

Recovers the original read (needed to tabulate methylation).

It allows us to align reads to nt database. 

## QuickStart

The commands:

```
python dragmap-meth.py buildhashtable -r ref.fa -o ref/
python  dragmap-meth.py dragmap -ht ref/ -r1 t_R1.fastq.gz -r2 t_R2.fastq.gz |samtools view -bS - -o dragmap-meth.bam
```

will create `dragmap-meth.bam`. 
To align single end-reads, specify only 1 file: `-r1 some_read.fastq.gz`

## Installation

```
conda env create -n env4dragmap-meth --file environment.yaml python=3
```

### Dependencies

```
Python 3.10.1
dragmap
```

## Acknowledgement

Special thanks to [bwa-meth](https://github.com/brentp/bwa-meth) because part of the codes were adapted from bwa-meth. 

