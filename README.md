# dragmap_meth (dragmap_meth.py)

Alignment of BS-Seq reads using [dragmap](https://github.com/Illumina/DRAGMAP). 

## Intro

This works for single-end reads and for **paired-end reads from the
directional protocol** (most common).

## Installation

### Install `dragmap-meth` via `pip`

`dragmap-meth` is available on https://pypi.org/project/dragmap-meth/

```
pip install dragmap-meth
```

## Install `dragmap-meth` via `conda`

```
conda env create -n env4dragmap-meth --file environment.yaml python=3
git clone https://github.com/xie186/dragmap-meth.git
cd dragmap-meth/
```

## QuickStart

The commands:

```
python dragmap-meth.py buildhashtable -r ref.fa -o ref/
python  dragmap-meth.py dragmap -ht ref/ -r1 t_R1.fastq.gz -r2 t_R2.fastq.gz |samtools view -bS - -o dragmap-meth.bam
```

will create `dragmap-meth.bam`. 
To align single end-reads, specify only 1 file: `-r1 some_read.fastq.gz`

### Dependencies

```
Python 3.10.1
dragmap
```

## Reference:

Introducing DRAGMAP, the new genome mapper in DRAGEN-GATK: https://gatk.broadinstitute.org/hc/en-us/articles/4410953761563-Introducing-DRAGMAP-the-new-genome-mapper-in-DRAGEN-GATK

Demystifying the versions of GRCh38/hg38 Reference Genomes, how they are used in DRAGEN™ and their impact on accuracy: https://www.illumina.com/science/genomics-research/articles/dragen-demystifying-reference-genomes.html

https://github.com/DavidStreid/bio_docker/tree/main/dragmap

https://github.com/HudsonAlpha/CSL_public_benchmark

## Acknowledgement

Special thanks to [bwa-meth](https://github.com/brentp/bwa-meth). Part of the codes were adapted from bwa-meth. 

