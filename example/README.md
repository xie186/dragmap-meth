These commands can be run from this directory once `dragmap-meth` is installed

1. Build hash table using refence genome file
2. Align the Reads using dragmap 

```Shell
python ../dragmap-meth.py buildhashtable -r ref.fa -o ref/
python  ../dragmap-meth.py dragmap -ht ref/ -r1 t_R1.fastq.gz -r2 t_R2.fastq.gz |samtools view -bS - -o dragmap-meth.bam
```

Then check the alignments:

```Shell
samtools flagstat dragmap-meth.bam

91537 + 1731 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
91463 + 1731 mapped (99.92%:100.00%)
91537 + 1731 paired in sequencing
45713 + 1035 read1
45824 + 696 read2
90398 + 0 properly paired (98.76%:0.00%)
91401 + 1721 with itself and mate mapped
62 + 10 singletons (0.07%:0.58%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
From here, it is recommended to use [PileOMeth](https://github.com/dpryan79/PileOMeth) for extraction and tabulation of the methylation.

