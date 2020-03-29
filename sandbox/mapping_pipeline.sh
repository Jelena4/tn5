#!/bin/bash

######                                                    ######
######		   PROCESSING OF TN5 DATA		  ######
######                                                    ######


### Arguments:
 
# $1 prefix of the raw read1 fastq file (up to .fastq.gz) (/data/nobel_backup/DATA/RAW/hic_tn5_2019_04/tn5_Nv_20dpf_S5_R1_001_trimmed.fastq.gz)
# $2 prefix of the raw read2 fastq file (up to .fastq.gz)
# $3 path to index file, only the prefix (/data/genomes/v_TRAPARSED_iter5)
# $4 
# $5 Genome size (approximate) 
# $6
# $7

### What the script does
# Trimmas read1, maps read1 onto the genome in single end mode; maps both reads in paired end mode
# Converts to BAM, filters MAPQ 2, removes duplicates with picard 

### Trimming
# The adapters we used are dedicated reads so no need for trimming them.
# TSO barcode sequence () needs to be trimmed from read1.

#cutadapt -a sequence-o $1_trimmed.fastq.gz $1.fastq.gz

### Mapping
bowtie2 -p 8 -x $3 -U $1 2> $1.mappingstats
samtools view -b $1.bam $1.sam

