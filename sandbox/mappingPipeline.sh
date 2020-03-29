#!/bin/bash
  
######                                                    ######
######             PROCESSING OF TN5 DATA                 ######
######                                                    ######


### Arguments:

# $1 prefix of the raw read1 fastq file (up to .fastq.gz)
# $2 prefix of the raw read2 fastq file (up to .fastq.gz)
# $3 path to index file, only the prefix (e.g. /data/genomes/v_TRAPARSED_iter5) 
# $4 sample prefix (name for writing files)
# $5
# $6

### What the script does
# Trims read1, maps read1 onto the genome in single end mode; maps both reads in paired end mode
# Converts to BAM, filters MAPQ 2 (multimapers), removes duplicates from paired end with picard 
# Separates forward and reverse strands from read1

### Trimming
# The adapters we used are dedicated reads so no need for trimming them.
# TSO barcode sequence (8 nt long) needs to be trimmed from read1.
seqtk trimfq -b 8 fastq/$1.fastq.gz | gzip - > fastq/$4.R1.trimmed.fastq.gz

### Mapping
# Read1 single end mode
bowtie2 -D 200 -R 3 -N 1 -L 20 -i S,1,0.50 --gbar 1 --no-unal -p 8 -x $3 -U fastq/$4.R1.trimmed.fastq.gz -S $4.R1.sam 2> $4.R1.mappingstats
# Paired end mode 
bowtie2 -D 200 -R 3 -N 1 -L 20 -i S,1,0.50 -p 8 -x $3 -1 fastq/$4.R1.trimmed.fastq.gz -2 fastq/$2.fastq.gz -S $4.PE.sam 2> $4.PE.mappingstats

### Sam to Bam and filtering
samtools view -bq 2 -@ 6 -o $4.R1.bam $4.R1.sam
samtools view -bq 2 -@ 6 -o $4.PE.bam $4.PE.sam
rm $4.R1.sam
rm $4.PE.sam

### Separating forward and reverse strand for read1
#samtools view -b -F 0x10 $4.R1.bam > $4.R1.forward.bam
#samtools view -b -f 0x10 $4.R1.bam > $4.R1.reverse.bam
samtools view -@ 6 $4.R1.forward.bam > $4.R1.forward.sam
samtools view -@ 6 $4.R1.reverse.bam > $4.R1.reverse.sam

### Removing duplicates from PE
java -XX:ParallelGCThreads=6 -jar /data/personal_folders/jscepanovic/picard.jar MarkDuplicates I=$4.PE.bam O=$4.rmdup.bam M=$4.marked.dup.metrics.txt REMOVE_DUPLICATES=true
samtools view -@ 6 $4.rmdup.bam > $4.rmdup.sam
