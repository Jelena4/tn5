
rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(Rsamtools)

################## Parameters to be changed #####################

organism = 'Nv'
sample = '20dpf'
tn5Dir = '/data/personal_folders/jscepanovic/tn5'
genomeIndex = '/data/genomes/v_TRAPARSED_iter5'

#################################################################

setwd(tn5Dir)
prefix = 'tn5_' %>% 
  paste0(organism) %>% 
  paste(sample, sep = '_')
fastq.files = list.files('fastq')
i = fastq.files %>% 
  str_split('_', simplify = T) %>% 
  as.data.frame() %>% 
  select(2) == organism
fastq.files = fastq.files[i] %>% str_subset('_R')
read1.file = fastq.files %>% str_subset('R1')
read2.file = fastq.files %>% str_subset('R2')
read1.prefix = (read1.file %>% str_split('.fastq.gz') %>% unlist)[1]
read2.prefix = (read2.file %>% str_split('.fastq.gz') %>% unlist)[1]

############# Running the mappingPipeline.sh script ############# 

command = "screen -d -m -S " %>% 
  paste0(prefix) %>% 
  paste(' /home/jscepanovic/git/tn5/sandbox/mappingPipeline.sh ') %>% 
  paste0(read1.prefix) %>% 
  paste(read2.prefix, sep = ' ') %>% 
  paste(genomeIndex, sep = ' ') %>% 
  paste(prefix, sep = ' ')

system(command, wait = TRUE)

############ Removing dupcicates from read 1 file ###############

# when the command is done:

source('/home/jscepanovic/git/tn5/sandbox/rmdupUsingPE.R')

############### Counting read starts with bbmap #################

### forward strand
strand = 'forward'
prefix.strand = prefix %>% paste0('.R1.') %>% paste0(strand)
srt = prefix.strand %>% paste0('.filtR.srt.sam')

command = 'samtools sort -@ 6 -o ' %>% paste0(srt) %>% paste(fileName)
system(command, wait = T)

command = '~/bbmap/pileup.sh in=' %>% paste0(srt) %>% paste('basecov=', sep = ' ') %>% paste0(prefix.strand) %>% paste0('.startCount.txt addfromreads=t startcov=t 32bit=t')
system(command, wait = T)


### reverse strand
strand = 'reverse'
prefix.strand = prefix %>% paste0('.R1.') %>% paste0(strand)
srt = prefix.strand %>% paste0('.filtR.srt.sam')

command = 'samtools sort -@ 6 -o ' %>% paste0(srt) %>% paste(fileName)
system(command, wait = T)

command = '~/bbmap/pileup.sh in=' %>% paste0(srt) %>% paste('basecov=', sep = ' ') %>% paste0(prefix.strand) %>% paste0('.startCount.txt addfromreads=t startcov=t 32bit=t')
system(command, wait = T)

######################## Peak calling ##########################

### forward strand
strand = 'forward'
prefix.strand = prefix %>% paste0('.R1.') %>% paste0(strand)
fileName = prefix.strand %>% paste0('.startCount.txt')
readStartCount = fread(fileName, header = T, stringsAsFactors = F, data.table = F)
names(readStartCount) = c('Chr', 'Position', 'Coverage')
idx = which(readStartCount[,3] > 5)
idx_up = idx + 1
idx_down = idx - 1
peaks = (readStartCount[idx,3] >  readStartCount[idx_up,3]) & (readStartCount[idx,3] >  readStartCount[idx_up,3])
idx_peaks = idx[peaks]
pos = readStartCount[idx_peaks,]
outFile = prefix.strand %>% paste0('.peaks.txt')
write.table(pos, outFile, sep='\t', quote=F, row.names=F, col.names = T)
bed = pos[,1:2] %>% cbind(pos[,2]+1)
outBed = prefix.strand %>% paste0('.peaks.bed')
write.table(bed, outBed, sep='\t', quote=F, row.names=F, col.names = F)

### reverse strand
strand = 'reverse'
prefix.strand = prefix %>% paste0('.R1.') %>% paste0(strand)
fileName = prefix.strand %>% paste0('.startCount.txt')
readStartCount = fread(fileName, header = T, stringsAsFactors = F, data.table = F)
names(readStartCount) = c('Chr', 'Position', 'Coverage')
idx = which(readStartCount[,3] > 5)
idx_up = idx + 1
idx_down = idx - 1
peaks = (readStartCount[idx,3] >  readStartCount[idx_up,3]) & (readStartCount[idx,3] >  readStartCount[idx_up,3])
idx_peaks = idx[peaks]
pos = readStartCount[idx_peaks,]
outFile = prefix.strand %>% paste0('.peaks.txt')
write.table(pos, outFile, sep='\t', quote=F, row.names=F, col.names = T)
bed = pos[,1:2] %>% cbind(pos[,2]+1)
outBed = prefix.strand %>% paste0('.peaks.bed')
write.table(bed, outBed, sep='\t', quote=F, row.names=F, col.names = F)


