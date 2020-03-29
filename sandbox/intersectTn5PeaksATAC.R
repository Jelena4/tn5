
rm(list = ls())

library(data.table)
library(GenomicRanges)

tn5Dir = '/data/personal_folders/jscepanovic/tn5/tn5PeakCalling/'
atacDir = '/data/personal_folders/jscepanovic/SK/ATAC-seq/'

setwd(tn5Dir)

organism = 'Sk'
sample = 'Gill'
strand = 'reverse'
params = 'param3'

tn5 = file.path(tn5Dir, 'tn5_' %>% 
                  paste0(organism) %>% 
                  paste(sample, sep = "_") %>% 
                  paste(strand, sep = '_') %>% 
                  paste0('_peaks.bed'))
tn5 = fread(tn5)

idrFile = "idrValues_" %>% paste0(tolower(organism)) %>% paste(params, sep = '_') %>% paste0('.txt')
idr = file.path(atacDir, idrFile)
atac = fread(idr)
names(atac)[1:3] = c('chr', 'start', 'end')

tn5.ranges = GRanges(seqnames = tn5$Chr, 
                     ranges = IRanges(start = tn5$Start, end = tn5$End))
atac.ranges = GRanges(seqnames = atac$chr, 
                      ranges = IRanges(start = atac$start, end = atac$end))

olaps = findOverlaps(tn5.ranges, atac.ranges)

olaps.df = olaps %>% as.data.frame()


tn5.in.atac = tn5[olaps.df$queryHits,]


tn5.in.atac.file.name = 'tn5_' %>% 
  paste0(organism) %>% 
  paste(sample, sep = "_") %>% 
  paste(strand, sep = '_') %>% 
  paste0('_peaks') %>% 
  paste0('_overlapping_') %>%  
  paste0(idrFile %>% strsplit('[.]') %>% lapply(function(x) x[1])) %>% 
  paste0('.bed')
write.table(tn5.in.atac, tn5.in.atac.file.name, row.names = F, col.names = F, sep = '\t', quote = F)


if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')

bed = paste0(tn5Dir, tn5.in.atac.file.name)
key = 'tn5_' %>% 
  paste0(organism) %>% 
  paste(sample, sep = '_') %>% 
  paste(strand, sep = '_') %>% 
  paste0('_peaks') %>% 
  paste0('_overlapping_atac_bulk')

bed
key
command = 'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(bed) %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')

JBdir
command



########################################################################################################################
#                                       FOR MERGED BULK AND CELL POPULATION ATAC                                       #
########################################################################################################################

#atacElavDir = '/data/personal_folders/jscepanovic/NV/ATAC_seq/elav_2017/elav_pos/'

# idrElavFile = 'idrValues.txt'
# idrElav = file.path(atacElavDir, idrElavFile)
# atacElav = fread(idrElav)

#atacMerge = rbind(atac, atacElav)
#names(atacMerge)[1:3] = c('chr', 'start', 'end')
# atacMerge.ranges = GRanges(seqnames = atacMerge$chr, 
#                            ranges = IRanges(start = atacMerge$start, end = atacMerge$end))
#olapsMerge = findOverlaps(tn5.ranges, atacMerge.ranges)
#olapsMerge.df = olapsMerge %>% as.data.frame()
#tn5.in.atacMerge = tn5[olapsMerge.df$queryHits,]
# tn5.in.atacMerge.file.name = 'tn5_' %>% 
#   paste0(organism) %>% 
#   paste(sample, sep = "_") %>% 
#   paste(strand, sep = '_') %>% 
#   paste0('_peaks') %>% 
#   paste0('_overlapping_') %>%  
#   paste0(idrFile %>% strsplit('[.]') %>% lapply(function(x) x[1])) %>% 
#   paste0('_merge_elav_idr') %>% 
#   paste0('.bed')
# write.table(tn5.in.atacMerge, tn5.in.atacMerge.file.name, row.names = F, col.names = T, sep = '\t', quote = F)

#bedMerge = paste0(tn5Dir, tn5.in.atacMerge.file.name)
# keyMerge = 'tn5_' %>% 
#   paste0(organism) %>% 
#   paste(sample, sep = '_') %>% 
#   paste(strand, sep = '_') %>% 
#   paste0('_peaks') %>% 
#   paste0('_overlapping_atac_bulk_elav_merge')
#bedMerge
#keyMerge
# if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
# if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
#JBdir
# commandMerge = 'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(bedMerge) %>% paste0(' --key ') %>% paste0(keyMerge) %>% paste0(' --trackLabel ') %>% paste0(keyMerge) %>% paste0(' --trackType HTMLFeatures')
# commandMerge

