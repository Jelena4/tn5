
############## Reading the gene models file #################

if(organism == 'Nv') models = '/data/personal_folders/jscepanovic/NV/nveGenes.vienna130208.gff'
if(organism == 'Sk') models = '/data/personal_folders/jscepanovic/SK/sacco_NCBI_valid_sorted_annot.higlass.annotation' # ?

models = fread(models)
models = models %>% mutate(name = models %>% select(V9) %>% as.list %>% unlist %>%  strsplit('=') %>% lapply(function(x) x[2]) %>% unlist%>% str_match_all('NVE[[:digit:]]+\\b') %>% unlist)
models = models %>% mutate(center = (V4 + V5) / 2)
utrs = models %>% filter(V3 == 'UTR')
cds = models %>% filter(V3 == 'CDS')

cds.starts = cds %>% filter(V7 == '+') %>% group_by(name) %>% dplyr::slice(which.min(V4)) %>% select(name, V4, V7)
cds.ends = cds %>% filter(V7 == '+') %>% group_by(name) %>% dplyr::slice(which.max(V5)) %>% select(name, V4, V7)
cds.starts = cds.starts %>% full_join(cds %>% filter(V7 == '-') %>% group_by(name) %>% dplyr::slice(which.max(V4)) %>% select(name, V4, V7))
cds.ends = cds.ends %>% full_join(cds %>% filter(V7 == '-') %>% group_by(name) %>% dplyr::slice(which.min(V5)) %>% select(name, V4, V7))

cds.coords = left_join(cds.starts, cds.ends, by = 'name') %>% select(c(1,2,4,3))
names(cds.coords) = c('name', 'start', 'end', 'strand')
cds.coords = cds.coords %>% mutate(center = (start + end)/2) %>% mutate(length = abs(end - start)) %>% as.data.frame()


#######################   UTRs   ########################

names = (cds.coords %>% select(name) %>% unique() %>% as.list())$name
utrs = utrs %>% filter(utrs$name %in% names)
type = vector(length = nrow(utrs))

for(row in 1:nrow(utrs))
{
      if(utrs[row, 'V7'] == '+')
      {
        if(utrs[row, 'center'] < cds.coords[match(utrs[row,'name'], cds.coords[,'name']),'center']) {type[row] = '5'}
        else {type[row] = '3'}
      }
      else
      {
        if(utrs[row, 'center'] > cds.coords[match(utrs[row,'name'], cds.coords[,'name']),'center']) {type[row] = '5'}
        else{type[row] = '3'}
      }
  if(row %% 100 == 0) print(row)
}

utrs = cbind(utrs, type)

utrs.5prime = utrs %>% filter(type == '5')
utrs.3prime = utrs %>% filter(type == '3')
utrs.5prime.ranges = GRanges(seqnames = utrs.5prime$V1, strand = utrs.5prime$V7, ranges = IRanges(start = utrs.5prime$V4, end = utrs.5prime$V5))
utrs.3prime.ranges = GRanges(seqnames = utrs.3prime$V1, strand = utrs.3prime$V7, ranges = IRanges(start = utrs.3prime$V4, end = utrs.3prime$V5))


mRNAs = models %>% filter(V3 == 'mRNA')
mRNAs.ranges = GRanges(seqnames = mRNAs$V1, strand = mRNAs$V7, ranges = IRanges(start = mRNAs$V4, end = mRNAs$V5))



######################################################
#                   FORWARD STRAND                   #
######################################################

# Reading the tn5 prime peaks file

strand = 'forward'
prefix.strand = prefix %>% paste0('.R1.') %>% paste0(strand)

bed = prefix.strand %>% paste0('.peaks.bed')
bed = fread(bed)
bed = bed %>% mutate(area = bed$V2 %>% as.character %>% str_sub(1, -3))

if(strand == 'forward') bed$strand = '+'
if(strand == 'reverse') bed$strand = '-'

bed.ranges = GRanges(seqnames = bed$V1, strand = bed$strand, ranges = IRanges(start = bed$V2, end = bed$V3))

############################################################################################################
############################################################################################################

# Finding UTR TSSs - directly in previously mapped UTRs
olaps = findOverlaps(bed.ranges, utrs.5prime.ranges, ignore.strand = F)
peaks = bed[as.data.frame(olaps)$queryHits,]
peaks = peaks %>% mutate(name = utrs.5prime[as.data.frame(olaps)$subjectHits, 'name'])

# Getting maximum coverage peak for each UTR group (since there is often multiple adjacent peaks per a single UTR)
top.peaks = peaks %>%
  group_by(V1, name) %>%
  dplyr::slice(which.max(V4)) %>% 
  ungroup

if(organism == 'Nv') assigned.peaks = top.peaks %>% mutate(type = 'directly at NVE UTRs')
if(organism == 'Sk') assigned.peaks = top.peaks %>% mutate(type = 'directly at NCBI UTRs')

assigned.peaks = assigned.peaks %>% ungroup

# Removing UTR TSSs from bed
bed = anti_join(bed, peaks)
bed.ranges = GRanges(seqnames = bed$V1, strand = bed$strand, ranges = IRanges(start = bed$V2, end = bed$V3))

############################################################################################################
############################################################################################################
############################################################################################################

findTSSsWithinRange = function(range, bed.ranges, utrs.5prime.ranges, bed, assigned.peaks){
  
  olaps <<- findOverlaps(bed.ranges, utrs.5prime.ranges, ignore.strand = F, select = 'first', maxgap = range)
  peaks <<- bed[!(olaps %>% is.na),] %>% mutate(name = (utrs.5prime[olaps,] %>% na.omit)$name)
  
  # Getting maximum coverage peak for each closely spaced peak group
  top.peaks <<- peaks %>% group_by(V1, area, name) %>% dplyr::slice(which.max(V4)) %>% ungroup
  
  assigned.peaks <<- assigned.peaks %>% full_join(top.peaks %>% mutate(type = 'within ' %>% paste0(range) %>% paste0('bp of NVE UTRs')))
  
  # Removing UTR TSSs from bed
  bed <<- anti_join(bed, peaks)
  bed.ranges <<- GRanges(seqnames = bed$V1, strand = bed$strand, ranges = IRanges(start = bed$V2, end = bed$V3))
  
}

############################################################################################################
############################################################################################################
############################################################################################################

range = 2000
findTSSsWithinRange(range = range, bed.ranges = bed.ranges, utrs.5prime.ranges = utrs.5prime.ranges, bed = bed, assigned.peaks = assigned.peaks)

assigned.peaks
bed %>% nrow

assigned.peaks %>% nrow
assigned.peaks %>% select(V1,V2,V3) %>% n_distinct

if(strand == 'forward') assigned.peaks.forward = assigned.peaks
if(strand == 'reverse') assigned.peaks.reverse = assigned.peaks

assigned.peaks = rbind(assigned.peaks.forward %>% ungroup, assigned.peaks.reverse %>% ungroup)

name =  prefix %>% paste0('.assigned.peaks.bed')
name.browser =  prefix %>% paste0('.assigned.peaks.for.browser.bed')
write.table(assigned.peaks, name, quote = F, row.names = F, col.names = F, sep = '\t')
write.table(assigned.peaks %>% select(V1,V2,V3, name), name.browser, quote = F, row.names = F, col.names = F, sep = '\t')

if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix %>% paste0('.assigned.peaks')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(name.browser, sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')






assigned.peaks


# assigned.ranges = GRanges(seqnames = assigned.peaks$V1, strand = assigned.peaks$strand, ranges = IRanges(start = assigned.peaks$V2, end = assigned.peaks$V3))
# utrs.3prime.ranges = GRanges(seqnames = utrs.3prime$V1, strand = utrs.3prime$V7, ranges = IRanges(start = utrs.3prime$V4, end = utrs.3prime$V5))
# 
# olaps = findOverlaps(assigned.ranges, utrs.3prime.ranges, ignore.strand = F) # assigned peaks that are in 3'-UTRs
# 
# olaps %>% as.data.frame() %>% nrow
# 
# assigned.peaks[as.data.frame(olaps)$queryHits,] %>% head
# utrs.3prime[as.data.frame(olaps)$subjectHits,] %>% head


# remove those? many are real (ovelapping or bad gene models?)


# 
# olaps = findOverlaps(bed.ranges, mRNAs.ranges, ignore.strand = F) # directly in previously mapped UTRs
# peaks = bed[as.data.frame(olaps)$queryHits,]
# peaks = mRNA.peaks %>% mutate(gene.name = mRNAs[as.data.frame(olaps)$subjectHits,'name'])
# 
