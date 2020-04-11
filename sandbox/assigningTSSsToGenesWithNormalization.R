
############## Reading the gene models file #################

if(organism == 'Nv') models = '/data/personal_folders/jscepanovic/NV/nveGenes.vienna130208.gff'
if(organism == 'Sk') models = '/data/personal_folders/jscepanovic/SK/sacco_NCBI_valid_sorted_annot.higlass.annotation' # ?

models = fread(models)
models = models %>% mutate(name = models %>% select(V9) %>% as.list %>% unlist %>%  strsplit('=') %>% lapply(function(x) x[2]) %>% unlist%>% str_match_all('NVE[[:digit:]]+\\b') %>% unlist)
models = models %>% mutate(center = (V4 + V5) / 2)
utrs = models %>% filter(V3 == 'UTR')
cds = models %>% filter(V3 == 'CDS')

cds.starts = cds %>% group_by(name) %>% dplyr::slice(which.min(V4)) %>% select(name, V4, V7)
cds.ends = cds %>% group_by(name) %>% dplyr::slice(which.max(V5)) %>% select(name, V5, V7)

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
mRNAs = mRNAs %>% mutate(length = abs(V4 - V5))
mRNAs.ranges = GRanges(seqnames = mRNAs$V1, strand = mRNAs$V7, ranges = IRanges(start = mRNAs$V4, end = mRNAs$V5))



#############################################################
#                Getting coverage of each gene              #
#############################################################


################# Getting tn5 bedGraph ####################

strand = 'forward'
bw = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.filtR.bw')
bg = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.filtR.bedGraph')
cmd = 'bigWigToBedGraph ' %>% paste0(bw) %>% paste(bg, sep = ' ')
system(cmd)

strand = 'reverse'
bw = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.filtR.bw')
bg = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.filtR.bedGraph')
cmd = 'bigWigToBedGraph ' %>% paste0(bw) %>% paste(bg, sep = ' ')
system(cmd)


################# Reading tn5 bedGraph ####################

######################################################
#                   FORWARD STRAND                   #
######################################################

strand = 'forward'
bg = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.filtR.bedGraph')
bg = fread(bg)
if(strand == 'forward') bg$strand = '+'
if(strand == 'reverse') bg$strand = '-'
bg.ranges = GRanges(seqnames = bg$V1, strand = bg$strand, ranges = IRanges(start = bg$V2, end = bg$V3))

# Assigning coverage to genes

olaps = findOverlaps(mRNAs.ranges, bg.ranges, ignore.strand = F) %>% as.data.frame()
mRNAs[olaps$queryHits,] %>% select(name) %>% n_distinct # number of genes with non-zero tn5 prime coverage on the given strand
gene.cvrg.strand = mRNAs[olaps$queryHits,] %>% cbind(bg[olaps$subjectHits,'V4'])
names(gene.cvrg.strand)[13] = 'cvrg'
gene.cvrg = gene.cvrg.strand %>% group_by(name) %>% summarise(sum.cvrg = sum(cvrg))

######################################################
#                   REVERSE STRAND                   #
######################################################

# Reading the tn5 prime bedGraph file

strand = 'reverse'
bg = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.filtR.bedGraph')
bg = fread(bg)
if(strand == 'forward') bg$strand = '+'
if(strand == 'reverse') bg$strand = '-'
bg.ranges = GRanges(seqnames = bg$V1, strand = bg$strand, ranges = IRanges(start = bg$V2, end = bg$V3))

# Assigning coverage to genes

olaps = findOverlaps(mRNAs.ranges, bg.ranges, ignore.strand = F) %>% as.data.frame()
mRNAs[olaps$queryHits,] %>% select(name) %>% n_distinct # number of genes with non-zero tn5 prime coverage on the given strand
gene.cvrg.strand = mRNAs[olaps$queryHits,] %>% cbind(bg[olaps$subjectHits,'V4'])
names(gene.cvrg.strand)[13] = 'cvrg'
gene.cvrg = gene.cvrg %>% rbind(gene.cvrg.strand %>% group_by(name) %>% summarise(sum.cvrg = sum(cvrg)))

gene.cvrg %>% select(sum.cvrg) %>% summary
gene.cvrg %>% filter(sum.cvrg==0)

gene.cvrg = gene.cvrg %>% mutate(length =(mRNAs[match(gene.cvrg$name, mRNAs$name),] %>% select(length) %>% as.list)$length)
gene.cvrg = gene.cvrg %>% mutate(avg.cvrg = sum.cvrg/length)

############################################################################################################
#                                           Peak calling with macs2                                        #
############################################################################################################

strand = 'forward'
prefix.strand = prefix %>% paste(strand, sep = '.')
bg = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.filtR.bedGraph')
out = prefix %>% paste(strand, sep = '.') %>% paste0('.bed')

system('macs2 bdgpeakcall -c 10 -l 15 -i ' %>% paste0(bg) %>% paste0(' -o ') %>% paste0(out))

peaks = fread(out, header = F)
peaks = peaks[,c(1,2,3,5)]
out = prefix %>% paste(strand, sep = '.') %>% paste0('.forbrowser.bed')
write.table(peaks, out, sep = '\t', quote = F, row.names = F, col.names = F)

if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix.strand %>% paste0('.macs2peaks')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(out, sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')


strand = 'reverse'
prefix.strand = prefix %>% paste(strand, sep = '.')
bg = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.filtR.bedGraph')
out = prefix %>% paste(strand, sep = '.') %>% paste0('.bed')

system('macs2 bdgpeakcall -c 10 -l 15 -i ' %>% paste0(bg) %>% paste0(' -o ') %>% paste0(out))

peaks = fread(out, header = F)
peaks = peaks[,c(1,2,3,5)]
out = prefix %>% paste(strand, sep = '.') %>% paste0('.forbrowser.bed')
write.table(peaks, out, sep = '\t', quote = F, row.names = F, col.names = F)

if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix.strand %>% paste0('.macs2peaks')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(out, sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')


############################################################################################################
############################################################################################################

strand = 'forward'

# Reading the tn5 prime macs2 peaks file

prefix.strand = prefix %>% paste(strand, sep = '.')
macs = prefix.strand %>% paste0('.bed')

macs = fread(macs, header = F)
macs = macs %>% select(V1, V2, V3, V5) %>% mutate(coverage = V5/10) %>% data.table()

if(strand == 'forward') macs$strand = '+'
if(strand == 'reverse') macs$strand = '-'

macs.ranges = GRanges(seqnames = macs$V1, strand = macs$strand, ranges = IRanges(start = macs$V2, end = macs$V3))


# Reading the tn5 prime bbmap peaks file

bbmap = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.peaks.bed')
bbmap = fread(bbmap)
bbmap = bbmap %>% mutate(area = bbmap$V2 %>% as.character %>% str_sub(1, -3))

if(strand == 'forward') bbmap$strand = '+'
if(strand == 'reverse') bbmap$strand = '-'

bbmap.ranges = GRanges(seqnames = bbmap$V1, strand = bbmap$strand, ranges = IRanges(start = bbmap$V2, end = bbmap$V3))

olaps = findOverlaps(macs.ranges, bbmap.ranges, maxgap = 3) %>% as.data.frame()
allpeaks = macs[olaps$queryHits,] %>% cbind(bbmap[olaps$subjectHits,])
allpeaks = allpeaks[,c(1,2,3,5,6,8:11)]
names(allpeaks) = c('chr', 'start.macs', 'end.macs', 'max.coverage', 'strand', 'start.bbmap', 'end.bbmap', 'n.reads', 'area')

allpeaks$keep = allpeaks$n.reads > allpeaks$max.coverage / 6

goodpeaks = allpeaks %>% filter(keep == T)

goodpeaks = goodpeaks %>% mutate(size.macs = abs(end.macs - start.macs)) %>% data.table



write.table(goodpeaks, prefix.strand %>% paste0('.peaks.threshold.adjusted.bed'), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(goodpeaks[,c(1,6,7)], prefix.strand %>% paste0('.peaks.threshold.adjusted.forbrowser.bed'), row.names = F, col.names = F, quote = F, sep = '\t')


if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix.strand %>% paste0('.peaks.adjusted')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(prefix.strand %>% paste0('.peaks.threshold.adjusted.forbrowser.bed'), sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')

############################################################################################################
############################################################################################################

strand = 'reverse'

# Reading the tn5 prime macs2 peaks file

prefix.strand = prefix %>% paste(strand, sep = '.')
macs = prefix.strand %>% paste0('.bed')

macs = fread(macs, header = F)
macs = macs %>% select(V1, V2, V3, V5) %>% mutate(coverage = V5/10) %>% data.table()

if(strand == 'forward') macs$strand = '+'
if(strand == 'reverse') macs$strand = '-'

macs.ranges = GRanges(seqnames = macs$V1, strand = macs$strand, ranges = IRanges(start = macs$V2, end = macs$V3))


# Reading the tn5 prime bbmap peaks file

bbmap = prefix %>% paste0('.R1.') %>% paste0(strand) %>% paste0('.peaks.bed')
bbmap = fread(bbmap)
bbmap = bbmap %>% mutate(area = bbmap$V2 %>% as.character %>% str_sub(1, -3))

if(strand == 'forward') bbmap$strand = '+'
if(strand == 'reverse') bbmap$strand = '-'

bbmap.ranges = GRanges(seqnames = bbmap$V1, strand = bbmap$strand, ranges = IRanges(start = bbmap$V2, end = bbmap$V3))

olaps = findOverlaps(macs.ranges, bbmap.ranges, maxgap = 3) %>% as.data.frame()
allpeaks = macs[olaps$queryHits,] %>% cbind(bbmap[olaps$subjectHits,])
allpeaks = allpeaks[,c(1,2,3,5,6,8:11)]
names(allpeaks) = c('chr', 'start.macs', 'end.macs', 'max.coverage', 'strand', 'start.bbmap', 'end.bbmap', 'n.reads', 'area')

allpeaks$keep = allpeaks$n.reads > allpeaks$max.coverage / 6

goodpeaks = allpeaks %>% filter(keep == T)

goodpeaks = goodpeaks %>% mutate(size.macs = abs(end.macs - start.macs)) %>% data.table



write.table(goodpeaks, prefix.strand %>% paste0('.peaks.threshold.adjusted.bed'), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(goodpeaks[,c(1,6,7)], prefix.strand %>% paste0('.peaks.threshold.adjusted.forbrowser.bed'), row.names = F, col.names = F, quote = F, sep = '\t')


if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix.strand %>% paste0('.peaks.adjusted')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(prefix.strand %>% paste0('.peaks.threshold.adjusted.forbrowser.bed'), sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')

############################################################################################################
############################################################################################################

#                                              FORWARD STRAND                                              #

############################################################################################################
############################################################################################################

# Reading the peaks file to be assigned

strand = 'forward'
prefix.strand = prefix %>% paste(strand, sep = '.')
peaks = fread(prefix.strand %>% paste0('.peaks.threshold.adjusted.bed'))
names(peaks) = c('chr', 'start.macs', 'end.macs', 'max.coverage', 'strand', 'start.bbmap', 'end.bbmap', 'n.reads', 'area', 'keep', 'size.macs')
peaks.ranges = GRanges(seqnames = peaks$chr, strand = peaks$strand, ranges = IRanges(start = peaks$start.bbmap, end = peaks$end.bbmap))

# Finding UTR TSSs - directly in previously mapped UTRs
olaps = findOverlaps(peaks.ranges, utrs.5prime.ranges, ignore.strand = F)
peaks.olap = peaks[as.data.frame(olaps)$queryHits,]
peaks.olap = peaks.olap %>% mutate(name = utrs.5prime[as.data.frame(olaps)$subjectHits, 'name'])

# Getting maximum coverage peak for each UTR group (since there is often multiple adjacent peaks per a single UTR)
top.peaks = peaks.olap %>%
  group_by(chr, name) %>%
  dplyr::slice(which.max(n.reads)) %>% 
  ungroup

if(organism == 'Nv') assigned.peaks = top.peaks %>% mutate(type = 'directly at NVE UTRs')
if(organism == 'Sk') assigned.peaks = top.peaks %>% mutate(type = 'directly at NCBI UTRs')

assigned.peaks = assigned.peaks %>% ungroup

# Removing UTR TSSs from bed
peaks = anti_join(peaks, peaks.olap)
peaks.ranges = GRanges(seqnames = peaks$chr, strand = peaks$strand, ranges = IRanges(start = peaks$start.bbmap, end = peaks$end.bbmap))

############################################################################################################

# Finding all peaks within gene body


olaps = findOverlaps(peaks.ranges, mRNAs.ranges, ignore.strand = F)
peaks.olap = peaks[as.data.frame(olaps)$queryHits,]
peaks.olap = peaks.olap %>% mutate(name = mRNAs[as.data.frame(olaps)$subjectHits, 'name'])

# Getting maximum coverage peak for each UTR group (since there is often multiple adjacent peaks per a single UTR)
top.peaks = peaks.olap %>%
  group_by(chr, name, area) %>%
  dplyr::slice(which.max(n.reads)) %>% 
  ungroup

assigned.peaks = assigned.peaks %>% full_join(top.peaks %>% mutate(type = 'within gene body')) %>% ungroup

# Removing UTR TSSs from bed
peaks = anti_join(peaks, peaks.olap)
peaks.ranges = GRanges(seqnames = peaks$chr, strand = peaks$strand, ranges = IRanges(start = peaks$start.bbmap, end = peaks$end.bbmap))



############################################################################################################

# Finding peaks within 2000bp of the 5'-UTR

findTSSsWithinRange = function(range, peaks.ranges, utrs.5prime.ranges, peaks, assigned.peaks){
  
  olaps <<- findOverlaps(peaks.ranges, utrs.5prime.ranges, ignore.strand = F, select = 'first', maxgap = range)
  peaks.olap <<- peaks[!(olaps %>% is.na),] %>% mutate(name = (utrs.5prime[olaps,] %>% na.omit)$name)
  
  # Getting maximum coverage peak for each closely spaced peak group
  top.peaks <<- peaks.olap %>% group_by(chr, area, name) %>% dplyr::slice(which.max(n.reads)) %>% ungroup
  
  assigned.peaks <<- assigned.peaks %>% full_join(top.peaks %>% mutate(type = 'within ' %>% paste0(range) %>% paste0('bp of NVE UTRs')))
  
  # Removing UTR TSSs from bed
  peaks <<- anti_join(peaks, peaks.olap)
  peaks.ranges <<- GRanges(seqnames = peaks$chr, strand = peaks$strand, ranges = IRanges(start = peaks$start.bbmap, end = peaks$end.bbmap))
  
}

range = 2000
findTSSsWithinRange(range = range, peaks.ranges = peaks.ranges, utrs.5prime.ranges = utrs.5prime.ranges, peaks = peaks, assigned.peaks = assigned.peaks)


if(strand == 'forward') assigned.peaks.forward = assigned.peaks
if(strand == 'reverse') assigned.peaks.reverse = assigned.peaks

############################################################################################################
############################################################################################################

#                                              REVERSE STRAND                                              #

############################################################################################################
############################################################################################################

# Reading the peaks file to be assigned

strand = 'reverse'
prefix.strand = prefix %>% paste(strand, sep = '.')
peaks = fread(prefix.strand %>% paste0('.peaks.threshold.adjusted.bed'))
names(peaks) = c('chr', 'start.macs', 'end.macs', 'max.coverage', 'strand', 'start.bbmap', 'end.bbmap', 'n.reads', 'area', 'keep', 'size.macs')
peaks.ranges = GRanges(seqnames = peaks$chr, strand = peaks$strand, ranges = IRanges(start = peaks$start.bbmap, end = peaks$end.bbmap))

# Finding UTR TSSs - directly in previously mapped UTRs
olaps = findOverlaps(peaks.ranges, utrs.5prime.ranges, ignore.strand = F)
peaks.olap = peaks[as.data.frame(olaps)$queryHits,]
peaks.olap = peaks.olap %>% mutate(name = utrs.5prime[as.data.frame(olaps)$subjectHits, 'name'])

# Getting maximum coverage peak for each UTR group (since there is often multiple adjacent peaks per a single UTR)
top.peaks = peaks.olap %>%
  group_by(chr, name) %>%
  dplyr::slice(which.max(n.reads)) %>% 
  ungroup

if(organism == 'Nv') assigned.peaks = top.peaks %>% mutate(type = 'directly at NVE UTRs')
if(organism == 'Sk') assigned.peaks = top.peaks %>% mutate(type = 'directly at NCBI UTRs')

assigned.peaks = assigned.peaks %>% ungroup

# Removing UTR TSSs from bed
peaks = anti_join(peaks, peaks.olap)
peaks.ranges = GRanges(seqnames = peaks$chr, strand = peaks$strand, ranges = IRanges(start = peaks$start.bbmap, end = peaks$end.bbmap))

############################################################################################################

# Finding all peaks within gene body


olaps = findOverlaps(peaks.ranges, mRNAs.ranges, ignore.strand = F)
peaks.olap = peaks[as.data.frame(olaps)$queryHits,]
peaks.olap = peaks.olap %>% mutate(name = mRNAs[as.data.frame(olaps)$subjectHits, 'name'])

# Getting maximum coverage peak for each UTR group (since there is often multiple adjacent peaks per a single UTR)
top.peaks = peaks.olap %>%
  group_by(chr, name, area) %>%
  dplyr::slice(which.max(n.reads)) %>% 
  ungroup

assigned.peaks = assigned.peaks %>% full_join(top.peaks %>% mutate(type = 'within gene body')) %>% ungroup

# Removing UTR TSSs from bed
peaks = anti_join(peaks, peaks.olap)
peaks.ranges = GRanges(seqnames = peaks$chr, strand = peaks$strand, ranges = IRanges(start = peaks$start.bbmap, end = peaks$end.bbmap))



############################################################################################################

# Finding peaks within 2000bp of the 5'-UTR

range = 2000
findTSSsWithinRange(range = range, peaks.ranges = peaks.ranges, utrs.5prime.ranges = utrs.5prime.ranges, peaks = peaks, assigned.peaks = assigned.peaks)


if(strand == 'forward') assigned.peaks.forward = assigned.peaks
if(strand == 'reverse') assigned.peaks.reverse = assigned.peaks

assigned.peaks = rbind(assigned.peaks.forward %>% ungroup, assigned.peaks.reverse %>% ungroup)

file = prefix %>% paste0('.assigned.peaks.for.browser.bed')
write.table(assigned.peaks, prefix.strand %>% paste0('.assigned.peaks.bed'), row.names = F, col.names = T, quote = F, sep = '\t')
write.table(assigned.peaks[,c(1,6,7,12)], file, row.names = F, col.names = F, quote = F, sep = '\t')


if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix %>% paste0('.assigned.peaks')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(file, sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')


############################################################################################################
############################################################################################################

# Normalizing using gene coverage

assigned.peaks = assigned.peaks %>% left_join(gene.cvrg, by = 'name')
filtered.peaks = assigned.peaks %>% filter(sum.cvrg / n.reads < 300)
filtered.peaks = filtered.peaks %>% full_join(assigned.peaks %>% filter(n.reads > 200))

write.table(filtered.peaks, prefix %>% paste0('.peaks.threshold.adjusted.normalized.bed'), row.names = F, col.names = T, quote = F, sep = '\t')
write.table(filtered.peaks[,c(1,6,7,12)], prefix %>% paste0('.peaks.threshold.adjusted.normalized.forbrowser.bed'), row.names = F, col.names = F, quote = F, sep = '\t')

if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix %>% paste0('.peaks.threshold.normalized')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(prefix %>% paste0('.peaks.threshold.adjusted.normalized.forbrowser.bed'), sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')


############################################################################################################
############################################################################################################


peaks = prefix %>% paste0('.peaks.threshold.adjusted.normalized.bed')
peaks = fread(peaks)

########## Removing very closely spaced peaks (E.g. NVE4175)

peaks = peaks %>% mutate(area = str_sub(start.bbmap, 1, -3)) %>% mutate(subarea = str_sub(start.bbmap, -2, -2))

peaks = peaks %>% mutate(subarea2 = subarea %>% as.integer()%/% 2)

peaks = peaks %>% group_by(chr, name, subarea2) %>% dplyr::slice(which.max(n.reads)) %>% ungroup


write.table(peaks, prefix %>% paste0('.peaks.threshold.adjusted.normalized.filt.bed'), row.names = F, col.names = T, quote = F, sep = '\t')
write.table(peaks[,c(1,6,7,12)], prefix %>% paste0('.peaks.threshold.adjusted.normalized.filt.forbrowser.bed'), row.names = F, col.names = F, quote = F, sep = '\t')

if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix %>% paste0('.peaks.threshold.normalized.filt')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(prefix %>% paste0('.peaks.threshold.adjusted.normalized.filt.forbrowser.bed'), sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')



########## Removing 3'-UTR peaks

peaks.ranges = GRanges(seqnames = peaks$chr, strand = peaks$strand, ranges = IRanges(start = peaks$start.bbmap, end = peaks$end.bbmap))
utrs.3prime.ranges
olaps = findOverlaps(peaks.ranges, utrs.3prime.ranges, ignore.strand = F, maxgap = 30) %>% as.data.frame()
filtout = peaks[olaps$queryHits,] %>% cbind(utrs.3prime[olaps$subjectHits,'name'])
names(filtout)[19] = 'name.3prime.UTR'
filtout = filtout %>% filter(name == name.3prime.UTR) %>% data.table
filtout = filtout %>% select(-name.3prime.UTR)
peaks = dplyr::setdiff(peaks, filtout)

write.table(peaks, prefix %>% paste0('.peaks.threshold.adjusted.normalized.filt.bed'), row.names = F, col.names = T, quote = F, sep = '\t')
write.table(peaks[,c(1,6,7,12)], prefix %>% paste0('.peaks.threshold.adjusted.normalized.filt.forbrowser.bed'), row.names = F, col.names = F, quote = F, sep = '\t')

if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix %>% paste0('.peaks.threshold.normalized.filt')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(prefix %>% paste0('.peaks.threshold.adjusted.normalized.filt.forbrowser.bed'), sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')


########## Removing downstream peaks (rlative to gene body)
peaks = prefix %>% paste0('.peaks.threshold.adjusted.normalized.filt.bed')
peaks = fread(peaks)

peaks = peaks %>% left_join(cds.coords, by = 'name') %>% data.table %>% select(-length.y,-strand.y)

names(peaks)[5] = 'strand'

a = peaks %>% filter(strand == '+' & start.bbmap < end)
b = peaks %>% filter(strand == '-' & start.bbmap > start)

peaks %>% nrow

peaks = full_join(a, b) %>% data.table

write.table(peaks, prefix %>% paste0('.peaks.final.bed'), row.names = F, col.names = T, quote = F, sep = '\t')
write.table(peaks[,c(1,6,7,12)], prefix %>% paste0('.peaks.final.forborowser.bed'), row.names = F, col.names = F, quote = F, sep = '\t')

if(organism == 'Nv') JBdir = ('/home/lab/website/JBrowse/JBrowse_NV_experimental_d/')
if(organism == 'Sk') JBdir = ('/home/lab/website/JBrowse/JBrowse_SK_d/')
key = prefix %>% paste0('.peaks.final')
'cd ' %>% paste0(JBdir)
'sudo bin/flatfile-to-json.pl --bed ' %>% paste0(tn5Dir) %>% paste(prefix %>% paste0('.peaks.final.forborowser.bed'), sep = '/') %>% paste0(' --key ') %>% paste0(key) %>% paste0(' --trackLabel ') %>% paste0(key) %>% paste0(' --trackType HTMLFeatures')

peaks %>% select(name) %>% n_distinct


