setwd('/home/jscepanovic')
library(data.table)
readStartCount = fread('tn5_Nv_20dpf_R1_Q2_forward_startCount', header=T, stringsAsFactors = F, data.table = F)
names(readStartCount) = c('Chr', 'Position', 'Coverage')
idx = which(readStartCount[,3] > 5)
idx_up = idx + 1
idx_down = idx - 1
peaks = (readStartCount[idx,3] >  readStartCount[idx_up,3]) & (readStartCount[idx,3] >  readStartCount[idx_up,3])
idx_peaks = idx[peaks]
#pos_reverse = readStartCount[idx_peaks,]
pos_forward = readStartCount[idx_peaks,]

write.table(pos_forward, '/home/jscepanovic/tn5_Nv_20dpf_forward_peaks.tsv', sep='\t', quote=F, row.names=F, col.names = T)

tab = fread('/home/jscepanovic/tn5_Nv_20dpf_reverse_peaks.tsv', sep='\t', fill=T, stringsAsFactors = F)
bed = tab[,1:2]
bed = cbind(bed, bed$Position)
names(bed) = c('Chr', 'Start', 'End')
write.table(bed, '/home/jscepanovic/tn5_Nv_20dpf_reverse_peaks.bed', sep='\t', quote=F, row.names=F, col.names = T)
