library(data.table)
library(dplyr)
bampe = fread('/home/jscepanovic/tn5_Nv_20dpf_Q2_srt_picardrmdup.sam', fill=T)
r1 = fread('/home/jscepanovic/tn5_Nv_20dpf_S5_R1_Q2_sort_reverse.sam', fill=T)
nrow(r1)
header = r1[1:17,]
r1 = r1[18:nrow(r1),]
nrow(r1)
sum(r1$V1 %in% bampe$V1)
r1.filt = r1 %>% dplyr::filter(r1$V1 %in% bampe$V1)
r1.filt.h = rbind(header,r1.filt)
write.table(r1.filt.h, '/home/jscepanovic/tn5_Nv_20dpf_R1_Q2_reverse_filtR.sam', sep='\t', quote=F, row.names = F, col.names = F)

r1 = fread('/home/jscepanovic/tn5_Nv_20dpf_S5_R1_Q2_sort_forward.sam', fill=T)
nrow(r1)
header = r1[1:17,]
r1 = r1[18:nrow(r1),]
nrow(r1)
sum(r1$V1 %in% bampe$V1)
r1.filt = r1 %>% dplyr::filter(r1$V1 %in% bampe$V1)
r1.filt.h = rbind(header,r1.filt)
write.table(r1.filt.h, '/home/jscepanovic/tn5_Nv_20dpf_R1_Q2_forward_filtR.sam', sep='\t', quote=F, row.names = F, col.names = F)