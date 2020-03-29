
rm(list = ls())

library(data.table)
library(dplyr)

organism = 'Nv'
sample = '20dpf'
strand = 'forward'

tn5Dir = '/data/personal_folders/jscepanovic/tn5'
setwd(tn5Dir)

bampe = file.path(tn5Dir, 'tn5_' %>% 
                    paste0(organism) %>% 
                    paste(sample, sep = "_") %>% 
                    paste0('.PE.rmdup.sam'))

bampe = fread(bampe, fill = T)

r1 = file.path(tn5Dir, 'tn5_' %>% 
                 paste0(organism) %>% 
                 paste(sample, sep = "_") %>% 
                 paste0('.R1.') %>% 
                 paste0(strand) %>% 
                 paste0('.sam'))

r1 = fread(r1, fill = T)

nrow(r1)
header = r1[1:17,]
r1 = r1[18:nrow(r1),]
nrow(r1)

sum(r1$V1 %in% bampe$V1)
r1.filt = r1 %>% dplyr::filter(r1$V1 %in% bampe$V1)
r1.filt.h = rbind(header,r1.filt)
fileName = 'tn5_' %>% 
  paste0(organism) %>% 
  paste(sample, sep = "_") %>% 
  paste0('.R1.') %>% 
  paste0(strand) %>% 
  paste0('.filtR.sam')
write.table(r1.filt.h, fileName, sep='\t', quote=F, row.names = F, col.names = F)


