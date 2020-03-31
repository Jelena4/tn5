

setwd(tn5Dir)

bampe = file.path(tn5Dir, prefix %>% 
                    paste0('.PE.rmdup.sam'))
bampe = fread(bampe, fill = T)

bampe.header.length = sum((bampe[1:100,] %>% select(V1) %>% as.list)$V1 %>% str_sub(1,1) == '@')
bampe.header = bampe[1:bampe.header.length,]
bampe = bampe[(bampe.header.length + 1) : nrow(bampe),]

######### Forward strand ######### 

strand = 'forward'
r1 = file.path(tn5Dir, prefix %>% 
                 paste0('.R1.') %>% 
                 paste0(strand) %>% 
                 paste0('.sam'))
r1 = fread(r1, fill = T)

r1.header.length = sum((r1[1:100,] %>% select(V1) %>% as.list)$V1 %>% str_sub(1,1) == '@')
r1.header = r1[1:r1.header.length,]
r1 = r1[(r1.header.length + 1) : nrow(r1),]


sum(r1$V1 %in% bampe$V1)
r1.filt = r1 %>% dplyr::filter(r1$V1 %in% bampe$V1)
r1.filt = rbind(r1.header,r1.filt)
fileName = prefix %>% 
  paste0('.R1.') %>% 
  paste0(strand) %>% 
  paste0('.filtR.sam')
write.table(r1.filt, fileName, sep='\t', quote=F, row.names = F, col.names = F)

######### Reverse strand ######### 

strand = 'reverse'
r1 = file.path(tn5Dir, prefix %>% 
                 paste0('.R1.') %>% 
                 paste0(strand) %>% 
                 paste0('.sam'))
r1 = fread(r1, fill = T)

r1.header.length = sum((r1[1:100,] %>% select(V1) %>% as.list)$V1 %>% str_sub(1,1) == '@')
r1.header = r1[1:r1.header.length,]
r1 = r1[(r1.header.length + 1) : nrow(r1),]


sum(r1$V1 %in% bampe$V1)
r1.filt = r1 %>% dplyr::filter(r1$V1 %in% bampe$V1)
r1.filt = rbind(r1.header,r1.filt)
fileName = prefix %>% 
  paste0('.R1.') %>% 
  paste0(strand) %>% 
  paste0('.filtR.sam')
write.table(r1.filt, fileName, sep='\t', quote=F, row.names = F, col.names = F)


