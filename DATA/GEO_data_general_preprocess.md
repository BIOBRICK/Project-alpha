# GEO Matrix General Preprocess

## Probe id & Gene id

Most of GEO data was collected through Gene Chip tech, thus the row names of the expression will be a little bit difficult to understand. So I have to transform all these probe id into familiar gene name.

Of course, if the gene id is symbol id, that will be good.

```R
setwd("Your_Work_Path")

#If you didn't 	install the GEOquery package, use the script below:
#BiocManager::install('GEOquery') 
library(GEOquery)
library(dplyr)
library(tidyr)

options(stringsAsFactors = F)

GPL6244 <-getGEO('GPL6244',destdir =".")
GPL6244_anno <- Table(GPL6244)

probe2symbol_df <- GPL6244_anno %>% 
  select(ID,gene_assignment) %>% 
  filter(gene_assignment != "---") %>% 
  separate(gene_assignment,c("drop","symbol"),sep="//") %>% 
  select(-drop)

exprSet<-read.table(file = "GSE42872_series_matrix.txt",header = T)
names(exprSet)[1] <- names(probe2symbol_df)[1]
exprSet$probe_id <- as.character(exprSet$probe_id)

exprSet <- exprSet %>% 
  
  inner_join(probe2symbol_df,by="probe_id") %>% #合并探针的信息
  
  select(-probe_id) %>% #去掉多余信息
  
  select(symbol, everything()) %>% #重新排列，
  
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  
  select(-rowMean) %>% #反向选择去除rowMean这一列
  
  tibble::column_to_rownames(colnames(.)[1]) # 把第一列变成行名并删除

write.csv(exprSet, file = 'GSE42872_expression_matrix.csv')

rm(list = ls())

```

**This script isn't my original creation, I get it from internet and thanks to [biotrainee](<http://www.biotrainee.com/>)  very much for this !**

## Ensembl id & Symbol id

