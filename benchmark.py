library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyverse)
library(GEOquery)
library(data.table)
library(ggthemes)


ST_DBiT_10um = read.table('/ST-DBiT-10um-GSM10.txt',header = T,row.names = 1)
ST_10X_55um = read.csv('/CBS-Visium.csv',row.names = 1)
pbmc.data <- Read10X(data.dir = "/ST_scope_0.8um/")
ST_scope_0.8um <- CreateSeuratObject(counts = pbmc.data, project = "ST_scope_0.8um")
write.csv(ST_scope_0.8um@meta.data,'/ST_scope_0.8um/meta.csv')
Convert('/ST_Stereo_0.7um/Mouse_brain.h5ad', "h5seurat",
        overwrite = TRUE,assay = "RNA")
ST_Stereo_0.7um <- LoadH5Seurat("/ST_Stereo_0.7um/Mouse_brain.h5seurat")
ST_Stereo_0.7um
ST_Stereo_0.7um@meta.data$type = 'ST_Stereo_0.7um'
write.csv(ST_Stereo_0.7um@meta.data,'/ST_Stereo_0.7um/meta.csv')

mat <- read.table(file = "/GSE104323_10X_expression_data_V2.tab.gz" , header = T, row.names=1, sep="", as.is=T)
colnames(mat) = substr(colnames(mat),2,nchar(colnames(mat))-1)
meta <- read.table("/GSE104323_metadata_barcodes_24185cells.txt.gz", header=T, sep="\t", as.is=T)#,check.names = F
meta_f = meta[1:24185,]
rownames(meta_f) <- meta_f[,1]
rownames(meta_f) = substr(rownames(meta_f),1,nchar(rownames(meta_f))-1)
meta_f$celltype = meta_f$characteristics..cell.cluster

rat <- CreateSeuratObject(counts = mat, project = "rat",meta.data = meta_f)
rat[["percent.mt"]] <- PercentageFeatureSet(rat, pattern = "^MT-")
rat$orig.ident = rat$characteristics..age
rat
table(rat$celltype)
saveRDS(rat,'/benchmark_data/scRNA/rat.rds')

rat <- readRDS('/benchmark_data/scRNA/rat.rds')

scRNA = rat@meta.data

STus = Read10X(data.dir = "/matrix/")
STus <- CreateSeuratObject(counts = STus, project = "STus")
STus = STus@meta.data

ST_DBiT_10um$type = 'ST_DBiT_10um'
ST_10X_55um$type = 'ST_10X_55um'
ST_scope_0.8um$type = 'ST_scope_0.8um'
ST_Stereo_0.7um$type = 'ST_Stereo_0.7um'
scRNA$type = 'scRNA'
STus$type = 'ST_10um'
ST_Stereo_0.7um$nCount_RNA = ST_Stereo_0.7um$nFeature_count

total = rbind(ST_DBiT_10um[,c('type','nCount_RNA')],
     ST_10X_55um[,c('type','nCount_RNA')],
     ST_scope_0.8um[,c('type','nCount_RNA')],
     ST_Stereo_0.7um[,c('type','nCount_RNA')],
    scRNA[,c('type','nCount_RNA')],
     STus[,c('type','nCount_RNA')])

write.csv(total,'/benchmark_data/total.csv')
table(total$type)

total = read.csv('/benchmark_data/total.csv')
total = subset(total,total$type != 'scRNA')

p = ggplot(total, aes(x=type, y=log2(nCount_RNA),fill = type)) + 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(outlier.colour = NA)+ theme_base()+scale_y_continuous(limits= c(0, 15))+
    scale_fill_manual(values = c('#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#AEC7E8FF'))#'#1F77B4FF',
p
ggsave(p,filename = '/benchmark_log2_nooutlier.pdf',width = 16,height = 13)

p = ggplot(total, aes(x=type, y=log10(nCount_RNA),fill = type)) + 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(outlier.colour = NA)+ theme_base()+scale_y_continuous(limits= c(0, 5))+
    scale_fill_manual(values = c('#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#AEC7E8FF'))#'#1F77B4FF',
p
ggsave(p,filename = '/benchmark_log10_nooutlier.pdf',width = 16,height = 13)

p = ggplot(total, aes(x=type, y=nCount_RNA,fill = type)) + 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(outlier.colour = NA)+ theme_base()+scale_y_continuous(limits= c(0, 20000))+
    scale_fill_manual(values = c('#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#AEC7E8FF'))#'#1F77B4FF',
p
ggsave(p,filename = '/benchmark_nooutlier.pdf',width = 16,height = 13)