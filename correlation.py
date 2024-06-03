library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggthemes)

Sample1 <- Read10X(data.dir = "/Sample1/matrix/")
Sample1 <- CreateSeuratObject(counts = Sample1)
Sample1
Sample2 <- Read10X(data.dir = "/Sample2/matrix/")
Sample2 <- CreateSeuratObject(counts = Sample2)
Sample2

barcode = read.csv('/barcode.csv')

tmp1 = Sample1@assays$RNA@counts[,barcode$barcode]
tmp2 = Sample2@assays$RNA@counts[,barcode$barcode]

table(rownames(tmp1) == rownames(tmp2))

a=log10((rowMeans(tmp1))+1)
b=log10((rowMeans(tmp2))+1)
data=as.data.frame(cbind(a,b))
colnames(data)=c("UMI_Sample1","UMI_Sample2")
head(data)

cor(a,b,method="pearson")
p = ggplot(data,aes(x=UMI_Sample1,y=UMI_Sample2))+ geom_point(size=0.5,shape=15,color= '#31C2C7')+#geom_smooth(method=lm)+
    theme_classic()+stat_cor(data=data, method = "pearson")+
    xlab('Sample1 log10 total UMls per gene')+
    ylab('Sample2 log10 total UMls per gene')
p    

ggsave(p,file = 'correlation.pdf')