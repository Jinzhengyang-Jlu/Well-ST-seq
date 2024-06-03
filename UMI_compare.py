library(Seurat)
library(patchwork)
library(ggplot2)
library(plyr)
library(dplyr)
library(harmony)
library(car)
library(ggpubr) 
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(ggsci)
library(Rgraphviz)
library(ggforce)
library(ggstatsplot)
library(ggrepel)
library(corrplot)
library(pheatmap)
library(readxl)
set.seed(2022)

10um <- Read10X(data.dir = "/S10/matrix")
10um <- CreateSeuratObject(counts = 10um, project = "10um", min.cells = 0, min.features = 0)
10um
20um <- Read10X(data.dir = "/S20/matrix")
20um <- CreateSeuratObject(counts = 20um, project = "20um", min.cells = 0, min.features = 0)
20um

dataset3_15um <- read_excel('/Dataset3_seqFISH_mmc6.xlsx', sheet=1,col_names = TRUE)
rownames(dataset3_15um) <- dataset3_15um$`'Tal1'`
dataset3_15um <- dataset3_15um[,-1]
dataset3_15um <- CreateSeuratObject(counts = dataset3_15um, project = "dataset3_15um", min.cells = 0, min.features = 0)
dataset3_15um
dataset43_5um <- read.delim('/Dataset43_Hippocampus/merged_counts.tsv')
rownames(dataset43_5um) <- dataset43_5um$X
dataset43_5um <- dataset43_5um[,-1]
dataset43_5um <- t(dataset43_5um)
dataset43_5um <- CreateSeuratObject(counts = dataset43_5um, project = "dataset43_5um", min.cells = 0, min.features = 0)
dataset43_5um
dataset41_10um <- read.table(gzfile("/Dataset41_Slide-seqV2_SCP815/other//Puck_200115_08.digital_expression.txt.gz"),sep="\t")
colnames(dataset41_10um) <- dataset41_10um[1,]
rownames(dataset41_10um) <- dataset41_10um[,1]
dataset41_10um <- dataset41_10um[-1,-1]
dataset41_10um <- CreateSeuratObject(counts = dataset41_10um, project = "dataset41_10um", min.cells = 0, min.features = 0)
dataset41_10um
dataset33_12um <- Read10X(data.dir = "/Dataset33/st")
dataset33_12um <- CreateSeuratObject(counts = dataset33_12um, project = "dataset33_10um", min.cells = 0, min.features = 0)
dataset33_12um

dataset43_5um@meta.data$orig.ident <- 'dataset43_5um'
dataset41_10um@meta.data$orig.ident <- 'Slide-seqV2_10um'
dataset33_12um@meta.data$orig.ident <- 'Visium_10um'
dataset3_15um@meta.data$orig.ident <- 'seqFish_15um'

a <- rbind(10um@meta.data, 20um@meta.data,dataset43_5um@meta.data,dataset41_10um@meta.data,dataset33_12um@meta.data,dataset3_15um@meta.data)
a
table(a$orig.ident)

mean(subset(a, orig.ident=='10um')$nCount_RNA)
mean(subset(a, orig.ident=='20um')$nCount_RNA)

p <- ggplot(a, aes(orig.ident, nCount_RNA))+
  geom_boxplot(aes(fill = orig.ident), notch = FALSE, outlier.alpha  = 1) +
  scale_fill_brewer(palette = "Set2") +
  #geom_jitter(aes(fill=orig.ident),width =0.2,shape = 21,size=0, alpha =0.5)++ 
  theme_classic() +
  ggtitle('') +
  labs(x = "", y = 'Number of UMI') +
  theme(axis.text.x=element_text(size=15,face = "bold",angle=90,hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=15,face = "bold"),
        axis.title.y=element_text(size=20,face = "bold")
       )
p

library(scales)
p <- ggplot(a, aes(orig.ident, nCount_RNA))+
  geom_boxplot(aes(fill = orig.ident), notch = FALSE, outlier.alpha  = 1) +
  scale_fill_brewer(palette = "Set2") +
  ylim(0,30000)+
  scale_y_continuous(breaks = c(1000, 3000, 5000,10000,20000,30000),
                    labels = c("1k", "3k", "5k",'10k','20k','30k')) +
  #geom_jitter(aes(fill=orig.ident),width =0.2,shape = 21,size=0, alpha =0.5)++ 
  theme_classic() +
  ggtitle('') +
  labs(x = "", y = 'Number of UMI') +
  #ylim(0,30000)+
  theme(axis.text.x=element_text(size=15,face = "bold",angle=90,hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=7,face = "bold"),
        axis.title.y=element_text(size=23,face = "bold")
       )+ NoLegend()
p

p <- ggplot(a, aes(orig.ident, nFeature_RNA))+
  geom_boxplot(aes(fill = orig.ident), notch = FALSE, outlier.alpha  = 1) +
  scale_fill_brewer(palette = "Set2") +
  #ylim(0,30000)+
  scale_y_continuous(breaks = c(1000,2000, 3000, 5000,10000,20000,30000),
                    labels = c("1k",'2k', "3k", "5k",'10k','20k','30k')) +
  #geom_jitter(aes(fill=orig.ident),width =0.2,shape = 21,size=0, alpha =0.5)++ 
  theme_classic() +
  ggtitle('') +
  labs(x = "", y = 'Detected Genes per Spot') +
  #ylim(0,30000)+
  theme(axis.text.x=element_text(size=15,face = "bold",angle=90,hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=7,face = "bold"),
        axis.title.y=element_text(size=23,face = "bold")
       )+ NoLegend()
p

library(tidyverse)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)

bam_file <- "/S10.featureCounts.query.bam"
bam <- BamFile(bam_file)
reads <- readGAlignments(bam, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE)))
reads

umi <- reads@seqnames
gene_counts <- as.data.frame(table(umi))
gene_counts

gene_sum <- cumsum(gene_counts$Freq)
gene_sum

plot(gene_sum, type="l", col = 'red',
     xlab="Gene", ylab="# of UMIs per pixel")

data <- read.csv('/10um/counts.csv',row.names = 1)
gene_list <- colnames(data)
gene_list

bam <- readGAlignments('/S10.featureCounts.query.bam')

gene_count <- c()
for (gene in gene_list){
  #获取基因片段
  gene_granges <- GRanges(seqnames=c(gene), IRanges(start=1, end=1000))
  #读取基因片段
  gene_alignments <- getAlignments(bam, gene_granges)
  #计算基因累计数量
  gene_count[gene] <- length(unique(gene_alignments$qname))
}

umi_count <- c()
for (umi in 1:max(gene_count)){
  umi_count[umi] <- length(which(gene_count >= umi))
}

library(ggplot2)
library(reshape2)

bamfile <- read.table("/S10.depth", header=FALSE)
colnames(bamfile) <- c("Chromosome","Position","Depth")
bamfile$Chromosome <- factor(bamfile$Chromosome, levels=unique(bamfile$Chromosome))

ggplot(bamfile, aes(x=Position, y=Depth)) +
  geom_line(color="blue", size=0.3) +
  geom_point(color="red", size=0.3) +
  #facet_wrap(~ Chromosome, scales="free") +
  theme_bw()

args <- commandArgs(trailingOnly=TRUE)
depth <- read.table('/S10.depth',header=FALSE,stringsAsFactors=FALSE)
# 计算reads数目
read_num <- sum(depth[,3])
# 计算深度
depth_num <- seq(1,max(depth[,3]),by=1)
# 计算基因数目
gene_num <- c()

sum(depth[,3])
for(i in depth_num){
  gene_num[i] <- sum(depth[,3]>=i)
}
ggplot() +
  geom_point(aes(x=depth_num,y=gene_num),size=0.5,colour='red') +
  geom_line(aes(x=depth_num,y=gene_num),size=0.5,colour='red') +
  ggtitle(paste('Detected gene number with reads saturation, total reads:',read_num)) +
  xlab('Depth') +
  ylab('Gene number')

library(org.Mm.eg.db)
library(ggplot2)
 k=keys(org.Mm.eg.db,keytype = "ENSEMBL")
list=select(org.Mm.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
dim(list)
head(list,5)

data <- read.delim('/S10_ENSEMBL.txt', header=F)
data <- as.data.frame(data)
dim(data)
head(data)

length(unique(data$V1))
s <- substring(data$V1, 6)
length(s)
data$ENSEMBL <- s
head(data)

ID <- data$ENSEMBL
ID_list=list[match(ID,list[,"ENSEMBL"]),]
ID_list

length(unique(ID_list$SYMBOL))

dd <- ID_list[complete.cases(ID_list),]
dim(dd)

length(unique(dd$SYMBOL))
data$SYMBOL <- ID_list$SYMBOL
dim(data)
head(data)

seq(0,88007597,by=1000000)

dd = data[sample(nrow(data), 0), ]
length(unique(dd$SYMBOL))
dat <- as.data.frame(seq(0,88007597,by=1000000))
colnames(dat) <- 'Reads'
dat$genes <- 0
for(i in 1:nrow(dat)){
    dd <- data[sample(nrow(data), dat$Reads[i]), ]
    dat$genes[i] <- length(unique(dd$SYMBOL))-1    
}
dat$genes[1] <- 0
head(dat)

dat$Reads_10um <- dat$Reads / 1000000

dat_10um <- dat
ggplot(dat_10um, aes(x = Reads_10um, y = genes)) + geom_point(col = '#9AC9DBFF') + geom_line(col = '#9AC9DBFF')

dd <- data[sample(nrow(data), 5e+06), ]
length(unique(dd$SYMBOL))-1  

data <- read.delim('/S20_ENSEMBL.txt', header=F)
data <- as.data.frame(data)
dim(data)
head(data)

length(unique(data$V1))
s <- substring(data$V1, 6)
length(s)
data$ENSEMBL <- s
head(data)

ID <- data$ENSEMBL
ID_list=list[match(ID,list[,"ENSEMBL"]),]
ID_list

length(unique(ID_list$SYMBOL))
data$SYMBOL <- ID_list$SYMBOL
dim(data)
head(data)

dat <- as.data.frame(seq(0,nrow(data),by=1000000))
colnames(dat) <- 'Reads'
dat$genes <- 0
for(i in 1:nrow(dat)){
    dd <- data[sample(nrow(data), dat$Reads[i]), ]
    dat$genes[i] <- length(unique(dd$SYMBOL))-1    
}
dat$genes[1] <- 0
head(dat)

dat_20um <- dat
dat_20um$Reads_20um <- dat_20um$Reads / 1000000
ggplot(dat_20um, aes(x = Reads_20um, y = genes)) + geom_point(col = '#DBBEAC') + geom_line(col = '#DBBEAC')

dat_10um
dat_20um
dat_10um$Read <- dat_10um$Reads_10um
dat_10um$label <- '10um'
dat_20um$Read <- dat_20um$Reads_20um
dat_20um$label <- '20um'
write.csv(dat_10um, '/S10_saturation_curve.csv')
write.csv(dat_20um, '/S20_saturation_curve.csv')
dat_10um <- read.csv('/S10_saturation_curve.csv', header = T, row.names = 1)
dat_20um <- read.csv('/S20_saturation_curve.csv', header = T, row.names = 1)

dat <- rbind(dat_10um[,-3],dat_20um[,-3])
dim(dat)
head(dat)

themes = theme_classic()+
    theme(axis.title.x=element_text(vjust=0, size=24,face = "plain"),
      axis.title.y=element_text(vjust=1, size=24,face = "plain"),
      axis.text.x=element_text(vjust=0, size=20,face = "plain"),
      axis.text.y=element_text(vjust=0, size=20,face = "plain"),
      plot.title = element_text(size=20,hjust = 0.5),
      legend.key.size = unit(1, 'cm'),
      legend.title = element_blank(), 
      legend.text = element_text(size=20),
      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
colors_manual <- c('#FF8884FF','#2878B5FF','#9AC9DBFF','#6B705C','#DBBEAC','#7F6DA9','#B5BB41','#4D98AB')
dat$label <- factor(dat$label, levels = c('10um','20um'))
ggplot(dat, aes(x = Read, y = genes),cols = label) + geom_point(aes(colour = label)) + geom_line(aes(colour = label))+
  xlab('Reads(M)')+
  ylab('# of gene detected')+
  scale_color_manual(values=colors_manual)+
  themes+
  theme(legend.position = c(0.8, 0.2))
ggsave('/Saturation curve.pdf', width = 11, height = 9)

dat1 <- read.delim('/accumulation_curves/S10_UMI_gene.txt', sep = '')
dat2 <- read.delim('/accumulation_curves/S20_UMI_gene.txt', sep = '')
dat1 <- subset(dat1, UMI != 0)
dat2 <- subset(dat2, UMI != 0)

cumi<-c()
total <- sum(dat1$UMI)
cumi[1]<-dat1$UMI[1]/total
for(j in 2:nrow(dat1)){
    cumi[j]=dat1$UMI[j]/total+cumi[j-1]
}
cgene <- c()
for(i in 1:nrow(dat1)){
    cgene[i]=i/nrow(dat1)
}
dat1 <- data.frame(cgene,cumi)
dat1$label <- '10um'
tail(dat1)

cumi<-c()
total <- sum(dat2$UMI)
cumi[1]<-dat2$UMI[1]/total
for(j in 2:nrow(dat2)){
    cumi[j]=dat2$UMI[j]/total+cumi[j-1]
}
cgene <- c()
for(i in 1:nrow(dat2)){
    cgene[i]=i/nrow(dat2)
}
dat2 <- data.frame(cgene,cumi)
dat2$label <- '20um'
dat <- rbind(dat1,dat2)
tail(dat)

dat$cgene <- dat$cgene *100
dat$cumi <- dat$cumi *100

dat$label <- factor(dat$label, levels = c('10um','20um'))
ggplot(dat, aes(x = cgene, y = cumi),cols = label) + geom_line(aes(colour = label))+
  xlab('Genes(%)')+
  ylab('UMIs(%)')+
  scale_color_manual(values=colors_manual)+
  themes+
  theme(legend.position = c(0.8, 0.2))#+
#geom_vline(aes(xintercept=20), colour="black", linetype="dashed")+
#geom_vline(aes(xintercept=40), colour="black", linetype="dashed")
ggsave('/accum_saturation.pdf', width = 11, height = 9)

dat1 <- read.delim('/accumulation_curves/S10_UMI_gene.txt', sep = '')
dat2 <- read.delim('/accumulation_curves/S20_UMI_gene.txt', sep = '')
dat1 <- subset(dat1, UMI != 0)
dat2 <- subset(dat2, UMI != 0)

ref <- read.delim('/ref_10x_mm10/genes/hgTables.txt', header = F)
ref <- ref[,1:4]
ref

list <- read.delim('/ref_10x_mm10/genes/gencode.vM10.metadata.MGI', header=F)
colnames(list) <- c('ENSEMBL','SYMBOL')
list

ID_list=list[match(ref$V4,list[,"ENSEMBL"]),]
#ID_list <- ID_list[complete.cases(ID_list),]
dim(ID_list)
head(ID_list)

ref$SYMBOL <- ID_list$SYMBOL
head(ref)

dat=ref[match(dat1$gene,ref[,"SYMBOL"]),]
dat$bp <- dat$V3 - dat$V2
dat <- dat[,5:6]
dat$UMI <- dat1$UMI
dat <- dat[complete.cases(dat),]
dat

min(dat$bp)
max(dat$bp)
s <- (max(dat$bp)-min(dat$bp))/10

dat$group<-0
for(i in 1:nrow(dat)){
    if(dat$bp[i] %in% seq(60,60+s)){
        dat$group[i] <- 'bin1'
        }
    else if(dat$bp[i] %in% seq(60+s,60+2*s)){
        dat$group[i] <- 'bin2'
        }
    else if(dat$bp[i] %in% seq(60+2*s,60+3*s)){
        dat$group[i] <- 'bin3'
        }
    else if(dat$bp[i] %in% seq(60+3*s,60+4*s)){
        dat$group[i] <- 'bin4'
        }
    else if(dat$bp[i] %in% seq(60+4*s,60+5*s)){
        dat$group[i] <- 'bin5'
        }
    else if(dat$bp[i] %in% seq(60+5*s,60+6*s)){
        dat$group[i] <- 'bin6'
        }
    else if(dat$bp[i] %in% seq(60+6*s,60+7*s)){
        dat$group[i] <- 'bin7'
        }
    else if(dat$bp[i] %in% seq(60+7*s,60+8*s)){
        dat$group[i] <- 'bin8'
        }
    else if(dat$bp[i] %in% seq(60+8*s,60+9*s)){
        dat$group[i] <- 'bin9'
        }
    else{
        dat$group[i] <- 'bin10'
        }
}
dat
table(dat$group)

library("tidyverse")
dat$logUMI <- log10(dat$UMI+1)
ggplot(dat, aes(x=group,y=logUMI,fill=group)) +
   geom_boxplot() + geom_jitter(width = 0.1,alpha = 0.2)