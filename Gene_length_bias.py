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
library("tidyverse")
set.seed(2022)

dat1 <- read.delim('/S10_UMI_gene.txt', sep = '')
dat2 <- read.delim('/S20_UMI_gene.txt', sep = '')
dat1 <- subset(dat1, UMI != 0)
dat2 <- subset(dat2, UMI != 0)

ref <- read.delim('/ref_10x_mm10/genes/hgTables.txt', header = F)
ref <- ref[,1:4]

reflist <- read.delim('/ref_10x_mm10/genes/gencode.vM10.metadata.MGI', header=F)
colnames(reflist) <- c('ENSEMBL','SYMBOL')
#reflist

ID_list=reflist[match(ref$V4,reflist[,"ENSEMBL"]),]
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

# named
dat$group<-0
for(i in 1:nrow(dat)){
    if(dat$bp[i] %in% seq(0,1564)){
        dat$group[i] <- 'bin1'
        }
    else if(dat$bp[i] %in% seq(1564,2218)){
        dat$group[i] <- 'bin2'
        }
    else if(dat$bp[i] %in% seq(2218,2799)){
        dat$group[i] <- 'bin3'
        }
    else if(dat$bp[i] %in% seq(2799,3394)){
        dat$group[i] <- 'bin4'
        }
    else if(dat$bp[i] %in% seq(3394,4044)){
        dat$group[i] <- 'bin5'
        }
    else if(dat$bp[i] %in% seq(4044,4743)){
        dat$group[i] <- 'bin6'
        }
    else if(dat$bp[i] %in% seq(4743,5658)){
        dat$group[i] <- 'bin7'
        }
    else if(dat$bp[i] %in% seq(5658,6946)){
        dat$group[i] <- 'bin8'
        }
    else if(dat$bp[i] %in% seq(6946,9102)){
        dat$group[i] <- 'bin9'
        }
    else{
        dat$group[i] <- 'bin10'
        }
}
#dat

# named
dat$range<-0
for(i in 1:nrow(dat)){
    if(dat$bp[i] %in% seq(0,1564)){
        dat$range[i] <- '0-1564'
        }
    else if(dat$bp[i] %in% seq(1564,2218)){
        dat$range[i] <- '1564-2218'
        }
    else if(dat$bp[i] %in% seq(2218,2799)){
        dat$range[i] <- '2218-2799'
        }
    else if(dat$bp[i] %in% seq(2799,3394)){
        dat$range[i] <- '2799-3394'
        }
    else if(dat$bp[i] %in% seq(3394,4044)){
        dat$range[i] <- '3394-4044'
        }
    else if(dat$bp[i] %in% seq(4044,4743)){
        dat$range[i] <- '4044-4743'
        }
    else if(dat$bp[i] %in% seq(4743,5658)){
        dat$range[i] <- '4743-5658'
        }
    else if(dat$bp[i] %in% seq(5658,6946)){
        dat$range[i] <- '5658-6946'
        }
    else if(dat$bp[i] %in% seq(6946,9102)){
        dat$range[i] <- '6946-9102'
        }
    else{
        dat$range[i] <- '9102-1642745'
        }
}
#dat

dat$logUMI <- log10(dat$UMI+1)
dat$group <- factor(dat$group,levels = c('bin1' ,'bin2' , 'bin3' , 'bin4' , 'bin5' , 'bin6' , 'bin7' , 'bin8' , 'bin9','bin10' ))
p = ggplot(dat, aes(x=range,y=logUMI)) +
   geom_boxplot(fill = 'grey')  +
    scale_y_continuous(name = "logUMI")+
    scale_x_discrete(name="Gene length")+
      theme_classic()+  
theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = 'black', lineend = 'round'),
        axis.text.x = element_text(size = 10, color = 'black',angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.text = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 10, color = 'black'),
        axis.title.x = element_text(size = 10, color = 'black', face="bold"),
        axis.title.y = element_text(size = 10, color = 'black', face="bold"))+
theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p
ggsave(p,filename = '/accumulation_curves/S10.pdf', width = 11, height = 9)

dat=ref[match(dat2$gene,ref[,"SYMBOL"]),]
dat$bp <- dat$V3 - dat$V2
dat <- dat[,5:6]
dat$UMI <- dat2$UMI
dat <- dat[complete.cases(dat),]

max(dat$bp)
min(dat$bp)

# named
dat$group<-0
for(i in 1:nrow(dat)){
    if(dat$bp[i] %in% seq(0,1564)){
        dat$group[i] <- 'bin1'
        }
    else if(dat$bp[i] %in% seq(1564,2218)){
        dat$group[i] <- 'bin2'
        }
    else if(dat$bp[i] %in% seq(2218,2799)){
        dat$group[i] <- 'bin3'
        }
    else if(dat$bp[i] %in% seq(2799,3394)){
        dat$group[i] <- 'bin4'
        }
    else if(dat$bp[i] %in% seq(3394,4044)){
        dat$group[i] <- 'bin5'
        }
    else if(dat$bp[i] %in% seq(4044,4743)){
        dat$group[i] <- 'bin6'
        }
    else if(dat$bp[i] %in% seq(4743,5658)){
        dat$group[i] <- 'bin7'
        }
    else if(dat$bp[i] %in% seq(5658,6946)){
        dat$group[i] <- 'bin8'
        }
    else if(dat$bp[i] %in% seq(6946,9102)){
        dat$group[i] <- 'bin9'
        }
    else{
        dat$group[i] <- 'bin10'
        }
}
#dat

# named
dat$range<-0
for(i in 1:nrow(dat)){
    if(dat$bp[i] %in% seq(0,1564)){
        dat$range[i] <- '0-1564'
        }
    else if(dat$bp[i] %in% seq(1564,2218)){
        dat$range[i] <- '1564-2218'
        }
    else if(dat$bp[i] %in% seq(2218,2799)){
        dat$range[i] <- '2218-2799'
        }
    else if(dat$bp[i] %in% seq(2799,3394)){
        dat$range[i] <- '2799-3394'
        }
    else if(dat$bp[i] %in% seq(3394,4044)){
        dat$range[i] <- '3394-4044'
        }
    else if(dat$bp[i] %in% seq(4044,4743)){
        dat$range[i] <- '4044-4743'
        }
    else if(dat$bp[i] %in% seq(4743,5658)){
        dat$range[i] <- '4743-5658'
        }
    else if(dat$bp[i] %in% seq(5658,6946)){
        dat$range[i] <- '5658-6946'
        }
    else if(dat$bp[i] %in% seq(6946,9102)){
        dat$range[i] <- '6946-9102'
        }
    else{
        dat$range[i] <- '9102-1642745'
        }
}
#dat

dat$logUMI <- log10(dat$UMI+1)
table(dat$group)
dat$group <- factor(dat$group,levels = c('bin1' ,'bin2' , 'bin3' , 'bin4' , 'bin5' , 'bin6' , 'bin7' , 'bin8' , 'bin9','bin10' ))

table(dat$range)

p = ggplot(dat, aes(x=range,y=logUMI)) +
   geom_boxplot(fill = 'grey')  +
    scale_y_continuous(name = "logUMI")+
    scale_x_discrete(name="Gene length")+
theme_classic()+
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = 'black', lineend = 'round'),
        legend.position = 'right',
        axis.text.x = element_text(size = 15, color = 'black',angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'),
        axis.title.x = element_text(size = 15, color = 'black', face="bold"),
        axis.title.y = element_text(size = 15, color = 'black', face="bold"))+
theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p

p = ggplot(dat, aes(x=range,y=logUMI)) +
   geom_boxplot(fill = 'grey')  +
    scale_y_continuous(name = "logUMI")+
    scale_x_discrete(name="Gene length")+
      theme_classic()+  
theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = 'black', lineend = 'round'),
        axis.text.x = element_text(size = 10, color = 'black',angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.text = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 10, color = 'black'),
        axis.title.x = element_text(size = 10, color = 'black', face="bold"),
        axis.title.y = element_text(size = 10, color = 'black', face="bold"))+
theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p
ggsave(p,filename = '/accumulation_curves/S20.pdf', width = 11, height = 9)