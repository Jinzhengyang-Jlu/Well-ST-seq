library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyverse)
library(GEOquery)
library(data.table)
library(ggthemes)

total = read.csv('/benchmark_data/total.csv')
table(total$type)

df = subset(total,total$type %in% c('scRNA','ST_10um') )
write.csv(df,'/scRNA_merge/UMIcompare/df.csv')

p = ggplot(df, aes(x=type, y=nCount_RNA,fill = type)) + 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(outlier.colour = NA)+ theme_base()+scale_y_continuous(limits= c(0, 13000))+
    scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF'))
p
ggsave(p,filename = '/scRNA_merge/UMIcompare/UMIcompare_nooutlier.pdf',width = 16,height = 13)

p = ggplot(df, aes(x=type, y=log2(nCount_RNA),fill = type)) + 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(outlier.colour = NA)+ theme_base()+scale_y_continuous(limits= c(8, 15))+
    scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF'))
p
ggsave(p,filename = '/scRNA_merge/UMIcompare/UMIcompare_log2_nooutlier.pdf',width = 16,height = 13)

p = ggplot(df, aes(x=type, y=log10(nCount_RNA),fill = type)) + 
    stat_boxplot(geom = "errorbar")+
    geom_boxplot(outlier.colour = NA)+ theme_base()+scale_y_continuous(limits= c(2,5))+
    scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF'))
p
ggsave(p,filename = '/scRNA_merge/UMIcompare/UMIcompare_log10_nooutlier.pdf',width = 16,height = 13)

df = readRDS('/scRNA_merge/mouse_data/mouse.RDS')
STus = Read10X(data.dir = "/S10/matrix/")
STus <- CreateSeuratObject(counts = STus, project = "STus")
ST10um = STus

length(rownames(df@assays$RNA@counts))
length(rownames(STus@assays$RNA@counts))

library(ggVennDiagram)
library(ggplot2)
mydata<-list(scRNA = rownames(df@assays$RNA@counts),ST_10um = rownames(STus@assays$RNA@counts))

# method2
x = mydata
venn <- Venn(x)
data <- process_data(venn)
p = ggplot() +
  geom_sf(aes(fill = count), 
          data = venn_region(data)) +
  geom_sf(color="grey", 
          size = 1, 
          data = venn_setedge(data), 
          show.legend = FALSE) +
  scale_fill_gradient(low="#1F77B4FF",high = "#FF7F0EFF",name = "gene count")+#'#1F77B4FF','#FF7F0EFF' low="white",high = "#b9292b",
  geom_sf_text(aes(label = name), 
               data = venn_setlabel(data),
               size = 5) +
  geom_sf_label(aes(label = count), 
                data = venn_region(data),
                size = 4) +
  theme_void()
p
ggsave(p,file="/scRNA_merge/figures/Venn.pdf",width = 12,height = 10)

library(VennDiagram)

length(mydata$scRNA)
length(mydata$ST_10um)

###指定分组，这里是两组mRNA_target 和de_mRNA

setwd('/scRNA_merge/figures/')
venn_list <- list(scRNA = rownames(df@assays$RNA@counts),ST_10um = rownames(STus@assays$RNA@counts))

###画图

venn.diagram(venn_list, #分组

filename = 'venn_refine.png', #输出文件名

imagetype = 'png', #输出文件类型

euler.d = TRUE,

scaled = F, #两个圆圈一样大小 scaled = T时按照数量比例来画圈的大小，如果数量相差很大建议画相同大小的圈

fill=c('#1F77B4FF','#FF7F0EFF'), #圈的填充颜色

alpha=c(1,1), ##圆圈颜色的透明度，建议设置不同的，因为论文打印灰度图时可以看出不同

cex=2, ##区域内部数字的字体大小，即个数的字体

cat.cex = 1.5, ##分类名称字体大小，即组别的字体

cat.pos = c(-5,5), ##分类名称在圆的位置，默认正上方从c(0,0)，通过角度进行调整

margin = 0.1); #Number giving the amount of whitespace around the diagram in grid units
length(intersect(mydata$scRNA,mydata$ST_10um))

venn.plot <- draw.pairwise.venn(  area1 = 27933,  #区域1的数 
area2 = 31253,   #区域2的数 
cross.area = 25707,  #重叠的个数 
category = c( "ST_10um","scRNA"),#分类命名
fill = c('#1F77B4FF','#FF7F0EFF'),#1 2 区域分别的填充颜色 #'#1F77B4FF','#FF7F0EFF'
lty = "blank",  #1 2 区域的边框线类型 
cex = 2,        #1 2 区域内部数字的字体大小 
cat.cex = 2,    # 分类名称的字体大小 
cat.dist = 0.14,   #分类名称距离边的距离 实际调整 
#cat.just = list(c(-1, -1), c(1, 1)),  #分类名称的位置  ，圈内或者圈外
scaled = F,
alpha=c(0.6,0.6),
ext.line.lty = "dashed" )  #外部线为虚线);
pdf("venn-pdf.pdf")
grid.draw(venn.plot)
dev.off()