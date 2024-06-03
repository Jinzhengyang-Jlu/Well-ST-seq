library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyverse)
library(GEOquery)
library(data.table)

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
saveRDS(rat,'benchmark_data/scRNA/rat.rds')

rat <- readRDS('/benchmark_data/scRNA/rat.rds')
df = rat
colnames(rat@meta.data)
table(rat$orig.ident)
table(rat$celltype)

df@active.ident <- as.factor(df$orig.ident)
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size = 0)

#进行SCT标准化：
df <- SCTransform(df, vars.to.regress = "percent.mt", verbose = FALSE)

#使用Seurat PCA对数据进行降维：
df <- RunPCA(df, features = VariableFeatures(object = df))
#之后进行UMAP聚类, SCT标准化一般选前30个PC进行聚类：
df <- RunUMAP(df, dims = 1:30, verbose=FALSE)
df <- FindNeighbors(df, dims = 1:30, berbose=FALSE)
df <- FindClusters(df, verbose = FALSE)

saveRDS(df,'/scRNA_merge/mouse_data/mouse.RDS')

p = DimPlot(df, label = T,group.by = 'celltype')
ggsave(p,file = '/scRNA_merge/mouse_data/UMAP_celltype.pdf',width = 10,height = 7)

col50 = c('#E78AC3','#E6AB02','#7FC97F','#FFED6F','#1B9E77','#FBB4AE','#BEAED4','#984EA3','#FDB462','#DECBE4','#A6D854','#CAB2D6','#8DA0CB','#B3CDE3','#FFF2AE','#E5D8BD','#1F78B4','#FB9A99','#E6F5C9','#CBD5E8','#B3E2CD','#FDCDAC','#6A3D9A','#8DD3C7','#F781BF','#A6CEE3','#377EB8','#FDC086','#E5C494','#4DAF4A','#FFFF99','#FF7F00','#FFD92F','#BC80BD','#E7298A','#66C2A5','#FB8072','#CCEBC5','#FFFFCC','#FCCDE5','#F1E2CC','#7570B3','#FDDAEC','#FDBF6F','#FC8D62','#FED9A6','#66A61E','#A65628','#80B1D3','#FFFFB3','#D95F02')
col20 = c('#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#8C564BFF','#E377C2FF','#7F7F7FFF','#BCBD22FF','#17BECFFF','#AEC7E8FF','#FFBB78FF','#98DF8AFF','#FF9896FF','#C5B0D5FF','#C49C94FF','#F7B6D2FF','#C7C7C7FF', '#DBDB8DFF','#9EDAE5FF','#4A8EBD','#7B539E','#1BA624','#EB852F','#F2E4B6','#A3DA7F','#B39B7C','#68A39A','#A45D28','#F9BCC3','#70CFD5','#E72D85','#BAABCF','#E1CA86','#BA84B3','#B0D1DC')

p = DimPlot(df, label = T,group.by = 'celltype',cols = col20)
p
ggsave(p,filename = '/scRNA_merge/figures/sc_celltype_Umap.pdf',width = 11,height = 8)

p = DimPlot(df, label = T,group.by = 'seurat_clusters',cols = col50)
p
ggsave(p,filename = '/scRNA_merge/figures/sc_cluster_umap.pdf',width = 9,height = 8)

#to h5ad
library(Seurat)
library(SeuratData)
library(SeuratDisk)

setwd('/scRNA_merge/mouse_data')
df = readRDS('/scRNA_merge/mouse_data/mouse.RDS')

mat <- read.table(file = "/GSE104323_10X_expression_data_V2.tab.gz" , header = T, row.names=1, sep="", as.is=T)
df@assays$RNA@counts = as(as.matrix(mat), "dgCMatrix")
DefaultAssay(df) <- 'RNA'#change to save unnormalized data to make cell2location work
df@assays$RNA@data = df@assays$RNA@counts
df <- DietSeurat(
    df,
    count = TRUE,
    data = TRUE,
    scale.data = FALSE,
    #assays = 'RNA',
    dimreducs = c('pca', 'umap'),
    graphs = NULL
    )
#df@assays$RNA@data = df@assays$RNA@counts
SaveH5Seurat(df, filename = "mouseunlog2.h5Seurat",overwrite=T)
Convert("mouseunlog2.h5Seurat", dest = "h5ad")

#merge
library(harmony)
sc = readRDS('/scRNA_merge/mouse_data/mouse.RDS')
st = readRDS('/scRNA_merge/s10/s10.rds')
sc$data = 'scRNA_mouse'
st$data = 'ST_10um'

scSTorg = merge(sc,st)
scSTorg <- NormalizeData(scSTorg, normalization.method = "LogNormalize", scale.factor = 10000)
scSTorg <- FindVariableFeatures(scSTorg, selection.method = "vst", nfeatures = 2000)
scSTorg <- ScaleData(scSTorg, features = rownames(scSTorg))
scSTorg <- RunPCA(scSTorg, features = VariableFeatures(object = scSTorg))

gc()
st$orig.ident = 'ST_10um'

# 首先，根据meta信息中不同的样本，分组信息，或测序技术对Seurat对象进行分割，构建不同的数据集；
split_seurat <- SplitObject(scSTorg, split.by = "orig.ident")
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt"))
}

# 寻找高变基因
integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = 3000)

# 准备好SCT列表对象进行聚合
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)

# 寻找最佳伙伴 —— 需要一定时间运行
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)
# 跨条件聚合，将这些识别好的anchors传递给IntegrateData函数，整合后的数据返一个Seurat对象，该对象中将包含一个新的Assay（integrated），里面存储了整合后表达矩阵，原始的表达矩阵存储在RNA这个Assay中
df <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

# 运行PCA
df <- RunPCA(object = df)

# 聚类分群
df <- df %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.1) %>%
  identity()
DimPlot(df, label = T)
saveRDS(df,'/scRNA_merge/ST_to_sc/seurat.rds')

df = readRDS('/scRNA_merge/ST_to_sc/seurat.rds')
df@meta.data[which(df@meta.data$data == 'scRNA_mouse'),]$data = 'scRNA'
df$data_type = df$data
p = DimPlot(df, label = F,group.by = 'data_type',cols = c('#1F77B4FF','#FF7F0EFF'))
p
ggsave(p,filename = '/scRNA_merge/ST_to_sc/data_type.pdf')

p2 = DimPlot(df, label = T,group.by = 'celltype')#,split.by = 'data'
p2
ggsave(p2,filename = '/scRNA_merge/ST_to_sc/celltype.pdf',width = 10,height = 8)
p3 = DimPlot(df, label = T,group.by = 'seurat_clusters')#,split.by = 'data'
p3
ggsave(p3,filename = '/scRNA_merge/ST_to_sc/seurat_clusters.pdf')#,width = 9,height = 8
p4 = DimPlot(df, label = T,group.by = 'celltype',split.by = 'data')#
p4
ggsave(p4,filename = '/scRNA_merge/ST_to_sc/celltype_split.pdf',width = 14,height = 8)