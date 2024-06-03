require(clusterProfiler)
require(org.Hs.eg.db)
library(BiocGenerics)
library(AnnotationHub)
library(clusterProfiler)
#library(org.Hs.eg.db)
library(org.Mm.eg.db) ##加载小鼠
keytypes(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library("ggsci")

group_markers <- read.csv('/Anno_DEGs_all.csv' ,header = T)
result = list()
for(i in colnames(group_markers)){
    mklist = bitr(group_markers[,i], fromType="SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
    ego <- enrichGO(mklist$ENTREZID,OrgDb = "org.Mm.eg.db",
        ont = "BP",
        readable = T)
    result[[i]] <- ego@result #%>% slice_min(p.adjust, n = 10)#[1:10,]
    #result[[i]] <- result[[i]][1:10,]
    result[[i]]$type <- i
}

result = do.call(rbind,result)
write.csv(result,'/GOtable.csv')

go_enrich_df <- result[,c('ID','Description','Count','type')]
go_enrich_df$number <- rev(factor(1:nrow(go_enrich_df)))
head(go_enrich_df)

CPCOLS <- c('#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#8C564BFF','#E377C2FF','#7F7F7FFF','#BCBD22FF','#17BECFFF',
           '#AEC7E8FF','#FFBB78FF','#98DF8AFF','#FF9896FF','#C5B0D5FF','#C49C94FF','#F7B6D2FF','#C7C7C7FF', '#DBDB8DFF','#9EDAE5FF')
p <- ggplot(data=go_enrich_df, aes(x=number, y=Count, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  scale_fill_manual(values = CPCOLS) + theme_bw() +
  scale_x_discrete(labels=rev(go_enrich_df$Description))+ 
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
p
ggsave(p,file = '/GO_barPlot.pdf',height = 15,width = 13)


# compare GO
group_g = as.data.frame(group_markers[,1])
group_g$cluster = colnames(group_markers)[1]
colnames(group_g) = c('gene','cluster')

for(i in 2:length(colnames(group_markers))){
    temp = as.data.frame(group_markers[,i])
    temp$cluster = colnames(group_markers)[i]
    colnames(temp) = c('gene','cluster')
    group_g = rbind(group_g,temp)
}

#group_g <- group_markers
tmp <- bitr(group_g$gene, fromType="SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
de_gene_clusters=merge(tmp,group_g,by.x='SYMBOL',by.y='gene')
# Run full GO enrichment test
formula_res <- compareCluster(
        ENTREZID~cluster, 
        data=de_gene_clusters, 
        fun="enrichGO", 
        OrgDb="org.Mm.eg.db",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.01,
        qvalueCutoff  = 0.05
)

# Run GO enrichment test and merge terms 
# that are close to each other to remove result redundancy
lineage1_ego <- clusterProfiler::simplify(
        formula_res, 
        cutoff=0.7, 
        by="p.adjust", 
        select_fun=min
)
enrichplot::dotplot(lineage1_ego, showCategory=5)

saveRDS(lineage1_ego,'GO.rds')

lineage1_ego = readRDS('GO.rds')

p = enrichplot::dotplot(lineage1_ego,
                        color = "-log10(pvalue)",
                        label_format = 70,
                        showCategory=10)
p
ggsave(p,filename = 'Compare_GO.pdf',width = 15,height = 18)

write.csv(lineage1_ego@compareClusterResult, 
          file="clusters_GO_simplify.csv")
write.csv(formula_res@compareClusterResult, 
          file="clusters_GO_all.csv")

# refine
picked = read.csv('GO_picked.csv',header = 0)
picked_10um$V1
GOtable[which(GOtable$Description %in% picked$V1),]
GOtable = read.csv('clusters_GO_simplify.csv',row.names = 1)

GO = subset(GOtable,GOtable$Description %in% picked$V1)
test = readRDS('GO.rds')
test@compareClusterResult = GO

loc = match(picked$V1,GO$Description)
GO = GO[loc,]

p = dotplot(
  test,
  x = "Cluster",
  color = "-log10(pvalue)",
  showCategory = picked$V1,#picked$V1,#
  size = NULL,
  split = NULL,
  font.size = 10,
  title = "",
  label_format = 50
)+scale_size('Fraction',limits=c(0.01,0.2))+#
theme(axis.text.x = element_text(angle = 45, hjust = 1))#+scale_x_continuous(labels = unique(GO_10$Cluster))
p
ggsave(p,file = 'GO_picked_dotplot.pdf',width = 8,height = 6)

