library(scatterpie)
library(ggplot2)

meta = read.csv('/scRNA_merge/pieplot/metadata.csv',row.names = 1)
head(meta)

celltype = c('')

meta[,celltype]

xylab = read.csv('/barcode_x_y.csv',row.names = 1)
head(xylab)

newx = rep(seq(2,24,by = 2),c(rep(300,12)))
newy = rep(rep(rev(seq(2,24,by = 2)),c(rep(5,12))),60)

xylab$newx = newx
xylab$newy = newy

head(xylab)

df <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(df) <- c('')
df[1,] = c(rep(0,16))
df

celltype = c('')

for(x in seq(2,24,by = 2)){
    for(y in seq(2,24,by = 2)){
        barcode = xylab[which(xylab$newx==x & xylab$newy==y),]$barcode
        df = rbind(df,c(x,y,as.data.frame(colSums(meta[barcode,celltype]))[,1]))
    }
}
df = df[-1,]

head(df)

col20 = c('#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#8C564BFF','#E377C2FF','#7F7F7FFF','#BCBD22FF','#17BECFFF','#AEC7E8FF','#FFBB78FF','#98DF8AFF','#FF9896FF','#C5B0D5FF','#C49C94FF','#F7B6D2FF','#C7C7C7FF', '#DBDB8DFF','#9EDAE5FF')

p = ggplot()+
  geom_scatterpie(data=df,
                  aes(x,y,r=0.9),
                  cols = c(''))+
  coord_equal()+
  theme_void()+
  theme(legend.position = "bottom",
       panel.border = element_blank())+
  scale_fill_manual(values = col20[1:14])
p

ggsave(p,filename = '/scRNA_merge/pieplot/piePlot.pdf')