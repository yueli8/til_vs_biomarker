library(pheatmap)
setwd("~/geo/gse162498")
#pheatmap
mat<-read.table("genes11_pheatmap.txt", head=TRUE, row.names = 1, sep="\t", quote="")
x<-as.matrix(mat)
#x<-t(x)
pheatmap(log((x+1),2),cellwidth=25, cellheight=20, cluster_cols=T,cluster_rows = T, 
        color=colorRampPalette(c("green", "black", "red"))(100),fontsize=15)



