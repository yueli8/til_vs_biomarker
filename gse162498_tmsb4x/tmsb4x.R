library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(pheatmap)

setwd("~/gse162498")
hms_cluster_id<-readRDS("hms_cluster_id_test_correct.rds")

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2","TMSB4X")
DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
a<-DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
a1<-a$data
write.table(a1,"a1_heatmap01")

RidgePlot(hms_cluster_id, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2","TMSB4X"),cols = c("green3","orangered"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

setwd("~/gse162498")
cd38<-read.table("cd38",header=TRUE)
expr_groups_cd38<-group_by(cd38,Identity)
mean_cd38<-summarise(expr_groups_cd38,expr_mean=mean(Expression))
write.csv(mean_cd38,"mean_cd38.csv")

cd69<-read.table("cd69",header=TRUE)
expr_groups_cd69<-group_by(cd69,Identity)
mean_cd69<-summarise(expr_groups_cd69,expr_mean=mean(Expression))
write.csv(mean_cd69,"mean_cd69.csv")

cd8a<-read.table("cd8a",header=TRUE)
expr_groups_cd8a<-group_by(cd8a,Identity)
mean_cd8a<-summarise(expr_groups_cd8a,expr_mean=mean(Expression))
write.csv(mean_cd8a,"mean_cd8a.csv")

cd8b<-read.table("cd8b",header=TRUE)
expr_groups_cd8b<-group_by(cd8b,Identity)
mean_cd8b<-summarise(expr_groups_cd8b,expr_mean=mean(Expression))
write.csv(mean_cd8b,"mean_cd8b.csv")

entpd1<-read.table("entpd1",header=TRUE)
expr_groups_entpd1<-group_by(entpd1,Identity)
mean_entpd1<-summarise(expr_groups_entpd1,expr_mean=mean(Expression))
write.csv(mean_entpd1,"mean_entpd1.csv")

gzma<-read.table("gzma",header=TRUE)
expr_groups_gzma<-group_by(gzma,Identity)
mean_gzma<-summarise(expr_groups_gzma,expr_mean=mean(Expression))
write.csv(mean_gzma,"mean_gzma.csv")

gzmh<-read.table("gzmh",header=TRUE)
expr_groups_gzmh<-group_by(gzmh,Identity)
mean_gzmh<-summarise(expr_groups_gzmh,expr_mean=mean(Expression))
write.csv(mean_gzmh,"mean_gzmh.csv")

myo1f<-read.table("myo1f",header=TRUE)
expr_groups_myo1f<-group_by(myo1f,Identity)
mean_myo1f<-summarise(expr_groups_myo1f,expr_mean=mean(Expression))
write.csv(mean_myo1f,"mean_myo1f.csv")

syne1<-read.table("syne1",header=TRUE)
expr_groups_syne1<-group_by(syne1,Identity)
mean_syne1<-summarise(expr_groups_syne1,expr_mean=mean(Expression))
write.csv(mean_syne1,"mean_syne1.csv")

tsc22d3<-read.table("tsc22d3",header=TRUE)
expr_groups_tsc22d3<-group_by(tsc22d3,Identity)
mean_tsc22d3<-summarise(expr_groups_tsc22d3,expr_mean=mean(Expression))
write.csv(mean_tsc22d3,"mean_tsc22d3.csv")

xcl2<-read.table("xcl2",header=TRUE)
expr_groups_xcl2<-group_by(xcl2,Identity)
mean_xcl2<-summarise(expr_groups_xcl2,expr_mean=mean(Expression))
write.csv(mean_xcl2,"mean_xcl2.csv")

tmsb4x<-read.table("tmsb4x",header=TRUE)
expr_groups_tmsb4x<-group_by(tmsb4x,Identity)
mean_tmsb4x<-summarise(expr_groups_tmsb4x,expr_mean=mean(Expression))
write.csv(mean_tmsb4x,"mean_tmsb4x.csv")

mat<-read.table("tmp02", head=TRUE, row.names = 1, sep="\t", quote="")
x<-as.matrix(mat)
#x<-t(x)
pheatmap(log((x+1),2),cellwidth=25, cellheight=20, cluster_cols=T,cluster_rows = T, 
         color=colorRampPalette(c("green", "black", "red"))(100),fontsize=15)
