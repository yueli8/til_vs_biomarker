library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(clustree)
setwd("~/gse162498")
#blood samples
a1 <- Read10X(data.dir = "~/gse162498/P57_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P57_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P57_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P57_B.rds")

a1 <- Read10X(data.dir = "~/gse162498/P58_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P58_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P58_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P58_B.rds")

a1 <- Read10X(data.dir = "~/gse162498/P60_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P60_B.rds")

a1 <- Read10X(data.dir = "~/gse162498/P60_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P60_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P61_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P61_B.rds")

a1 <- Read10X(data.dir = "~/gse162498/P61_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P61_T.rds")

P57_B<-readRDS(file="P57_B.rds")
P57_T<-readRDS(file="P57_T.rds")
P58_B<-readRDS(file="P58_B.rds")
P58_T<-readRDS(file="P58_T.rds")
P60_B<-readRDS(file="P60_B.rds")
P60_T<-readRDS(file="P60_T.rds")
P61_B<-readRDS(file="P61_B.rds")
P61_T<-readRDS(file="P61_T.rds")

P57_B<-RenameCells(P57_B,add.cell.id="P57_B",for.merge=T)
P57_B@meta.data$tech<-"Blood"
P57_B@meta.data$celltype<-"Blood_P57"

P57_T<-RenameCells(P57_T,add.cell.id="P57_T",for.merge=T)
P57_T@meta.data$tech<-"Tumor"
P57_T@meta.data$celltype<-"Tumor_P57"

P58_B<-RenameCells(P58_B,add.cell.id="P58_B",for.merge=T)
P58_B@meta.data$tech<-"Blood"
P58_B@meta.data$celltype<-"Blood_P58"

P58_T<-RenameCells(P58_T,add.cell.id="P58_T",for.merge=T)
P58_T@meta.data$tech<-"Tumor"
P58_T@meta.data$celltype<-"Tumor_P58"

P60_B<-RenameCells(P60_B,add.cell.id="P60_B",for.merge=T)
P60_B@meta.data$tech<-"Blood"
P60_B@meta.data$celltype<-"Blood_P60"

P60_T<-RenameCells(P60_T,add.cell.id="P60_T",for.merge=T)
P60_T@meta.data$tech<-"Tumor"
P60_T@meta.data$celltype<-"Tumor_P60"

P61_B<-RenameCells(P61_B,add.cell.id="P61_B",for.merge=T)
P61_B@meta.data$tech<-"Blood"
P61_B@meta.data$celltype<-"Blood_P61"

P61_T<-RenameCells(P61_T,add.cell.id="P61_T",for.merge=T)
P61_T@meta.data$tech<-"Tumor"
P61_T@meta.data$celltype<-"Tumor_P61"

P5758_T<-merge(P57_T,P58_T)
P6061_T<-merge(P60_T,P61_T)
T57586061<-merge(P5758_T,P6061_T)

P5758_B<-merge(P57_B,P58_B)
P6061_B<-merge(P60_B,P61_B)
B57586061<-merge(P5758_B,P6061_B)

library(Seurat)
library(ggplot2)

saveRDS(B57586061, file="B57586061_before_integrate.rds")
hms<-B57586061

hms<-readRDS(file="B57586061_before_integrate.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
#reference.list <- pancreas.list[c("Blood_P57","Tumor_P57","Blood_P58","Tumor_P58", 
#"Blood_P60","Tumor_P60", "Blood_P61", "Tumor_P61","Juxta_P60","Juxta_P61")]
reference.list <- pancreas.list[c("Blood_P57","Blood_P58", "Blood_P60", "Blood_P61" )]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
#saveRDS(pancreas.integrated, file = "hms_after_integrated01.rds")
saveRDS(pancreas.integrated, file = "B57586061_after_integrated.rds")

hms_individual_integrated<-readRDS(file="B57586061_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1

#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster,label = TRUE, reduction = "umap")
saveRDS(hms_cluster, file = "Blood_cluster_test.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.001, 
                                logfc.threshold = 0.001
)
write.table(scRNA.markers,file="blood_Markers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="blood_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="blood_genes.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="Blood_cluster_test.rds")
DimPlot(hms_cluster, label=TRUE,reduction = "umap")

new.cluster.ids <- c("naive_cd4", "memory_cd4","nkt", "nkt","Effect_memory_cd4", "effect_cd8","effect_cd8", 
                     "treg","nkt", "naive_cd8", "nkt", "nkt","effect_cd8", "effect_cd8","nkt", "nkt",
                     "Effect_memory_cd4","Pre_exhausted_cd8","monocyte", "Effect_memory_cd4","epithelial",
                     "monocyte","B", "secretory","cycling_NK","megakaryocyte") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "blood_cluster_id_test.rds")

hms_cluster_id<-readRDS("blood_cluster_id_test.rds")

setwd("~/gse162498/CytoTRACE01")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="blood_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4,subsamplesize = 1000)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p

#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="naive_cd4", "1"="memory_cd4","2"="nkt", "3"="nkt","4"="Effect_memory_cd4","5"= "effect_cd8","6"="effect_cd8", 
                                                 "7"="treg","8"="nkt", "9"="naive_cd8","10"= "nkt", "11"="nkt","12"="effect_cd8", "13"="effect_cd8","14"="nkt","15" ="nkt",
                                                "16"= "Effect_memory_cd4","17"="Pre_exhausted_cd8","18"="monocyte","19"= "Effect_memory_cd4","20"="epithelial",
                                                  "21"="monocyte","22"="B","23"= "secretory","24"="cycling_NK","25"="megakaryocyte") 
colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)
plotCytoTRACE(result01, phenotype = phe)

#tumor samples
library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(clustree)
setwd("~/gse162498")

a1 <- Read10X(data.dir = "~/gse162498/P57_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P57_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P58_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P58_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P60_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P60_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P61_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P61_T.rds")

P57_T<-readRDS(file="P57_T.rds")
P58_T<-readRDS(file="P58_T.rds")
P60_T<-readRDS(file="P60_T.rds")
P61_T<-readRDS(file="P61_T.rds")

P57_T<-RenameCells(P57_T,add.cell.id="P57_T",for.merge=T)
P57_T@meta.data$tech<-"Tumor"
P57_T@meta.data$celltype<-"Tumor_P57"

P58_T<-RenameCells(P58_T,add.cell.id="P58_T",for.merge=T)
P58_T@meta.data$tech<-"Tumor"
P58_T@meta.data$celltype<-"Tumor_P58"

P60_T<-RenameCells(P60_T,add.cell.id="P60_T",for.merge=T)
P60_T@meta.data$tech<-"Tumor"
P60_T@meta.data$celltype<-"Tumor_P60"

P61_T<-RenameCells(P61_T,add.cell.id="P61_T",for.merge=T)
P61_T@meta.data$tech<-"Tumor"
P61_T@meta.data$celltype<-"Tumor_P61"

P5758_T<-merge(P57_T,P58_T)
P6061_T<-merge(P60_T,P61_T)
T57586061<-merge(P5758_T,P6061_T)

library(Seurat)
library(ggplot2)

saveRDS(T57586061, file="T57586061_before_integrate.rds")
hms<-T57586061

hms<-readRDS(file="T57586061_before_integrate.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
#reference.list <- pancreas.list[c("Blood_P57","Tumor_P57","Blood_P58","Tumor_P58", 
#"Blood_P60","Tumor_P60", "Blood_P61", "Tumor_P61","Juxta_P60","Juxta_P61")]
reference.list <- pancreas.list[c("Tumor_P57","Tumor_P58", "Tumor_P60", "Tumor_P61" )]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
#saveRDS(pancreas.integrated, file = "hms_after_integrated01.rds")
saveRDS(pancreas.integrated, file = "T57586061_after_integrated.rds")

hms_individual_integrated<-readRDS(file="T57586061_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1

#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster,label = TRUE, reduction = "umap")
saveRDS(hms_cluster, file = "Tumor_cluster_test.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.001, 
                                logfc.threshold = 0.001
)
write.table(scRNA.markers,file="tumor_Markers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="tumor_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="tumor_genes.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="Tumor_cluster_test.rds")
DimPlot(hms_cluster, label=TRUE,reduction = "umap")
new.cluster.ids <- c("naive_cd4","nkt", "epithelial","treg", "exhausted_cd8","exhausted_cd4","effect_memory_cd8", "exhausted_cd8",
                   "nkt","macrophage", "epithelial","epithelial","cycling NK","naive_cd4", "secretory", "macrophage",
                   "exhausted_cd8", "macrophage", "epithelial","pre_exhauste_cd8","B","macrophage","macrophage","B",
                    "treg","secretory")

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "tumor_cluster_id_test.rds")

hms_cluster_id<-readRDS("tumor_cluster_id_test.rds")

setwd("~/gse162498/CytoTRACE01")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="tumor_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4,subsamplesize = 1000)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p

#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
"0"="naive_cd4","1"="nkt","2"= "epithelial","3"="treg","4"= "exhausted_cd8","5"="exhausted_cd4","6"="effect_memory_cd8",
"7"= "exhausted_cd8", "8"="nkt","9"="macrophage", "10"="epithelial","11"="epithelial","12"="cycling NK","13"="naive_cd4", 
"14"="secretory", "15"="macrophage","16"="exhausted_cd8","17"= "macrophage","18"= "epithelial","19"="pre_exhauste_cd8",
"20"="B","21"="macrophage","22"="macrophage","23"="B","24"="treg","25"="secretory" ) 
colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)
plotCytoTRACE(result01, phenotype = phe)

