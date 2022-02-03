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

setwd("~/geo/gse162498")

a1 <- Read10X(data.dir = "~/geo/gse162498/P34_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P34_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P35_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P35_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P42_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P42_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P43_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P43_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P46_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P46_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P47_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P47_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P55_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P55_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P57_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P57_B.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P57_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P57_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P58_Blood_filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P58_B.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P58_Tumor_filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P58_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P60_Blood_filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P60_B.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P60_Juxta_ffiltered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P60_J.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P60_Tumor_filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P60_T.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P61_Blood_filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P61_B.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P61_Juxta_filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P61_J.rds")

a1 <- Read10X(data.dir = "~/geo/gse162498/P61_Tumor_filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "P61_T.rds")

P34_T<-readRDS(file="P34_T.rds")
P35_T<-readRDS(file="P35_T.rds")
P42_T<-readRDS(file="P42_T.rds")
P43_T<-readRDS(file="P43_T.rds")
P46_T<-readRDS(file="P46_T.rds")
P47_T<-readRDS(file="P47_T.rds")
P55_T<-readRDS(file="P55_T.rds")
P57_B<-readRDS(file="P57_B.rds")
P57_T<-readRDS(file="P57_T.rds")
P58_B<-readRDS(file="P58_B.rds")
P58_T<-readRDS(file="P58_T.rds")
P60_B<-readRDS(file="P60_B.rds")
P60_J<-readRDS(file="P60_J.rds")
P60_T<-readRDS(file="P60_T.rds")
P61_B<-readRDS(file="P61_B.rds")
P61_J<-readRDS(file="P61_J.rds")
P61_T<-readRDS(file="P61_T.rds")

P34_T<-RenameCells(P34_T,add.cell.id="P34_T",for.merge=T)
P34_T@meta.data$tech<-"Tumor"
P34_T@meta.data$celltype<-"Tumor_P34"

P35_T<-RenameCells(P35_T,add.cell.id="P35_T",for.merge=T)
P35_T@meta.data$tech<-"Tumor"
P35_T@meta.data$celltype<-"Tumor_P35"

P42_T<-RenameCells(P42_T,add.cell.id="P42_T",for.merge=T)
P42_T@meta.data$tech<-"Tumor"
P42_T@meta.data$celltype<-"Tumor_P42"

P43_T<-RenameCells(P43_T,add.cell.id="P43_T",for.merge=T)
P43_T@meta.data$tech<-"Tumor"
P43_T@meta.data$celltype<-"Tumor_P43"

P46_T<-RenameCells(P46_T,add.cell.id="P46_T",for.merge=T)
P46_T@meta.data$tech<-"Tumor"
P46_T@meta.data$celltype<-"Tumor_P46"

P47_T<-RenameCells(P47_T,add.cell.id="P47_T",for.merge=T)
P47_T@meta.data$tech<-"Tumor"
P47_T@meta.data$celltype<-"Tumor_P47"

P55_T<-RenameCells(P55_T,add.cell.id="P55_T",for.merge=T)
P55_T@meta.data$tech<-"Tumor"
P55_T@meta.data$celltype<-"Tumor_P55"

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

P60_J<-RenameCells(P60_J,add.cell.id="P60_J",for.merge=T)
P60_J@meta.data$tech<-"Juxta"
P60_J@meta.data$celltype<-"Juxta_P60"

P60_T<-RenameCells(P60_T,add.cell.id="P60_T",for.merge=T)
P60_T@meta.data$tech<-"Tumor"
P60_T@meta.data$celltype<-"Tumor_P60"

P61_B<-RenameCells(P61_B,add.cell.id="P61_B",for.merge=T)
P61_B@meta.data$tech<-"Blood"
P61_B@meta.data$celltype<-"Blood_P61"

P61_J<-RenameCells(P61_J,add.cell.id="P61_J",for.merge=T)
P61_J@meta.data$tech<-"Juxta"
P61_J@meta.data$celltype<-"Juxta_P61"

P61_T<-RenameCells(P61_T,add.cell.id="P61_T",for.merge=T)
P61_T@meta.data$tech<-"Tumor"
P61_T@meta.data$celltype<-"Tumor_P61"

P57<-merge(P57_B,P57_T)
P58<-merge(P58_B,P58_T)
P60<-merge(P60_B,P60_T)
P61<-merge(P61_B,P61_T)
P5758<-merge(P57,P58)
P6061<-merge(P60,P61)
P6061_J <- merge(P60_J, P61_J)
P57586061<-merge(P5758,P6061)
P57586061J<-merge(P57586061,P6061_J)

saveRDS(P57586061J, file="P57586061J_before_integrate.rds")
saveRDS(P57586061, file="P57586061_before_integrate.rds")
hms01<-P57586061J
hms<-P57586061

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
#reference.list <- pancreas.list[c("Blood_P57","Tumor_P57","Blood_P58","Tumor_P58", "Blood_P60","Tumor_P60", "Blood_P61", "Tumor_P61","Juxta_P60","Juxta_P61")]
reference.list <- pancreas.list[c("Blood_P57","Tumor_P57","Blood_P58","Tumor_P58",
                                  "Blood_P60","Tumor_P60", "Blood_P61", "Tumor_P61")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")
saveRDS(pancreas.integrated, file = "P57586061_after_integrated.rds")

hms_individual_integrated<-readRDS(file="P57586061_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1
#find how many 15cluster
ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:10)
hms_cluster <- FindClusters( hms_neighbor, resolution = 0.5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:10)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test.rds")

hms_cluster<-readRDS(file="hms_cluster_test.rds")
cluster0.markers <- FindMarkers(hms_cluster, ident.1=0, min.pcr=0.25)
head(cluster0.markers, n=10)
cluster1.markers <- FindMarkers(hms_cluster, ident.1=1, min.pcr=0.25)
head(cluster1.markers, n=10)
cluster2.markers <- FindMarkers(hms_cluster, ident.1=2, min.pcr=0.25)
head(cluster2.markers, n=10)
cluster3.markers <- FindMarkers(hms_cluster, ident.1=3, min.pcr=0.25)
head(cluster3.markers, n=10)
cluster4.markers <- FindMarkers(hms_cluster, ident.1=4, min.pcr=0.25)
head(cluster4.markers, n=10)
cluster5.markers <- FindMarkers(hms_cluster, ident.1=5, min.pcr=0.25)
head(cluster5.markers, n=10)
cluster6.markers <- FindMarkers(hms_cluster, ident.1=6, min.pcr=0.25)
head(cluster6.markers, n=10)
cluster7.markers <- FindMarkers(hms_cluster, ident.1=7, min.pcr=0.25)
head(cluster7.markers, n=10)
cluster8.markers <- FindMarkers(hms_cluster, ident.1=8, min.pcr=0.25)
head(cluster8.markers, n=10)
cluster9.markers <- FindMarkers(hms_cluster, ident.1=9, min.pcr=0.25)
head(cluster9.markers, n=10)
cluster10.markers <- FindMarkers(hms_cluster, ident.1=10, min.pcr=0.25)
head(cluster10.markers, n=10)
cluster11.markers <- FindMarkers(hms_cluster, ident.1=11, min.pcr=0.25)
head(cluster11.markers, n=10)

#new.cluster.ids <- c("Memory CD4+ T", "Natural Killer","CD8+ T", "Naive CD4+ T","Memory CD4+ T",
#                     "B", "Astrocyte","Dendritic Cell",
#                     "Dendritic Cell", "Dendritic Cell","B") 


new.cluster.ids <- c("Memory CD4+ T", "Memory CD4+ T", "Naive CD4+ T","Memory CD4+ T",
                     "Natural Killer T", "CD8+ T","B","Dendritic Cell","Astrocyte",
                     "Natural Killer T", "Dendritic Cell","B") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")
  
#hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
Natural_Killer_T <-subset(hms_cluster_id, idents=c('Natural Killer T'))
DimPlot(Natural_Killer_T, reduction = "umap")
saveRDS(Natural_Killer_T, file="Natural_Killer_T.rds")

Naive_CD4_T<-subset(hms_cluster_id, idents=c('Naive CD4+ T'))
DimPlot(Naive_CD4_T, reduction = "umap")
saveRDS(Naive_CD4_T, file="Naive_CD4_T.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Astrocyte<-subset(hms_cluster_id, idents=c('Astrocyte'))
DimPlot(Astrocyte, reduction = "umap")
saveRDS(Astrocyte, file="Astrocyte.rds")

Dendritic_cell<-subset(hms_cluster_id, idents=c('Dendritic Cell'))
DimPlot(Dendritic_cell, reduction = "umap")
saveRDS(Dendritic_cell, file="Dendritic_cell.rds")

CD8_T<-subset(hms_cluster_id, idents=c('CD8+ T'))
DimPlot(CD8_T, reduction = "umap")
saveRDS(CD8_T, file="CD8_T.rds")

Memory_CD4_T<-subset(hms_cluster_id, idents=c('Memory CD4+ T'))
DimPlot(Memory_CD4_T, reduction = "umap")
saveRDS(Memory_CD4_T, file="Memory_CD4_T.rds")

#input each cluster
Naive_CD4_T<-readRDS("Naive_CD4_T.rds")
Natural_Killer_T<-readRDS("Natural_Killer_T.rds")
Astrocyte<-readRDS("Astrocyte.rds")
B_cell<-readRDS("B.rds")
Memory_CD4_T<-readRDS("Memory_CD4_T.rds")
Dendritic_cell<-readRDS("Dendritic_cell.rds")
CD8_T<-readRDS("CD8_T.rds")
hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

#DimPlot
DimPlot(Naive_CD4_T, reduction = "umap", split.by = "tech")
DimPlot(Natural_Killer_T, reduction = "umap", split.by = "tech")
DimPlot(Astrocyte, reduction = "umap", split.by = "tech")
DimPlot(B_cell, reduction = "umap", split.by = "tech")
DimPlot(Memory_CD4_T, reduction = "umap", split.by = "tech")
DimPlot(Dendritic_cell, reduction = "umap", split.by = "tech")
DimPlot(CD8_T, reduction = "umap", split.by = "tech")

RidgePlot(CD8_T, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"),cols = c("green3","orangered"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5,label = FALSE)

genes11_heatmap<-DotPlot(hms_cluster_id,features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()
genes11_heatmap<-genes11_heatmap$data
write.csv(genes11_heatmap,"genes11_heatmap")

#expression level in each patients
CD8_T<-readRDS("CD8_T.rds")
a<-DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1")

