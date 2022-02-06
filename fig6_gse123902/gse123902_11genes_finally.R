#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("scran")

library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
setwd("~/geo/gse123902")

#input, CreateSeuratObject, filter, save
a<-"GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
#构建Seurat对象，这里会有个初筛，保证所有基因在至少3个细胞中表达（0.1%细胞数），保证每个细胞至少能检测到200个基因。
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
#用subset函数，质控：筛选检测到基因数目超过2500或低于200的细胞，单个细胞中线粒体基因数目占比超过>5%
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#数据标准化：默认使用数据标准化方法是LogNormalize, 每个细胞总的表达量都标准化到10000，然后log取对数
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#变化基因鉴定：鉴定在细胞间表达高度变化的基因，后续研究需要集中于这部分基因，首先计算每一个基因的均值和方差，并且直接模拟其关系。默认返回2000个基因。
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#数据缩放 线性转换缩放数据，ScaleData()函数可以实现此功能。最终每个基因均值为0，方差为1。结果存放于pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
#设置参数features是因为ScaleData默认处理前面鉴定的差异基因。这一步怎么做都不会影响到后续pca和聚类，但是会影响做热图
pbmc <- ScaleData(pbmc, features = all.genes)
#移除影响方差的因素
#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
#对缩放后的数据进行PCA分析，默认使用前面鉴定表达变化大的基因。使用features参数可以重新定义数据集。
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p62.rds")

a<-"GSM3516663_MSK_LX661_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p2", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p63.rds")

a<-"GSM3516664_MSK_LX666_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "m1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "m64.rds")

a<-"GSM3516665_MSK_LX675_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p3", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p65.rds")


a<-"GSM3516666_MSK_LX675_NORMAL_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "n1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "n66.rds")


a<-"GSM3516667_MSK_LX676_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p4", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p67.rds")


a<-"GSM3516668_MSK_LX255B_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "m2", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "m68.rds")


a<-"GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p5", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p69.rds")


a<-"GSM3516670_MSK_LX680_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p6", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p70.rds")

a<-"GSM3516671_MSK_LX681_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "m3", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "m71.rds")



a<-"GSM3516672_MSK_LX682_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p7", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p72.rds")

a<-"GSM3516673_MSK_LX682_NORMAL_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "n2", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "n73.rds")


a<-"GSM3516674_MSK_LX684_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p8", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "p74.rds")


a<-"GSM3516675_MSK_LX684_NORMAL_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "n3", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "n75.rds")


a<-"GSM3516676_MSK_LX685_NORMAL_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "n4", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "n76.rds")

a<-"GSM3516677_MSK_LX699_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "m4", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "m77.rds")


a<-"GSM3516678_MSK_LX701_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "m5", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "m78.rds")


#input readRDS
m64<-readRDS(file="m64.rds")
m68<-readRDS(file="m68.rds")
m71<-readRDS(file="m71.rds")
m77<-readRDS(file="m77.rds")
m78<-readRDS(file="m78.rds")

n66<-readRDS(file="n66.rds")
n73<-readRDS(file="n73.rds")
n75<-readRDS(file="n75.rds")
n76<-readRDS(file="n76.rds")

p62<-readRDS(file="p62.rds")
p63<-readRDS(file="p63.rds")
p65<-readRDS(file="p65.rds")
p67<-readRDS(file="p67.rds")
p69<-readRDS(file="p69.rds")
p70<-readRDS(file="p70.rds")
p72<-readRDS(file="p72.rds")
p74<-readRDS(file="p74.rds")

#setup tech and celltype
m64<-RenameCells(m64,add.cell.id="m64",for.merge=T)
m64@meta.data$tech<-"meta"
m64@meta.data$celltype<-"meta_64"

m68<-RenameCells(m68,add.cell.id="m68",for.merge=T)
m68@meta.data$tech<-"meta"
m68@meta.data$celltype<-"meta_68"

m71<-RenameCells(m71,add.cell.id="m71",for.merge=T)
m71@meta.data$tech<-"meta"
m71@meta.data$celltype<-"meta_71"

m77<-RenameCells(m71,add.cell.id="m77",for.merge=T)
m77@meta.data$tech<-"meta"
m77@meta.data$celltype<-"meta_77"

m78<-RenameCells(m71,add.cell.id="m78",for.merge=T)
m78@meta.data$tech<-"meta"
m78@meta.data$celltype<-"meta_78"

n66<-RenameCells(m64,add.cell.id="n66",for.merge=T)
n66@meta.data$tech<-"nontumor"
n66@meta.data$celltype<-"nontumor_66"

n73<-RenameCells(m64,add.cell.id="n73",for.merge=T)
n73@meta.data$tech<-"nontumor"
n73@meta.data$celltype<-"nontumor_73"

n75<-RenameCells(m64,add.cell.id="n75",for.merge=T)
n75@meta.data$tech<-"nontumor"
n75@meta.data$celltype<-"nontumor_75"

n76<-RenameCells(m64,add.cell.id="n76",for.merge=T)
n76@meta.data$tech<-"nontumor"
n76@meta.data$celltype<-"nontumor_76"

p62<-RenameCells(p62,add.cell.id="p62",for.merge=T)
p62@meta.data$tech<-"primary"
p62@meta.data$celltype<-"primary_62"

p63<-RenameCells(p62,add.cell.id="p63",for.merge=T)
p63@meta.data$tech<-"primary"
p63@meta.data$celltype<-"primary_63"

p65<-RenameCells(p62,add.cell.id="p65",for.merge=T)
p65@meta.data$tech<-"primary"
p65@meta.data$celltype<-"primary_65"

p67<-RenameCells(p62,add.cell.id="p67",for.merge=T)
p67@meta.data$tech<-"primary"
p67@meta.data$celltype<-"primary_67"

p69<-RenameCells(p62,add.cell.id="p69",for.merge=T)
p69@meta.data$tech<-"primary"
p69@meta.data$celltype<-"primary_69"

p70<-RenameCells(p62,add.cell.id="p70",for.merge=T)
p70@meta.data$tech<-"primary"
p70@meta.data$celltype<-"primary_70"

p72<-RenameCells(p62,add.cell.id="p72",for.merge=T)
p72@meta.data$tech<-"primary"
p72@meta.data$celltype<-"primary_72"

p74<-RenameCells(p62,add.cell.id="p74",for.merge=T)
p74@meta.data$tech<-"primary"
p74@meta.data$celltype<-"primary_74"


#merge
m64_68<-merge(m64,m68)
m71_77<-merge(m71,m77)
m64_68_71_77<-merge(m64_68,m71_77)
m64_68_71_77_78<-merge(m64_68_71_77,m78)

n66_73<-merge(n66,n73)
n75_76<-merge(n75,n76)
n66_73_75_76<-merge(n66_73,n75_76)

p62_63<-merge(p62,p63)
p65_67<-merge(p65,p67)
p69_70<-merge(p69,p70)
p72_74<-merge(p72,p74)

p62_63_65_67<-merge(p62_63,p65_67)
p69_70_72_74<-merge(p69_70,p72_74)
p62_63_65_67_69_70_72_74<-merge(p62_63_65_67,p69_70_72_74)


m_n<-merge(m64_68_71_77_78,n66_73_75_76)
mnp<-merge(m_n,p62_63_65_67_69_70_72_74)
saveRDS(mnp,file="mnp_17_before_integrate.rds")

mnp<-readRDS(file="mnp_17_before_integrate.rds")

#before integrate
mnp[["percent.mt"]] <- PercentageFeatureSet(mnp, pattern = "^Mt-")
VlnPlot(mnp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = mnp, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
plot_grid(p1,p2)

#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}

reference.list <- pancreas.list[c("meta_64","meta_68","meta_71","meta_77","meta_78","nontumor_66","nontumor_73","nontumor_75","nontumor_76",
                                  "primary_62","primary_63","primary_65","primary_67","primary_69","primary_70","primary_72","primary_74")]

#setup k.anchor and k.filter correctly
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20,k.anchor=5,k.filter=30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "mnp_after_integrated.rds")

hms_individual_integrated<-readRDS(file="mnp_after_integrated.rds")

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
cluster12.markers <- FindMarkers(hms_cluster, ident.1=12, min.pcr=0.25)
head(cluster12.markers, n=10)
cluster13.markers <- FindMarkers(hms_cluster, ident.1=13, min.pcr=0.25)
head(cluster13.markers, n=10)
cluster14.markers <- FindMarkers(hms_cluster, ident.1=14, min.pcr=0.25)
head(cluster14.markers, n=10)
cluster15.markers <- FindMarkers(hms_cluster, ident.1=15, min.pcr=0.25)
head(cluster15.markers, n=10)


new.cluster.ids <- c("Naive CD8+ T", "CD8+ T","Natural Killer T","Naive CD8+ T","Mast Cell",
                     "Natural Killer T","Exhausted CD4/8+ T","Progenitor","Exhausted CD8+ T",
                     "Progenitor","Progenitor","Microglial","Dendritic Cell","B","Progenitor","Microglial") 

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

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Dendritic_cell<-subset(hms_cluster_id, idents=c('Dendritic Cell'))
DimPlot(Dendritic_cell, reduction = "umap")
saveRDS(Dendritic_cell, file="Dendritic_cell.rds")

CD8_T<-subset(hms_cluster_id, idents=c('CD8+ T'))
DimPlot(CD8_T, reduction = "umap")
saveRDS(CD8_T, file="CD8_T.rds")

Mast<-subset(hms_cluster_id, idents=c('Mast Cell'))
DimPlot(Mast, reduction = "umap")
saveRDS(Mast, file="Mast.rds")

Exhausted_CD48<-subset(hms_cluster_id, idents=c('Exhausted CD4/8+ T'))
DimPlot(Exhausted_CD48, reduction = "umap")
saveRDS(Exhausted_CD48, file="Exhausted_CD48.rds")

Exhausted_CD8<-subset(hms_cluster_id, idents=c('Exhausted CD8+ T'))
DimPlot(Exhausted_CD8, reduction = "umap")
saveRDS(Exhausted_CD8, file="Exhausted_CD8.rds")

Naive_CD8<-subset(hms_cluster_id, idents=c('Naive CD8+ T'))
DimPlot(Naive_CD8, reduction = "umap")
saveRDS(Naive_CD8, file="Naive_CD8.rds")


Progenitor<-subset(hms_cluster_id, idents=c('Progenitor'))
DimPlot(Progenitor, reduction = "umap")
saveRDS(Progenitor, file="Progenitor.rds")

Microglial<-subset(hms_cluster_id, idents=c('Microglial'))
DimPlot(Microglial, reduction = "umap")
saveRDS(Microglial, file="Microglial.rds")



#input each cluster
Natural_Killer_T<-readRDS("Natural_Killer_T.rds")
B_cell<-readRDS("B.rds")
Dendritic_cell<-readRDS("Dendritic_cell.rds")
CD8_T<-readRDS("CD8_T.rds")
Progenitor<-readRDS("Progenitor.rds")
Mast<-readRDS("Mast.rds")
Microglial<-readRDS("Microglial.rds")
Exhausted_CD8<-readRDS("Exhausted_CD8.rds")
Exhausted_CD48<-readRDS("Exhausted_CD48.rds")
Naive_CD8<-readRDS("Naive_CD8.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

#DimPlot
DimPlot(B_cell, reduction = "umap", split.by = "tech")
DimPlot(Natural_Killer_T, reduction = "umap", split.by = "tech")
DimPlot(Dendritic_cell, reduction = "umap", split.by = "tech")
DimPlot(CD8_T, reduction = "umap", split.by = "tech")
DimPlot(Progenitor, reduction = "umap", split.by = "tech")
DimPlot(Mast, reduction = "umap", split.by = "tech")
DimPlot(Exhausted_CD8, reduction = "umap", split.by = "tech")
DimPlot(Microglial, reduction = "umap", split.by = "tech")
DimPlot(Exhausted_CD48, reduction = "umap", split.by = "tech")
DimPlot(Naive_CD8, reduction = "umap", split.by = "tech")

RidgePlot(CD8_T, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"),cols = c("green3","orangered","red"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
DoHeatmap(subset(Exhausted_CD8,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
DoHeatmap(subset(Exhausted_CD48,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
DoHeatmap(subset(Naive_CD8,downsample=50000),features = markers.to.plot,size=5,label = FALSE)


genes11_heatmap<-DotPlot(hms_cluster_id,features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()
genes11_heatmap<-genes11_heatmap$data
write.csv(genes11_heatmap,"genes11_heatmap")

#expression level in each patients
setwd("~/geo/gse123902")
CD8_T<-readRDS("CD8_T.rds")
a<-DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_cd8T")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")
a<-DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_hms")

ex_CD8_T<-readRDS("Exhausted_CD8.rds")
a<-DoHeatmap(subset(ex_CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_ex_cd8T")

ex_CD48_T<-readRDS("Exhausted_CD48.rds")
a<-DoHeatmap(subset(ex_CD48_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_ex_cd48T")

Naive_CD8<-readRDS("Naive_CD8.rds")
a<-DoHeatmap(subset(Naive_CD8,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_Naive_CD8")