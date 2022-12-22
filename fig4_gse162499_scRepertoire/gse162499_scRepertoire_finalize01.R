library(Seurat)
library(scRepertoire)
library(devtools)
#devtools::install_github("ncborcherding/scRepertoire@dev")
#manual: https://ncborcherding.github.io/vignettes/vignette.html
#not install from bioconductor,it is from github. 
setwd("~/geo/gse162499/scRepertoire_finalize")
P57_T <- read.csv("GSM4952972_P57_Tumor_filtered_contig_annotations.csv")
P57_B <- read.csv("GSM4952973_P57_Blood_filtered_contig_annotations.csv")
P58_T <- read.csv("GSM4952974_P58_Tumor_filtered_contig_annotations.csv")
P58_B <- read.csv("GSM4952975_P58_Blood_filtered_contig_annotations.csv")
P60_T <- read.csv("GSM4952976_P60_Tumor_filtered_contig_annotations.csv")
P60_B <- read.csv("GSM4952978_P60_Blood_filtered_contig_annotations.csv")
P61_T <- read.csv("GSM4952979_P61_Tumor_filtered_contig_annotations.csv")
P61_B <- read.csv("GSM4952981_P61_Blood_filtered_contig_annotations.csv")
contig_list <- list(P57_T,P57_B,P58_T,P58_B,P60_T,P60_B,P61_T,P61_B)
combined <- combineTCR(contig_list, 
                       samples = c("P57", "P57", "P58", "P58", "P60","P60","P61","P61"), 
                       ID = c("T", "B", "T", "B", "T","B","T","B"), cells ="T-AB")
quantContig(combined, cloneCall="gene+nt", scale = T)
quantContig(combined, cloneCall="gene+nt", chain = "TRA")
quantContig(combined, cloneCall="gene+nt", chain = "TRB")
quantContig(combined, cloneCall="gene+nt", chain = "TRA", scale = TRUE)
quantContig(combined, cloneCall="gene+nt", scale = T, chain = "both")
quantContig_output <- quantContig(combined, cloneCall="gene+nt", 
                                  scale = T, exportTable = T)
quantContig_output

abundanceContig(combined, cloneCall = "gene", scale = F)
#abundanceContig(combined, cloneCall = "gene", exportTable = T)
lengthContig(combined, cloneCall="aa", chain = "both") 
lengthContig(combined, cloneCall="nt", chain = "TRA") 

compareClonotypes(combined, numbers = 5, samples = c("P57_B", "P57_T"), 
                  cloneCall="aa", graph = "alluvial")

vizGenes(combined, gene = "V", 
         chain = "TRA", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)
vizGenes(combined, gene = "V", 
         chain = "TRB", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)
vizGenes(combined[c(1,3,5,7)], 
         gene = "V", 
         chain = "TRB", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")
vizGenes(combined[c(2,4,6,8)], 
         gene = "V", 
         chain = "TRB", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")
clonalHomeostasis(combined, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))


clonalHomeostasis(combined, cloneCall = "aa")
clonalProportion(combined, cloneCall = "gene",
                 split = c(10, 100, 1000, 10000, 30000, 1e+05)) 
clonalProportion(combined, cloneCall = "nt") 
clonalOverlap(combined, cloneCall = "gene+nt", 
              method = "morisita")
clonesizeDistribution(combined, cloneCall = "gene+nt", 
                      method="ward.D2")
clonalDiversity(combined, cloneCall = "gene", group.by = "sample", 
                x.axis = "ID", n.boots = 100)

scatterClonotype(combined, cloneCall ="gene", 
                 x.axis = "P57_T", 
                 y.axis = "P57_B",
                 dot.size = "total",
                 graph = "proportion")

#tumor and blood samples
Tumor <- read.csv("Tumor_total.csv")
Blood <- read.csv("Blood_total.csv")
contig_list <- list(Tumor,Blood)
combined <- combineTCR(contig_list, 
                       samples = c("Tumor", "Blood"), 
                       ID = c("T", "B"), cells ="T-AB")
quantContig(combined, cloneCall="gene+nt", scale = T)
quantContig(combined, cloneCall="gene+nt", chain = "TRA")
quantContig(combined, cloneCall="gene+nt", chain = "TRB")
quantContig(combined, cloneCall="gene+nt", chain = "TRA", scale = TRUE)
quantContig(combined, cloneCall="gene+nt", scale = T, chain = "both")
quantContig_output <- quantContig(combined, cloneCall="gene+nt", 
                                  scale = T, exportTable = T)
quantContig_output

abundanceContig(combined, cloneCall = "gene", scale = F)
#abundanceContig(combined, cloneCall = "gene", exportTable = T)
lengthContig(combined, cloneCall="aa", chain = "both") 
lengthContig(combined, cloneCall="nt", chain = "TRA") 
compareClonotypes(combined, numbers = 5, samples = c("P57_B", "P57_T"), 
                  cloneCall="aa", graph = "alluvial")

vizGenes(combined, gene = "V", 
         chain = "TRB", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)

clonalHomeostasis(combined, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))
clonalHomeostasis(combined, cloneCall = "aa")
clonalProportion(combined, cloneCall = "gene",
                 split = c(10, 100, 1000, 10000, 30000, 1e+05)) 
clonalProportion(combined, cloneCall = "nt") 
clonalOverlap(combined, cloneCall = "gene+nt", 
              method = "morisita")
clonesizeDistribution(combined, cloneCall = "gene+nt", 
                      method="ward.D2")
clonalDiversity(combined, cloneCall = "gene", group.by = "sample", 
                x.axis = "ID", n.boots = 100)
#add seurat
seurat <- readRDS(file="hms_cluster_id_test_correct.rds")
DimPlot(seurat, label = T) + NoLegend()
table(Idents(seurat))

setwd("~/geo/gse162499/scRepertoire_finalize")
P57_T <- read.csv("GSM4952972_P57_Tumor_filtered_contig_annotations.csv")
P57_B <- read.csv("GSM4952973_P57_Blood_filtered_contig_annotations.csv")
P58_T <- read.csv("GSM4952974_P58_Tumor_filtered_contig_annotations.csv")
P58_B <- read.csv("GSM4952975_P58_Blood_filtered_contig_annotations.csv")
P60_T <- read.csv("GSM4952976_P60_Tumor_filtered_contig_annotations.csv")
P60_B <- read.csv("GSM4952978_P60_Blood_filtered_contig_annotations.csv")
P61_T <- read.csv("GSM4952979_P61_Tumor_filtered_contig_annotations.csv")
P61_B <- read.csv("GSM4952981_P61_Blood_filtered_contig_annotations.csv")
contig_list <- list(P57_T,P57_B,P58_T,P58_B,P60_T,P60_B,P61_T,P61_B)
combined <- combineTCR(contig_list, 
                       samples = c("P57", "P57", "P58", "P58", "P60","P60","P61","P61"), 
                       ID = c("T", "B", "T", "B", "T","B","T","B"), cells ="T-AB")


seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", group.by = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

names(seurat@meta.data)
head(seurat@meta.data)

DimPlot(seurat, group.by = "tech") +
  scale_color_manual(values=colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

DimPlot(seurat, group.by = "celltype") +
  scale_color_manual(values=colorblind_vector(8), na.value="grey") + 
  theme(plot.title = element_blank())

slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
DimPlot(seurat, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())


clonalOverlay(seurat, reduction = "umap", 
              freq.cutpoint = 30, bins = 10, facet = "tech") + 
  guides(color = "none")


clonalOverlay(seurat, reduction = "umap", 
              freq.cutpoint = 30, bins = 10, facet = "celltype") + 
  guides(color = "none")

library(ggraph)

saveRDS(seurat, file = "seurat.rds")


seurat <- highlightClonotypes(seurat,  cloneCall= "aa", sequence = c("CALSGVTSYDKVIF_CASSLLGGGNNEQFF","CAASRNAGNMLTF_CASSISGTGEIGEAFF",
                                                                     "CAEGALGSYIPTF;CAFYSSASKIIF_CSGGTSGYEQYF","CALLGGINTGNQFYF_CASISWDRETGHPLHF","CALSEAGAAGNKLTF_CASSPPLGDTEAFF",
                                                                     "CALSEVNQAGTALIF_CASSDAGVSTNEKLFF","CAMREGNTGGFKTIF_CASSQEDRVEETQYF"))

DimPlot(seurat, group.by = "highlight", pt.size = 0.5)
theme(plot.title = element_blank())


occupiedscRepertoire(seurat, x.axis = "ident")

library(circlize)
library(scales)

circles <- getCirclize(seurat, 
                       group.by = "ident")

#Just assigning the normal colors to each cluster
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)

#Graphing the chord diagram
circlize::chordDiagram(circles,
                       self.link = 1, 
                       grid.col = grid.cols)


#Just assigning the normal colors to each cluster
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)

#Graphing the chord diagram
circlize::chordDiagram(circles, self.link = 1, grid.col = grid.cols)

sub_combined <- clusterTCR(combined[[2]], 
                           chain = "TRA", 
                           sequence = "aa", 
                           threshold = 0.85, 
                           group.by = NULL)

seurat <- clusterTCR(seurat,
                     chain = "TRB",
                     group.by = "celltype", 
                     sequence = "aa", 
                     threshold = 0.85)

DimPlot(seurat, group.by = "TRB_cluster") + 
  scale_color_manual(values = colorblind_vector(length(unique(seurat@meta.data[,"TRB_cluster"])))) + 
  NoLegend()



StartracDiversity(seurat, type = "celltype", 
                  sample = "celltype", by = "overall")

library(Seurat)
library(scRepertoire)
library(devtools)
setwd("~/geo/gse162499/scRepertoire_finalize")
P57_T <- read.csv("GSM4952972_P57_Tumor_filtered_contig_annotations.csv")
P57_B <- read.csv("GSM4952973_P57_Blood_filtered_contig_annotations.csv")
P58_T <- read.csv("GSM4952974_P58_Tumor_filtered_contig_annotations.csv")
P58_B <- read.csv("GSM4952975_P58_Blood_filtered_contig_annotations.csv")
P60_T <- read.csv("GSM4952976_P60_Tumor_filtered_contig_annotations.csv")
P60_B <- read.csv("GSM4952978_P60_Blood_filtered_contig_annotations.csv")
P61_T <- read.csv("GSM4952979_P61_Tumor_filtered_contig_annotations.csv")
P61_B <- read.csv("GSM4952981_P61_Blood_filtered_contig_annotations.csv")
contig_list <- list(P57_T,P57_B,P58_T,P58_B,P60_T,P60_B,P61_T,P61_B)
combined <- combineTCR(contig_list, 
                       samples = c("P57", "P57", "P58", "P58", "P60","P60","P61","P61"), 
                       ID = c("T", "B", "T", "B", "T","B","T","B"), cells ="T-AB")
seurat <- readRDS(file="hms_cluster_id_test_correct.rds")
DimPlot(seurat, label = T) + NoLegend()
table(Idents(seurat))
seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", group.by = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
DimPlot(seurat, group.by = "tech") +
  scale_color_manual(values=colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())
slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
DimPlot(seurat, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())
clonalOverlay(seurat, reduction = "umap", 
              freq.cutpoint = 30, bins = 10, facet = "tech") + 
  guides(color = "none")


