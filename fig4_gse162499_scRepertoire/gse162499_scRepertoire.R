#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("scRepertoire",force = TRUE)


#library(devtools)
#devtools::install_github("ncborcherding/scRepertoire@dev")

library(Seurat)
library(scRepertoire)
setwd("~/geo/gse162499")
P57_T <- read.csv("GSM4952972_P57_Tumor_filtered_contig_annotations.csv")
P57_B <- read.csv("GSM4952973_P57_Blood_filtered_contig_annotations.csv")
P58_T <- read.csv("GSM4952974_P58_Tumor_filtered_contig_annotations.csv")
P58_B <- read.csv("GSM4952975_P58_Blood_filtered_contig_annotations.csv")
P60_T <- read.csv("GSM4952976_P60_Tumor_filtered_contig_annotations.csv")
P60_B <- read.csv("GSM4952978_P60_Blood_filtered_contig_annotations.csv")
P61_T <- read.csv("GSM4952979_P61_Tumor_filtered_contig_annotations.csv")
P61_B <- read.csv("GSM4952981_P61_Blood_filtered_contig_annotations.csv")

contig_list <- list(P57_T,P57_B,P58_T,P58_B,P60_T,P60_B,P61_T,P61_B)

#data("contig_list") #the data built into scRepertoire
#head(contig_list[[1]])

combined <- combineTCR(contig_list, 
                       samples = c("P57", "P57", "P58", "P58", "P60","P60","P61","P61"), 
                       ID = c("T", "B", "T", "B", "T", "B","T","B"), cells ="T-AB")


quantContig(combined, cloneCall="gene+nt", scale = T)
quantContig(combined, cloneCall="gene+nt", chain = "TRA")
quantContig(combined, cloneCall="gene+nt", chain = "TRB")
quantContig(combined, cloneCall="gene+nt", chain = "TRA", scale = TRUE)
quantContig(combined, cloneCall="gene+nt", scale = T, chain = "both")
quantContig_output <- quantContig(combined, cloneCall="gene+nt", 
                                  scale = T, exportTable = T)
quantContig_output
abundanceContig(combined, cloneCall = "gene", scale = F)
abundanceContig(combined, cloneCall = "gene", exportTable = T)
lengthContig(combined, cloneCall="aa", chain = "both") 
lengthContig(combined, cloneCall="nt", chain = "TRA") 
compareClonotypes(combined, numbers = 5, samples = c("P57_B", "P57_T"), 
                     cloneCall="aa", graph = "alluvial")

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

seurat<-readRDS(file="P57586061_after_integrated.rds")
     