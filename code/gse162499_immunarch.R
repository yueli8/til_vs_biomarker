library(immunarch)
setwd("~/geo/gse162499")
file_path="~/geo/gse162499/immunarch"
immdata<-repLoad(file_path)
repExplore(immdata$data, "lens") %>% vis() #input data
ensembl_list<-read.table("tmp01",header=FALSE)
a<-immdata$meta#meta data
write.csv(a,"a")
b<-top(immdata$data[[1]])
write.csv(b,"b")
#each file
P57_tumor<-immdata$data$P57_Tumor
write.csv(P57_tumor,"P57_tumor")
P57_blood<-immdata$data$P57_Blood
write.csv(P57_blood,"P57_blood")
P58_tumor<-immdata$data$P58_Tumor
write.csv(P58_tumor,"P58_tumor")
P58_blood<-immdata$data$P58_Blood
write.csv(P58_blood,"P58_blood")
P60_tumor<-immdata$data$P60_Tumor
write.csv(P60_tumor,"P60_tumor")
P60_blood<-immdata$data$P60_Blood
write.csv(P60_blood,"P60_blood")
P61_tumor<-immdata$data$P61_Tumor
write.csv(P61_tumor,"P61_tumor")
P61_blood<-immdata$data$P61_Blood
write.csv(P61_blood,"P61_blood")

repExplore(immdata$data, "lens") %>% vis() 
p1 <- repOverlap(immdata$data) %>% vis()
p2 <- repDiversity(immdata$data) %>% vis()

target <- c("CVVSADNDYKLSF;CASSRPGAQENGAFF", "CALSGVTSYDKVIF;CASSLLGGGNNEQFF", "CALSEAGAAGNKLTF;CASSPPLGDTEAFF", 
            "CCALEGPNAGKSTF;CAISEFGNEAFF", "CAVRGETSGSRLTF;CASSLGPSYEQYF","CIVRVETGANNLFF;CASSWTDLPNSPLHF",
            "CASSLLGGGNNEQFF","CALSEVNQAGTALIF;CASSDAGVSTNEKLFF","CAVAAGNKLTF;CASSNLDRTGQETQYF",
            "CAMRLRRSYKLIF;CASSPWGNTGELFF",
            "CIVRVTGNQFYF;CASSFLGGNQPQHF","CAVLRDDKIIF;CASSSRREPSYNEQFF",
            "CATVLGTGKLIF;CASSELNQPQHF","CAASRNAGNMLTF;CASSISGTGEIGEAFF","CAVRLTTPPYGQNFVF;CASCPGQGLSPLHF",
  "CALLGGINTGNQFYF;CASISWDRETGHPLHF","CAEGALGSYIPTF;CAFYSSASKIIF;CSGGTSGYEQYF",
            "CAMREGNTGGFKTIF;CASSQEDRVEETQYF","CASSWTDLPNSPLHF","CASSFLGGNQPQHF",
 "CSAREPGGYGYTF","CAYRSAWGPGNQFYF;CASSPLVGSSYNEQFF;CASSKGTGLAKNIQYF","CAVGESGGSNYKLTF;CASRQRDSYEQYF",
             "CAMREPYSGAGSYQLTF;CASSHEPGGGPNEQFF","CLTPGGFKTIF;CASSLGGLYEQYF",
              "CAVSGEHTDKLIF;CASSLWATFNYGYTF",
            "CASSSRREPSYNEQFF","CAFYSSASKIIF;CAEGALGSYIPTF;CSGGTSGYEQYF","CIVDTFHLKYNTDKLIF;CAIILWIGLNTEAFF",
            "CASSLGGLYEQYF")
p3 <- trackClonotypes(immdata$data, target, .col = "aa") %>% vis()
p1
p2
p3

target <- c("CVVSADNDYKLSF;CASSRPGAQENGAFF", "CALSGVTSYDKVIF;CASSLLGGGNNEQFF", "CALSEAGAAGNKLTF;CASSPPLGDTEAFF", 
            "CCALEGPNAGKSTF;CAISEFGNEAFF", "CAVRGETSGSRLTF;CASSLGPSYEQYF")
p4 <- trackClonotypes(immdata$data, target, .col = "aa") %>% vis()
p4

target <- c("CVVSADNDYKLSF;CASSRPGAQENGAFF", "CALSGVTSYDKVIF;CASSLLGGGNNEQFF", "CALSEAGAAGNKLTF;CASSPPLGDTEAFF", 
            "CCALEGPNAGKSTF;CAISEFGNEAFF", "CAVRGETSGSRLTF;CASSLGPSYEQYF","CIVRVETGANNLFF;CASSWTDLPNSPLHF",
            "CASSLLGGGNNEQFF","CALSEVNQAGTALIF;CASSDAGVSTNEKLFF","CAVAAGNKLTF;CASSNLDRTGQETQYF",
            "CAMRLRRSYKLIF;CASSPWGNTGELFF")
p5 <- trackClonotypes(immdata$data, target, .col = "aa") %>% vis()
p5

target <- c("CVVSADNDYKLSF;CASSRPGAQENGAFF", "CALSGVTSYDKVIF;CASSLLGGGNNEQFF", "CALSEAGAAGNKLTF;CASSPPLGDTEAFF", 
            "CCALEGPNAGKSTF;CAISEFGNEAFF", "CAVRGETSGSRLTF;CASSLGPSYEQYF","CIVRVETGANNLFF;CASSWTDLPNSPLHF",
            "CASSLLGGGNNEQFF","CALSEVNQAGTALIF;CASSDAGVSTNEKLFF","CAVAAGNKLTF;CASSNLDRTGQETQYF",
            "CAMRLRRSYKLIF;CASSPWGNTGELFF","CIVRVTGNQFYF;CASSFLGGNQPQHF","CAVLRDDKIIF;CASSSRREPSYNEQFF",
            "CATVLGTGKLIF;CASSELNQPQHF","CAASRNAGNMLTF;CASSISGTGEIGEAFF","CAVRLTTPPYGQNFVF;CASCPGQGLSPLHF")
p3 <- trackClonotypes(immdata$data, target, .col = "aa") %>% vis()
p3

target <- c("CVVSADNDYKLSF;CASSRPGAQENGAFF", "CALSGVTSYDKVIF;CASSLLGGGNNEQFF", "CALSEAGAAGNKLTF;CASSPPLGDTEAFF", 
            "CCALEGPNAGKSTF;CAISEFGNEAFF", "CAVRGETSGSRLTF;CASSLGPSYEQYF","CIVRVETGANNLFF;CASSWTDLPNSPLHF",
            "CASSLLGGGNNEQFF","CALSEVNQAGTALIF;CASSDAGVSTNEKLFF","CAVAAGNKLTF;CASSNLDRTGQETQYF",
            "CAMRLRRSYKLIF;CASSPWGNTGELFF",
            "CIVRVTGNQFYF;CASSFLGGNQPQHF","CAVLRDDKIIF;CASSSRREPSYNEQFF",
            "CATVLGTGKLIF;CASSELNQPQHF","CAASRNAGNMLTF;CASSISGTGEIGEAFF","CAVRLTTPPYGQNFVF;CASCPGQGLSPLHF",
            "CALLGGINTGNQFYF;CASISWDRETGHPLHF","CAEGALGSYIPTF;CAFYSSASKIIF;CSGGTSGYEQYF",
            "CAMREGNTGGFKTIF;CASSQEDRVEETQYF","CASSWTDLPNSPLHF","CASSFLGGNQPQHF")
p3 <- trackClonotypes(immdata$data, target, .col = "aa") %>% vis()
p3

repExplore(immdata$data, "lens") %>% vis() #length

exp_vol <- repExplore(immdata$data, .method = "volume")
p1 <- vis(exp_vol, .by = c("Status"), .meta = immdata$meta)
p1
p2 <- vis(exp_vol, .by = c("Sex"), .meta = immdata$meta)
p2
#smoke not work
#p3 <- vis(exp_vol, .by = c("Smoke"), .meta = immdata$meta)
#p3

p3 <- vis(exp_vol, .by = c("Smoke"), .meta = immdata$meta)
p3
p4 <- vis(exp_vol, .by = c("Stage"), .meta = immdata$meta)
p4
p5 <- vis(exp_vol, .by = c("Status", "Sex"), .meta = immdata$meta)
p5
p6 <- vis(exp_vol, .by = c("Status", "Smoke"), .meta = immdata$meta)
p6
p7 <- vis(exp_vol, .by = c("Status", "Stage"), .meta = immdata$meta)
p7
#the same as the above P1
exp_vol <- repExplore(immdata$data, .method = "volume")
by_vec <- c( "Tumor", "Tumor", "Tumor", "Tumor","Blood", "Blood", "Blood", "Blood")
p1 <- vis(exp_vol, .by = by_vec)
p1
#length
exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
p1 <- vis(exp_len)
fixVis(p1)

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
head(exp_len)
p1 <- vis(exp_len)
p1

exp_cnt <- repExplore(immdata$data, .method = "count")
head(exp_cnt)
p2 <- vis(exp_cnt)
p2

exp_vol <- repExplore(immdata$data, .method = "volume")
head(exp_vol)
p3 <- vis(exp_vol)
p3

p4 <- vis(exp_len, .by = "Status", .meta = immdata$meta)
p5 <- vis(exp_cnt, .by = "Sex", .meta = immdata$meta)
p6 <- vis(exp_vol, .by = c("Status", "Sex"), .meta = immdata$meta)
p4
p5
p6

imm_pr  <-  repClonality ( immdata$data , .method  =  "clonal.prop" )
imm_pr
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top
imm_rare <- repClonality(immdata$data, .method = "rare")
imm_rare
imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
imm_hom
vis(imm_top) 
imm_pr  <-  repClonality ( immdata$data , .method  =  "clonal.prop" )
imm_pr
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top
imm_rare <- repClonality(immdata$data, .method = "rare")
imm_rare
imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)

vis(imm_top)
vis(imm_top, .by = "Status", .meta = immdata$meta)
vis(imm_rare) 
vis(imm_rare, .by = "Status", .meta = immdata$meta)
vis(imm_hom)
vis(imm_hom, .by = "Status", .meta = immdata$meta)

imm_ov1 <- repOverlap(immdata$data, .method = "public", .verbose = F)
head(imm_ov1)

imm_ov2 <- repOverlap(immdata$data, .method = "morisita", .verbose = F)
head(imm_ov2)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 4)
p1
p2

vis(imm_ov1, "heatmap2")
p1  <-  vis ( imm_ov2 , .text.size  =  4 , .signif.digits  =  1 )
p1

repOverlapAnalysis(imm_ov1, "tsne") %>% vis()
#repOverlapAnalysis(imm_ov1, "mds+kmeans") %>% vis()  #error
#构建公共克隆型库
pr.nt <- pubRep(immdata$data, "nt", .verbose = F)
pr.nt
pr.aav <- pubRep(immdata$data, "aa+v", .verbose = F)
pr.aav
pr.aav.cod <- pubRep(immdata$data, "aa+v", .coding = T)
pr <- pubRep(immdata$data, "aa+v", .coding = T, .verbose = F)
head(pr)
pr1 <- pubRepFilter(pr, immdata$meta, .by = c(Status = "Tumor"))
pr2 <- pubRepFilter(pr, immdata$meta, .by = c(Status = "Blood"))
pr3 <- pubRepApply(pr1, pr2)
p <- ggplot() +
  geom_jitter(aes(x = "Treatment", y = Result), data = pr3)
p

#基因使用计算: no number,so no value
gene_stats()
imm_gu <- geneUsage(immdata$data, "hs.trbv")
imm_gu
imm_gu <- geneUsage(immdata$data[c(1, 2)], "hs.trbv", .norm = T)
vis(imm_gu)
imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T)
vis(imm_gu, .by = "Status", .meta = immdata$meta)
vis(imm_gu, .grid = T)

#基因使用分析
imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T)
head(imm_gu)
imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = F)
head(imm_gu_js)
imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)
head(imm_gu_cor)
p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 3)
p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 3)

p1
p2
imm_gu_js [ is.na ( imm_gu_js )] <-  0
vis ( geneUsageAnalysis ( imm_gu , "cosine+hclust" , .verbose  =  F ))
vis(geneUsageAnalysis(imm_gu, "js+dbscan", .verbose = F))
imm_cl_pca <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = F)
imm_cl_mds <- geneUsageAnalysis(imm_gu, "js+mds+kmeans", .verbose = F)
imm_cl_tsne <- geneUsageAnalysis(imm_gu, "js+tsne+kmeans", .perp = .01, .verbose = F)
## Perplexity should be lower than K!
p1 <- vis(imm_cl_pca, .plot = "clust")
p2 <- vis(imm_cl_mds, .plot = "clust")
p3 <- vis(imm_cl_tsne, .plot = "clust")
p1
p2
p3

imm_cl_pca2 <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .k = 3, .verbose = F)
vis(imm_cl_pca2)

p1 <- vis(spectratype(immdata$data[[1]], .quant = "id", .col = "nt"))
p2 <- vis(spectratype(immdata$data[[1]], .quant = "count", .col = "aa+v"))
p1+p2
p2

data(immdata)

# Compute statistics and visualise them
# Chao1 diversity measure
div_chao <- repDiversity(immdata$data, "chao1")
head(div_chao)
div_hill <- repDiversity(immdata$data, "hill")
head(div_hill)
div_d50 <- repDiversity(immdata$data, "d50")
head(div_d50)
div_div <- repDiversity(immdata$data, "div")
head(div_div)
p1 <- vis(div_chao)
p2 <- vis(div_chao, .by = c("Status", "Sex"), .meta = immdata$meta)
p3 <- vis(div_hill, .by = c("Status", "Sex"), .meta = immdata$meta)

p4 <- vis(div_d50)
p5 <- vis(div_d50, .by = "Status", .meta = immdata$meta)
p6 <- vis(div_div)
p1 + p2
p3 + p6
p4 + p5

imm_raref <- repDiversity(immdata$data, "raref", .verbose = F)
head(imm_raref)
p1 <- vis(imm_raref)
p2 <- vis(imm_raref, .by = "Status", .meta = immdata$meta)
p1 + p2
repDiversity(immdata$data, "raref", .verbose = F) %>% vis(.log = TRUE)

tc1 <- trackClonotypes(immdata$data, list(1, 10), .col = "aa+v")
head(tc1)
p1  <-  vis ( tc1 )
p1

file_path="~/geo/gse162499/immunarch_blood"
immdata<-repLoad(file_path)
tc2  <-  trackClonotypes ( immdata$data , list (1,10), .col  =  "nt" )
head(tc2)
p2  <-  vis ( tc2 )
p2
tc2  <-  trackClonotypes ( immdata$data , list (1,10), .col  =  "aa+v" )
head(tc2)
p2  <-  vis ( tc2 )
p2

file_path="~/geo/gse162499/immunarch_tumor"
immdata<-repLoad(file_path)
tc2  <-  trackClonotypes ( immdata$data , list (1,10), .col  =  "nt" )
head(tc2)
p2  <-  vis ( tc2 )
p2
tc2  <-  trackClonotypes ( immdata$data , list (1,10), .col  =  "aa+v" )
head(tc2)
p2  <-  vis ( tc2 )
p2

file_path="~/geo/gse162499/immunarch"
immdata<-repLoad(file_path)
target <- immdata$data[[1]] %>%
  select(CDR3.aa, V.name) %>%
  head(10)
tc <- trackClonotypes(immdata$data, target)
vis(tc)

tc <- trackClonotypes(immdata$data, target, .col = "aa")
vis(tc, .plot = "smooth")
vis(tc, .plot = "area")
vis(tc, .plot = "line")

# Passing indices
names(immdata$data)[c(1, 3, 5,7)] # check sample names
## [1] "A2-i129" "A2-i133" "A4-i191"
vis(tc, .order = c(1, 3, 5,7))


immdata$meta$Timepoint <- sample(1:length(immdata$data))
immdata$meta

vdjdb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb")
vdjdb
#我们可以通过设置.species，.chain以及.pathology等参数进行过滤筛选出需要的信息
vdjdb  =  dbLoad ( "https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz" , "vdjdb" , .species  =  "HomoSapiens" , .chain  =  "TRB" , .pathology  =  "CMV" )
vdjdb
vdjdb_st = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/SearchTable-2019-10-17%2012_36_11.989.tsv.gz", "vdjdb-search", .species = "HomoSapiens", .chain = "TRB", .pathology = "CMV")
vdjdb_st
mcpas = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human", .chain = "TRB", .pathology = "Cytomegalovirus (CMV)")
mcpas
tbadb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/TBAdb.xlsx", "tbadb", .species = "Homo Sapiens", .chain = c("TRB", "TRA-TRB"), .pathology = "CMV")
tbadb

dbAnnotate(immdata$data, vdjdb, "CDR3.aa", "cdr3")
dbAnnotate ( immdata $ data , mcpas , c ( "CDR3.aa" , "V.name" ), c ( "CDR3.beta.aa" , "TRBV" ))










repClonality(immdata$data, "homeo") %>% vis() 
repOverlap(immdata$data) %>% vis()  
geneUsage(immdata$data[[1]]) %>% vis()
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta) 