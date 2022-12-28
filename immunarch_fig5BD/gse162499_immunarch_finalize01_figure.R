#install.packages(c("devtools","pkgload"))
#devtools::install_github("immunomind/immunarch",ref="dev")#sometime not work, because of internet.
#work on gnv at home computer,mm:1.
library(immunarch)
setwd("~/gse162499")
file_path="~/gse162499/immunarch_data/immunarch"
immdata<-repLoad(file_path)
repExplore(immdata$data, "lens") %>% vis() #input data
#ensembl_list<-read.table("tmp01",header=FALSE)
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



#finalize figure
target <- c("CVVSADNDYKLSF;CASSRPGAQENGAFF", "CALSGVTSYDKVIF;CASSLLGGGNNEQFF", "CALSEAGAAGNKLTF;CASSPPLGDTEAFF", 
            "CCALEGPNAGKSTF;CAISEFGNEAFF", "CAVRGETSGSRLTF;CASSLGPSYEQYF","CIVRVETGANNLFF;CASSWTDLPNSPLHF",
            "CASSLLGGGNNEQFF","CALSEVNQAGTALIF;CASSDAGVSTNEKLFF","CAVAAGNKLTF;CASSNLDRTGQETQYF",
            "CAMRLRRSYKLIF;CASSPWGNTGELFF",
            "CIVRVTGNQFYF;CASSFLGGNQPQHF","CAVLRDDKIIF;CASSSRREPSYNEQFF",
            "CATVLGTGKLIF;CASSELNQPQHF","CAASRNAGNMLTF;CASSISGTGEIGEAFF","CAVRLTTPPYGQNFVF;CASCPGQGLSPLHF",
            "CALLGGINTGNQFYF;CASISWDRETGHPLHF","CAEGALGSYIPTF;CAFYSSASKIIF;CSGGTSGYEQYF",
            "CAMREGNTGGFKTIF;CASSQEDRVEETQYF","CASSWTDLPNSPLHF","CASSFLGGNQPQHF",
            "CSAREPGGYGYTF","CAYRSAWGPGNQFYF;CASSPLVGSSYNEQFF;CASSKGTGLAKNIQYF","CAVGESGGSNYKLTF;CASRQRDSYEQYF",
            "CAMREPYSGAGSYQLTF;CASSHEPGGGPNEQFF","CLTPGGFKTIF;CASSLGGLYEQYF")
p3 <- trackClonotypes(immdata$data, target, .col = "aa") %>% vis()
p3

repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta) 

imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)

vis(imm_hom)
vis(imm_hom, .by = "Status", .meta = immdata$meta)

repClonality(immdata$data, "homeo") %>% vis() 
repOverlap(immdata$data) %>% vis()  
geneUsage(immdata$data[[1]]) %>% vis()
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta) 

kmers  <-  getKmers ( immdata $ data [[ 1 ]], 5 )
kp  <-  kmer_profile ( kmers , "self" )
p1  <-  vis ( kp )
p2  <-  vis ( kp , .plot  =  "seq" )

p1 + p2