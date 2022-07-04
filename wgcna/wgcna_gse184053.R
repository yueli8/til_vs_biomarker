#install.packages("BiocManager")
#BiocManager::install("WGCNA")

setwd("~/gse184053/wgcna_gse184053")
##导入数据##
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

samples=read.csv('Sam_info.txt',sep = '\t',row.names = 1)
expro=read.csv('9samples.txt',sep = '\t',row.names = 1)
dim(expro)

##筛选方差前25%的基因
##这一步是为了减少运算量，因为一个测序数据可能会有好几万个探针，而可能其中很多基因在各个样本中的表达情况并没有什么太大变化，
#为了减少运算量，这里我们筛选方差前25%的基因。

#当然了，以下几种情况，你可以忽略此步，转而执行下列代码
#情况1：所用的数据是比较老的芯片数据，探针数量较少；
#情况2：你的电脑足够强大，不必减少运算；
#情况3：你第一步导入的数据是差异基因的数据，已经作过初筛，比如这一文章的套路（PMID：27133569）（个人不建议这么做）。

#下面两段程序二选一， 第一段筛选方差前25%的基因，第二段全部保留
m.vars=apply(expro,1,var)
expro.upper=expro[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
dim(expro.upper)
datExpr=as.data.frame(t(expro.upper));
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#datExpr= as.data.frame(t(expro));
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)

##样本聚类检查离群值##选择一个β值建立临近矩阵根据连接度使我们的基因分布符合无尺度网络
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
#save(datExpr, file = "hasm_yueli.RData")

##软阈值筛选##软阈值的筛选原则是使构建的网络更符合无标度网络特征。
#目的是为了帮助用户选择一个合适的软阈值进行网络构建。

# number must be very large, then can find the threshold
#设置软阈值调参范围，powers是数组，包括1，2，...10,12,14，...,20
powers = c(c(1:10), seq(from = 5, to=100, by=5))
# 网络拓扑分析
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 绘图
sizeGrWindow(9, 5)	# 图片的宽度和高度
# 1行2列排列
par(mfrow = c(1,2));# 一页多图，一页被分为一行，两列
cex1 = 0.9;
# 无标度拓扑拟合指数与软阈值的函数(左图)，下面的会用就行
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

#select line 这条线对应于h的R^2截止点
abline(h=0.7,col="red")#change h=0.00
## Mean Connectivity与软阈值的函数(右图)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#软阈值是WGCNA的算法中非常重要的一个环节，简单的说硬阈值是一种一刀切的算法，
#比如高考分数>500分能上一本，低于500就不行，软阈值的话切起来比较柔和一些，会考虑你学校怎么样，平时成绩怎么样之类。
#网盘中的数据跑起来其实是不太好的，没有合适的软阈值，这根线是要划到0.9的。
#这一步是为了算powers的值。一般来说，powers取红线（0.9）左右的数字都是可以的。
#如果你天秤座特征比较明显，你也可以运行下列代码，让程序推荐你一个值（本案例中返回值是NA，所以后面为了让程序能够进行下去，选了powers=28）。

sft$powerEstimate

#这个powers的值在后面的代码中会一直用到，所以你在跑别的数据的时候一定要更改powers的数值。

##一步法网络构建：One-step network construction and module detection##
#把所有的基因分为不同的基因模块
#deepSplit 参数调整划分模块的敏感度，值越大，越敏感，得到的模块就越多，默认是2；
#minModuleSize 参数设置最小模块的基因数，值越小，小的模块就会被保留下来；
#mergeCutHeight 设置合并相似性模块的距离，值越小，就越不容易被合并，保留下来的模块就越多
#这时的net里就已经包含了基因分类为模块的相关信息了

net = blockwiseModules(datExpr, power = 28, maxBlockSize = 6000,# 表达矩阵，软阈值
                       TOMType = "unsigned", minModuleSize = 30,# 数据为无符号类型，最小模块大小为30
                       reassignThreshold = 0, mergeCutHeight = 0.25,#mergeCutHeight合并模块的阈值，越大模块越少
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)

#结果是每个模块中包含的基因数量。一般来说，结果包含十几个到二十几个模块是比较正常的，此外一个模块中的基因数量不宜过多，
#像我们这个结果里模块1的基因数量达到了2651，这个就有点太多了，主要是因为我们powers=14，软阈值太低了导致的。
#所以说上述软阈值的筛选可以对我们的模块分析起到微调的作用。

##绘画结果展示## open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting# 将标签转换为颜色
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath# 绘制树状图和模块颜色图
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#由于我们的软阈值比较低，所以这一结果中几乎没有grey模块，grey模块中的基因是共表达分析时没有被接受的基因，
#可以理解为一群散兵游勇。当然如果分析结果中grey模块中的基因数量比较多也是不太好的，表示样本中的基因共表达趋势不明显,
#不同特征的样本之间差异性不大，或者组内基因表达一致性比较差。

##结果保存###在这里插入图片描述
#保存模块赋值和模块特征基因信息，以供后续分析。
#moduleLabels代表每个基因对应的模块所属颜色（略去了基因名，但顺序同上图一样)
moduleLabels = net$colors#moduleLabels代表每个基因对应的模块序号（在基因名下面）
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs;#MEs是以基因模块为单位，各个样本在这个模块中的表达量（别管是怎么来的）
geneTree = net$dendrograms[[1]];
#save(MEs, moduleLabels, moduleColors, geneTree,
    # file = "hasm_yueli01.RData")

#这一步就是保存上面跑出来的结果了，同时哪个模块有多少基因一目了然。

##表型与模块相关性##
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.8,   
               xColorWidth =0.01,
               yColorWidth =0.02, 
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1.3, zlim = c(-1.5,1.5)
               , main = paste("Module-trait relationships"))

#cex.lab可以更改X轴Y轴label字体的大小，cex.text可以更改热图中字体的大小，colors可以改变颜色。
#样本特征和共表达模块的相关性热图中，grey模块中的相关性应该很小，如果你与样本特征相关性最显著的模块是grey模块,
#那肯定是有问题的，毕竟grey模块中的基因是一群散兵游勇，它们的表达在各个样本中杂乱无章，根本说明不了问题。

#导出网络到Cytoscape#### Recalculate topological overlap if needed

TOM = TOMsimilarityFromExpr(datExpr, power = 28)

# Read in the annotation file
# annot = read.csv(file = "GeneAnnotation.csv");
# Select modules需要修改，选择需要导出的模块颜色
#modules = c("darkred");
#modules = c("blue");
modules = c("brown");
table(net$colors)
table(moduleColors)


#Select module probes选择模块探测
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule]

#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("3_edges", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("3_nodes", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);
#模块导出后可以用cytoscape构建网络



#如果模块包含的基因太多，网络太复杂，还可以进行筛选，比如：
#filter 里面包含了前八十个基因
nTop = 20;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)


filter <- modTOM[top, top]
write.table(filter, "3_20.txt")
#排序前面八十个基因
top.genes.logical = cbind(modProbes,top)
top.genes.sort = as.matrix (top.genes.logical[order(-top, modProbes)])


#下面的都没有用

#altNodeNames = modGenes,nodeAttr = moduleColors[inModule]); 

##可视化基因网络##
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
##运行时间长
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 3);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function 
#sizeGrWindow(9,9)

TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


#随便选取1000个基因来可视化
nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];

# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
#sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#这里是随机选取1000个基因来可视化模块内基因的相关性，你也可以多取一点，不过取太多容易报错，也没有必要。
#像结果中天青色和蓝色两个模块的共表达聚类结果还是不错的。
#此处画的是根据基因间表达量进行聚类所得到的各模块间的相关性图

MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
sizeGrWindow(7, 6) 
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

#这个是分析共表达模块之间的相关性分析。
