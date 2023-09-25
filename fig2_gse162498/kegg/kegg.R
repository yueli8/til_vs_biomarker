library(ggplot2)
setwd("~/gse162498/kegg")
a<-read.table("Terminally_Exhausted_CD8", head=TRUE, sep="\t")
#fill=padj fill颜色填充，使用fdr
p <- ggplot(data=a,aes(x=Term,y=Count,fill=FDR))
#coord_flip()颠倒坐标轴
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(color='gray'),
                 axis.text.y=element_text(color="black",face='bold',size=20),axis.text.x = element_text(color="black",face='bold',size=20))
#ylim(0,30) 更改横坐标的范围这里坐标轴颠倒了，虽然看起来是x轴，但其实是y轴
p3 <- p2 + ylim(0,35) + scale_fill_gradient(low="red",high="blue") 
p4 <- p3 + scale_x_discrete(limits=rev(a[,1])) +labs(x="",y="",face="bold",title="KEGG")
p4
