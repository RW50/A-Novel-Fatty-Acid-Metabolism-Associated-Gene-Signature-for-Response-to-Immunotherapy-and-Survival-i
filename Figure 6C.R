#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("digest")
#install.packages("GOplot")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)

pvalueFilter=0.05      #p值过滤条件
qvalueFilter=0.05      #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

setwd("C:\\biowolf\\fatty\\38.KEGG")                    #设置工作目录
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]      #删除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#通路富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=8, height=7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format=50, color = colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=8, height=7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", label_format=50, color = colorSel)
dev.off()

#获取通路的信息
kegg=data.frame(Category="ALL", ID = KEGG$ID, Term=KEGG$Description, Genes = gsub("/", ", ", KEGG$geneID), adj_pval = KEGG$p.adjust)
#读取基因的差异情况
genelist <- data.frame(ID = rt$gene, logFC = rt$logFC)
row.names(genelist)=genelist[,1]
#设置圈图参数
circ <- circle_dat(kegg, genelist)
termNum =8         #设置显示通路的数目
termNum=ifelse(nrow(kegg)<termNum,nrow(kegg),termNum)
geneNum=300        #显示基因的数目
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#绘制通路的圈图
chord <- chord_dat(circ, genelist[1:geneNum,], kegg$Term[1:termNum])
pdf(file="KEGGcircos.pdf", width=11, height=11)
GOChord(chord, 
        space = 0.001,           #基因之间的距离
        gene.order = 'logFC',    #基因的排序方式，按照logFC值对基因排序
        gene.space = 0.25,       #基因名称与圆圈的相对距离
        gene.size = 5,           #基因名称字体的大小 
        border.size = 0.1,       #线条粗细
        process.label = 6)       #通路名称字体的大小
dev.off()