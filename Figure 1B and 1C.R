#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)      #引用包
setwd("C:\\biowolf\\fatty\\17.maftools")    #设置工作目录

#读取单因素的结果文件，获取预后基因的列表
geneRT=read.table("TCGA.uniCox.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)

#绘制瀑布图
pdf(file="oncoplot.pdf", width=8, height=7)
maf=read.maf(maf="input.maf")
oncoplot(maf=maf, genes=gene)
dev.off()

#绘制共突变图形
pdf(file="cooccur.pdf", width=6, height=6)
somaticInteractions(maf=maf, genes=gene, pvalue = c(0.01, 0.05), showSum =F)
dev.off()
