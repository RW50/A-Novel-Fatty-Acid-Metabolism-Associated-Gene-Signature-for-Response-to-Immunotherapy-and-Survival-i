#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")


#引用包
library(limma)
library(ggplot2)

expFile="FAMexp.txt"         #表达数据文件
riskFile="risk.TCGA.txt"     #风险文件
setwd("C:\\biowolf\\fatty\\19.PCA")     #设置工作目录

######绘制模型基因的PCA图######
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
data=risk[,3:(ncol(risk)-2)]      #提取模型基因表达量
Risk=as.vector(risk[,"Risk"])     #提取病人的风险信息
#对模型基因进行PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
PCA = data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Risk)
#绘制模型基因的PCA图
pdf(file="PCA.modelGene.pdf", width=5.5, height=4.5)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Risk)) +
  scale_colour_manual(name="Risk",  values =c("red", "blue"))+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


######绘制脂肪酸代谢基因的PCA图######
#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#删除正常样品
type=sapply(strsplit(colnames(data),"\\-"),"[",4)
type=sapply(strsplit(type,""),"[",1)
type=gsub("2","1",type)
data=t(data[,type==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))

#数据合并
sameSample=intersect(rownames(data),rownames(risk))
data=data[sameSample,]
Risk=risk[sameSample,"Risk"]

#PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
PCA = data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Risk)

#绘制脂肪酸代谢基因的PCA图
pdf(file="PCA.FAMgene.pdf", width=5.5, height=4.5)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Risk)) +
  scale_colour_manual(name="Risk",  values =c("red", "blue"))+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()
