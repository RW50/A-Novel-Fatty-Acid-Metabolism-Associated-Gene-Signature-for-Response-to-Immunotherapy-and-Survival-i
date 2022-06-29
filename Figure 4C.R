#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("pheatmap")


#引用包
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

expFile="symbol.txt"          #表达数据文件
riskFile="risk.TCGA.txt"      #风险文件
gmtFile="c2.cp.kegg.v7.4.symbols.gmt"     #基因集文件
setwd("C:\\biowolf\\fatty\\32.GSVA")      #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#GSVA分析
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)

#删除正常样品
data=t(gsvaResult)
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#读取risk文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"Risk",drop=F]
gsvarisk=cbind(data, risk)

#样品分组
con=gsvarisk[gsvarisk$Risk=="low",]
treat=gsvarisk[gsvarisk$Risk=="high",]
data=rbind(con, treat)

#对通路进行差异分析
Type=as.vector(data$Risk)
Type=factor(Type, levels=c("low", "high"))
ann=data[,ncol(data),drop=F]
data=t(data[,-ncol(data),drop=F])
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
fit=lmFit(data, design)
cont.matrix=makeContrasts(high-low, levels=design)
fit2=contrasts.fit(fit, cont.matrix)
fit2=eBayes(fit2)

#输出所有通路的差异情况
allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#输出差异显著的通路
diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < 0.05)), ]
diffSigOut=rbind(id=colnames(diffSig), diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

#定义热图注释的颜色
ann_colors=list()
bioCol=c("blue", "red")
names(bioCol)=c("low", "high")
ann_colors[["Risk"]]=bioCol

#绘制差异通路热图
termNum=50     #设置显示通路的数目
diffTermName=as.vector(rownames(diffSig))
diffLength=length(diffTermName)
if(diffLength<termNum){termNum=diffLength}
hmGene=diffTermName[1:termNum]
hmExp=data[hmGene,]
pdf(file="heatmap.pdf", width=10, height=6)
pheatmap(hmExp, 
         annotation=ann,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
         cluster_cols=F,
         show_colnames = F,
         gaps_col=as.vector(cumsum(table(Type))),
         scale="row",
         fontsize = 7,
         fontsize_row=6,
         fontsize_col=7)
dev.off()