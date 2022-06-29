#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)

expFile="geneExp.txt"       #表达数据文件
cliFile="clinical.txt"      #临床数据文件
setwd("C:\\biowolf\\Gene\\13.cliCor")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

#删掉正常样品
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

#合并数据
samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)

#临床相关性分析，输出图形结果
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c(gene, clinical)]
  colnames(data)=c(gene, "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  #设置比较组
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #绘制箱线图
  boxplot=ggboxplot(data, x="clinical", y=gene, fill="clinical",
                    xlab=clinical,
                    ylab=paste(gene, " expression"),
                    legend.title=clinical)+ 
    stat_compare_means(comparisons = my_comparisons)
  #输出图片
  pdf(file=paste0("clinicalCor_", clinical, ".pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
}
