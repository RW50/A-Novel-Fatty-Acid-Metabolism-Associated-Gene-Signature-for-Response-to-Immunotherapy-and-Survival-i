#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)

pFilter=0.001                #pvalue的过滤条件
expFile="symbol.txt"         #表达数据文件
riskFile="risk.TCGA.txt"     #风险文件
setwd("C:\\biowolf\\fatty\\34.pRRophetic")     #设置工作目录

#获取药物列表
data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
allDrugs=unique(drugData2016$Drug.name)

#读取表达输入文件,并对数据进行处理
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

#读取风险输入文件
riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
riskRT$riskScore[riskRT$riskScore>quantile(riskRT$riskScore,0.99)]=quantile(riskRT$riskScore,0.99)

for(drug in allDrugs){
  #预测药物敏感性
  possibleError=tryCatch(
    {senstivity=pRRopheticPredict(data, drug, selection=1, dataset = "cgp2016")},
    error=function(e) e)
  if(inherits(possibleError, "error")){next}
  senstivity=senstivity[senstivity!="NaN"]
  senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
  
  #将风险文件与药物敏感性的结果进行合并
  sameSample=intersect(row.names(riskRT), names(senstivity))
  risk=riskRT[sameSample, c("riskScore","Risk"),drop=F]
  senstivity=senstivity[sameSample]
  rt=cbind(risk, senstivity)
  
  #设置比较组
  rt$Risk=factor(rt$Risk, levels=c("low", "high"))
  type=levels(factor(rt[,"Risk"]))
  comp=combn(type, 2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #获取高低风险组差异的pvalue
  test=wilcox.test(senstivity~Risk, data=rt)
  diffPvalue=test$p.value
  #获取相关性分析pvalue
  x=as.numeric(rt[,"riskScore"])
  y=as.numeric(rt[,"senstivity"])
  corT=cor.test(x, y, method="spearman")
  corPvalue=corT$p.value
  
  if((diffPvalue<pFilter) & (corPvalue<pFilter)){
    #绘制箱线图
    boxplot=ggboxplot(rt, x="Risk", y="senstivity", fill="Risk",
                      xlab="Risk",
                      ylab=paste0(drug, " senstivity (IC50)"),
                      legend.title="Risk",
                      palette=c("#0066FF","#FF0000")
    )+ 
      stat_compare_means(comparisons=my_comparisons)
    pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
    print(boxplot)
    dev.off()
    #绘制相关性图形
    df1=as.data.frame(cbind(x,y))
    p1=ggplot(df1, aes(x, y)) + 
      xlab("Risk score") + ylab(paste0(drug, " senstivity (IC50)"))+
      geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    #相关性图形
    pdf(file=paste0("Cor.", drug, ".pdf"), width=5, height=4.6)
    print(p1)
    dev.off()
  }
}