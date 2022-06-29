#install.packages("ggpubr")


#引用包
library(ggpubr)
riskFile="risk.TCGA.txt"     #风险文件
mutFile="mutMatrix.txt"      #突变的矩阵文件
mutGene="TP53 or PLK1"                #选择突变分组基因
setwd("C:\\biowolf\\fatty\\33.mutRisk")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)

#读取突变数据文件
mut=read.table(mutFile, header=T, sep="\t", check.names=F, row.names=1)
mut=t(mut[mutGene,,drop=F])
colnames(mut)=c("Type")

#合并数据
sameSample=intersect(row.names(mut), row.names(risk))
mut=mut[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(as.data.frame(risk), as.data.frame(mut))
data$Type=paste0(mutGene, " " , data$Type)

#设置比较组
data$Type=factor(data$Type, levels=c(paste0(mutGene, " Wild"), paste0(mutGene, " Mutation")) )
group=levels(factor(data$Type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
data1=data[,c("riskScore", "Type")]
boxplot=ggboxplot(data1, x="Type", y="riskScore", fill="Type",
                  xlab="",
                  ylab="Risk score",
                  legend.title="",
                  palette=c("green2", "red2") )+ 
  stat_compare_means(comparisons = my_comparisons)
#输出图片
pdf(file=paste0(mutGene, ".pdf"), width=5, height=4.5)
print(boxplot)
dev.off()
