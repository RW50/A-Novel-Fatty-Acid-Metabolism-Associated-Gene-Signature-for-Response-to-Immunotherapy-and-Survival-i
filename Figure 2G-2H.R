#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)
setwd("C:\\biowolf\\fatty\\20.survival")     #设置工作目录

#定义生存分析的函数
bioSurvival=function(inputFile=null, outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #比较高低风险组的生存差异，得到差异显著性pvalue
  diff=survdiff(Surv(futime, fustat) ~ Risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     ylab="Overall survival",
                     break.time.by = 1,
                     palette=c("red", "blue"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, width=6.5, height=5.5, onefile=FALSE)
  print(surPlot)
  dev.off()
}

#调用函数,绘制生存曲线
bioSurvival(inputFile="risk.TCGA.txt", outFile="survival.TCGA.pdf")
bioSurvival(inputFile="risk.GEO.txt", outFile="survival.GEO.pdf")
