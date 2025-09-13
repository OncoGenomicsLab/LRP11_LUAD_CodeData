######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

expFile="geneExp.txt"          #表达数据文件
scoreFile="TMEscores.txt"      #肿瘤微环境的打分文件
setwd("C:\\biowolf\\Gene\\24.TMEvioplot")     #设置工作目录

#读取表达数据文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(exp)[1]
exp=exp[exp$Type=="Tumor",1,drop=F]
exp$Type=ifelse(exp[,gene]>median(exp[,gene]), "High", "Low")

#读取肿瘤微环境打分文件，并对数据进行整理
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]

#样品取交集
sameSample=intersect(row.names(exp), row.names(score))
exp=exp[sameSample,"Type",drop=F]
score=score[sameSample,,drop=F]
rt=cbind(score, exp)
rt$Type=factor(rt$Type, levels=c("Low", "High"))

#将合并后的数据转换为ggplot2的输入文件
data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "scoreType", "Score")

#绘制小提琴图
p=ggviolin(data, x="scoreType", y="Score", fill = "Type",
	     xlab="",
	     ylab="TME score",
	     legend.title=gene,
	     add = "boxplot", add.params = list(color="white"),
	     palette = c("blue","red"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出图形
pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

