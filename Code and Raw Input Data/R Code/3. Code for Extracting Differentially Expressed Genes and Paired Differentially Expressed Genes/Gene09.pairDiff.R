######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
expFile="geneExp.txt"      #表达数据文件
setwd("C:\\biowolf\\Gene\\09.pairDiff")      #设置工作目录

#读取表达输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
geneName=colnames(rt)[1]

#提取配对样品的数据
normalData=rt[rt$Type=="Normal",1,drop=F]
normalData=as.matrix(normalData)
rownames(normalData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(normalData))
normalData=avereps(normalData)
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
tumorData=avereps(tumorData)
sameSample=intersect(row.names(normalData), row.names(tumorData))
data=cbind(normalData[sameSample,,drop=F], tumorData[sameSample,,drop=F])
colnames(data)=c("Normal", "Tumor")
data=as.data.frame(data)

#绘制配对差异分析的图形
pdf(file="pairDiff.pdf", width=5, height=4.5)
ggpaired(data, cond1="Normal", cond2="Tumor", fill="condition",
	xlab="", ylab=paste0(geneName, " expression"),
	legend.title="Type",
	palette=c("blue","red"))+
    #stat_compare_means(paired = TRUE, label = "p.format", label.x = 1.35)
    stat_compare_means(paired = TRUE, symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",label.x = 1.35)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

