######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)

gene="VCAN"                 #基因名称
expFile="exp.txt"           #表达数据文件
clikFile="clinical.txt"     #临床数据文件
setwd("C:\\biowolf\\panImmune\\24.IMvigor")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#如果数据没有取log2，可以进行log2处理
qx=as.numeric(quantile(data, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ((qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	data[data<0]=0
	data=log2(data+1)}
#data=normalizeBetweenArrays(data)

#提取目标基因表达量
exp=t(data[gene,,drop=F])

#读取临床数据文件
clinical=read.table(clikFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(clinical)[1]

#合并数据
sameSample=intersect(row.names(exp), row.names(clinical))
exp=exp[sameSample,,drop=F]
clinical=clinical[sameSample,,drop=F]
data=cbind(as.data.frame(exp), as.data.frame(clinical))

#设置比较组
group=levels(factor(data[,cliName]))
data[,cliName]=factor(data[,cliName], levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
boxplot=ggboxplot(data, x=cliName, y=gene, fill=cliName,
			      xlab="",
			      ylab=paste0(gene, " expression"),
			      legend.title="",
			      palette = c("blue", "red") )+ 
	stat_compare_means(comparisons = my_comparisons)	

#输出箱线图
pdf(file=paste0(gene, ".boxplot.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

