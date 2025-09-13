######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#BiocManager::install("limma")


#引用包
library(limma)
library(ComplexHeatmap)
expFile="geneExp.txt"       #表达数据文件
cliFile="clinical.txt"      #临床数据文件
setwd("C:\\biowolf\\Gene\\14.cliHeatmap")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

#删掉正常样品
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#根据目标基因表达量对样品进行分组
Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
Type=factor(Type, levels=c("Low","High"))
data=cbind(as.data.frame(data), Type)
data=data[order(data[,gene]),] 

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

#合并数据
samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,"Type",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)

#对临床性状进行循环，观察临床性状在高低表达组之间是否具有差异，得到显著性标记
sigVec=c(gene)
for(clinical in colnames(rt[,2:ncol(rt)])){
	data=rt[c("Type", clinical)]
	colnames(data)=c("Type", "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	tableStat=table(data)
	stat=chisq.test(tableStat)
	pvalue=stat$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(clinical, Sig))
}
colnames(rt)=sigVec

#定义热图注释的颜色
#rt=rt[apply(rt,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
colorList[[gene]]=c("Low"="blue", "High"="red")
j=0
for(cli in colnames(rt[,2:ncol(rt)])){
	cliLength=length(levels(factor(rt[,cli])))
	cliCol=bioCol[(j+1):(j+cliLength)]
	j=j+cliLength
	names(cliCol)=levels(factor(rt[,cli]))
	cliCol["unknow"]="grey75"
	colorList[[cli]]=cliCol
}

#绘制热图
ha=HeatmapAnnotation(df=rt, col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha)

#输出热图
pdf(file="heatmap.pdf", width=7, height=5)
draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

