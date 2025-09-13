######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("survival")
#install.packages("survminer")
#install.packages("forestplot")


#引用包
library(limma)
library(survival)
library(survminer)
library(forestplot)

pFilter=0.05                #KM方法pvalue过滤条件
expFile="geneExp.txt"       #表达数据文件
cliFile="Survival_SupplementalTable_S1_20171025_xena_sp"      #临床数据文件
setwd("C:\\biowolf\\panImmune\\12.DFS")         #设置工作目录

############绘制森林图函数############
bioForest=function(coxFile=null, forestFile=null, forestCol=null, titleName=null){
    #读取输入文件
	rt=read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	data=as.matrix(rt)
	HR=data[,1:3]
	hr=sprintf("%.3f",HR[,"HR"])
	hrLow=sprintf("%.3f",HR[,"HR.95L"])
	hrHigh=sprintf("%.3f",HR[,"HR.95H"])
	pVal=data[,"pvalue"]
	pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
	#定义颜色
	clrs <- fpColors(box=forestCol,line="darkblue", summary="royalblue")      #定义颜色
	#定义图形的文字
	tabletext <- 
	  list(c(NA, rownames(HR)),
	       append("pvalue", pVal),
	       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   #定义图片文字
	#绘制森林图
	pdf(file=forestFile, width=9, height=6, onefile = FALSE)
	forestplot(tabletext, 
	           rbind(rep(NA, 3), HR),
	           col=clrs,
	           graphwidth=unit(50, "mm"),
	           xlog=T,
	           lwd.ci=4,
	           boxsize=0.6,
	           title=titleName,
	           xlab="Hazard ratio",
	           txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))
	           )
	dev.off()
}
############绘制森林图函数############

#读取表达数据文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]
gene=colnames(exp)[1]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("DFI.time","DFI","DSS.time","DSS","PFI.time","PFI")]

#数据合并
sameSample=intersect(row.names(cli), row.names(exp))
data=cbind(cli[sameSample,,drop=F], exp[sameSample,,drop=F])

#定义森林图颜色
cols=c("#EEEE00", "#00DD00", "#FF00FF")
project=c("DFS", "DSS", "PFS")
showNames=c("Disease free survival", "Disease specific survival", "Progression free survival")

#对DFS、DSS和PFS分别进行循环，得到生存分析结果
for(i in 1:3){
	#提取生存分析的数据
	k=2*i-1
	rt=data[,c(k,k+1,7:ncol(data))]
	rt=na.omit(rt)
	colnames(rt)[1:3]=c("futime","fustat","gene")
	rt$futime=rt$futime/365
	
	#对肿瘤类型进行循环，进行KM分析和COX分析
	outTab=data.frame()
	for(cancerType in levels(factor(rt[,"CancerType"]))){
		rt1=rt[(rt[,"CancerType"]==cancerType),]
		if(nrow(rt1)<3){next}
		#根据基因表达量对样品进行分组，比较高低表达组的生存差异
		group=ifelse(rt1[,"gene"]>median(rt1[,"gene"]), "high", "low")
		diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
		pValue=1-pchisq(diff$chisq,df=1)
		if(pValue<pFilter){
			if(pValue<0.001){
				pValue="p<0.001"
			}else{
				pValue=paste0("p=",sprintf("%.03f",pValue))
			}
			fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
			#绘制生存曲线
			surPlot=ggsurvplot(fit, 
					    data=rt1,
					    title=paste0("Cancer: ",cancerType),
					    pval=pValue,
					    pval.size=6,
					    conf.int=F,
					    legend.title=paste0(gene," levels"),
					    legend.labs=c("high","low"),
					    font.legend=12,
					    fontsize=4,
					    xlab="Time(years)",
					    ylab=showNames[i],
					    break.time.by=2,
					    palette=c("red","blue"),
					    risk.table=TRUE,
					    risk.table.title="",
					    risk.table.height=.25)
			#输出生存曲线
			pdf(file=paste0(project[i],"_",cancerType,".pdf"), onefile=FALSE, width=6, height=5)
			print(surPlot)
			dev.off()
		}
		
		#COX分析
		cox=coxph(Surv(futime, fustat) ~ gene, data = rt1)
		coxSummary=summary(cox)
	    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
		outTab=rbind(outTab,
		             cbind(cancer=cancerType,
		                   HR=coxSummary$conf.int[,"exp(coef)"],
		                   HR.95L=coxSummary$conf.int[,"lower .95"],
		                   HR.95H=coxSummary$conf.int[,"upper .95"],
				           pvalue=coxP) )
	}
	
	#输出COX结果并且绘制森林图
	write.table(outTab,file=paste0(project[i],".cox.txt"),sep="\t",row.names=F,quote=F)
	bioForest(coxFile=paste0(project[i],".cox.txt"),
	          forestFile=paste0(project[i],".forest.pdf"),
	          forestCol=cols[i],
	          titleName=showNames[i])
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

