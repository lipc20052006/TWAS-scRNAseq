####1.library packages####
#install.packages("devtools");
#install_github("cheuerde/cpgen", ref = "master", build_vignettes=FALSE)
#install.packages("psych");install.packages("reshape2")
#install.packages("VennDiagram") ##下载并加载包
library("tidyverse");library("ggplot2")
library(UpSetR);library(VennDiagram)
library("reshape2");library("psych")
library(devtools);library(cpgen)
library(data.table);library(dplyr)
library("readxl")
library("writexl")
library("PlantNGSTools")
####2.test####
rm(list=ls())

####3.read phe exp PCA####
getwd()
###3.1 读入表达数据和表型，处理成每个基因一列的矩阵，并确定表达和表型材料排序完全一致###
####3.1 表型数据处理####

phe=read.table("./phe.LRL.txt",sep="\t",header=T)
Lines=read.table("./Lines357.txt")[2:358,]###
##Lines=phe$SampleID
phe=left_join(data.frame(SampleID=Lines),phe,by="SampleID")#
trait="LRL"
#phenotypic distribute


p.dist=ggplot(phe, aes(x=LRL)) +geom_density(fill="darkgreen",colour="white")+
  theme_classic()+ylab("Density")+xlab("LRL")
ggsave("3-plot/1-phenotypic distributeLRL.tiff",plot=p.dist,
       width=12,height=6,units="cm",dpi=300)
ggsave("3-plot/1-phenotypic distributeLRL.pdf",plot=p.dist,
       width=12,height=6,units="cm",dpi=300) 

####3.2 exp matix####
exp=fread("./tpm.txt",sep="\t",header = T)
if(!identical(phe$SampleID,exp$SampleID) ){
  print("Sort ID") 
exp=left_join(data.frame(SampleID=Lines),exp,by="SampleID")}
identical(phe$SampleID,exp$SampleID) 
row.names(exp)=exp$SampleID;exp$SampleID=NULL 
identical(phe$SampleID,row.names(exp))

#3.3 read kinship####===========================================
kin=read.table("1-Kinship.txt",sep="\t",header=F,skip = 1)
identical(phe$SampleID,kin$V1)
kin=as.matrix(kin[,-1])
####4.TWAS####
####4.1 运行TWAS####=========================================
m=as.matrix(exp)
y=phe$LRL
twas.p=cGWAS.emmax(y,m,verbose=TRUE,A = as.matrix(kin))

write.table(twas.p,"2-twas/1-TWAS.csv",sep=",",row.names=F,col.names=T)

####4.2  Calculate correlation coefficient and P-Value######=========================================

zm.cor.r=matrix(nrow=ncol(exp),ncol=length(trait));zm.cor.p=matrix(nrow=ncol(exp),ncol=length(trait))
colnames(zm.cor.p)=colnames(zm.cor.r)=trait
rownames(zm.cor.p)=rownames(zm.cor.r)=colnames(exp)
for (i in 1:length(trait)){
  print(trait[i])
  cor.r=vector(mode="numeric",length=0)
  cor.p=vector(mode="numeric",length=0)
  t=phe[,trait[i]]
  cor.r=as.numeric(corr.test(t,exp.t)$r)
  cor.p=as.numeric(corr.test(t,exp.t)$p)
  zm.cor.r[,i]=cor.r;  zm.cor.p[,i]=cor.p
}

write.table(zm.cor.p,"2-twas/gene.trait.cor.p.csv",sep=",")
write.table(zm.cor.r,"2-twas/gene.trait.cor.r.csv",sep=",")

####5.FDR####=========================================================
twas.fdr=as.data.frame(twas.p)
for (i in 1:length(trait)){
  #i=1
  print(trait[i])
  twas.fdr[,i+1]=p.adjust(twas.fdr[,i+1],method ="fdr",n=length(twas.fdr[,i+1]))
}

write.table(twas.fdr,"2-twas/1-TWAS.fdr.csv",sep=",",row.names=F,col.names=T)









