library(Seurat)  
library(reticulate)  
library(tidyverse)   
library(patchwork)
source("3-Bulk反卷积/Cibersort_function.R")
library(data.table)
library("psych")
library("ggpubr")
rm(list=ls())
sc=readRDS("../../04_Cluster/All_sample_combined.ck.rds")
rm(immune.combined)
exp <- read.table("./tpm.txt",header = T)
sc_exp <- AverageExpression(sc,assays = "RNA",group.by = "seurat_clusters")
sc_exp <- as.data.frame(sc_exp$RNA)
colnames(sc_exp) <- paste0("cluster",colnames(sc_exp))

sc.markers=readRDS("sc.markers.rds")
top50 <- sc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

sc_exp50 <- sc_exp[row.names(sc_exp)%in%top50$gene,]

write.table(x = sc_exp50,"./sc_exp50.tsv",quote=F,sep = "\t")

####2.反卷积####
CIBERSORT_res50 <- CIBERSORT('./sc_exp50.tsv',"./tpm.w.txt", perm = 1000, QN = T)  #perm置换次数=1000，QN分位数归一化=TRUE

CIBERSORT_res50=as.data.frame(CIBERSORT_res50)


write.table(CIBERSORT_res50,"./cibersort378Bulk50.tsv",sep="\t",quote = F,row.names=F)


##LRL-----
scBulk=fread("3-Bulk反卷积/cibersort378Bulk50.tsv")

phe=read.table("./0-LRL.txt",sep="\t",header=T)

Lines=read.table("./Lines357.txt")[2:358,]##
phe.sc=left_join(phe,scBulk[,1:22],by="SampleID")
cor.r=corr.test(phe.sc[,2],phe.sc[,3:23])$r

cor.p=corr.test(phe.sc[,2],phe.sc[,3:23])$p

cor.rp=as.data.frame(cbind(t(cor.r),t(cor.p)))
names(cor.rp)=c("r","p.value")


p50=cor.rp50 %>% filter(p.value<0.05) %>% arrange(r) %>% 
  mutate(cluster=row.names(.)) %>%
  mutate(cluster=factor(cluster,levels = .$cluster)) %>% 
  ggplot()+
  geom_point(aes(cluster,r,size=abs(r),color=r))+
  scale_size(range = c(3, 8))+
  scale_color_gradient(low="green",high="red")+
  geom_segment( aes(x=cluster, xend=cluster, y=0, yend=r),color="grey",
                size=1.5,linetype=1)+#使用reorder()排序变量
  geom_text(aes(cluster,r,label =sprintf("%.2f",r)), color = "black", size = 3)+ 
  scale_y_continuous(trans = "reverse")+
  geom_hline(yintercept=0)+
  xlab("")+ ylab("")+ coord_flip()+ 
  theme_classic()+theme(axis.text=element_text(color="black"))
p50
ggsave("./deconvolution.corTOP50cLRL.pdf",plot = p50,width=15,height=9,units="cm",dpi=600)  
ggsave("./deconvolution.corTOP50cLRL.png",plot = p50,width=15,height=9,units="cm",dpi=600)  
write.table(cor.rp,"./deconvolution.cor.LRL.csv",sep=",",quote = F)
