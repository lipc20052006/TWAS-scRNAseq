rm(list=ls())
library(Seurat)  
library(reticulate)  
library(tidyverse)   
library(patchwork)
#devtools::install_github("sajuukLyu/ggunchull", type = "source")
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
library(readxl)
library("writexl")
library(scales)
library(viridis)
#install.packages("nord")
library(nord)

getwd()
sc=readRDS("../../04_Cluster/All_sample_combined.ck.rds")
raw.ident=Idents(sc) 
sc@meta.data$raw.ident=Idents(sc)
head(sc@meta.data)
table(sc$orig.ident)
Idents(sc)="orig.idents" 
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 2)
VlnPlot(sc, features = c("nFeature_RNA"), ncol = 2)

Idents(sc)=raw.ident
head(Idents(sc))
####1.ump####
pdf('ump.pdf', width=6,height=4)
DimPlot(sc, reduction = "umap",label=T,pt.size=1)
dev.off()

DimPlot(sc, reduction = "umap",label=T)

####2.marker####
sc.markers <- FindAllMarkers(sc, 
                                only.pos = TRUE,  
                                min.pct = 0.25, 
                                logfc.threshold = 0.25
)

marker=read.csv("./marker.csv",header = T) %>% pull(Gene)
DotPlot(sc, features =marker ) + RotatedAxis()

pdf('dotplot0-20.pdf', width=15,height=10)
DotPlot(sc, features =marker ) +ggplot2:::coord_flip()
dev.off()

sc@active.ident=factor(sc@active.ident,
                       levels =c(0,3,9,19,4,5,8,12,14,6,16,7,1,10,11,13,15,17,18,20,2) )
pdf('dotplotSort.pdf', width=15,height=10)
DotPlot(sc, features =marker ) +ggplot2:::coord_flip()
dev.off()

