library(Seurat)
#library(SeuratData)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(magrittr)
library(AUCell)
#BiocManager::install("AUCell")
library(patchwork)
library(ggplot2)
library(readxl)
Sys.setenv(LANGUAGE = "en") #
options(stringsAsFactors = FALSE) #
rm(list=ls())

LRL.gene=read_xlsx("./6-twas.fdr.05.647LRL.gene.xlsx")##
LRL.gene=LRL.gene %>% pull(gene)
geneSets <- list(LRL=LRL.gene)


sc=readRDS("../../04_Cluster/All_sample_combined.ck.rds")
# Build gene expression rankings for each cell
sc@assays$RNA@data
sc@assays[["RNA"]]@data
cells_rankings <- AUCell_buildRankings(sc@assays$RNA@data, nCores=10)  #

# Quantiles for the number of genes detected by cell:
# (Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
# Calculates the 'AUC' for each gene-set in each cell
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores = 10, 
                            aucMaxRank = nrow(cells_rankings)*0.05) 
# 
pdf("./4-AuCell/1-cells_assignmentLRL.pdf")
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE, nCores=10)
dev.off()
# 提取UMAP坐标和AUC值，在UMAP图上展示AUC值
geneSet.name <- "LRL"
AUC_Exp <- as.numeric(getAUC(cells_AUC)[geneSet.name, ])
sc$AUC <- AUC_Exp
plot.df<- data.frame(sc@meta.data, sc@reductions$umap@cell.embeddings)

p <- ggplot() + 
  geom_point(data=plot.df, aes(x=UMAP_1,y=UMAP_2,colour=AUC), size =1) +
  geom_point(data=subset(plot.df,AUC>0.05), aes(x=UMAP_1,y=UMAP_2,colour=AUC), size =1.5) +
  scale_colour_gradient2(low = "#E8EAF6", mid="#E8EAF6",high = "red",midpoint=0.07) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text= element_text(colour= 'black',size=14),
        axis.title= element_text(size = 14),
        axis.line= element_line(colour= 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"), 
        aspect.ratio = 1)
p

ggsave(filename = "./4-AuCell/2-UMAP_AUC_value1LRL.pdf", plot = p1, height = 6, width = 6)
