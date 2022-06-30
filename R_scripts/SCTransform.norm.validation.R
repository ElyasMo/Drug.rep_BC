

library(ggplot2)
library(SingleR)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(stringr)
library(slinky)
library(ggthemr)
library(tidyverse)




# Step1. Data loading
## Loading and clustering \nthe data


Normal <- readRDS("Directory/Normal.rds")


Normal <- lapply(Normal, function(X){
  X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = "^MT-");X
})


# filtering the data based on number of features and percentage of mitochondrial genes
nfeat <- c(6000, 6000, 6000, 6000, 4000,7000,6000, 6000, 6000, 4000,2000,4000,7000)
mt <- c(10,15,40,7.5,20,15,20,20,20,20,5,20,30)

for (num in c(1:13)) {
  Normal[[num]] <- subset(Normal[[num]], subset = nFeature_RNA > 300 & nFeature_RNA < nfeat[num] &
                            percent.mt < mt[num])
}

# Merging the objects
Normal<- Reduce(merge,Normal )
# merge(x = Normal$Normal1, y = c(Normal[-1]))


Normal <- SCTransform(Normal, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.3), verbose = FALSE)

#Choosing the proper resolution
Normal <- SetIdent(Normal, value = Normal@meta.data$SCT_snn_res.0.1)


BC <- readRDS("Directory/BC7_new.rds")

BC <- lapply(BC, function(X){
  X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = "^MT-");X
}) 


BC <-merge(x = BC$CID3586, y = c(BC$CID3838,BC$CID3921,BC$CID3941,BC$CID3946,BC$CID3948,BC$CID3963,BC$CID4040, 
                                 BC$CID4066,BC$CID4067,BC$CID4290A,BC$CID4398,BC$CID44041,BC$CID4461,BC$CID4463, 
                                 BC$CID4465,BC$CID4471,BC$CID4495,BC$CID44971,BC$CID44991,BC$CID4513,BC$CID4515, 
                                 BC$CID45171,BC$CID4523,BC$CID4530N,BC$CID4535 ))



Malignant <- SplitObject(BC, split.by = "celltype_major")
Malignant <- Malignant$`Cancer Epithelial`

Malignant$Type <- "Malignant"



protein_normalized <- read_csv("Directory/protein_normalized.csv")
protein_normalized <- as.data.frame(protein_normalized[!duplicated(protein_normalized$Gene_Symbol),])

protein_normalized <- protein_normalized[complete.cases(protein_normalized$Gene_Symbol), ]
rownames(protein_normalized) <- protein_normalized$Gene_Symbol
protein.table <- protein_normalized[,-c(1:6, 20:22)]
protein.table <- protein.table[,-c(14,17,24,26)]
protein.table[is.na(protein.table)] = 0

protein_seurat <- CreateSeuratObject(protein.table, min.cells = 3, min.genes = 200, project = "Protein" )
protein_seurat@assays$RNA@data <- protein_seurat@assays$RNA@counts
protein_seurat <- FindVariableFeatures(protein_seurat, nfeatures = 6000)
protein_seurat <- ScaleData(protein_seurat, verbose = FALSE)
protein_seurat <- RunPCA(protein_seurat, verbose = FALSE,npcs = 15,approx=FALSE)
protein_seurat <- RunUMAP(protein_seurat, dims = 1:15, verbose = FALSE,n.neighbors =24 )
protein_seurat$subtype <- c(rep("HER2",6), rep("ER",7), rep("TNBC",11))
protein_seurat <- SetIdent(protein_seurat, value = protein_seurat$subtype)


path="Directory"
setwd(path)
save(Normal,Malignant,protein_seurat, file = "SCT_Inputs_step1.RData" ) 



## UMAPs



DimPlot(protein_seurat, raster = FALSE)+ggtitle("Proteomics data")+theme(title =element_text(size=16, face='bold'),
                                                                         legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                         axis.title=element_text(size=14,face="bold"))



DimPlot(Normal, raster = FALSE)+ggtitle("Normal cells from breast tissue")+theme(title =element_text(size=16, face='bold'),
                                                                                 legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                                 axis.title=element_text(size=14,face="bold"))


###############################
## Annotation of Normal cells #
###############################

path="Directory"
setwd(path)
load(file = "SCT_Inputs_step1.RData" )

#Changing the Seurat object to sSinglecell object
Normal.sce <- as.SingleCellExperiment(Normal,assay = "SCT")

#Determining the path to the annotation libraries
wpath="Directory"

#Load references
blue <- readRDS(paste(wpath,"singleR.BlueprintEncodeData.rds",sep="/"))
hpca <- readRDS(paste(wpath,"singleR.HumanPrimaryCellAtlasData.rds",sep="/"))


#Combining to libraries for a better annotation
ref <- list(BP=blue,HPCA=hpca)

# Disinguishing the main and fine libraries (main has the major annotations and fine includes the cell subtypes as well)
main.labels <- list(blue$label.main, hpca$label.main)
fine.labels <-  list(blue$label.fine, hpca$label.fine)

#Annotationg the cells. It is time consuming so I saved it in .rds files for the sake of time and computational process in code repetition.
# pred.Normal.main <- SingleR(test=Normal.sce, ref=ref, assay.type.test=1, labels=main.labels)
# pred.Normal.fine <- SingleR(test=Normal.sce, ref=ref, assay.type.test=1, labels=fine.labels)
# saveRDS(pred.Normal.main,"../New/SCT/pred.Normal.main.rds")
# saveRDS(pred.Normal.fine,"../New/SCT/pred.Normal.fine.rds")

#Loading the annotations
pred.Normal.main <- readRDS("Directory/pred.Normal.main.rds")
pred.Normal.fine <- readRDS("Diretory/pred.Normal.fine.rds")

#Adding the annotations to the objects 
Normal.celllabels2 <- cbind(celllabelsb=pred.Normal.main$labels)
Normal.celllabels1 <- cbind(celllabels=pred.Normal.fine$labels)
rownames(Normal.celllabels1) <- rownames(Normal.celllabels2) <- colnames(Normal)
Normal <- AddMetaData(Normal, as.data.frame(cbind(Normal.celllabels1,Normal.celllabels2)))

## Evaluate annotation based on scores
# plotScoreHeatmap(pred.Normal.main ,max.labels=20, fontsize=7,cluster_cols = TRUE)

# par(mfrow=c(1,1), mar=c(2.5,7,1,.5))
# graphics::boxplot(pred.Normal.fine$scores, xlab="prediction.score", ylab="", main="Prediction scores for cell labels", boxwex=.8,
#                   cex=.5, cex.main=.8, cex.lab=.8, cex.axis=.6, las=1, horizontal=T, varwidth=T, col=rainbow(36))

#Tables of cell annotations
tab.main <- table(Normal$celllabelsb, Normal$SCT_snn_res.0.1)
tab.fine <- table(Normal$celllabels, Normal$SCT_snn_res.0.1)

knitr::kable(tab.main, caption = "Normal-main")
knitr::kable(tab.fine, caption = "Normal-fine")



## label the annotations

#Adding new cell labels
cluster.ids.Normal <- c("Fibroblasts","Epithelial cells", "Keratinocytes",  "Epithelial cells","Endothelial cells","Adipocytes", "T Cells", "Monocytes", "Fibroblasts")

names(cluster.ids.Normal) <- levels(Normal)
Normal <- RenameIdents(Normal, cluster.ids.Normal)


# Extracting the Epithelial cells to be compared with the epithelial malignant cells
Normal$annotation <- Normal@active.ident
Normal$subtype <- "Normal"

Normal.Epithelial <- SplitObject(Normal, split.by = "annotation")
Normal.Epithelial <- Normal.Epithelial$`Epithelial cells`

Normal.Epithelial$Type <- "Normal"


#Ploting to see the new annotated cells
DimPlot(Normal, raster = FALSE)+ggtitle("Cell type annotation using \nSCT normalization")+theme(title =element_text(size=16, face='bold'),
                                                                                                legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                                                axis.title=element_text(size=14,face="bold"))


DimPlot(Normal.Epithelial, raster = FALSE)+ggtitle("Normal epithelial cells from breast tissue")+theme(title =element_text(size=16, face='bold'),
                                                                                                       legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                                                       axis.title=element_text(size=14,face="bold"))


## All objects

Malignant <- SCTransform(Malignant, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

Malignant <- SetIdent(Malignant, value = Malignant$Type)


Epithelial.list <- list(ctrl = Normal.Epithelial, stim = Malignant)
features <- SelectIntegrationFeatures(object.list = Epithelial.list, nfeatures = 3000)
Epithelial.list <- PrepSCTIntegration(object.list = Epithelial.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = Epithelial.list, normalization.method = "SCT",
                                  anchor.features = features)

combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells:

combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindClusters(combined.sct, resolution = 0.1)

Idents(combined.sct) <- "Type"


combined.sct <- PrepSCTFindMarkers(combined.sct)


path="Directory"
setwd(path)
save(Normal.Epithelial,Malignant,protein_seurat, combined.sct, file = "Directory/SCT_Inputs_step2-DataLoading.RData" )


DimPlot(combined.sct, raster = FALSE, group.by = "subtype")+ggtitle("scRNAseq data integration using \nSCT normalization")+theme(title =element_text(size=16, face='bold'),
                                                                                                                                 legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                                                                                axis.title=element_text(size=14,face="bold"))



################################## 
## Differential gene expression ##
##################################

path="Directory"
setwd(path)
load("Directory/SCT_Inputs_step2-DataLoading.RData")

library(future)
plan()
# change the current plan to access parallelization
plan("multiprocess", workers = 8)
plan()

options(future.globals.maxSize = 8000 * 1024^2)

scMalvsNorm.SCT <- FindMarkers(combined.sct, group.by = "Type", assay = "SCT", ident.1 = "Malignant", ident.2 = "Normal",
                               verbose = FALSE,logfc.threshold =FALSE,densify = TRUE)

scTNBCvsNorm.SCT <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "TNBC", ident.2 = "Normal",
                                verbose = FALSE, group.by = "subtype",logfc.threshold =FALSE)

scHER2vsNorm.SCT <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "HER2+", ident.2 = "Normal",
                                verbose = FALSE, group.by = "subtype",logfc.threshold =FALSE,densify = TRUE)

scERvsNorm.SCT <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "ER+", ident.2 = "Normal",
                              verbose = FALSE, group.by = "subtype",logfc.threshold =FALSE,densify = TRUE)


path="Directory"
setwd(path)
save(scMalvsNorm.SCT,scTNBCvsNorm.SCT,scHER2vsNorm.SCT, scERvsNorm.SCT, file = "scDEGs.SCT.Rdata")
head(scMalvsNorm, n = 15)


## Investigating the status of differentialy expressed genes un LogNormalization method in SCTransform normalized data
ACPs <- c("BST2",
          "GAPDH",
          "H2AFJ",
          "SCGB2A1",
          "HMGN2",
          "S100A8",
          "S100A7",
          "S100A9",
          "ZG16B",
          "HMGB1"
)

load("Directory/scDEGs.SCT.Rdata")

## DEGs: Ploting the volcanoplot for comparing Different BC subtypes vs normal epithelial cells

p1 <- EnhancedVolcano(scMalvsNorm.SCT,
                      lab = rownames(scMalvsNorm.SCT),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      selectLab = c(rownames(subset(na.omit(scMalvsNorm.SCT[ACPs,]), avg_log2FC>0.5 & p_val_adj<0.05))),
                      title = 'Malignant vs Normal',
                      subtitle = "Differential expression",
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 1,
                      pointSize = 3.0,
                      labSize = 3.5,
                      legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                     'p-value & Log (base 2) FC'),
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      colAlpha = 2/5,
                      legendPosition = 'right',
                      legendLabSize = 14,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      lengthConnectors = unit(0.02, "npc"),
                      col=c('lightgrey', 'lightgrey', 'lightgrey', 'red3',"green"),
                      max.overlaps = Inf,
                      maxoverlapsConnectors = Inf,
                      typeConnectors = "closed",
                      xlim = c(-2.5,2.5),
                      caption = NULL
)+theme(legend.position = "none")



ACPs.up.MalvsNorm <- rownames(subset(na.omit(scMalvsNorm.SCT[ACPs,]), avg_log2FC>0.5 & p_val_adj<0.05))


p2 <- EnhancedVolcano(scTNBCvsNorm.SCT,
                      lab = rownames(scTNBCvsNorm.SCT),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      selectLab = ACPs.up.MalvsNorm,
                      title = 'TNBC vs Normal',
                      subtitle = "Differential expression",
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 0,
                      pointSize = 3.0,
                      labSize = 3.5,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      colAlpha = 2/5,
                      legendPosition = NULL,
                      legendLabSize = NULL,
                      legendIconSize = NULL,
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      lengthConnectors = unit(0.02, "npc"),
                      col=c('lightgrey', 'lightgrey', 'lightgrey', 'lightgrey',"lightgrey"),
                      max.overlaps = Inf,
                      maxoverlapsConnectors = Inf,
                      typeConnectors = "closed",
                      xlim = c(-2.5,2.5),
                      caption = NULL
)+theme(legend.position = "none")

p3 <- EnhancedVolcano(scHER2vsNorm.SCT,
                      lab = rownames(scHER2vsNorm.SCT),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      selectLab = ACPs.up.MalvsNorm,
                      title = 'HER2 vs Normal',
                      subtitle = "Differential expression",
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 0,
                      pointSize = 3.0,
                      labSize = 3.5,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      colAlpha = 2/5,
                      legendPosition = NULL,
                      legendLabSize = NULL,
                      legendIconSize = NULL,
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      lengthConnectors = unit(0.02, "npc"),
                      col=c('lightgrey', 'lightgrey', 'lightgrey', 'lightgrey',"lightgrey"),
                      max.overlaps = Inf,
                      maxoverlapsConnectors = Inf,
                      typeConnectors = "closed",
                      xlim = c(-2.5,2.5),
                      caption = NULL
)+theme(legend.position = "none")

p4 <- EnhancedVolcano(scERvsNorm.SCT,
                      lab = rownames(scERvsNorm.SCT),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      selectLab = ACPs.up.MalvsNorm,
                      title = 'ER vs Normal',
                      subtitle = "Differential expression",
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 0,
                      pointSize = 3.0,
                      labSize = 3.5,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      colAlpha = 2/5,
                      legendPosition = NULL,
                      legendLabSize = NULL,
                      legendIconSize = NULL,
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      lengthConnectors = unit(0.02, "npc"),
                      col=c('lightgrey', 'lightgrey', 'lightgrey', 'lightgrey',"lightgrey"),
                      max.overlaps = Inf,
                      maxoverlapsConnectors = Inf,
                      typeConnectors = "closed",
                      xlim = c(-2.5,2.5),
                      caption = NULL
)+theme(legend.position = "none")



pdf("SCT.DEGs.pdf", width=10, height=8)
par(mar = c(6, 6, 6, 6));
ggarrange(p1, p2, p3,p4,
          ncol = 2, nrow = 2)
dev.off()

rownames(subset(na.omit(scMalvsNorm.SCT[ACPs$ACPs,]), avg_log2FC>0.3 & p_val_adj<0.05))
rownames(subset(na.omit(scTNBCvsNorm.SCT[ACPs$ACPs,]), avg_log2FC>0 & p_val_adj<0.05))
rownames(subset(na.omit(scHER2vsNorm.SCT[ACPs$ACPs,]), avg_log2FC>0 & p_val_adj<0.05))
rownames(subset(na.omit(scERvsNorm.SCT[ACPs$ACPs,]), avg_log2FC>0 & p_val_adj<0.05))
