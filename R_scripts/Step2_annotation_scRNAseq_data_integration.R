
# Step2 Annotation_Data integration
## Normal cell annotation

############################
## Normal cell annotation ##
############################


# Loading the data from step 1
path="path"
setwd(path)
load(file = "Inputs_step1.RData" )


#Changing the Seurat object to sSinglecell object
Normal.sce <- as.SingleCellExperiment(Normal,assay = "RNA")

#Determining the path to the annotation libraries
wpath="wpath"

#Load references
blue <- readRDS(paste(wpath,"singleR.BlueprintEncodeData.rds",sep="/"))
hpca <- readRDS(paste(wpath,"singleR.HumanPrimaryCellAtlasData.rds",sep="/"))


#Combining to libraries for a better annotation
ref <- list(BP=blue,HPCA=hpca)

# Disinguishing the main and fine libraries (main has the major annotations and fine includes the cell subtypes as well)
main.labels <- list(blue$label.main, hpca$label.main)
fine.labels <-  list(blue$label.fine, hpca$label.fine)

# Annotationg the cells. It is time consuming so I saved it in .rds files for the sake of time and computational intensity.
# pred.Normal.main <- SingleR(test=Normal.sce, ref=ref, assay.type.test=1, labels=main.labels)
# pred.Normal.fine <- SingleR(test=Normal.sce, ref=ref, assay.type.test=1, labels=fine.labels)
# saveRDS(pred.Normal.main,"pred.Normal.main.rds")
# saveRDS(pred.Normal.fine,"pred.Normal.fine.rds")

refpath="refpath"

#Loading the annotations
pred.Normal.main <- readRDS(paste0(refpath,"pred.Normal.main.rds"))
pred.Normal.fine <- readRDS(paste0(refpath,"pred.Normal.fine.rds"))

#Adding the annotations to the objects 
Normal.celllabels2 <- cbind(celllabelsb=pred.Normal.main$labels)
Normal.celllabels1 <- cbind(celllabels=pred.Normal.fine$labels)
rownames(Normal.celllabels1) <- rownames(Normal.celllabels2) <- colnames(Normal)
Normal <- AddMetaData(Normal, as.data.frame(cbind(Normal.celllabels1,Normal.celllabels2)))

## Evaluate annotation based on scores

plotScoreHeatmap(pred.Normal.main ,max.labels=10, fontsize=10,cluster_cols = FALSE, scores.use = 0,calls.use = 0)

plotScoreHeatmap(pred.Normal.main ,max.labels=10, fontsize=12,cluster_cols = FALSE, scores.use = 1,calls.use = 1)

plotScoreHeatmap(pred.Normal.main ,max.labels=10, fontsize=12,cluster_cols = FALSE, scores.use = 2,calls.use = 2)


tab <- table(cluster=Normal$RNA_snn_res.0.1, label=pred.Normal.main$labels) 
pheatmap::pheatmap(log10(tab+20),fontsize = 12) # using a larger pseudo-count for smoothing. 

# par(mfrow=c(1,1), mar=c(2.5,7,1,.5))
# graphics::boxplot(pred.Normal.fine$scores, xlab="prediction.score", ylab="", main="Prediction scores for cell labels", boxwex=.8,
#                   cex=.5, cex.main=.8, cex.lab=.8, cex.axis=.6, las=1, horizontal=T, varwidth=T, col=rainbow(36))

#Tables of cell annotations
tab.main <- table(Normal$celllabelsb, Normal$RNA_snn_res.0.1)
tab.fine <- table(Normal$celllabels, Normal$RNA_snn_res.0.1)

knitr::kable(tab.main, caption = "Normal-main")
knitr::kable(tab.fine, caption = "Normal-fine")

library(bluster)
pairwiseRand(Normal.sce$RNA_snn_res.0.1, pred.Normal.main$labels, mode="index")


plotDeltaDistribution(pred.Normal.main,ncol = 8)



## label the annotations


###########################
## label the annotations ##
###########################

#Adding new cell labels
cluster.ids.Normal <- c("Epithelial cells","Fibroblasts", "Keratinocytes",  "Endothelial cells","Epithelial cells","Adipocytes\n-Tissue Stem cells", "CD4+/CD8+ T cells", "Monocytes", "Endothelial cells","Epithelial cells", "B cells")

names(cluster.ids.Normal) <- levels(Normal)
Normal <- RenameIdents(Normal, cluster.ids.Normal)


# Extracting the Epithelial cells to be compared with the epithelial malignant cells
Normal$annotation <- Normal@active.ident
Normal$subtype <- "Normal"
Normal.Epithelial <- SplitObject(Normal, split.by = "annotation")
Normal.Epithelial <- Normal.Epithelial$`Epithelial cells`
Normal.Epithelial <- RunPCA(Normal.Epithelial, verbose = FALSE,npcs = 15)
Normal.Epithelial <- RunUMAP(Normal.Epithelial, dims = 1:15, verbose = FALSE,n.neighbors =24 )
Normal.Epithelial$Type <- "Normal"




#Ploting to see the new annotated cells
DimPlot(Normal, raster = FALSE)+ggtitle("Normal cells from breast tissue")+theme(title =element_text(size=16, face='bold'),
                                                                                 legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                                 axis.title=element_text(size=14,face="bold"))


DimPlot(Normal.Epithelial, raster = FALSE)+ggtitle("Normal epithelial cells from breast tissue")+theme(title =element_text(size=16, face='bold'),
                                                                                                       legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                                                       axis.title=element_text(size=14,face="bold"))


## scRNAseq Data integration

###############################
## scRNAseq Data integration ##
###############################

Malignant$Type <- "Malignant"
Normal.Epithelial$Type <- "Normal"

all_objects <- merge(Malignant, Normal.Epithelial)
all_objects <- SplitObject(all_objects, split.by = "Type")
# normalize and identify variable features for each dataset independently
all_objects <- lapply(X = all_objects, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst")
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = all_objects, nfeatures = 10000)

anchors <- FindIntegrationAnchors(object.list = all_objects, anchor.features = features)
# this command creates an 'integrated' data assay
all.int <- IntegrateData(anchorset = anchors)


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(all.int) <- "integrated"

# Run the standard workflow for visualization and clustering
all.int <- ScaleData(all.int, verbose = FALSE)
all.int <- RunPCA(all.int, npcs = 30, verbose = FALSE)
all.int <- RunUMAP(all.int, reduction = "pca", dims = 1:30)
all.int <- FindNeighbors(all.int, reduction = "pca", dims = 1:30)
all.int <- FindClusters(all.int, resolution = 0.5)

#Setting the sybtypes as default identity
all.int <- SetIdent(all.int, value = all.int$subtype)



DimPlot(all.int, reduction = "umap", group.by = "Type")
DimPlot(all.int, reduction = "umap", label = FALSE, repel = TRUE)+ggtitle("Malignant and normal \nintegrated epithelial cells")+theme(legend.position = "top")


path="path"
setwd(path)
save(Normal.Epithelial,Malignant, all.int, file = "Step2 Annotation_Data integration.RData" )
