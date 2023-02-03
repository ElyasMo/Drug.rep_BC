---
title: "Visium.data.analysis.Rmd"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#Spatial

```{r,cache.path=}
library(Seurat)

st.mtx <- function(path,name){
 x1 <- Read10X(paste0(path,"filtered_feature_bc_matrix/"),gene.column = 1)
x2 <- Read10X_Image(paste0(path,"spatial/"))
x1 <- CreateSeuratObject(counts = x1,assay = "Spatial",project = name)
return(x1)
}


pre_process <- function(x){
x <- NormalizeData(x,verbose = F)
x <- FindVariableFeatures(x,nfeatures = 6000,verbose = F)
 x<- ScaleData(x,verbose = F)
x <- RunPCA(x,npcs = 30,verbose = F)
x <- RunUMAP(x,dims = 1:30,verbose = F)
x <- FindNeighbors(x, dims = 1:30, verbose = FALSE)
x <- FindClusters(x, verbose = FALSE,resolution = c(0.05,0.08,0.1,0.2))
return(x)
}


dir1 <- "./Visium/CID4290/" 
dir2 <- "./Visium/CID4465/" 
dir3 <- "./Visium/CID44971/" 
dir4 <- "./Visium/CID4535/" 

CID4290 <- st.mtx(dir1,"CID4290")
CID4465 <- st.mtx(dir2,"CID4465")
CID44971 <- st.mtx(dir3,"CID44971")
CID4535 <- st.mtx(dir4,"CID4535")

CID4290@images$CID4290 <- Read10X_Image(paste0(dir1,"spatial/"),filter.matrix = T)
CID4290@images$CID4290@assay <- "Spatial"
CID4290@images$CID4290@key <- "cid4290_"

CID4465@images$CID4465 <- Read10X_Image(paste0(dir2,"spatial/"))
CID4465@images$CID4465@assay <- "Spatial"
CID4465@images$CID4465@key <- "cid4465_"

CID44971@images$CID44971 <- Read10X_Image(paste0(dir3,"spatial/"))
CID44971@images$CID44971@assay <- "Spatial"
CID44971@images$CID44971@key <- "cid44971_"

CID4535@images$CID4535 <- Read10X_Image(paste0(dir4,"spatial/"))
CID4535@images$CID4535@assay <- "Spatial"
CID4535@images$CID4535@key <- "cid4535_"


CID4290@images$CID4290@coordinates <- CID4290@images$CID4290@coordinates[na.omit(match(names(as.data.frame(CID4290@assays$Spatial@counts)),row.names(CID4290@images$CID4290@coordinates))),]
CID4465@images$CID4465@coordinates <- CID4465@images$CID4465@coordinates[na.omit(match(names(as.data.frame(CID4465@assays$Spatial@counts)),row.names(CID4465@images$CID4465@coordinates))),]
CID44971@images$CID44971@coordinates <- CID44971@images$CID44971@coordinates[na.omit(match(names(as.data.frame(CID44971@assays$Spatial@counts)),row.names(CID44971@images$CID44971@coordinates))),]
CID4535@images$CID4535@coordinates <- CID4535@images$CID4535@coordinates[na.omit(match(names(as.data.frame(CID4535@assays$Spatial@counts)),row.names(CID4535@images$CID4535@coordinates))),]

objects <- Reduce(merge,c(CID4290,CID4465,CID44971,CID4535))
 objects <- pre_process(objects)
 
saveRDS(objects,"./Visium/Objects.rds")
```

```{r}
objects <- readRDS("./Visium/Objects.rds")
 
DimPlot(objects)+ggtitle("Distinguishment of  Visium profiles in patients")
```

```{r}
# Loading the data from step 1
path="./New/"
setwd(path)
load(file = "Inputs_step1.RData" )

BC <- SplitObject(BC, split.by = "Patient")
BC <- Reduce(merge, c(BC$CID4290A, BC$CID4465,BC$CID44971,BC$CID4535))

save(BC,St.Objects,file = "./Visium/Visium.set.Rdata")
```

```{r}
load(file = "./Visium/Visium.set.Rdata")

library(dplyr)

Label.transfer.sct <- function(sc,ST){
sc <- SCTransform(sc, ncells = 3000, verbose = FALSE,assay = "RNA") %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
# After subsetting, we renormalize cortex
ST <- SCTransform(ST, assay = "Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)

anchors <- FindTransferAnchors(reference = sc, query = ST, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = sc$celltype_major, prediction.assay = FALSE,
    weight.reduction = ST[["pca"]], dims = 1:30)
SP1.predictions.assay1 <- TransferData(anchorset = anchors, refdata = sc$celltype_major, prediction.assay = TRUE,
    weight.reduction = ST[["pca"]], dims = 1:30)

ST[["prediction.labels"]] <- SP1.predictions.assay1
ST[["predictions"]] <- predictions.assay$predicted.id
return(ST)
}
```

```{r,fig.height=10}

sc.ST.int.sct <- Label.transfer.sct(BC,St.Objects)
sc.ST.int.sct$predictions <- factor(sc.ST.int.sct$predictions,levels = c("Cancer Epithelial","Normal Epithelial","Endothelial", "CAFs", "PVL",   "Plasmablasts", "Myeloid", "T-cells","B-cells" ))

pdf("./Visium1.pdf", width=5, height=12)
par(mar = c(6, 6, 6, 6));
SpatialDimPlot(sc.ST.int.sct,group.by = "predictions",cols = c("Cancer Epithelial"="#E76BF3","Normal Epithelial"="#F8766D","Endothelial"="#00BF7D", "CAFs"="#00B0F6", "PVL"="#A3A500",   "Plasmablasts"="#D0DF07", "Myeloid"="#05E1E1", "T-cells"="#E5E591","B-cells"= "#AB32B4"),ncol = 1,pt.size.factor = 1.8)&guides(fill=guide_legend(title="Cell type"))&theme(legend.text = element_text(size = 12),legend.title = element_text(size=14))
dev.off()
  
# SpatialDimPlot(sc.ST.int,group.by = "predictions",images = c("CID4290","CID4465","CID44971","CID4535"),ncol = 2)

saveRDS(sc.ST.int.sct,"./Visium/ST.list.final.rds")

```

```{r,fig.width=20}
sc.ST.int.sct <- readRDS("./Visium/ST.list.final.rds")
gene.list <- c("BST2","GAPDH","H2AFJ","S100A8","S100A9","SCGB2A1","HMGN2")

sc.ST.int.sct <- SetIdent(sc.ST.int.sct,value = "prediction")

sc.ST.int.sct$predictions <- factor(sc.ST.int.sct$predictions,levels = c("Cancer Epithelial","Normal Epithelial","Endothelial", "CAFs", "PVL",   "Plasmablasts", "Myeloid", "T-cells","B-cells" ))

sc.ST.int.sct <- AddModuleScore(sc.ST.int.sct,features = c("BST2","GAPDH","H2AFJ","S100A8","S100A9","SCGB2A1","HMGN2"),name = c("BST2","GAPDH","H2AFJ","S100A8","S100A9","SCGB2A1","HMGN2"),assay = "SCT")


P1 <- VlnPlot(sc.ST.int.sct,features = names(sc.ST.int.sct@meta.data)[13],group.by = "predictions",pt.size = 0,cols = c("#E76BF3","#F8766D","#00BF7D","#00B0F6","#A3A500","#D0DF07","#05E1E1","#E5E591","#AB32B4"))+ylim(0,1.5)+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+ylab("Expression scores")+theme(legend.position = "none")+rotate_x_text(55)+xlab("")+theme(axis.text.x = element_text(size=0))

P2 <- VlnPlot(sc.ST.int.sct,features = names(sc.ST.int.sct@meta.data)[14],group.by = "predictions",pt.size = 0,cols = c("#E76BF3","#F8766D","#00BF7D","#00B0F6","#A3A500","#D0DF07","#05E1E1","#E5E591","#AB32B4"))+ylim(0,3)+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+ylab("Expression scores")+theme(legend.position = "none")+rotate_x_text(55)+xlab("")+theme(axis.text.x = element_text(size=0))

P3 <- VlnPlot(sc.ST.int.sct,features = names(sc.ST.int.sct@meta.data)[15],group.by = "predictions",pt.size = 0,cols = c("#E76BF3","#F8766D","#00BF7D","#00B0F6","#A3A500","#D0DF07","#05E1E1","#E5E591","#AB32B4"))+ylim(0,2.5)+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+ylab("Expression scores")+theme(legend.position = "none")+rotate_x_text(55)+xlab("")+theme(axis.text.x = element_text(size=0))

P4 <- VlnPlot(sc.ST.int.sct,features = names(sc.ST.int.sct@meta.data)[16],group.by = "predictions",pt.size = 0,cols = c("#E76BF3","#F8766D","#00BF7D","#00B0F6","#A3A500","#D0DF07","#05E1E1","#E5E591","#AB32B4"))+ylim(0,3.5)+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+ylab("Expression scores")+theme(legend.position = "none")+rotate_x_text(55)+xlab("")+theme(axis.text.x = element_text(size=0))

P5 <- VlnPlot(sc.ST.int.sct,features = names(sc.ST.int.sct@meta.data)[17],group.by = "predictions",pt.size = 0,cols = c("Cancer Epithelial"="#E76BF3","Normal Epithelial"="#F8766D","Endothelial"="#00BF7D", "CAFs"="#00B0F6", "PVL"="#A3A500",   "Plasmablasts"="#D0DF07", "Myeloid"="#05E1E1", "T-cells"="#E5E591","B-cells"= "#AB32B4"))+ylim(0,3.5)+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+ylab("Expression scores")+theme(legend.position = "none")+rotate_x_text(55)+xlab("")+theme(axis.text.x = element_text(size=0))

P6 <- VlnPlot(sc.ST.int.sct,features = names(sc.ST.int.sct@meta.data)[18],group.by = "predictions",pt.size = 0,cols = c("Cancer Epithelial"="#E76BF3","Normal Epithelial"="#F8766D","Endothelial"="#00BF7D", "CAFs"="#00B0F6", "PVL"="#A3A500",   "Plasmablasts"="#D0DF07", "Myeloid"="#05E1E1", "T-cells"="#E5E591","B-cells"= "#AB32B4"))+ylim(0,3.5)+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+ylab("Expression scores")+theme(legend.position = "none")+rotate_x_text(55)+xlab("")+theme(axis.text.x = element_text(size=0))

P7 <- VlnPlot(sc.ST.int.sct,features = names(sc.ST.int.sct@meta.data)[19],group.by = "predictions",pt.size = 0,cols = c("#E76BF3","#F8766D","#00BF7D","#00B0F6","#A3A500","#D0DF07","#05E1E1","#E5E591","#AB32B4"))+ylim(0,1.5)+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+ylab("Expression scores")+theme(legend.position = "none")+rotate_x_text(55)+xlab("")+theme(axis.text.x = element_text(size=0))

P8 <- VlnPlot(sc.ST.int.sct,features = names(sc.ST.int.sct@meta.data)[19],group.by = "predictions",pt.size = 0,cols = c("#E76BF3","#F8766D","#00BF7D","#00B0F6","#A3A500","#D0DF07","#05E1E1","#E5E591","#AB32B4"))+ylim(0,1.5)+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))+rotate_x_text(55)+xlab("")+theme(axis.text.x = element_text(size=0),
                                                                                                           legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=16))


library(ggpubr)
leg <- get_legend(P8)
P8 <- as_ggplot(leg)

pdf("../New/PIC/vln.pdf", width=6.5, height=10.5)
par(mar = c(6, 6, 6, 6));
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,ncol = 2,nrow = 4)
dev.off()

wilcox.test(as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Cancer Epithelial"))$BST21,
            as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Normal Epithelial"))$BST21)

wilcox.test(as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Cancer Epithelial"))$GAPDH2,
            as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Normal Epithelial"))$GAPDH2)

wilcox.test(as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Cancer Epithelial"))$H2AFJ3,
            as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Normal Epithelial"))$H2AFJ3)

wilcox.test(as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Cancer Epithelial"))$S100A84,
            as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Normal Epithelial"))$S100A84)

wilcox.test(as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Cancer Epithelial"))$S100A95,
            as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Normal Epithelial"))$S100A95)

wilcox.test(as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Cancer Epithelial"))$SCGB2A16,
            as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Normal Epithelial"))$SCGB2A16)

wilcox.test(as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Cancer Epithelial"))$HMGN27,
            as.data.frame(subset(sc.ST.int.sct@meta.data,predictions=="Normal Epithelial"))$HMGN27)
```

```{r}

pre_process <- function(ST){
ST <- FindVariableFeatures(ST,nfeatures=4000)
ST <- ScaleData(ST)
ST <- RunPCA(ST,npcs = 30)
ST <- RunUMAP(ST,dims = 1:30)
return(ST)
}


DimPlot(sc.ST.int.sct,group.by = "predictions",split.by = "orig.ident")

ST3.DEGS <- FindMarkers(sc.ST.int.sct,ident.1 = "Cancer Epithelial",ident.2 = "Normal Epithelial",group.by = "predictions",logfc.threshold = 0)

saveRDS(ST3.DEGS,"./Visium/DEGs.rds")

library(EnhancedVolcano)

p1 <- EnhancedVolcano(ST3.DEGS,
    lab = rownames(ST3.DEGS),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    xlab = bquote(~Log[2]~ 'fold change'),
      selectLab = c("BST2","GAPDH","H2AFJ","S100A8","S100A9","SCGB2A1","HMGN2"),
    title = 'Malignant vs Normal epithelial cells',
    subtitle = "Differential expression",
    pCutoff = 0.05,
    FCcutoff = 0.3,
        cutoffLineType = 'twodash',
    cutoffLineWidth = 1,
    pointSize = 3.0,
    labSize = 4,
    legendLabels=c('Not sig.','Log (base 2) FC','p-value',
      'p-value & Log (base 2) FC'),
        labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 1/5,
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
    xlim = c(-4,4),
    caption = NULL
)+theme(legend.position = "none")
```
