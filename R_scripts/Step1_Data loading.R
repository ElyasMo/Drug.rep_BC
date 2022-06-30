library(Seurat)
library(ggplot2)
library(SingleR)
library(WGCNA)
library(DESeq2)
library(EnhancedVolcano)
library(stringr)
library(slinky)
library(ggthemr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)
library(Hmisc)
library(tidyverse)
library(ggpubr)
library(ggfortify)




# Step1. Data loading
## Loading and clustering \nthe data

#############################
## Reading the Normal data ##
#############################
Directory1 <- "Directory1"
Directory2 <- "Directory2"
Directory3 <- "Directory3"

# Reading the Normal data
Normal <- readRDS(paste0(Directory1,"Normal.rds"))

# Assigning the percentage of expression of mitochondrial geneses
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

# Normalization
Normal <- NormalizeData(Normal)

#Further pre-processing steps
Normal <- FindVariableFeatures(Normal, selection.method = "vst", nfeatures = 8000)
Normal <- ScaleData(Normal, verbose = FALSE)

#Clustering
Normal <- RunPCA(Normal, verbose = FALSE)
Normal <- RunUMAP(Normal, dims = 1:30, verbose = FALSE)
Normal <- FindNeighbors(Normal, dims = 1:30, verbose = FALSE)
Normal <- FindClusters(Normal, verbose = FALSE,resolution = 0.1)

#Choosing the proper resolution
Normal <- SetIdent(Normal, value = Normal@meta.data$RNA_snn_res.0.1)

# Adding the normal status to the metadata 
Normal$Type <- "Normal"


#########################
## Reading the BC data ##
#########################

# Reading the data
BC <- readRDS(paste0(Directory2,"BC7_new.rds"))

# Assigning the percentage of expression of mitochondrial geneses
BC <- lapply(BC, function(X){
  X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = "^MT-");X
}) 

# Merging the datasets from patients which the Malignanat and Nonmalignant cells could be distinguished in their Epithelial cells
BC <-merge(x = BC$CID3586, y = c(BC$CID3838,BC$CID3921,BC$CID3941,BC$CID3946,BC$CID3948,BC$CID3963,BC$CID4040, 
                                 BC$CID4066,BC$CID4067,BC$CID4290A,BC$CID4398,BC$CID44041,BC$CID4461,BC$CID4463, 
                                 BC$CID4465,BC$CID4471,BC$CID4495,BC$CID44971,BC$CID44991,BC$CID4513,BC$CID4515, 
                                 BC$CID45171,BC$CID4523,BC$CID4530N,BC$CID4535 ))


# Data normalization
BC <- NormalizeData(BC)
##Finidng variable features
BC <- FindVariableFeatures(BC, selection.method = "vst",nfeatures = 10000)
##Scaling the data and regressing out the MT genes,
BC<- ScaleData(BC, assay = "RNA", verbose = FALSE)

# Required steps forclustering
BC <- RunPCA(BC, verbose = FALSE)
BC <- RunUMAP(BC, dims = 1:30, verbose = FALSE)



# Extracting the Malignant epitelial cells
Malignant <- SplitObject(BC, split.by = "celltype_major")
Malignant <- Malignant$`Cancer Epithelial`


# # Required steps for clustering
Malignant <- RunPCA(Malignant, verbose = FALSE)
Malignant <- RunUMAP(Malignant, dims = 1:30, verbose = FALSE)
Malignant <- FindNeighbors(Malignant, dims = 1:30, verbose = FALSE)
Malignant <- FindClusters(Malignant,resolution = 0.1, verbose = FALSE)

# Adding the Malignant status to the metadata 
Malignant$Type <- "Malignant"


#################################
## Reading the Proteomics data ##
#################################

# Reading the data
protein_normalized <- read_csv(paste0(Directory3,"protein_normalized.csv"))
protein_normalized <- as.data.frame(protein_normalized[!duplicated(protein_normalized$Gene_Symbol),])

#Removing the outliers after hierarchical clustering the proteomics profiles
protein_normalized <- protein_normalized[complete.cases(protein_normalized$Gene_Symbol), ]
rownames(protein_normalized) <- protein_normalized$Gene_Symbol
protein.table <- protein_normalized[,-c(1:6, 20:22)]
protein.table <- protein.table[,-c(14,17,24,26)]
protein.table[is.na(protein.table)] = 0

# Preprocessing
protein_seurat <- CreateSeuratObject(protein.table, min.cells = 3, min.genes = 200, project = "Protein" )
protein_seurat@assays$RNA@data <- protein_seurat@assays$RNA@counts
protein_seurat <- FindVariableFeatures(protein_seurat, nfeatures = 6000)
protein_seurat <- ScaleData(protein_seurat, verbose = FALSE)
protein_seurat <- RunPCA(protein_seurat, verbose = FALSE,npcs = 15,approx=FALSE)
protein_seurat <- RunUMAP(protein_seurat, dims = 1:15, verbose = FALSE,n.neighbors =24 )
protein_seurat$subtype <- c(rep("HER2",6), rep("ER",7), rep("TNBC",11))
protein_seurat <- SetIdent(protein_seurat, value = protein_seurat$subtype)


## UMAPs

DimPlot(protein_seurat, raster = FALSE)+ggtitle("Proteomics data")+theme(title =element_text(size=16, face='bold'),
                                                                         legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                         axis.title=element_text(size=14,face="bold"))

p1 <-DimPlot(BC, group.by = "celltype_major", split.by = "subtype", raster = FALSE)+ggtitle("Cells from Malignant breast tissue")+theme(title =element_text(size=16, face='bold'),
                                                                                                                                        legend.text = element_text(size=15),  legend.position = "top",                                                                                       axis.text=element_text(size=14),
                                                                                                                                        axis.title=element_text(size=14,face="bold"))

DimPlot(Normal, raster = FALSE)+ggtitle("Normal cells from breast tissue")+theme(title =element_text(size=16, face='bold'),
                                                                                 legend.text = element_text(size=15),                                                                                         axis.text=element_text(size=14),
                                                                                 axis.title=element_text(size=14,face="bold"))



path="path"
setwd(path)
save(BC,Normal,Malignant,protein_seurat, file = "Inputs_step1.RData" ) 

