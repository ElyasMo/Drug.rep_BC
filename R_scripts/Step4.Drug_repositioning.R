

#Step4. Drug repositioning
# Checking the available BC cell lines in LINCS L1000 database

Cells <- data.frame(cells= c("BT20", "HS578T","MDAMB231", "SKBR3", "MCF7"),
                    type=c("TNBC", "TNBC", "TNBC", "HER2", "ER"))

# Loading the metadata
gene_info <- read.delim("Directory/GSE92742_Broad_LINCS_gene_info.txt", header = TRUE, sep = "\t", dec = ".")

meta_info <- read.delim("Directory/GSE92742_Broad_LINCS_inst_info.txt", header = TRUE, sep = "\t", dec = ".")

meta_info.cell <- read.delim("Directory/GSE92742_Broad_LINCS_cell_info-2.txt", header = TRUE, sep = "\t", dec = ".")

meta_info.cell <- subset(meta_info.cell, primary_site=="breast")

#Reading in the data
setwd("Directory")
getwd()

key <- "144f828115e01872dbcf07a048c0032e"

gctx <- "Directory/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"

info <- "Directory/GSE92742_Broad_LINCS_inst_info.txt"

sl <- Slinky(key, gctx, info)


# Extracting the required information from LINCS L1000 database for different BC cell lines

col.ix <- base::which(slinky::metadata(sl)$cell_id == Cells$cells[1]  & slinky::metadata(sl)$pert_dose=="10"  & slinky::metadata(sl)$pert_time == "24" 
                      & slinky::metadata(sl)$pert_iname != "DMSO" )

data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
BT20_pert <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")
rownames(BT20_pert) <- BT20_pert$pr_gene_symbol
names(BT20_pert[, grep("LJP", names(BT20_pert))])


col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[1] & slinky::metadata(sl)$pert_iname == "DMSO" & slinky::metadata(sl)$pert_time == "24")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
BT20_DMSO <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")
rownames(BT20_DMSO) <- BT20_DMSO$pr_gene_symbol
names(BT20_DMSO[, grep("LJP", names(BT20_DMSO))])


col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[2] & slinky::metadata(sl)$pert_dose=="10"  & slinky::metadata(sl)$pert_time == "24" 
                      & slinky::metadata(sl)$pert_iname != "DMSO")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
HS578T_pert <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")

col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[2] & slinky::metadata(sl)$pert_time == "24" & slinky::metadata(sl)$pert_iname == "DMSO")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
HS578T_DMSO <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")


col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[3] & slinky::metadata(sl)$pert_dose=="10"  & slinky::metadata(sl)$pert_time == "24" 
                      & slinky::metadata(sl)$pert_iname != "DMSO")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
MDAM_pert <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")

col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[3] & slinky::metadata(sl)$pert_time == "24" & slinky::metadata(sl)$pert_iname == "DMSO")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
MDAM_DMSO <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")


col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[4] & slinky::metadata(sl)$pert_dose=="10"  & slinky::metadata(sl)$pert_time == "24" 
                      & slinky::metadata(sl)$pert_iname != "DMSO")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
SKBR3_pert <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")

col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[4] & slinky::metadata(sl)$pert_time == "24" & slinky::metadata(sl)$pert_iname == "DMSO")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
SKBR3_DMSO <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")


col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[5] & slinky::metadata(sl)$pert_dose=="10"  & slinky::metadata(sl)$pert_time == "24" 
                      & slinky::metadata(sl)$pert_iname != "DMSO")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
MCF7_pert <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")
rownames(MCF7_pert) <- MCF7_pert$pr_gene_symbol
MCF7_pert <- MCF7_pert[,-c(1:2)]

col.ix <- base::which(slinky::metadata(sl)$cell_id ==Cells$cells[5] & slinky::metadata(sl)$pert_time == "24" & slinky::metadata(sl)$pert_iname == "DMSO")
data <- as.data.frame(readGCTX(sl[, col.ix]))
data$col <- rownames(data)
MCF7_DMSO <- merge(gene_info[,1:2], data, by.x="pr_gene_id", by.y= "col")
rownames(MCF7_DMSO) <- MCF7_DMSO$pr_gene_symbol
MCF7_DMSO <- MCF7_DMSO[,-c(1:2)]

perts <- data.frame(Celltypes=Cells$cells, Pert.numbers= c(length(BT20_pert), length(HS578T_pert),length(MDAM_pert),length(SKBR3_pert),length(MCF7_pert)  ))

# Investigating the number of perturbagens for each BC cell line in LINCS L1000 database
ggplot(data=perts, aes(x=Celltypes, y=Pert.numbers)) +
  geom_bar(stat="identity", position=position_dodge(),fill = c("#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352"))
+theme(title =element_text(size=16, face='bold', color = "black"),
     legend.text = element_text(size=15, colour = "black"), 
 axis.text=element_text(size=14, colour = "black"))+ rotate_x_text(55)+
  ggtitle("Number of available perturbagens \nfor each cancer cell line")+theme_light()+rotate_x_text(45)

# saveRDS(MCF7_DMSO, "../PD/New/MCF7_DMSO.rds")
# saveRDS(MCF7_pert, "../PD/New/MCF7_pert.rds")
## Reading the saved control and perturbed expression profiles for MCF7 BC cell line
MCF7_DMSO <- readRDS("../New/MCF7_DMSO.rds")
MCF7_pert <- readRDS("../New/MCF7_pert.rds")

# Unifying the metadata
Mean_ctl <- data.frame(matrix(ncol = 0, nrow = 12328))
for (i in 1:length(unique(str_extract(names(MCF7_DMSO), "(.*?_.*?_.*?)_..")))) {
  Mean_ctl <- cbind(Mean_ctl,as.data.frame(rowMeans(MCF7_DMSO[, grep(unique(str_extract(names(MCF7_DMSO), "(.*?_.*?_.*?)_.."))[i], names(MCF7_DMSO))])))
}

colnames(Mean_ctl) <- unique(str_extract(names(MCF7_DMSO), "(.*?_.*?_.*?)_.."))

# Calculating the LogFoldChange for MCF7 BC cell line
Mean_ctl <- log2(Mean_ctl+1)
MCF7_pert <- log2(MCF7_pert+1)

patterns <- intersect(unique(str_extract(names(MCF7_DMSO), "(.*?_.*?_.*?)_..")), unique(str_extract(names(MCF7_pert), "(.*?_.*?_.*?)_..")))

MCF7_pert.sub <- data.frame(matrix(ncol = 0, nrow = 12328))
for (i in 1:length(patterns)){
  MCF7_pert.sub <- cbind(MCF7_pert.sub, as.data.frame(MCF7_pert[grep(patterns[i],names(MCF7_pert))], drop=FALSE))
}

subs.df <- function(df1, df2) {
  x <- data.frame(matrix(ncol = 0, nrow = length(rownames(df1))))
  for (i in 1:length(df1)) {
    x <- cbind(x, as.data.frame(df1[,i])-as.data.frame(df2))
  }
  return(x)
}


LFC_MCF7 <- data.frame(matrix(ncol = 0, nrow = 12328))
for (i in 1:length(patterns)) {
  LFC_MCF7 <- cbind( LFC_MCF7,(subs.df(as.data.frame(MCF7_pert.sub[,grep(patterns[i], names(MCF7_pert.sub))]), as.data.frame(Mean_ctl[,patterns[i]]))))
  print(i)
}

meta.perts <- read_csv("../New/GEO/GSE92742_Broad_LINCS_inst_info.csv")


names(LFC_MCF7) <- names(MCF7_pert.sub)
rownames(LFC_MCF7) <- rownames(MCF7_pert.sub)


# names(LFC_MCF7) <- na.omit(meta.perts[match(names(LFC_MCF7),meta.perts$inst_id ),4])
#Investigating the deregulation of candidate immunomodulatory peptides by different drugs
# Classifying the immunomodulatory genes to three categories based on their role on cancer progression i.e Up, Down, NoChange
pos.pert <- na.omit(LFC_MCF7[pos_all,])
pos.pert.Down <- as.data.frame(t(pos.pert[c(1,2,6,8),]))
pos.pert.Up <- as.data.frame(t(pos.pert[c(3:5),]))


pert.list <- c()
for (i in 1:length(names(pos.pert.Down))) {
  sort <- pos.pert.Down[order(pos.pert.Down[,i]),]
  pert.list <- c(pert.list, rownames(sort)[1:20])
}

pos.pert.Down.candidate <- data.frame(perts=pert.list,
                                      Genes= c(rep("BST2",20),rep("GAPDH",20),rep("S100A8",20),rep("S100A9",20)),
                                      Effect= c(rep("Down.BST2",20),rep("Down.GAPDH",20),rep("Down.S100A8",20),rep("Down.S100A9",20)))

table(pos.pert.Down.candidate$perts)

subset(pos.pert.Down.candidate, perts%in%c("CPC011_MCF7_24H_X4_B5_DUO52HI53LO:I24"  ,"CPC011_MCF7_24H_X4_B5_DUO52HI53LO:D24", "CPC011_MCF7_24H_X4_B5_DUO52HI53LO:E22"))

pos.pert.Down.candidate[c(41,44,51,64,66,77),3] <- "Down.S100A8/\n   Down.S100A9"

pert.list <- c()
for (i in 1:length(names(pos.pert.Up))) {
  sort <- pos.pert.Up[order(pos.pert.Up[,i],decreasing = T),]
  pert.list <- c(pert.list, rownames(sort)[1:20])
}

pos.pert.Up.candidate <- data.frame(perts=pert.list,
                                    Genes= c(rep("H2AFJ",20),rep("SCGB2A1",20),rep("HMGN2",20)),
                                    Effect= c(rep("Up.H2AFJ",20),rep("Up.SCGB2A1",20),rep("Up.HMGN2",20)))

table(pos.pert.Up.candidate$perts)
intersect(pos.pert.Up.candidate$perts, pos.pert.Down.candidate$perts)

subset(pos.pert.Up.candidate,perts=="CPC007_MCF7_24H_X2_B5_DUO52HI53LO:K19")
subset(pos.pert.Down.candidate,perts=="CPC007_MCF7_24H_X2_B5_DUO52HI53LO:K19")

pos.pert.Up.candidate[c(20),3] <- "Up.H2AFJ/\n   Down.S100A9"
pos.pert.Down.candidate[c(74),3] <- "Up.H2AFJ/\n   Down.S100A9"


pos.pert.Up.Down.candidate <- rbind(pos.pert.Down.candidate,pos.pert.Up.candidate)
pos.pert.Up.Down.candidate  <-  pos.pert.Up.Down.candidate[!duplicated(pos.pert.Up.Down.candidate$perts),]



pert.subset <- LFC_MCF7[match(unique(c(pos.pert.Down.candidate$Genes, pos.pert.Up.candidate$Genes,"B2M","HLA-B", "HLA-C")), 
                              rownames(LFC_MCF7)), match(unique(c(pos.pert.Down.candidate$perts, pos.pert.Up.candidate$perts)), names(LFC_MCF7))]


names(pert.subset) <- c(meta.perts$pert_iname[match(names(pert.subset), meta.perts$inst_id)])

names(pert.subset)[c(50,91,93)] <- c("PRISM001_MCF7_24H_X1_B7_DUO52HI53LO:G03","PRISM001_MCF7_24H_X1_B7_DUO52HI53LO:M05",
                                     "PCLB003_MCF7_24H_X1_B13:D07")

pos.pert.Up.Down.candidate.pub <- as.data.frame(pos.pert.Up.Down.candidate[!duplicated(pos.pert.Up.Down.candidate$perts),])

pos.pert.Up.Down.candidate.pub.final <- meta_info[match(pos.pert.Up.Down.candidate.pub$perts, meta_info$inst_id),]
rownames(pos.pert.Up.Down.candidate.pub.final) <- NULL
pos.pert.Up.Down.candidate.pub.final$Effect <- pos.pert.Up.Down.candidate.pub$Effect

write.csv(pos.pert.Up.Down.candidate.pub.final, "../New/PIC/Supplementary.DRUG.metadata.csv")
write.csv(pert.subset, "../New/PIC/Supplementary.perturbagen.LFC.csv")



# Reading the perturbagen matrix as a seurat object

seurat.perts <- CreateSeuratObject(pert.subset)
seurat.perts <- FindVariableFeatures(seurat.perts, nfeatures = 10)
seurat.perts <- ScaleData(seurat.perts, do.scale = F,do.center = F)
seurat.perts <- RunPCA(seurat.perts)
seurat.perts <- RunUMAP(seurat.perts, dims = 1:9)
seurat.perts$Repos.plan <- pos.pert.Up.Down.candidate$Effect

## Looking at the effect of candidate drugs on immunomodulatory peptides
pdf("../New/PIC/pert.heat.pdf", width=13, height=9)
par(mar = c(6, 6, 6, 6));
DoHeatmap(seurat.perts, features = rownames(pert.subset), group.by = "Repos.plan", size = 4)+
  ggtitle(paste("Drug repositioning to deregulate the immunomodulatory peptides"))+ theme(title =element_text(size=16, face='bold'),legend.position = "right",
                                             legend.text = element_text(size=15), 
                                             axis.text=element_text(size=14))
dev.off()

P <- RidgePlot(seurat.perts, features = rownames(pert.subset)[1:7], group.by = "Repos.plan", ncol = 4,combine = F) 
P <- lapply(X = P, FUN = function(x) x + theme(title =element_text(size=18, face='bold'),legend.position="non",
                                               legend.text = element_text(size=15), 
                                               axis.text=element_text(size=22),plot.title = element_text(size = 35)))



pdf("../New/PIC/pert.RidgePlot.4.pdf", width=30, height=18)
par(mar = c(6, 6, 6, 6));
CombinePlots(plots = P,ncol = 4,label_size=16)
dev.off()


P <- VlnPlot(seurat.perts, features = rownames(pert.subset)[1:7], group.by = "Repos.plan",combine = F,pt.size = 0)

P <- lapply(X = P, FUN = function(x) x + theme(title =element_text(size=18, face='bold'),legend.position="non",
                                               legend.text = element_text(size=15), 
                                               axis.text=element_text(size=22),plot.title = element_text(size = 35)))

pdf("../New/PIC/pert.vln.pdf", width=27, height=20)
par(mar = c(6, 6, 6, 6));
CombinePlots(plots = P,ncol = 3,label_size=16)
dev.off()


DotPlot(seurat.perts, features = rownames(pert.subset)[1:7], group.by = "Repos.plan") + RotatedAxis()+
  ggtitle("The effect of drug repositioning on \npercentage and average expression \nof genes")+theme(title =element_text(size=18, face='bold'),legend.position="right",
                                          legend.text = element_text(size=15), 
                                              axis.text=element_text(size=16),plot.title = element_text(size = 25))