
# Step3. Differential gene expression
## scRNAseq Data

path="Directory"
setwd(path)
load("Step2 Annotation_Data integration.RData")


# Differential gene expression between BC subtypes and normal epithelial cells to track the immunomodulatory peptides

all.int <- SetIdent(all.int, value = all.int$Type)
sc.MalvsNorm <- FindMarkers(all.int, ident.1 = "Malignant", ident.2 = "Normal" ,logfc.threshold = FALSE, assay = "RNA",densify = TRUE)
all.int <- SetIdent(all.int, value = all.int$subtype)

sc.TNBCvsNorm <- FindMarkers(all.int, ident.1 = "TNBC", ident.2 = "Normal", assay = "RNA",logfc.threshold = FALSE, densify = TRUE)
sc.HER2vsNorm <- FindMarkers(all.int, ident.1 = "HER2+", ident.2 = "Normal" , assay = "RNA",logfc.threshold = FALSE, densify = TRUE)
sc.ERvsNorm <- FindMarkers(all.int, ident.1 = "ER+", ident.2 = "Normal", assay = "RNA",logfc.threshold = FALSE, densify = TRUE)

save(sc.MalvsNorm,sc.TNBCvsNorm,sc.HER2vsNorm,sc.ERvsNorm, file = "../New/DEGs.RDATA")
load("../New/DEGs.RDATA")

#######################################################
# Loading the theraputic peptides from UDAMP database #
#######################################################
ACPs <- read_csv("Directory/UDAMP-Immu-modul.csv")
ACPs <- c(ACPs$`Gene name`, c("HNP1", "MENK", "BRD2", "BRD3", "BRD4", "BRDT", "HMGB1", "HAS2", "HAS1"))

#Volcanoplot to show the deregulation of immunomodulatory peptides between BC subtypes and nomral epithelial cells
p1<- EnhancedVolcano(sc.MalvsNorm,
                     lab = rownames(sc.MalvsNorm),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     selectLab = c(rownames(subset(na.omit(sc.MalvsNorm[ACPs,]), avg_log2FC>0.5 & p_val_adj<0.05))),
                     title = 'Malignant vs Normal',
                     subtitle = "Differential expression",
                     pCutoff = 0.05,
                     FCcutoff = 0.5,
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
                     xlim = c(-2.5,2.5),
                     caption = NULL
)+theme(legend.position = "none")



ACPs.up.MalvsNorm <- rownames(subset(na.omit(sc.MalvsNorm[ACPs,]), avg_log2FC>0.5 & p_val_adj<0.05))


p2<- EnhancedVolcano(sc.TNBCvsNorm,
                     lab = rownames(sc.TNBCvsNorm),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     selectLab = ACPs.up.MalvsNorm,
                     title = 'TNBC vs Normal',
                     subtitle = "Differential expression",
                     pCutoff = 0.05,
                     FCcutoff = 3,
                     cutoffLineType = 'twodash',
                     cutoffLineWidth = 0,
                     pointSize = 3.0,
                     labSize = 4,
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

p3<- EnhancedVolcano(sc.HER2vsNorm,
                     lab = rownames(sc.HER2vsNorm),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlab = bquote(~Log[2]~ 'fold change'),
                     selectLab = ACPs.up.MalvsNorm,
                     title = 'HER2 vs Normal',
                     subtitle = "Differential expression",
                     pCutoff = 0.05,
                     FCcutoff = 3,
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

p4 <- EnhancedVolcano(sc.ERvsNorm,
                      lab = rownames(sc.ERvsNorm),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      selectLab = ACPs.up.MalvsNorm,
                      title = 'ER vs Normal',
                      subtitle = "Differential expression",
                      pCutoff = 0.05,
                      FCcutoff = 3,
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


pdf("DEGs_Log.pdf", width=13, height=12)
par(mar = c(6, 6, 6, 6));
ggarrange(p1,p2, p3,p4,
          ncol = 2, nrow = 2)
dev.off()

rownames(subset(na.omit(sc.MalvsNorm[ACPs,]), avg_log2FC>0.5 & p_val_adj<0.05))
rownames(subset(na.omit(sc.TNBCvsNorm[ACPs,]), avg_log2FC>0 & p_val_adj<0.05))
rownames(subset(na.omit(sc.HER2vsNorm[ACPs,]), avg_log2FC>0 & p_val_adj<0.05))
rownames(subset(na.omit(sc.ERvsNorm[ACPs,]), avg_log2FC>0 & p_val_adj<0.05))



## ACPs

# Distribution of LFC for the expression of candidate peptides

pos_all <- c(rownames(subset(na.omit(sc.MalvsNorm[ACPs,]), avg_log2FC>0.5 & p_val_adj<0.05)))

box.LFC <- data.frame(avg=c(sc.MalvsNorm[pos_all,]$avg_log2FC,sc.TNBCvsNorm[pos_all,]$avg_log2FC,sc.HER2vsNorm[pos_all,]$avg_log2FC, sc.ERvsNorm[pos_all,]$avg_log2FC),
                      subtype=c(rep("MalvsNorm",10),rep("TNBCvsNorm",10),rep("HER2vsNorm",10),rep("ERvsNorm",10) ),
                      Gene=rep(pos_all,4),
                      Repositioning_strategy=rep(c("Down", "Down", "Up","Up","Up", "Down", "No change", "Down", "No change", "Down"),4) )

box.LFC <- na.omit(box.LFC)


pdf("../New/PIC/LFC.pdf", width=5, height=8)
par(mar = c(6, 6, 6, 6));
ggboxplot(box.LFC , x = "subtype", y = "avg"
          ,color = "subtype"
          ,ylab = ("LogFoldChange"), xlab = "")+rotate_x_text(45)+ggtitle(paste("Distribution of LFC for the \nexpression of candidate peptides"))+ theme(title =element_text(size=16, face='bold'),
                                                                                                                                                          legend.text = element_text(size=15), 
                                                                                                                                                          axis.text.y=element_text(size=14), axis.text.x=element_text(size=0), legend.position = "right")+facet_wrap(.~Repositioning_strategy,scales="free", ncol = 1)
dev.off()


###########################
## Gene-Gene correlation ##
###########################

# Gene Gene correlation between highli variable features across BC subtypes and normal epithelial cells

up.MalvsNorm <- rownames(subset(na.omit(sc.MalvsNorm), avg_log2FC>0.3 & p_val_adj<0.05))
int.data <- na.omit(as.data.frame(all.int@assays$integrated@data)[c(up.MalvsNorm),])

all.var.scaled <- all.int@assays$integrated@scale.data
dim(all.var.scaled)
vari <- rowVars(all.var.scaled)

quantile(vari,probs = seq(0, 1, 1/40))

all.var.scaled <- as.data.frame(all.var.scaled)
all.var.scaled$vari <- vari
all.var.scaled <- subset(all.var.scaled,vari>0.3)
all.var.scaled <- all.var.scaled[,-length(colnames(all.var.scaled))]

corr <- rcorr(as.matrix(t(all.var.scaled)), type = "pearson")


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

corr.out <- flattenCorrMatrix(corr$r, corr$P)
corr.out$adj.p <- p.adjust(corr.out$p, method = "fdr")
saveRDS(corr.out, "../New/all.int.corr.out.rds")
corr.out <- readRDS("all.int.corr.out.rds")


corr.sub <- subset(corr.out, adj.p < 0.05 & cor > 0 )
corr.sub <- subset(corr.sub,row%in% pos_all |column%in%pos_all)

corr.sub <- corr.sub[order(corr.sub$cor,decreasing = T),]

rownames(corr.sub) <- NULL

# Extracting the top 20 correlated genes with deregulated immunomodulatory peptides

cor.df <- as.data.frame(matrix(nrow = 0,ncol = 0))
for (i in pos_all) {
  x <- subset(corr.sub, row==i | column==i)
  x <- x[order(x$cor,decreasing = T),]
  cor.df <- rbind(cor.df, x[1:20,])
}

x <- list(
  BST2 <- subset(cor.df, row==pos_all[1] | column==pos_all[1]),
  GAPDH <- subset(cor.df, row==pos_all[2] | column==pos_all[2]),
  H2AFJ <- subset(cor.df, row==pos_all[3] | column==pos_all[3]),
  SCGB2A1 <- subset(cor.df, row==pos_all[4] | column==pos_all[4]),
  HMGN2 <- subset(cor.df, row==pos_all[5] | column==pos_all[5]),
  S100A8 <- subset(cor.df, row==pos_all[6] | column==pos_all[6]),
  S100A7 <- subset(cor.df, row==pos_all[7] | column==pos_all[7]),
  S100A9 <- subset(cor.df, row==pos_all[8] | column==pos_all[8]),
  ZG16B <- subset(cor.df, row==pos_all[9] | column==pos_all[9]),
  HMGB1 <- subset(cor.df, row==pos_all[10] | column==pos_all[10]))


# Top 10 enriched GO pathways for correlated genes with each immnomodulatory peptide

name <- list()
for (i in 1:10) {
  geneGO <- enrichGO(gene = c(x[[i]]$row, x[[i]]$column),
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.2,
                     qvalueCutoff = 0.2)
  
  name[[i]] <- print(barplot(geneGO,showCategory = 10, font.size = 14)+ggtitle(pos_all[i])+ theme(title =element_text(size=16, face='bold'),
                                                                                                  legend.text = element_text(size=15), 
                                                                                                  axis.text=element_text(size=14), legend.position = "right"))
}


figure.GO <- ggarrange(name[[1]],name[[2]], name[[5]], name[[6]], name[[7]], name[[8]], name[[9]], name[[10]],
                       ncol = 3, nrow = 3)

pdf("Directory/figure.GO.pdf", width=18, height=12)
par(mar = c(6, 6, 6, 6));
figure.GO
dev.off()
# cnetplot(geneGO, categorySize="pvalue", circular = TRUE,
#          colorEdge = TRUE, node_label= "gene")+
#   theme(legend.text=element_text(size=14),
#     legend.title=element_text(size=14))

# 
# 

# Top 10 enriched KEGG pathways for correlated genes with each immnomodulatory peptide

hsens=EnsDb.Hsapiens.v86
name <- list()
for (i in 1:10) {
  my.symbols <- c(x[[i]]$row, x[[i]]$column)
  enterz<- AnnotationDbi::select(hsens,
                                 keys = my.symbols,
                                 columns = c("ENTREZID", "SYMBOL", "GENEID"),
                                 keytype = "SYMBOL")
  enterz=enterz[complete.cases(enterz), ]
  
  
  
  geneKEGG <- enrichKEGG(gene         = enterz$ENTREZID,
                         organism = 'hsa',
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.15,
                         qvalueCutoff  = 0.15,
                         keyType = "ncbi-geneid")
  
  
  edox <- setReadable(geneKEGG, 'org.Hs.eg.db', 'ENTREZID')
  # 
  # 
  name[[i]] <- print(barplot(geneKEGG,showCategory = 10, font.size = 14)+ggtitle(pos_all[[i]]))+ theme(title =element_text(size=16, face='bold'),
                                                                                                       legend.text = element_text(size=15), 
                                                                                                       axis.text=element_text(size=14), legend.position = "right")
}


figure.KEGG <- ggarrange(name[[1]],name[[2]], name[[4]], name[[5]], name[[6]], name[[7]], name[[8]],
                         ncol = 3, nrow = 3)

pdf("Directory/figure.KEGG.pdf", width=18, height=13)
par(mar = c(6, 6, 6, 6));
figure.KEGG
dev.off()


DotPlot(all.int, assay = "RNA", features = pos_all)+rotate_x_text(55)+ggtitle("Gene expression of candidate peptides")
DotPlot(protein_seurat, features = pos_all)+rotate_x_text(55)+ggtitle("Protein expression of \ncandidate peptides")




# Average of gene and protein expression of deregulated immunomodulatory peptides in different BC subtypes
avg.gene <- as.data.frame(AverageExpression(all.int, features = pos_all, assays = "RNA"))
avg.protein <- as.data.frame(AverageExpression(protein_seurat,features = pos_all))



avg.gene.list <- c()
for (i in 1:length(names(avg.gene))) {
  avg.gene.list <- c(avg.gene.list,c(avg.gene[,i]))
}


avg.protein.list <- c()
for (i in 1:length(names(avg.protein))) {
  avg.protein.list <- c(avg.protein.list,c(avg.protein[,i]))
}



avg.gene.plot <- data.frame(avg=avg.gene.list,
                            Subtype=c(rep("HER2", 10), rep("ER", 10),rep("TNBC", 10), rep("Normal", 10)),
                            Mal=c(rep("Malignant",30), rep("Normal",10)),
                            Gene=rep(pos_all,4),
                            Repositioning_strategy=rep(c("Down", "Down", "Up","Up","Up", "Down", "No change", "Down", "No change", "Down"),4))


avg.protein.plot <- data.frame(avg=avg.protein.list,
                               Subtype=c(rep("HER2", 4), rep("ER", 4),rep("TNBC", 4)),
                               Gene=rep(rownames(avg.protein),3),
                               Repositioning_strategy=rep(c("Down","Down","Up","Down"),3))



ggboxplot(avg.gene.plot , x = "Subtype", y = "avg"
          ,color = "Subtype"
          ,ylab = ("Average expression"), xlab = "")+rotate_x_text(45)+facet_wrap(.~Repositioning_strategy,scales="free", ncol = 3)+ggtitle(paste("Average gene expression of candidate peptides"))+ theme(title =element_text(size=16, face='bold'),
                                                                                                                                                                                                           legend.text = element_text(size=15), 
                                                                                                                                                                                                           axis.text=element_text(size=14))


ggboxplot(avg.protein.plot , x = "Subtype", y = "avg"
          ,color = "Subtype"
          ,ylab = ("Average expression"), xlab = "")+rotate_x_text(45)+facet_wrap(.~Repositioning_strategy,scales="free", ncol = 2)+ggtitle(paste("Average protein expression of \ncandidate peptides"))+ theme(title =element_text(size=16, face='bold'),
                                                                                                                                                                                                                legend.text = element_text(size=15), 
                                                                                                                                                                                                                axis.text=element_text(size=14))



pdf("../New/PIC/avg.gene.sep.pdf", width=6, height=10)
par(mar = c(6, 6, 6, 6));
ggplot(data=avg.gene.plot, aes(x=Subtype, y=avg, fill= Gene)) +
  geom_bar(stat="identity" )+
  theme_minimal()+rotate_x_text(45)+ggtitle(paste("Average gene expression of \ncandidate peptides"))+facet_wrap(.~Repositioning_strategy,scales="free", ncol = 1)+ theme(title =element_text(size=16, face='bold'),
                                                                                                                                                                          legend.text = element_text(size=15), 
                                                                                                                                                                          axis.text=element_text(size=14),strip.text = element_text(size = 14))+ylab("Average gene expression")
dev.off()

pdf("Directory/avg.pro.sep.pdf", width=6, height=7)
par(mar = c(6, 6, 6, 6));
ggplot(data=avg.protein.plot, aes(x=Subtype, y=avg, fill= Gene)) +
  geom_bar(stat="identity" )+
  theme_minimal()+rotate_x_text(45)+ggtitle(paste("Average protein expression of \ncandidate peptides"))+facet_wrap(.~Repositioning_strategy,scales="free", ncol = 1)+ theme(title =element_text(size=16, face='bold'),
                                                                                                                                                                             legend.text = element_text(size=15), 
                                                                                                                                                                             axis.text=element_text(size=14),strip.text = element_text(size = 14))+ylab("Average protein expression")
dev.off()
