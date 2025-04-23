library(ggplot2)
library(Seurat)
library(dplyr)
library(slingshot)
library(viridis)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)

# Installation of rDBEC package from Github: https://github.com/s-shichino1989/TASSeq-paper)　is required for the reproduction of results.
library(rDBEC)

# Download the file "GSM8288568_matrix_inflection_demulti_DBEC_MiyakeWTA14.txt.gz" from NCBI GEO.
# Download file from this link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8288568
# Set the working directory to the folder containing the "GSM8288568_matrix_inflection_demulti_DBEC_MiyakeWTA14.txt.gz" file.
setwd("/media/owner/Data1/Lung ARDS scRNA-seq 2nd/MiyakeWTA14/MiyakeWTA14_results/Seurat")

#Distibution-based background subtration,hashtag/sampletag demultiplexing, and Seurat Object creation for 
# BD Rhapsody-derived WTA data was conducted by using rDBEC package as previously described 
# Shichino, S. et al. Commun Biol. (2022). doi:10.1038/s42003-022-03536-0
# Github: https://github.com/s-shichino1989/TASSeq-paper)

#=================Creation of  Seurat object and Quality Control===============================#
# read processed count data matrix 
tablelist = lapply(as.vector("matrix_inflection_demulti_DBEC_MiyakeWTA14.txt.gz"), tableread_fast_sparse)

mBC = CreateSeuratObject(tablelist[[1]], min.cells = 5)
View(mBC@meta.data)
table(mBC$orig.ident)

#Remove doublet sample from dataset
mBC <- subset(x = mBC, subset = (mBC$orig.ident =="A961")|(orig.ident =="A962"))
table(mBC$orig.ident)

mBC$log10GenesPerUMI <- log10(mBC$nFeature_RNA) / log10(mBC$nCount_RNA)
mBC$mitoRatio <- PercentageFeatureSet(object = mBC, pattern = "mt-")
mBC$mitoRatio <- mBC@meta.data$mitoRatio / 100

metadata = mBC@meta.data
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA, sample=orig.ident )
metadata$Mouse_origin <- NA
metadata$Mouse_origin[which(str_detect(metadata$sample, "A961"))] <- "WT"
metadata$Mouse_origin[which(str_detect(metadata$sample, "A962"))] <- "IL-4 KO"
mouse_levels = c("WT","IL-4 KO")
metadata$Mouse_origin= factor(x=metadata$Mouse_origin, levels = mouse_levels)
View(metadata)
mBC@meta.data = metadata


# QC visualization (can be skipped to QC filtering)
# Visualize the number of cell counts per Mouse_origin
metadata %>%   ggplot(aes(x=Mouse_origin, fill=Mouse_origin)) + geom_bar() +  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +ggtitle("NCells")
# Visualize the number UMIs/transcripts per cell
metadata %>%   ggplot(aes(color=Mouse_origin, x=nUMI, fill= Mouse_origin)) + 
  geom_density(alpha = 0.2) +   scale_x_log10() + theme_classic() +
  ylab("Cell density") +  geom_vline(xintercept = 500)
# Visualize the distribution of genes detected per cell via histogram
metadata %>% ggplot(aes(color=Mouse_origin, x=nGene, fill= Mouse_origin)) + 
  geom_density(alpha = 0.2) +   theme_classic() +  scale_x_log10() + geom_vline(xintercept = 300)
# Visualize the distribution of genes detected per cell via boxplot
metadata %>% ggplot(aes(x=Mouse_origin, y=log10(nGene), fill=Mouse_origin)) +   geom_boxplot() +   theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +  ggtitle("NCells vs NGenes")
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>%   ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +   geom_point() +   scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +  scale_x_log10() +   scale_y_log10() + 
  theme_classic() +  geom_vline(xintercept = 500) +  geom_hline(yintercept = 250) +  facet_wrap(~Mouse_origin)
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%ggplot(aes(x=log10GenesPerUMI, color = Mouse_origin, fill=Mouse_origin)) +
  geom_density(alpha = 0.2) +  theme_classic() +  geom_vline(xintercept = 0.70)

# QC Filtering
# Filter out low quality reads using selected thresholds - these will change with experiment
# Change filtering criteria for log10GenesPerUMI to 0.70 (originally 0.80)
filtered_seurat <- subset(x = mBC, subset= (nUMI >= 500) &  (nGene >= 250) &  (log10GenesPerUMI > 0.70) & (mitoRatio < 0.20))
View(filtered_seurat@meta.data)
# Output a logical vector for every gene on whether the more than zero counts per cell
counts <- GetAssayData(object = filtered_seurat, slot = "counts") # Extract counts
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

save(filtered_seurat, file="filtered_seurat_ARDS.rda")

#===========================Data normalization & PCA=========================================#
# Data normalization & PCA
load("filtered_seurat_ARDS.rda")

mBC = NormalizeData(object = filtered_seurat, scale.factor=1000000)
all.genes <- rownames(mBC)
mBC = ScaleData(object = mBC, vars.to.regress = c("nUMI"), features = all.genes, verbose = TRUE)
mBC = FindVariableFeatures(mBC, selection.method = "vst", nfeatures = 2000)
mBC <- RunPCA(mBC, npcs = 100, verbose = FALSE)
mBC = JackStraw(mBC, num.replicate = 100, dims = 100) 
mBC = ScoreJackStraw(mBC, dims = 1:100)
JackStrawPlot(mBC, dims = 1:100)
ElbowPlot(object = mBC)

#Clustering
tmp = as.data.frame(mBC@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))

mBC <- RunUMAP(mBC, reduction = "pca", dims = dims)
mBC <- FindNeighbors(mBC, reduction = "pca", dims = dims)
mBC <- FindClusters(mBC, resolution = 1.0)

#FindAllMarkers
markers <- FindAllMarkers(mBC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file="ARDS_lung_clustermarkergenes.csv", sep=",")


# Visualization
DimPlot(mBC, reduction = "umap", label=T)+coord_fixed()
DimPlot(mBC, reduction = "umap", label=T, split.by="Mouse_origin", ncol=1)+coord_fixed()
View(mBC@meta.data)
table(mBC$Mouse_origin)

# Add cell type metadata
metadata <- mBC@meta.data  
View(metadata)
metadata$cell_type <- "Neutrophil"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(0,9,10,14,16), collapse = "|")))] <- "MoMac"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(17), collapse = "|")))] <- "non-classical mono"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(4), collapse = "|")))] <- "AM"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(5,21,24), collapse = "|")))] <- "Endothelial cell"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(6), collapse = "|")))] <- "B cell"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(7,20), collapse = "|")))] <- "T cell"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(12), collapse = "|")))] <- "NK cell"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(13), collapse = "|")))] <- "Proliferating AM"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(15), collapse = "|")))] <- "pDC"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(18), collapse = "|")))] <- "Fibroblast"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(22), collapse = "|")))] <- "AT1/AT2"
metadata$cell_type[which(str_detect(metadata$seurat_clusters, paste(c(23), collapse = "|")))] <- "Basophil"
celltype_levels = c("Neutrophil","MoMac","AM","Proliferating AM","T cell","B cell","NK cell","non-classical mono","pDC",
                    "Endothelial cell","Fibroblast","AT1/AT2","Basophil")
metadata$cell_type= factor(x=metadata$cell_type, levels = celltype_levels)
mBC@meta.data <- metadata
View(metadata)

# Visualization of WT lung cells
Idents(mBC) <-"Mouse_origin"
mBC_WT= subset(mBC, idents ="WT")
Idents(mBC)<-"seurat_clusters"
Idents(mBC_WT)<-"seurat_clusters"

Idents(mBC_WT) <-"cell_type"
DimPlot(mBC_WT, reduction = "umap", label=F)+coord_fixed()
DimPlot(mBC, reduction = "umap", split.by = "Mouse_origin", label=F)+coord_fixed()

#FindAllMarkers
markers <- FindAllMarkers(mBC_WT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, file="ARDS_WT_lung_celltype_markergenes.csv", sep=",")
View(markers)
Baso_marker =filter(markers, cluster =="Basophil") %>% top_n(-10, p_val_adj)
View(Baso_marker)

celltype_levels2 = c("Basophil", "Neutrophil","MoMac","AM","Proliferating AM","T cell","B cell","NK cell","non-classical mono","pDC",
                    "Endothelial cell","Fibroblast","AT1/AT2")
metadata_WT=mBC_WT@meta.data
metadata_WT$cell_type= factor(x=metadata_WT$cell_type, levels = celltype_levels2)
mBC_WT@meta.data=metadata_WT
Idents(mBC_WT) <-"cell_type"

DotPlot(mBC_WT, features = rownames(Baso_marker))+RotatedAxis()+coord_flip()+  
  theme(axis.text.y = element_text(face="italic"))

metadata_WT=mBC_WT@meta.data
metadata_WT$cell_type= factor(x=metadata_WT$cell_type, levels = celltype_levels)
mBC_WT@meta.data=metadata_WT
Idents(mBC_WT) <-"cell_type"

#select 3 marker genes from each cell type
Marker = c("Hcar2","Il1b","S100a8","C1qc","C1qb","Apoe","Ear1","Ear2","Fabp1","Ube2c","Ccna2","Top2a","Cd3g","Cd3d","Cd3e",
           "Cd79a","Igkc","Ighd","Klra8","Gzma","Ncr1","Ccl22","Flt3","Itgae","Cldn5","Cxcl12","Cdh5","Mgp","Eln","Col1a2",
           "Cxcl15","Sftpb","Cbr2")
Marker = union (Marker, rownames(Baso_marker))
DotPlot(mBC_WT, features = Marker)+RotatedAxis()+coord_flip()+  
  theme(axis.text.y = element_text(face="italic"))

#calculate the number of cells and the proportion of each cluster among Neu
pie(table(Idents(mBC_WT)),clockwise = TRUE,
    init.angle = 90,border = TRUE)
cellnum_table = table(Idents(mBC_WT))
ggtable = as.data.frame(cellnum_table)
colnames(ggtable)=c("cell_type","Frequency")
ggtable
ggtable$cell_type <- factor(ggtable$cell_type, levels = celltype_levels3)
celltype_levels3 = c("Basophil","AT1/AT2","Fibroblast","Endothelial cell","pDC","non-classical mono","NK cell","B cell","T cell","Proliferating AM","AM","MoMac","Neutrophil")

ggplot(ggtable, aes(x = 1, y=Frequency,fill=cell_type)) +
  coord_polar(theta="y")+ geom_col()+ geom_col(color = "black")+ theme_void()+ theme(legend.position = "right")+
  scale_fill_manual(values = c("blue", "beige","darkgreen","royalblue","chocolate","yellowgreen","magenta","orange","brown","yellow","purple","red"))


#===========================Neutrophil re-clustering========================================#
Neu = subset(mBC, idents = "Neutrophil")

# Dotplot visualization
Genes = c("Atf5","Ddit3","Bcl2a1d","Bcl2a1b","Bcl2a1a","S100a9","S100a8","Itgam","Fpr2","Fpr1","Il1b","Il1a","Cxcl2")
Idents(Neu)<-"Mouse_origin"
# obtain count matrix by using Seurat dotplot function
DotPlot(Neu,features = Genes)+scale_y_discrete(position = "right")+coord_flip()+  
  theme(axis.text.y = element_text(face="italic"))

# Neutrophil re-clustering
Neu = NormalizeData(Neu, scale.factor = 1000000)
Neu = FindVariableFeatures(Neu, mean.function = ExpMean, 
                             dispersion.function = LogVMR, 
                             mean.cutoff = c(0.1,Inf),
                             dispersion.cutoff = c(0.5,Inf))
all.genes <- rownames(Neu)

Neu = ScaleData(Neu, vars.to.regress = c("nReads"), features = all.genes)
Neu = RunPCA(Neu, features = VariableFeatures(Neu), npcs = 100)
Neu = ProjectDim(Neu)
Neu = JackStraw(Neu, num.replicate = 100, dims = 100)
Neu = ScoreJackStraw(Neu, dims = 1:100)
JackStrawPlot(Neu, dims = 1:40)
ElbowPlot(object = Neu)
tmp = as.data.frame(Neu@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))
Neu = FindNeighbors(Neu, dims = dims)
Neu = FindClusters(Neu, resolution = 0.3)
Neu = RunUMAP(Neu, dims = dims)

# Add cluster names
metadata <- Neu@meta.data  
View(metadata)
metadata$cluster <-NA
metadata$cluster[which(str_detect(metadata$RNA_snn_res.0.3, "1"))] <- "Neu1"
metadata$cluster[which(str_detect(metadata$RNA_snn_res.0.3, "0"))] <- "Neu2"
metadata$cluster[which(str_detect(metadata$RNA_snn_res.0.3, "3"))] <- "Neu3"
metadata$cluster[which(str_detect(metadata$RNA_snn_res.0.3, "4"))] <- "Neu4"
metadata$cluster[which(str_detect(metadata$RNA_snn_res.0.3, "2"))] <- "Neu5"
cluster_levels = c("Neu1","Neu2","Neu3","Neu4","Neu5")
metadata$cluster= factor(x=metadata$cluster, levels = cluster_levels)
Neu@meta.data = metadata

#DimPlot visualization
Idents(Neu) <-"cluster"
cols = c('Neu2'="blue", "Neu1"="purple","Neu5"="#FF0000","Neu3"="#00B050","Neu4"="pink")
DimPlot(Neu, reduction = "umap", pt.size = 1, label=T, cols=cols)+coord_fixed()
DimPlot(Neu, reduction = "umap", split.by="Mouse_origin", pt.size = 1, ncol=1, label=T, cols=cols)+coord_fixed()

#calculate the number of cells and the proportion of each cluster among Neu
table(Idents(Neu))
cellnum_table = table(Idents(Neu), Neu$Mouse_origin)
head(cellnum_table)
ggtable = as.data.frame(cellnum_table)
head(ggtable)

colnames(ggtable)=c("cluster","mouse","Frequency")
ggplot(ggtable, aes(x=mouse, y=Frequency, fill=cluster))+
  geom_bar(stat = "identity", position = "fill", width=0.5, color="black")+ 
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x=element_text(size=30))+
  scale_fill_manual(values = c("purple", "blue", "#00B050", "pink","red"))+
  theme_classic()

save(Neu, dims, file="Neu_reclustering.rda")

#===========================Differential expression analysis between Neu1 and Neu5========================================#
#Differential gene expression analysis using FindMarkers function
DEGs <- FindMarkers(Neu, ident.1 = c("Neu5"), ident.2 = c("Neu1"), 
                    logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
View(DEGs)
DEGs$q_val <- p.adjust(DEGs$p_val, method = "BH")
DEGs$gene <- rownames(DEGs)
write.csv(DEGs, file="ARDS_Neu5_vs_Neu1.csv")

# Gene Set Enrichment Analysis by using clusterProfiler package
DEGlist = DEGs[,c("gene", "avg_log2FC")] %>% arrange(desc(avg_log2FC))
DEGgenes = DEGlist$gene
DEGgenes = bitr(DEGgenes, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
colnames(DEGgenes)[1]="gene" 
DEGlist= merge(DEGlist, DEGgenes, by="gene")
DEG = as.vector(DEGlist$avg_log2FC)
names(DEG) <- as.vector(DEGlist$ENTREZID)
DEG = sort(DEG, decreasing = T)

# GO GSEA
gseGO_diff <- gseGO(geneList     = DEG,
                    OrgDb        = org.Mm.eg.db,
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    verbose      = FALSE)

res_GO = as.data.frame(gseGO_diff)
gseGO_diff_simplified = simplify(gseGO_diff)
View(res_GO)
write.csv(res_GO, file="ARDS_Neu5_vs_Neu1_GSEA_GO.csv")


Interest = c("GO:0032611","GO:0045088","GO:0002221","GO:0002274")

GsID = Interest[1]
GsName = res_GO$Description[res_GO$ID==GsID]
Gspval=res_GO$qvalue[res_GO$ID==GsID]
GsNES=res_GO$NES[res_GO$ID==GsID]
g = gseaplot(gseGO_diff, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
g[[2]]=g[[2]]+
    annotate("text",x=15000,y=0.4,label=paste0("NES: ",round(GsNES,3)),size=6)+
    annotate("text",x=16500,y=0.35,label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
g
   
GsID = Interest[2]
GsName = res_GO$Description[res_GO$ID==GsID]
Gspval=res_GO$qvalue[res_GO$ID==GsID]
GsNES=res_GO$NES[res_GO$ID==GsID]
g = gseaplot(gseGO_diff, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
g[[2]]=g[[2]]+
  annotate("text",x=15000,y=0.4,label=paste0("NES: ",round(GsNES,3)),size=6)+
  annotate("text",x=16500,y=0.35,label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
g

GsID = Interest[3]
GsName = res_GO$Description[res_GO$ID==GsID]
Gspval=res_GO$qvalue[res_GO$ID==GsID]
GsNES=res_GO$NES[res_GO$ID==GsID]
g = gseaplot(gseGO_diff, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
g[[2]]=g[[2]]+
  annotate("text",x=15000,y=0.35,label=paste0("NES: ",round(GsNES,3)),size=6)+
  annotate("text",x=16500,y=0.3,label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
g

GsID = Interest[4]
GsName = res_GO$Description[res_GO$ID==GsID]
Gspval=res_GO$qvalue[res_GO$ID==GsID]
GsNES=res_GO$NES[res_GO$ID==GsID]
g = gseaplot(gseGO_diff, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
g[[2]]=g[[2]]+
  annotate("text",x=15000,y=0.35,label=paste0("NES: ",round(GsNES,3)),size=6)+
  annotate("text",x=16500,y=0.3,label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
g
#KEGG GSEA
gseGO_KEGG <- gseKEGG(geneList     = DEG,
                      organism        = "mmu",
                      keyType = "ncbi-geneid",
                      minGSSize = 2,
                      pvalueCutoff = 1,
                      use_internal_data = F, 
                      verbose      = FALSE)
res_KEGG = as.data.frame(gseGO_KEGG)
View(res_KEGG)

#KEGG GSEA visualization
GsID = "mmu04620"
GsName = res_KEGG$Description[res_KEGG$ID==GsID]
GsName = str_split(GsName, pattern = " - ")[[1]][1]
Gspval=res_KEGG$qvalue[res_KEGG$ID==GsID]
GsNES=res_KEGG$NES[res_KEGG$ID==GsID]
p=gseaplot(gseGO_KEGG, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]]=p[[2]] +
  annotate("text",x=15000,y=0.45,label=paste0("NES: ",round(GsNES,3)),size=6)+
  annotate("text",x=16500,y=0.4,label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p

GsID = "mmu04670"
GsName = res_KEGG$Description[res_KEGG$ID==GsID]
GsName = str_split(GsName, pattern = " - ")[[1]][1]
Gspval=res_KEGG$qvalue[res_KEGG$ID==GsID]
GsNES=res_KEGG$NES[res_KEGG$ID==GsID]
p=gseaplot(gseGO_KEGG, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]]=p[[2]] +
  annotate("text",x=15000,y=0.4,label=paste0("NES: ",round(GsNES,3)),size=6)+
  annotate("text",x=16500,y=0.35,label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p

GsID = "mmu04621"
GsName = res_KEGG$Description[res_KEGG$ID==GsID]
GsName = str_split(GsName, pattern = " - ")[[1]][1]
Gspval=res_KEGG$qvalue[res_KEGG$ID==GsID]
GsNES=res_KEGG$NES[res_KEGG$ID==GsID]
p=gseaplot(gseGO_KEGG, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]]=p[[2]] +
  annotate("text",x=15000,y=0.4,label=paste0("NES: ",round(GsNES,3)),size=6)+
  annotate("text",x=16500,y=0.35,label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p

GsID = "mmu04668"
GsName = res_KEGG$Description[res_KEGG$ID==GsID]
GsName = str_split(GsName, pattern = " - ")[[1]][1]
Gspval=res_KEGG$qvalue[res_KEGG$ID==GsID]
GsNES=res_KEGG$NES[res_KEGG$ID==GsID]
p=gseaplot(gseGO_KEGG, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]]=p[[2]] +
  annotate("text",x=15000,y=0.4,label=paste0("NES: ",round(GsNES,3)),size=6)+
  annotate("text",x=16500,y=0.35,label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p




