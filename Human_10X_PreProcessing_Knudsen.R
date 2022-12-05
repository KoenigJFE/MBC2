#conda activate velo3
setwd("/home/npkdk/drives/x/Research/NPKDK/10xCombined/SCT-update/feat3000_new_paper")
#source("/home/npkdk/drives/x/Research/NPKDK/10xCombined/SCT-update/feat3000_new_paper/10x_combi_SCT_integration_SCT-update_feat2000_02_from_scratch_re_new_paper.R")


library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(reshape2)
library(patchwork)
library(SingleR) #BiocManager::install("SingleR")
library(SingleCellExperiment) #BiocManager::install("SingleCellExperiment")
library(AnnotationDbi) #BiocManager::install("AnnotationDbi")
library(EnsDb.Hsapiens.v86) #BiocManager::install("EnsDb.Hsapiens.v86")
library(celldex) #BiocManager::install("celldex")
library(ensembldb) #BiocManager::install("ensembldb")
library(xlsx)
library(dplyr) 
library(flextable) #conda install -c conda-forge r-flextable
library(janitor) # fconda install -c conda-forge r-janitor
library(velocyto.R) #remotes::install_github("rrydbirk/velocyto.R") or #conda install -c bioconda r-velocyto.r
library(SeuratWrappers) #remotes::install_github('satijalab/seurat-wrappers')
library(Matrix)
theme_set(theme_cowplot())
library(viridis) #install.packages("viridis")
#


## 10x data

#10x single-cell 5' gene-expression data from 10x runs 5, 6, 7, documented in
#ELN104626, ELN105020 and ELN105157.

### Read individual data sets and filter


#setwd("/home/npkdk/drives/x/Research/NPKDK/10xCombined/DATA/cellranger-6.1.2-GRCh38-2020-A/5GEX/")
print("Load data 10x data")
dir <- "/home/npkdk/drives/x/Research/NPKDK/10xCombined/DATA/cellranger-6.1.2-GRCh38-2020-A/5GEX/"
print(dir)
filenames <- c("TT04_subj1_V2", 
    "TT04_subj1_V3", 
    "TT04_subj2_V2", 
    "TT04_subj2_V3",
    "TT04_subj3_V2", 
    "TT04_subj3_V3", 
    "TT04_subj4_V2", 
    "TT04_subj4_V3", 
    "TT04_subj5_V2", 
    "TT04_subj5_V3", 
    "TT04_subj6_V2", 
    "TT04_subj6_V3",
    "allergic_1", 
    "allergic_2", 
    "allergic_3", 
    "allergic_4",
    "non_allergic_1", 
    "non_allergic_2", 
    "non_allergic_3", 
    "non_allergic_4", 
    "non_allergic_5")

data.10x = list()

for (index in 1:length(filenames)) {
    fn <- filenames[index]
    print(fn)
    print(index)
    data.10x[[index]] <- Read10X(data.dir = paste(dir,fn,"/filtered_feature_bc_matrix", sep=""))    
}

print("CreateSeuratObject")
scrna.list = list()

for (i in 1:length(data.10x)) {
    print(i)
    scrna.list[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=10, min.features=200, project=filenames[i]);
    scrna.list[[i]][["SampleID"]] = filenames[i];
}

names(scrna.list) <- filenames

print("PercentageFeatureSet")
for (i in seq_len(length(scrna.list))) {
  scrna.list[[i]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[i]], pattern = "^MT-")
  scrna.list[[i]][["percent.ribo"]] <- PercentageFeatureSet(scrna.list[[i]], pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA")
}


print("Subset nFeature and percent.mt")
for (i in seq_len(length(scrna.list))) {
  scrna.list[[i]]  <- subset(scrna.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 7)
}

scrna.list$allergic_1[["Type"]] <- "Allergic"
scrna.list$allergic_2[["Type"]] <- "Allergic"
scrna.list$allergic_3[["Type"]] <- "Allergic"
scrna.list$allergic_4[["Type"]] <- "Allergic"

scrna.list$non_allergic_1[["Type"]] <- "Non-Allergic"
scrna.list$non_allergic_2[["Type"]] <- "Non-Allergic"
scrna.list$non_allergic_3[["Type"]] <- "Non-Allergic"
scrna.list$non_allergic_4[["Type"]] <- "Non-Allergic"
scrna.list$non_allergic_5[["Type"]] <- "Non-Allergic"

scrna.list$TT04_subj1_V2[["Type"]] <- "Allergic"
scrna.list$TT04_subj2_V2[["Type"]] <- "Allergic"
scrna.list$TT04_subj3_V2[["Type"]] <- "Allergic"
scrna.list$TT04_subj4_V2[["Type"]] <- "Allergic"
scrna.list$TT04_subj5_V2[["Type"]] <- "Allergic"
scrna.list$TT04_subj6_V2[["Type"]] <- "Allergic"

scrna.list$TT04_subj1_V3[["Type"]] <- "Allergic SLIT"
scrna.list$TT04_subj2_V3[["Type"]] <- "Allergic SLIT"
scrna.list$TT04_subj3_V3[["Type"]] <- "Allergic SLIT"
scrna.list$TT04_subj4_V3[["Type"]] <- "Allergic SLIT"
scrna.list$TT04_subj5_V3[["Type"]] <- "Allergic SLIT"
scrna.list$TT04_subj6_V3[["Type"]] <- "Allergic SLIT"


scrna.list$allergic_1[["Group"]] <- "Allergic"
scrna.list$allergic_2[["Group"]] <- "Allergic"
scrna.list$allergic_3[["Group"]] <- "Allergic"
scrna.list$allergic_4[["Group"]] <- "Allergic"

scrna.list$non_allergic_1[["Group"]] <- "Non-Allergic"
scrna.list$non_allergic_2[["Group"]] <- "Non-Allergic"
scrna.list$non_allergic_3[["Group"]] <- "Non-Allergic"
scrna.list$non_allergic_4[["Group"]] <- "Non-Allergic"
scrna.list$non_allergic_5[["Group"]] <- "Non-Allergic"

scrna.list$TT04_subj1_V2[["Group"]] <- "TT-04 baseline"
scrna.list$TT04_subj2_V2[["Group"]] <- "TT-04 baseline"
scrna.list$TT04_subj3_V2[["Group"]] <- "TT-04 baseline"
scrna.list$TT04_subj4_V2[["Group"]] <- "TT-04 baseline"
scrna.list$TT04_subj5_V2[["Group"]] <- "TT-04 baseline"
scrna.list$TT04_subj6_V2[["Group"]] <- "TT-04 baseline"

scrna.list$TT04_subj1_V3[["Group"]] <- "TT-04 1month"
scrna.list$TT04_subj2_V3[["Group"]] <- "TT-04 1month"
scrna.list$TT04_subj3_V3[["Group"]] <- "TT-04 1month"
scrna.list$TT04_subj4_V3[["Group"]] <- "TT-04 1month"
scrna.list$TT04_subj5_V3[["Group"]] <- "TT-04 1month"
scrna.list$TT04_subj6_V3[["Group"]] <- "TT-04 1month"


print("Celldex:NormalizeData")
scrna.list <- lapply(X = scrna.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
})
print("Celldex:as.SingleCellExperiment")
scrna.list.sce <- lapply(X = scrna.list, FUN = function(x) {
    x <- as.SingleCellExperiment(x)
})

print("Predict cell type with Celldex")
ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)

for (i in 1:length(scrna.list.sce)) {
    print(i)
    ens <- mapIds(EnsDb.Hsapiens.v86,keys = rownames(scrna.list.sce[[i]]),column = 'GENEID',keytype = 'SYMBOL')
    all(rownames(scrna.list.sce[[i]]) == names(ens))
    keep <- !is.na(ens)
    ens <- ens[keep]
    scrna.list.sce[[i]] <- scrna.list.sce[[i]][keep,]
    rownames(scrna.list.sce[[i]]) <- ens
    predictions <- SingleR(test=scrna.list.sce[[i]], assay.type.test=1, ref=ref.data, labels=ref.data$label.main)
    print(table(predictions$labels))
    scrna.list[[i]][["predictions"]] <- predictions$labels
}



print("Subset B cells")
scrna.list.Bcell = list()

for (i in 1:length(scrna.list)) {
    print(i)
    meta <- scrna.list[[i]]@meta.data
    meta$cellID <- rownames(meta)
    len <- length(meta$cellID)
    print("Total cells:")    
    print(len)
    Bcells <- meta %>% dplyr::filter(predictions=="B_cell") %>% dplyr::select(cellID)
    len <- nrow(Bcells)
    print("Predicted Bcells:")    
    print(len)
    scrna.list.Bcell[[i]] <- subset(scrna.list[[i]], cells = Bcells$cellID)
}



names(scrna.list.Bcell) <- filenames

print("Merge seurat files")
memB.big <- merge(x=scrna.list.Bcell[[1]], y=scrna.list.Bcell[2:length(scrna.list.Bcell)], add.cell.ids = c(1:21), project = "10xCombi")
sample.list <- SplitObject(memB.big, split.by = "SampleID")

# normalize data with SCTransform()
print("SCTransform")
sample.list <- lapply(X = sample.list, FUN = function(x) {
    x <- SCTransform(x, assay = 'RNA', method = "glmGamPoi", new.assay.name = 'SCT', verbose = TRUE)
})

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Perform cell cycle analysis
print("Perform cell cycle analysis")
sample.list <- lapply(X = sample.list, FUN = function(x) {
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)
})

# Normalize again but this time including also the cell cycle scores,
# notice the use of RNA assay to normalize the original count data

print("SCTransform w regression")
vars_to_regress <- c("nCount_RNA", "S.Score", "G2M.Score", "percent.mt")
#vars_to_regress <- c("nCount_RNA", "percent.mt")
sample.list <- lapply(X = sample.list, FUN = function(x) {
    x <- SCTransform(x, assay = 'RNA', method = "glmGamPoi", new.assay.name = 'SCT', vars.to.regress = vars_to_regress, verbose = TRUE)
})

print("Prep SCT")
features1 <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 1599)

features <- grep(pattern = "^IGHV|^IGK|^IGL|^HLA-D|^MT-", features1,value=TRUE, invert=TRUE) #Filter genes 
sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = features)
sample.list <- lapply(X = sample.list, FUN = RunPCA, features = features)


memB.anchors <- FindIntegrationAnchors(object.list = sample.list, 
    assay = NULL,
    reference = NULL, 
    anchor.features = features, 
    scale = TRUE, 
    normalization.method = c("SCT"),  
    sct.clip.range = NULL,  
    reduction = c("rpca"),  
    l2.norm = TRUE,  
    dims = 1:50,   
    k.anchor = 10,  
    k.filter = 200, 
    k.score = 30, 
    max.features = 200,  
    nn.method = "annoy",  
    n.trees = 50,  
    eps = 0,  
    verbose = TRUE)


print("IntegrateData")
memB.combined2 <- IntegrateData(anchorset = memB.anchors, normalization.method = "SCT",dims=1:50)
#saveRDS(memB.combined2, file="memB.combined2.SCT.2022-01-13.rds")

#Findclusters
print("PCA, tSNE, UMAP etc")
memB.combined2 <- RunPCA(memB.combined2, npcs = 50, verbose = FALSE)
memB.combined2 <- RunTSNE(memB.combined2, dims=1:50)
memB.combined2 <- RunUMAP(memB.combined2, reduction = "pca", dims = 1:50)
memB.combined2 <- FindNeighbors(memB.combined2, reduction = "pca", dims = 1:50)

resolutions <- seq(0, 1, 0.1)
memB.combined2 <- FindClusters(memB.combined2, resolution = resolutions, algorithm = 1)
res <- 0.6
Idents(object = memB.combined2) <- memB.combined2@meta.data[,paste0("integrated_snn_res.",res)]
Idents(memB.combined2) <- factor(x = Idents(memB.combined2), levels = sort(as.numeric(levels(memB.combined2))))

memB.combined2 <- subset(memB.combined2, idents = c(0:23))

DefaultAssay(memB.combined2) <- "RNA"
memB.big <- DietSeurat(memB.combined2, assays = "RNA")
sample.list <- SplitObject(memB.big, split.by = "SampleID")


print("SCTransform w regression")
#vars_to_regress <- c("nCount_RNA", "S.Score", "G2M.Score", "percent.mt")
vars_to_regress <- c("nCount_RNA", "percent.mt")
sample.list <- lapply(X = sample.list, FUN = function(x) {
    x <- SCTransform(x, assay = 'RNA', method = "glmGamPoi", new.assay.name = 'SCT', vars.to.regress = vars_to_regress, verbose = TRUE)
})

print("Prep SCT")
features1 <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 1599)

features <- grep(pattern = "^IGHV|^IGK|^IGL", features1,value=TRUE, invert=TRUE) #Filter genes 

sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = features)
sample.list <- lapply(X = sample.list, FUN = RunPCA, features = features)

memB.anchors <- FindIntegrationAnchors(object.list = sample.list, 
    assay = NULL,
    reference = NULL, 
    anchor.features = features, 
    scale = TRUE, 
    normalization.method = c("SCT"),  
    sct.clip.range = NULL,  
    reduction = c("rpca"),  
    l2.norm = TRUE,  
    dims = 1:50,   
    k.anchor = 10,  
    k.filter = 200, 
    k.score = 30, 
    max.features = 200,  
    nn.method = "annoy",  
    n.trees = 50,  
    eps = 0,  
    verbose = TRUE)


print("IntegrateData")
memB.combined <- IntegrateData(anchorset = memB.anchors, normalization.method = "SCT",dims=1:50)


#Findclusters
print("PCA, tSNE, UMAP etc")
memB.combined <- RunPCA(memB.combined, npcs = 50, verbose = FALSE)
memB.combined <- RunTSNE(memB.combined, dims=1:50)
memB.combined <- RunUMAP(memB.combined, reduction = "pca", dims = 1:50)
memB.combined <- FindNeighbors(memB.combined, reduction = "pca", dims = 1:50)


resolutions <- seq(0, 1, 0.1)
memB.combined <- FindClusters(memB.combined, resolution = resolutions, algorithm = 1)
ResolutionList <- grep("integrated_snn_res", colnames(memB.combined@meta.data), value = TRUE)

for (Resolution in ResolutionList){
    jpeg(paste0(Resolution, "_umap_scaled.jpeg"), width = 4000, height = 4000)
    g <- DimPlot(object = memB.combined, reduction = "umap", group.by = Resolution,pt.size = 2, label =TRUE,label.size = 25)& 
    theme(text = element_text(face = "bold",size=50))
    print(g)
    dev.off()
}


res <- 0.4

pdf("PRetty_umap_scaled2.pdf", width=10, height=7)
    g <- DimPlot(object = memB.combined, reduction = "umap", group.by = paste0("integrated_snn_res.",res), label =TRUE)
    print(g)
    dev.off()

Idents(object = memB.combined) <- memB.combined@meta.data[,paste0("integrated_snn_res.",res)]
Idents(memB.combined) <- factor(x = Idents(memB.combined), levels = sort(as.numeric(levels(memB.combined))))
memB.combined[["seurat_clusters"]] <- Idents(memB.combined)

# Normalize RNA data for visualization purposes
DefaultAssay(memB.combined) <- "RNA"
print("Scale all genes for plotting")
memB.combined <- NormalizeData(memB.combined, verbose = FALSE)
memB.combined <- ScaleData(memB.combined,features = rownames(memB.combined), vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = TRUE)
memB.combined <- FindVariableFeatures(memB.combined,verbose = FALSE)






########## UMAP per SampleID to ensure representation 

DefaultAssay(memB.combined) <- "RNA"
jpeg("UMAP_SampleID.jpg", width = 10000, height = 4000)
g <- DimPlot(object = memB.combined, reduction = "umap", group.by = "seurat_clusters", split.by ="SampleID",pt.size = 2, label =TRUE, combine=TRUE, ncol =7) & 
  theme(text = element_text(face = "bold",size=50))
print(g)
dev.off()




#####################################VDJ###############################################

#Read V(D)J annotation, filter for functional heavy chain and add constant region info to Seurat object meta data.
#source("/home/npkdk/drives/x/Research/NPKDK/10xCombined/VDJ.R")
IGH_info <- data.frame()

print("Read V(D)J annotation")
dir <- "/home/npkdk/drives/x/Research/NPKDK/10xCombined/DATA/cellranger-6.1.2-GRCh38-2020-A/VDJ/"
print(dir)
filenames <- c("TT04_subj1_V2", 
    "TT04_subj1_V3", 
    "TT04_subj2_V2", 
    "TT04_subj2_V3",
    "TT04_subj3_V2", 
    "TT04_subj3_V3", 
    "TT04_subj4_V2", 
    "TT04_subj4_V3", 
    "TT04_subj5_V2", 
    "TT04_subj5_V3", 
    "TT04_subj6_V2", 
    "TT04_subj6_V3",
    "allergic_1", 
    "allergic_2", 
    "allergic_3", 
    "allergic_4",
    "non_allergic_1", 
    "non_allergic_2", 
    "non_allergic_3", 
    "non_allergic_4", 
    "non_allergic_5")


IGH_info1 <- data.frame()
id_index <- 1

for ( index in c(1:21) ) {

    fn <- filenames[index]
    print(fn)
    print(index)
    print(id_index)
    vdj_annot1 <- read.csv(paste(dir,fn,"/all_contig_annotations.csv", sep="")) %>%              
                dplyr::filter(chain=="IGH" & productive=="true" & raw_consensus_id != "") %>% 
                mutate(Sample=fn, sample_index=id_index) %>% 
                dplyr::select(Sample, sample_index, barcode, c_gene, reads) %>% 
                separate(barcode, into=c("cellID", "rest"), sep="-",remove=T) %>%
                unite(col=barcode, sample_index, cellID,  sep="_") %>% 
                dplyr::select(Sample, barcode, c_gene,reads) %>%
                group_by(barcode)%>%
                top_n(1, wt=reads)  
    IGH_info1 <- rbind(IGH_info1, vdj_annot1)
    id_index <- id_index+1
}



IGH_info<- IGH_info1




DefaultAssay(memB.combined) <- "RNA"
tempsid <- as.character(gsub("-1","",colnames(memB.combined)))
memB.combined <- RenameCells(memB.combined, new.names=tempsid
)




################################# Highlight based on isotype ###############################

memB.combined[["barcode"]] <- colnames(memB.combined)
memB.combined[["Isotype"]] <- "None" 

meta <- memB.combined@meta.data
meta$barcode <- colnames(memB.combined)

IGHA1 <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHA1") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHA1] <- "IGHA1"

IGHA2 <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHA2") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHA2] <- "IGHA2"

IGHG1 <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHG1") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHG1] <- "IGHG1"

IGHG2 <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHG2") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHG2] <- "IGHG2"

IGHG3 <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHG3") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHG3] <- "IGHG3"

IGHG4 <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHG4") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHG4] <- "IGHG4"

IGHM <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHM") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHM] <- "IGHM"

IGHD <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHD") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHD] <- "IGHD"

IGHE <- unlist(IGH_info %>% dplyr::filter(c_gene=="IGHE") %>% dplyr::select(barcode), use.names=FALSE)
meta$Isotype[meta$barcode %in% IGHE] <- "IGHE"

#> table(IGH_info$c_gene)
#IGHA1 IGHA2  IGHD  IGHE IGHG1 IGHG2 IGHG3 IGHG4  IGHM  None
#27765  8464   159     8 27797  9676  4196  1016  1106  1092

#> table(meta$Isotype)
#IGHA1 IGHA2  IGHD  IGHE IGHG1 IGHG2 IGHG3 IGHG4  IGHM  None
#23547  7417   116     8 25529  8782  3872   932   891 19224



memB.combined@meta.data <- meta
