#Outputs ----

#Set Working Path ----
#This path should have both the 10X and SMARTSeq2 data sets.
path <- "E:/USER DATA/Koenig/20221112 Submission Copies/"

# Packages ----
library(Seurat)
library(dplyr) 
library(magrittr)
library(ggplot2)
library(patchwork)
library(viridis)
library(EnhancedVolcano)
library(data.table)
library(pheatmap)
library(purrr)

# Read in data ----
#10x dataset: switched MBCs from allergic and non-allergic participants.
mem <- readRDS(paste0(path,"mem.combined_an.rds"))

#SMARTSeq dataset: Sorted and grass-allergic
ss <- readRDS(paste0(path,"smart.combined.18-10-2022.rds"))


#Assign cluster names for DimPlot

fc <- as.data.table(read.csv("FigureClusters.csv")) #load in names

cl <- as.data.table(mem@meta.data$seurat_clusters) #get seurat clusters
colnames(cl) <- "Cluster" #name column "cluster"

cl <- cl %>%  #R was being annoying and wouldn't merge based on columns so I mutated.
  mutate(Cluster = as.numeric(Cluster)) 

cl$n <- rownames(cl) #I need the cells to be added in otherwise the merge function loses indexes

fc2 <- merge(cl, fc, by = "Cluster") %>%   # fc2 has the metadata, arranged by cell number
  arrange(by = as.numeric(n))

mem <- AddMetaData(mem, fc2$FigureClusters, col.name = "FigureCluster") #apply the metadata

mem <- SetIdent(mem, value = "FigureCluster") #set FigureCluster as the active metadata.

#write.csv(fc2, file = "FigureClusterFinal.csv")

# At this point, we have the data read in with the appropriate metadata active.


#Figure 1 ----

# Figure 1B----
Fig1B <- function(object = mem, ip = path){
  #create plot of selected clusters from scRNA-seq switched MBC atlas.
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #   eps/png file
  
  #Ensure that FigureCluster is the active metadata
  object <- SetIdent(object, value = "FigureCluster")
  
  #custom colours
  cols <- c('#E6E9F4', '#3D4CA0', '#F26B63', "#F9A02F", '#9F276E', '#E864A5','#ABC637')
  
  #order for the data to be plotted by DimPlot
  lev <- levels(object)
  lev <- lev[c(6,2,5,7,4,3,1)]
  
  #Dimplot
  DimPlot(object, reduction = "umap", label = F, pt.size = 1.2, order = lev)+
    scale_color_manual(values = cols)+
    theme(axis.title = element_blank()
          , axis.line=element_blank(),axis.ticks=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none")
  
  #Saving parameters
  ggsave(paste0(ip, "Figure 1/Figure1B.eps"), width = 12, height = 10)
  ggsave(paste0(ip, "Figure 1/Figure1B.png"), width = 12, height = 10)
}

Fig1B()

# Figure 1C ----
Fig1C <- function(object = mem, ip = path){
  #create UMAP plot colorized by MBC2 markers.
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #   4 eps files/ 4 png files
  
  #Feature plot loop for MBC2 markers.
  for(i in c("IGHE", "FCER2", "IL4R", "IL13RA1")){
    FeaturePlot(object = object, features = i, max.cutoff = "q95", pt.size = 1.5, reduction="umap", order = T)+
      scale_color_viridis(option = "C")+
      theme(plot.title = element_blank(), legend.position = 'none', axis.title = element_blank(),
            axis.line=element_blank(),axis.ticks=element_blank(),
            axis.text.x=element_blank(), axis.text.y=element_blank())
    
    #saving parameters
    ggsave(paste0(path,"Figure 1/Figure 1C/",i,".eps"), width = 6, height = 5)
    ggsave(paste0(path,"Figure 1/Figure 1C/",i,".png"), width = 6, height = 5)
  }
}
Fig1C()

# Figure 1D ----
Fig1D <- function(object = mem, ip = path){
  #Heatmap of top 10 upregulated DEG per cluster. 
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #  png file
  
  
  #Ensure that FigureCluster is the active metadata
  object <- SetIdent(object, value = "FigureCluster")
  
  #Find top upregulated DEGs for each cluster
  fam <- FindAllMarkers(object, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
  
  #Get top 10 upregulated DEGs for each cluster
  fam <- fam[fam$cluster != "NR",] %>%
    group_by(cluster)%>%
    slice_max(order_by = avg_log2FC, n = 10)
  
  #Retrieve average expression of each DEG per cluster.
  ae <- AverageExpression(object, assays = "RNA", features = unique(fam2$gene))
  
  #Create output heatmap.
  png(paste0(ip,"Figure 1/Figure 1D.png"), width = 19, height = 3, units = "in", res = 720)
  pheatmap(t(ae$RNA[,-1]), scale = "column", col = viridis(100),
           cellwidth = 25, cellheight = 25, angle_col = "45", legend = F,
           treeheight_col = 0, treeheight_row = 10, fontsize = 12)
  dev.off()
}
Fig1D()

# Figure 1F ----
Fig1F <- function(object = mem, ip = path){
  #create heatmap of Euclidean distances
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #   eps/png file
  
  #Ensure that FigureCluster is the active metadata
  object <- SetIdent(object, value = "FigureCluster")
  
  #Get average expression by cluster for all genes.
  ae <- AverageExpression(object, assay = "RNA", features = rownames(mem@assays$RNA))
  
  #Calculate Euclidean distances between clusters.
  euc.dist <- dist(t(ae$RNA[,-1]), method = "euclidean", diag = T, upper = T)
  
  #Output png
  png(paste0(path, "Figure 1/Figure 1F.png"), width = 300, height = 300)
  pheatmap(euc.dist, color = viridis(50), labels_row = rownames(as.matrix(euc.dist)), labels_col = colnames(as.matrix(euc.dist)), legend = F,
           fontsize = 11)
  dev.off()
  
  #Output eps
  postscript(paste0(path, "Figure 1/Figure 1F.eps"), width=4, height=4)
  pheatmap(euc.dist, color = viridis(50), labels_row = rownames(as.matrix(euc.dist)), labels_col = colnames(as.matrix(euc.dist)), legend = F,
           fontsize = 11)
  dev.off()
}

Fig1F()

# Figure 1I ----
Fig1I <- function(object = mem, ip = path){
  #Create plot of frequency of isotypes among clusters as determined by VDJ library
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #   eps/png file
  
  #Math for calculating frequencies
  fq <- count(object@meta.data, FigureCluster, Isotype)
  fq <- merge(fq, as.data.table(object@meta.data)[,.N, by= FigureCluster], by = "FigureCluster")
  fq$P <- fq$n/fq$N*100
  
  #Remove NR cluster
  fq <- fq[fq$FigureCluster!="NR",]
  
  #Reassign factor levels for ggplot representation
  fq$FigureCluster <- factor(fq$FigureCluster, levels = c("CD19hi CD11chi", "CD95+", "AREG+", "IGHE+ MBC2", "IGHE- MBC2", "IgM+ IgD+"))
  fq$Isotype <- factor(fq$Isotype, levels = c("IGHA2","IGHA1","IGHE","IGHG4","IGHG3","IGHG2","IGHG1", "IGHD", "IGHM", "None"))
  
  #Custom color palette
  cc <- c("#fdcd26", #yellow
          "#F9A02F", #orange
          "#c0df25", #green
          "#42BB8F", #seafoam
          "#2A7D7F", #teal
          "#6B7DBD", #periwinkle
          "#293691", #blue
          "#6B3A96", #violet
          "#442C56", #purple
          '#E6E9F4') #grey
  
  #Plot
  ggplot(fq, aes(x=P, y = FigureCluster, fill = Isotype))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = cc)+
    xlab("% of Total")+
    ylab("Cluster")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
  
  #Saving parameters
  ggsave(paste0(path, "Figure 1/Figure1I.png"), width = 6, height = 3)
  ggsave(paste0(path, "Figure 1/Figure1I.eps"), width = 6, height = 3)
  
}

Fig1I()

# Figure 1J ----
Fig1J<- function(object = mem, ip = path){
  #create UMAP plot colorized by MBC2 markers.
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #   eps/png file
  
  #Ensure that Isotype from the VDJ library is the active metadata
  object <- SetIdent(object, value = "Isotype")
  
  #order for the data to be plotted by DimPlot
  lev <- levels(object)
  lev <- lev[c(8,1:7,9,10)]
  
  #custom colours
  cols <- c('#E6E9F4','#E6E9F4','#E6E9F4','#E6E9F4','#E6E9F4','#E6E9F4','#E6E9F4','#E6E9F4','#E6E9F4','#42BB8F')
  
  #DimPlot
  DimPlot(object, reduction = "umap", label = F, pt.size = 1.2, order = lev)+
    scale_color_manual(values = cols)+
    theme(axis.title = element_blank(), axis.line=element_blank(),axis.ticks=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none")
  
  #saving parameters
  ggsave(paste0(ip, "Figure 1/Figure1J.eps"), width = 6, height = 5)
  ggsave(paste0(ip, "Figure 1/Figure1J.png"), width = 6, height = 5)
}

Fig1J(mem, path)
# Figure 2 ----
# Figure 2B RNA ----
Fig2BRNA <- function(object = mem, ip = path){
  # Make RNA heatmap of mean cluster expression to match mass cytometry data
  # Inputs:
  #   object - SeuratObject
  #   path - path to output images directory
  # Outputs:
  #   png file
  
  #Gene annotations for markers in Glass et al CyTOF experiment.
  f <- c("CD24", "CD1C", "IGHM", "CCR7", "CXCR4", "CR2", "CXCR5", "CD40", "CD79B", 
         "CD72", "IGHD", "LAIR1", "ITGAX", "CD22", "CXCR3", "CD38", "ENTPD1", "FCER2", 
         "CD9", "NT5E", "PTPRC", "FCGR2B", "MS4A1", "CD19", "CD27", "FAS")
  
  #Get average expression for the defined markers.
  ae1 <- AverageExpression(object, features = f, assay = "RNA")
  
  #below: surface stains for IgA, IgG, HLA-DR, and light chain used reagents that bind proteins that are encoded by multiple genes.
  #to approximate an equivalent to the surface stain, we took the average of expression of the multiple genes.
  IGHA <- colMeans(AverageExpression(object, features = c("IGHA1", "IGHA2"), assay = "RNA")$RNA) 
  IGHG <- colMeans(AverageExpression(object, features = c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), assay = "RNA")$RNA)
  HLA <- colMeans(AverageExpression(object, features = grep("HLA-DR", rownames(mem@assays$RNA@counts), ignore.case = T, value = T ), assay = "RNA")$RNA) 
  IGL <- colMeans(AverageExpression(object, features = c("IGKC", grep("IGLC", rownames(mem@assays$RNA@counts), ignore.case = T, value = T )), assay = "RNA")$RNA)
  
  #Assemble the data in order of the protein heatmap, with appropriate feature labels.
  ae2 <- rbind(ae1$RNA[1:17,], IGHA, ae1$RNA["FCER2",],IGHG, ae1$RNA[19:23,], HLA, ae1$RNA[24,],IGL, ae1$RNA[25:26,])
  rownames(ae2) <- c("CD24", "CD1C", "IGHM", "CCR7", "CXCR4", "CR2", "CXCR5", "CD40", "CD79B", "CD72", "IGHD", "LAIR1", "ITGAX", "CD22", "CXCR3",
                     "CD38", "ENTPD1", "IGHA", "FCER2", "IGHG", "CD9", "NT5E", "PTPRC", "FCGR2B", "MS4A1", "HLA-DR", "CD19", "IGL", "CD27", "FAS")
  ae2 <- ae2[,c(4,2,3,6,7,5)]
  
  #Export heatmap at path in png format
  #png(paste0(path, "Figure 2/Figure 2B_RNA.png"), width=11, height=3, units = "in", res = 720)
  postscript(paste0(ip, "Figure 2/ps.eps"), width=9.5, height=2.5)
  pheatmap(mat=t(ae2),
           scale = "column",
           color=inferno(n=20),
           border_color=NA,
           angle_col="45",
           cluster_cols = F,
           cluster_rows = F,
           legend = F)
  dev.off()
}

Fig2BRNA(object = mem, ip = path)

# Figure 2F ----
  #Make GO figure

#Figure 2 G----
markers <- c("HOPX", "JUN", "JUNB", "JUND", "FOS", "FOSB")
Fig2G <- function(object = mem, ip = path, features = markers){
  #Create Violin plots with white dot for mean, this is done for HOPX and AP-1 transcription factors FOS and JUN.
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #   markers - features to plot
  #Outputs:
  #   eps/png files
  
  #Specify idents to show
  id <- c("CD19hi CD11chi", "AREG+","CD95+", "IgM+ IgD+", "IGHE- MBC2", "IGHE+ MBC2")
  
  #Specify the order to show the idents
  
  object@meta.data$FigureCluster <- factor(object@meta.data$FigureCluster, c("AREG+","CD95+","IGHE- MBC2", "IGHE+ MBC2", "IgM+ IgD+", "CD19hi CD11chi", "NR"))
  
  #Ensure that the appropriate ident is active.
  object <- SetIdent(object, value = "FigureCluster")
  
  #Specify custom colors for the idents
  cols <- c('#F26B63', #Salmon
            '#E864A5', #Pink
            '#3D4CA0', #Blue
            '#ABC637', #Green
            "#F9A02F", #Orange
            '#9F276E') #Maroon
  
  #Loop across all the specified markers
  for(i in markers){
    #Create Violin plot using Seurat
    VlnPlot(object = object, features = i, pt.size = 0, idents = id, cols = cols)+
      stat_summary(fun=mean, geom="point", shape=23, size=2, color="black", fill="white")+
      theme( axis.title = element_blank(), axis.text.x = element_blank(), 
             plot.title = element_blank(), legend.position = "none")
    
    #Saving parameters for both png and EPS.
    ggsave(filename = paste0(ip,"Figure 2/Figure 2G/",i,".png"), width = 3.2, height = 2.4)
    ggsave(filename = paste0(ip,"Figure 2/Figure 2G/",i,".eps"), width = 3.2, height = 2.4)
  }
}
Fig2G(mem, path, markers)

#Figure 3 ----
#Figure 3G ----
Fig3FHPreProcess <- function(){
# Required to subset and appropriately label antigen-binding cells.
# This will be replaced following final meta data assembly.
# For now, need to separate antigen-specific cells with a metadata line that
# allows to compare Allergen vs RBD.
 
ss <- SetIdent(ss, value = "sorted") #Ensure sorted is the active metadata

ss2 <- subset(x = ss, idents = unique(ss@meta.data$sorted)[c(1,2,4,5,6)]) #Subset out based on what specificity they were sorted against.
  # NA = a supplemented dataset from Croote et al. that was used as a reference during analysis.

ss2@meta.data$Specific2[is.na(ss2@meta.data$Specific2)] = "RBD" #NA cells in the Specific2 line are from RBD-sorted experiment.
ss2 <- SetIdent(ss2, value = "Specific2") #Set this metadata as active.

ss3 <- subset(x = ss2, idents = unique(ss2@meta.data$Specific2)[c(1,4:6)]) #Select cells binding Birch (Bet), Oak (Que), Grass (Phl), and RBD.

ss3@meta.data$Specific[is.na(ss3@meta.data$Specific)] = "RBD" #the Specific metadata lists whether or not the cells bind allergen.
#used this metadata to simply assemble the NAs as RBD and the "YES" cells as Allergen binding.
ss3@meta.data$Specific[ss3@meta.data$Specific == "Yes"] = "Allergen"
ss3 <- SetIdent(ss3, value = "Specific")
}

ss2 <- Fig3GHPreProcess() #this can be used as the object for creating Isotype and cluster output.

Fig3F <- function(object = ss2, ip = path){
  # Make plot of cluster identity for Smartseq sorted cells
  # Inputs:
  #   object - SeuratObject
  #   path - path to output images directory
  # Outputs:
  #   png and eps file
  
  # Calculate percentage of each cluster by antigen specificity (Allergen vs RBD)
  ss.cl <- count(object@meta.data, Specific, predicted.cluster)
  ss.cl <- merge(ss.cl, count(ss2@meta.data, Specific), by = "Specific")
  ss.cl$P <- ss.cl$n.x/ss.cl$n.y*100
  
  # Apply factor to predicted cluster so it will appear in order on the plot
  ss.cl$predicted.cluster <- factor(ss.cl$predicted.cluster, sort(as.numeric(unique(ss.cl$predicted.cluster)), decreasing = T))
  
  # select custom colors aligned with the factors above.
  cols <- c("#ABC637","#aab4da", "#cdd3e9", "#E6E9F4", "#3D4CA0", "#E864A5")
  
  #Plot it
  ggplot(ss.cl, aes(P, Specific, fill = predicted.cluster))+
    geom_bar(position = "stack", stat = "identity")+
    scale_fill_manual(values = cols)+
    xlab("Percent")+
    theme_bw()+
    theme(axis.title.y = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  #Saving parameters
  ggsave(paste0(ip, "Figure 3/Figure3F.png"), height = 2, width = 3.75)
  ggsave(paste0(ip, "Figure 3/Figure3F.eps"), height = 2, width = 3.75)
}

Fig3F(object = ss2, ip = path)

Fig3H <- function(object = ss2, ip = path){
  # Make plot of isotype expression from RNA.
  # Inputs:
  #   object - SeuratObject
  #   path - path to output images directory
  # Outputs:
  #   png and eps file
  
  # Get expression of each isotype
  ss.iso <- object@assays$RNA@data[c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2"),]

  # Identify row with max expression and rename columns.
  dt <- as.data.table(cbind(ss2$Specific,rownames(ss.iso)[apply(ss.iso, which.max, MARGIN = 2)]))
  colnames(dt)<- c("Specific", "Isotype")

  # Do the math to get percentages of each isotype per specificity (allergen vs RBD).
  dt2 <- count(dt, Specific, Isotype)
  dt2 <- merge(dt2, dt[,.N,by = Specific], by = "Specific")
  dt2$P <- dt2$n/dt2$N*100
  
  # Provide factor levels to the isotypes so they come out correctly.
  dt2$Isotype <- factor(dt2$Isotype,levels = c("IGHA2","IGHA1","IGHG4","IGHG3","IGHG1", "IGHM"))

  # Custom colors aligned with factors above.
  cc <- c("#fdcd26", #yellow
        "#F9A02F", #orange
        "#42BB8F", #seafoam
        "#293691", #blue
        "#442C56") #purple

  # Plot it
  ggplot(dt2, aes(P, Specific, fill = Isotype))+
    geom_bar(position = "stack", stat = "identity")+
    scale_fill_manual(values = cc)+
    xlab("Percent")+
    theme_bw()+
    theme(axis.title.y = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  #Saving parameters
  ggsave(paste0(ip, "Figure 3/Figure3H.png"), height = 2, width = 3.75)
  ggsave(paste0(ip, "Figure 3/Figure3H.eps"), height = 2, width = 3.75)
}

Fig3H(object = ss2, ip = path)

#Figure 4 ----
#Figure 4B ----
#This code works but requires editing.
ccmd <- read.csv("E:/USER DATA/Koenig/Human/meta_ige_cocluster.csv")

mem@meta.data$ClonalRelative <- ccmd$cocluster

mem <- SetIdent(mem, value = "ClonalRelative")

m1 <- as.data.table(mem[,mem@meta.data$ClonalRelative == "no"]@reductions$umap@cell.embeddings)
m1$CR <- "no"

m2 <- as.data.table(mem[,mem@meta.data$ClonalRelative == "yes"]@reductions$umap@cell.embeddings)
m2$CR <- "yes"

m1 <- rbind(m1, m2)

ggplot(m1, aes(UMAP_1, UMAP_2, color = CR, size = CR))+
  geom_point()+
  scale_color_manual(values = c("#E6E9F4", "#293691"))+
  scale_size_manual(values = c(1.2,4))+
  theme_void()+
  theme(axis.title = element_blank()
        , axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none")

ggsave("E:/USER DATA/Koenig/20221112 Submission Copies/Figure 4/Figure 4B.png", width = 12, height = 10)
ggsave("E:/USER DATA/Koenig/20221112 Submission Copies/Figure 4/Figure 4B.eps", width = 12, height = 10)

#Figure 4D ----
mem@meta.data$Group <- ccmd$Group
mem@meta.data$ClonalRelative <- ccmd$cocluster

mem <- SetIdent(mem, value = "Group")

m1 <- SplitObject(mem, split.by = "Group")
bl1 <- as.data.table(m1$`TT-04 baseline`[,m1$`TT-04 baseline`@meta.data$ClonalRelative == "no"]@reductions$umap@cell.embeddings)
bl1$CR <- "no"

bl2 <- as.data.table(m1$`TT-04 baseline`[,m1$`TT-04 baseline`@meta.data$ClonalRelative == "yes"]@reductions$umap@cell.embeddings)
bl2$CR <- "yes"

bl1 <- rbind(bl1, bl2)

ggplot(bl1, aes(UMAP_1, UMAP_2, color = CR, size = CR))+
  geom_point()+
  scale_color_manual(values = c("#E6E9F4", "#293691"))+
  scale_size_manual(values = c(2,6))+
  theme_void()+
  theme(axis.title = element_blank()
        , axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none")

ggsave("E:/USER DATA/Koenig/20221112 Submission Copies/Figure 4/Figure 4C Baseline.png", width = 12, height = 10)
ggsave("E:/USER DATA/Koenig/20221112 Submission Copies/Figure 4/Figure 4C Baseline.eps", width = 12, height = 10)


om1 <- as.data.table(m1$`TT-04 1month`[,m1$`TT-04 1month`@meta.data$ClonalRelative == "no"]@reductions$umap@cell.embeddings)
om1$CR <- "no"

om2 <- as.data.table(m1$`TT-04 1month`[,m1$`TT-04 1month`@meta.data$ClonalRelative == "yes"]@reductions$umap@cell.embeddings)
om2$CR <- "yes"

om1 <- rbind(om1, om2)

ggplot(om1, aes(UMAP_1, UMAP_2, color = CR, size = CR))+
  geom_point()+
  scale_color_manual(values = c("#E6E9F4", "#293691"))+
  scale_size_manual(values = c(2,6))+
  theme_void()+
  theme(axis.title = element_blank()
        , axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none")

ggsave("E:/USER DATA/Koenig/20221112 Submission Copies/Figure 4/Figure 4C 1month.png", width = 12, height = 10)
ggsave("E:/USER DATA/Koenig/20221112 Submission Copies/Figure 4/Figure 4C 1month.eps", width = 12, height = 10)

# Figure 4E----

mem <- SetIdent(mem, value = "ClonalRelative")
m1 <- subset(mem, idents = "yes") 

#Math for calculating frequencies
m1 <- as.data.table(m1@meta.data)
m1 <- as.data.table(cbind(m1[,.N, by = "Isotype"],
                          m1[,.N, by = "Isotype"] %>%
                            .$N/nrow(m1)*100))
colnames(m1) <- c("Isotype", "N", "P")
m1$group <- "group"

#Reassign factor levels for ggplot representation
m1$Isotype <- factor(m1$Isotype, levels = c("IGHA2","IGHA1","IGHE","IGHG4","IGHG3","IGHG2","IGHG1", "None"))

#Custom color palette
cc <- c("#fdcd26", #yellow
        "#F9A02F", #orange
        "#c0df25", #green
        "#42BB8F", #seafoam
        "#2A7D7F", #teal
        "#6B7DBD", #periwinkle
        "#293691", #blue
        '#E6E9F4') #grey

#Plot
ggplot(m1, aes(x= group, y = P, fill = Isotype))+
  geom_bar(stat = "identity")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = cc)+
  ylab("% of Total")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

#Saving parameters
ggsave(paste0(path, "Figure 4/Figure4E.png"), width = 1.2, height = 3)
ggsave(paste0(path, "Figure 4/Figure4E.eps"), width = 1.2, height = 3)

# Figure 4F----

mem <- SetIdent(mem, value = "ClonalRelative")
m1 <- subset(mem, idents = "yes") 

#Math for calculating frequencies
m1 <- as.data.table(m1@meta.data)
m1 <- as.data.table(cbind(m1[,.N, by = "seurat_clusters"],
                          m1[,.N, by = "seurat_clusters"] %>%
                            .$N/nrow(m1)*100))%>%
  arrange(desc(seurat_clusters))
colnames(m1) <- c("Cluster", "N", "P")
m1$group <- "group"

m2 <- as.data.table(cbind(0:13,c("Other", "CD95+", "AREG+", "IGHE- MBC2", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "IGHE+ MBC2", "Other")))
colnames(m2) <- c("Cluster", "FigureCluster")

m1 <- merge(m1,m2, by = "Cluster")

#Custom color palette
cc <- c('#F26B63', #Samlon  AREG+
        '#E864A5', #Pink    CD95+
        '#3D4CA0', #Blue IGHE- MBC2
        '#ABC637', #Green IGHE+ MBC2
        '#E6E9F4') #Grey Other

#Plot
ggplot(m1, aes(x= group, y = P, fill = FigureCluster))+
  geom_bar(stat = "identity")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = cc)+
  ylab("% of Total")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

#Saving parameters
ggsave(paste0(path, "Figure 4/Figure4F.png"), width = 1.2, height = 3)
ggsave(paste0(path, "Figure 4/Figure4F.eps"), width = 1.2, height = 3)

#Figure 5 ----
#Figure 5I ----
# Read in `matrix.mtx`
counts <- readMM("filtered/matrix.mtx")
#counts <- readMM("X:\\Research\\IHO\\IP-4\\ELN110697_SmartSeq\\STAR_solo\\filtered\\matrix.mtx")

# Read in `genes.tsv`
genes <- read_tsv("filtered/features.tsv", col_names = FALSE)
#genes <- read_tsv("X:\\Research\\IHO\\IP-4\\ELN110697_SmartSeq\\STAR_solo\\filtered\\features.tsv", col_names = FALSE)
gene_ids <- genes$X2

# Read in `barcodes.tsv`
cell_ids <- read_tsv("filtered/barcodes.tsv", col_names = FALSE)$X1
#cell_ids <- read_tsv("X:\\Research\\IHO\\IP-4\\ELN110697_SmartSeq\\STAR_solo\\filtered\\barcodes.tsv", col_names = FALSE)$X1


# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids

# Create Seurat Object
smartseq <- CreateSeuratObject(counts = counts, project = "smartseq")
smartseq@meta.data[["cellid"]] <- colnames(smartseq)
smartseq[["group"]] <- "OVA+"

# Annotations. Figure this out before publication.
Unswitch_OVApos <- c("A08", "A09", "A11")
smartseq@meta.data[Unswitch_OVApos,"group"] <- "Unswitch_OVApos"
IgMpos_OVApos_MBC <- c("B03", "B04", "B05", "B08", "B09", "B11")
smartseq@meta.data[IgMpos_OVApos_MBC,"group"] <- "IgMpos_OVApos_MBC"
OVAneg_MBC <- c("E03", "F02", "F03", "F04", "F05", "G02", "G03", "G04", "G05", "H03")
smartseq@meta.data[OVAneg_MBC,"group"] <- "OVA-"

Idents(smartseq) <- "group"

# Cleaning
smartseq[["percent.mt"]] <- PercentageFeatureSet(smartseq, pattern = "^mt-")
smartseq[["SampleID"]] <- "smartseq"
smartseq[["ProjectID"]] <- "smartseq"
smartseq[["Visit"]] <- "None"

# Features and mitochondrial DNA filtering
FeatureScatter(smartseq, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 3)
plot1 <- FeatureScatter(smartseq, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 3)
plot2 <- FeatureScatter(smartseq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 3)
plot1 + plot2
smartseq <- subset(smartseq, subset = nFeature_RNA > 200 & percent.mt < 8)

#Normalize
smartseq <- NormalizeData(smartseq, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

#Subset OVA+ and OVA- 
a <- subset(smartseq, idents = c("OVA+","OVA-"))

#Select markers to make plots of
markers <- c("Fcer2a", "Il4ra", "Fcgr2b", "Hopx", "Cd40", "Ciita", "Il13ra1")

#Create and save Violin plots
for(i in markers){
  VlnPlot(a, i, pt.size = 0)+
    scale_fill_manual(values = c("#6B7DBD","#F9A02F"))+
    stat_summary(fun=median, geom="point", shape=23, size=2, color="black", fill="white")+
    theme(legend.position = "none", title = element_blank(), plot.title = element_blank(), axis.text.x = element_blank())
  ggsave(paste0("E:/USER DATA/Koenig/20221112 Submission Copies/Figure 5/Figure 5I/",i,".png"), width = 2, height = 2.5)
  ggsave(paste0("E:/USER DATA/Koenig/20221112 Submission Copies/Figure 5/Figure 5I/",i,".eps"), width = 2, height = 2.5)
}

# Figure S1A/S1B ----
FigS1AB <- function(object = mem, ip = path){
  
  #Ensure that seurat_clusters is the active metadata
  object <- SetIdent(object, value = "seurat_clusters") #set seurat_clusters as the active metadata.
  
  #DimPlot
  DimPlot(object = object, reduction = "umap", label = T, pt.size = 1, label.size = 8)+
    theme(axis.title = element_blank()
          , axis.line=element_blank(),axis.ticks=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none")
  ggsave(paste0(ip, "Figure S1/FigureS1A.eps"), width = 12, height = 10)
  ggsave(paste0(ip, "Figure S1/FigureS1A.png"), width = 12, height = 10)
  
  #Non-Allergic
  DimPlot(object = object[,object@meta.data$Type == "Non-Allergic"], reduction = "umap", label = T, pt.size = 1, label.size = 8)+
    theme(axis.title = element_blank()
          , axis.line=element_blank(),axis.ticks=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none")
  ggsave(paste0(ip, "Figure S1/FigureS1BNonAllergic.eps"), width = 12, height = 10)
  ggsave(paste0(ip, "Figure S1/FigureS1BNonAllergic.png"), width = 12, height = 10)
  
  #Allergic
  DimPlot(object = object[,object@meta.data$Type != "Non-Allergic"], reduction = "umap", label = T, pt.size = 1, label.size = 8)+
    theme(axis.title = element_blank()
          , axis.line=element_blank(),axis.ticks=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position = "none")
  ggsave(paste0(ip, "Figure S1/FigureS1BAllergic.eps"), width = 12, height = 10)
  ggsave(paste0(ip, "Figure S1/FigureS1BAllergic.png"), width = 12, height = 10)
}
FigS1ABDimPlots(object = mem)


#Table S1----
TableS1 <- function(object = mem, ip = path){
  
  #Ensure that seurat_clusters is the active metadata
  object <- SetIdent(object, value = "seurat_clusters")
  
  #Find all DEGs, positive and negative
  fam <- FindAllMarkers(object, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Write CSV file.
  write.csv(fam, file = paste0(ip,"Table S1/Table S1.csv"))
}
TableS1()

#THIS IS INCOMPLETE -----------
FigS2B <- function(object = mem, ip = path){
  #Heatmap of top 10 upregulated DEG per cluster. 
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #  png file
  
  
  #Ensure that FigureCluster is the active metadata
  object <- SetIdent(object, value = "seurat_clusters")
  
  #Find top upregulated DEGs for each cluster
  fam <- FindAllMarkers(object, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
  
  #Get top 10 upregulated DEGs for each cluster
  fam <- fam %>%
    group_by(cluster)%>%
    slice_max(order_by = avg_log2FC, n = 10)
  
  #Retrieve average expression of each DEG per cluster.
  ae <- AverageExpression(object, assays = "RNA", features = unique(fam$gene))
  
  #Create output heatmap.
  png(paste0(ip,"Figure S2/FigureS2B.png"), width = 5, height = 19, units = "in", res = 720)
  pheatmap(ae$RNA[,-1], scale = "row", col = viridis(100),
           cellwidth = 25, cellheight = 25, angle_col = "45", legend = F,
           treeheight_col = 0, treeheight_row = 10, fontsize = 12)
  
  dev.off()
}
FigS2B()


# Figure S2C ----
FigS2C <- function(object = mem, ip = path){
  #create heatmap of Euclidean distances
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #   eps/png file
  
  #Ensure that FigureCluster is the active metadata
  object <- SetIdent(object, value = "seurat_clusters")
  
  #Get average expression by cluster for all genes.
  ae <- AverageExpression(object, assay = "RNA", features = rownames(mem@assays$RNA))
  
  #Calculate Euclidean distances between clusters.
  euc.dist <- dist(t(ae$RNA[,-1]), method = "euclidean", diag = T, upper = T)
  
  #Output png
  png(paste0(ip, "Figure S2/Figure S2C.png"), width = 500, height = 500)
  pheatmap(euc.dist, color = viridis(50), labels_row = rownames(as.matrix(euc.dist)), labels_col = colnames(as.matrix(euc.dist)), legend = F,
           fontsize = 11)
  dev.off()
  
  #Output eps
  postscript(paste0(ip, "Figure S2/Figure S2C.eps"), width=8, height=8)
  pheatmap(euc.dist, color = viridis(50), labels_row = rownames(as.matrix(euc.dist)), labels_col = colnames(as.matrix(euc.dist)), legend = F,
           fontsize = 11)
  dev.off()
}

FigS2C()

# Figure S2D ----
FigS2D <- function(object = mem, ip = path){
  #Create plot of frequency of isotypes among clusters as determined by VDJ library
  #Inputs:
  #   object - SeuratObject
  #   ip - path for output images directory
  #Outputs:
  #   eps/png file
  
  #Math for calculating frequencies
  fq <- count(object@meta.data, seurat_clusters, Isotype)
  fq <- merge(fq, as.data.table(object@meta.data)[,.N, by= seurat_clusters], by = "seurat_clusters")
  fq$P <- fq$n/fq$N*100
  
  #Reassign factor levels for ggplot representation
  fq$Isotype <- factor(fq$Isotype, levels = c("IGHA2","IGHA1","IGHE","IGHG4","IGHG3","IGHG2","IGHG1", "IGHD", "IGHM", "None"))
  fq$seurat_clusters <- factor(fq$seurat_clusters, levels = 20:0)
  
  #Custom color palette
  cc <- c("#fdcd26", #yellow
          "#F9A02F", #orange
          "#c0df25", #green
          "#42BB8F", #seafoam
          "#2A7D7F", #teal
          "#6B7DBD", #periwinkle
          "#293691", #blue
          "#6B3A96", #violet
          "#442C56", #purple
          '#E6E9F4') #grey
  
  #Plot
  ggplot(fq, aes(x=P, y = seurat_clusters, fill = Isotype))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = cc)+
    xlab("% of Total")+
    ylab("Cluster")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
  
  #Saving parameters
  ggsave(paste0(ip, "Figure S2/FigureS2D.png"), width = 6, height = 3)
  ggsave(paste0(ip, "Figure S2/FigureS2D.eps"), width = 6, height = 3)
  
}

FigS2D()