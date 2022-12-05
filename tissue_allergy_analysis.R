####################################################################################################################
#
# Script: tissue_allergy_analysis.R
# Project: Allergy
# Author: David Glass
# Date: 10-10-22
#
# Purpose: Analysis of MBC2 population
#
#######################################################################################################################



###### LIBRARIES ######

require(uwot)
require(ggplot2)
require(pheatmap)
require(viridis)
require(nationalparkcolors)
require(magrittr)
require(data.table)



##### INPUTS #####

# path to head directory
path <- "~/Research/phd/PheB/"
# path to images directory
images.path <- paste0(path, "images_allergy/")
# path to csv file
csv.path <- paste0(path, "tables/tissue_final_dat.csv")

clusters <- c("Naive", "27- Memory", "RB+ 27+ Memory", "RB- Memory", "95+ Memory",
              "19++ 11c+ Memory", "MBC2", "Plasma", "Germinal Center", "39+ Tonsilar")
cluster.colors <- c("#376186", "#4CA8C9", "#A7738D", "#F0616E", "#4C0099",
                    "#990000", "#B0C551", "#FFADAF", "#55D6BE", "#861657") %>%
  setNames(clusters)
isotypes <- c("IgD", "IgMD", "IgM", "IgG", "IgA", "ND")
isotype.colors <- magma(6) %>%
  rev() %>%
  setNames(isotypes)
isotype.colors["ND"] <- "#D1D2D4"
isotype.channels <- isotypes[c(1,3,4,5)]
tissues <- c("PB", "BM", "T", "LN")
tissue.colors <- park_palette("Everglades", 5)[5:2] %>%
  setNames(tissues)
factors <- c("tissue", "donor", "run", "isotype", "cluster", "meta")


##### FUNCTIONS #####

getDensity <- function(x, y, ...) {
  # Returns the 2D density vector for two numeric vectors
  # Inputs:
  #   x, y - numeric vectors
  #   ... additional arguments passed onto kde2d (probably n, number of grid point in each direction)
  # Outputs:
  #   vector of density values
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


makeCD23Biaxials <- function(dt=mature.dat, ip=images.path) {
  # Makes CD23 Biaxials
  # Inputs:
  #   dt - data.table of mature cells
  #   ip - path to image directory
  # Outputs:
  #   png

  cd27.cutoff <- 0.19
  proportions <- dt[, .N, by=.(gate, CD27>cd27.cutoff)] %>%
    .[, P:=N*100/sum(N)] %>%
    .[gate=="MBC2"]
  
  dt[, density:=getDensity(CD27, CD23, n=100)]
  ggplot(dt, aes(CD27, CD23, color=density)) +
    geom_point(size=0.1) +
    theme_bw() +
    scale_color_viridis() +
    labs(title="Non-naive IgD- B cells") +
    theme(legend.position="none") +
    geom_hline(yintercept=0.32) +
    geom_vline(xintercept=cd27.cutoff) +
    annotate(geom="label",
             label.padding=unit(0.15, "lines"),
             x=0.065,
             y=1.16,
             label=proportions[CD27==F, paste0(round(P, 2), "%")]) +
    annotate(geom="label",
             label.padding=unit(0.15, "lines"),
             x=1.18,
             y=1.16,
             label=proportions[CD27==T, paste0(round(P, 2), "%")])
  ggsave(paste0(ip, "CD27_vs_CD23_biaxial.png"), width=3.5, height=3.5)
  
  rb.cutoff <- 0.5
  proportions <- dt[, .N, by=.(gate, CD45RB>rb.cutoff)] %>%
    .[, P:=N*100/sum(N)] %>%
    .[gate=="MBC2"]
  
  dt[, density:=getDensity(CD45RB, CD23, n=100)]
  ggplot(dt, aes(CD45RB, CD23, color=density)) +
    geom_point(size=0.1) +
    theme_bw() +
    scale_color_viridis() +
    labs(title="Non-naive IgD- B cells") +
    theme(legend.position="none") +
    geom_hline(yintercept=0.32) +
    geom_vline(xintercept=rb.cutoff) +
    annotate(geom="label",
             label.padding=unit(0.15, "lines"),
             x=0.065,
             y=1.16,
             label=proportions[CD45RB==F, paste0(round(P, 2), "%")]) +
    annotate(geom="label",
             label.padding=unit(0.15, "lines"),
             x=1.05,
             y=1.16,
             label=proportions[CD45RB==T, paste0(round(P, 2), "%")])
  ggsave(paste0(ip, "CD45RB_vs_CD23_biaxial.png"), width=3.5, height=3.5)
}


makeAllergyHeatmap <- function(dt=pb.dat, ip=images.path) {
  # Make heatmap of median metacluster expression
  # Inputs:
  #   dt - data.table of PB B cells
  #   ip - path to images directory
  # Outputs:
  #   eps file
  
  medians <- dt[, lapply(.SD, median), .SDcols=all.markers, by=meta]
  med.matrix <- as.matrix(medians, rownames="meta")
  
  postscript(paste0(ip, "heatmap.eps"), width=12, height=4)
  pheatmap(mat=med.matrix,
           color=inferno(n=20),
           border_color=NA,
           angle_col="45")
  dev.off()
}


makeAllergyIsotype <- function(dt=pb.dat, ip=images.path, ic=isotype.colors) {
  # Stacked bars of isotype usage by memory cell population
  # Inputs:
  #   dt - data.table of PB B cells
  #   ip - path to images directory
  #   ic - vector of isotype hex colors
  # Outputs:
  #   eps file
  
  #### By populations
  temp <- dt[!meta %in% c("Naive", "Plasma"), .N, by=.(meta, isotype)] %>%
    .[, P:=N*100/sum(N), by=.(meta)]
  ggplot(temp, aes(meta, P, fill=isotype)) +
    geom_col() +
    theme_bw() +
    scale_fill_manual(values=ic) +
    labs(x=NULL, y="% of population", title="MBC isotype usage") +
    theme(panel.grid.major.x=element_blank(),
          axis.text.x=element_text(hjust=1, angle=45))
  ggsave(paste0(ip, "isotype_bars.eps"), width=3.5, height=3)
}


makeAllergyUmap <- function(dt=pb.dat,
                            ic=isotype.colors,
                            cc=cluster.colors,
                            pa=paste0(images.path, "umap/"),
                            markers=all.markers,
                            iso.channels=isotype.channels) {
  # Generates umap subsampled by tissue and meta
  # Inputs:
  #   dt - data.table
  #   tc - named vector of tissue colors
  #   ic - named vector of isotype colors
  #   cc - named vector of cluster colors
  #   pa - path to images folder
  #   fa - vector of factor column names
  #   iso.channels - vector of isotype column names
  # Outputs:
  #   umap pngs
  if (!dir.exists(pa)) dir.create(pa)
  setorder(dt, CD23)
  
  set.seed(666)
  umap.out <- umap(dt[, setdiff(markers, isotype.channels), with=F], n_neighbors=15, min_dist=0.3)
  dt[, `:=`(umap.1=umap.out[,1], umap.2=umap.out[,2])]
  
  # meta
  ggplot(dt, aes(umap.1, umap.2)) +
    geom_point(size=0.2, aes(color=meta)) +
    theme_bw() +
    scale_color_manual(values=cc) +
    labs(x=NULL, y=NULL) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          legend.position="none")
  ggsave(paste0(pa, "umap_meta.png"), height=7, width=7)
  
  
  # isotype
  ggplot(dt, aes(umap.1, umap.2)) +
    geom_point(size=0.2, aes(color=isotype)) +
    theme_bw() +
    labs(x=NULL, y=NULL) +
    scale_color_manual(values=ic) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          legend.position="none")
  ggsave(paste0(pa, "umap_isotype.png"), width=7, height=7)
  
  
  # marker
  for (marker in markers) {
    dt[eval(parse(text=marker))>1, eval(quote(marker)):=1]
    ggplot(dt, aes(umap.1, umap.2)) +
      geom_point(size=0.2, aes(color=eval(parse(text=marker)))) +
      theme_bw() +
      scale_color_viridis(option="B", limits=c(0,1)) +
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.ticks=element_blank(),
            axis.text=element_blank(),
            legend.position="none")
    ggsave(paste0(pa, "umap_", marker, ".png"), width=7, height=7)
  }
}


makeTissueBars <- function(dt=dat, ip=images.path, tc=tissue.colors) {
  # Makes bars of proportion of MBC2 cells by donor/tissue
  # Inputs:
  #   dt - data.table of all tissues
  #   ip - path to images directory
  #   tc - vector of tissue color hex codes
  # Outputs:
  #   eps files
  
  temp <- dt[, .N, by=.(donor, tissue, meta=="MBC2")] %>%
    .[, P:=N*100/sum(N), by=.(donor, tissue)] %>%
    .[meta==T] %>%
    .[order(P)]
  temp[, donor:=factor(donor, donor)]
  ggplot(temp, aes(donor, P, fill=tissue)) +
    geom_col() +
    theme_bw() +
    scale_fill_manual(values=tc) +
    labs(x="Donor", y="% of B cells", title="MBC2") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major.x=element_blank())
  ggsave(paste0(ip, "tissue_bars.eps"), width=3.5, height=2.5)
}


makeTissueHeatmap <- function(dt=dat, ip=images.path) {
  # Make heatmap of MBC2 by tissue
  # Inputs:
  #   dt - data.table
  #   ip - path to images directory
  # Outputs:
  #   eps
  medians <- dt[meta=="MBC2", lapply(.SD, median), .SDcols=all.markers, by=.(tissue)]
  med.matrix <- as.matrix(medians, rownames="tissue")
  
  postscript(paste0(ip, "tissue_heatmap.eps"), width=12, height=2.7)
  pheatmap(mat=med.matrix,
           color=inferno(n=20),
           border_color=NA,
           angle_col="45")
  dev.off()
}


makeLegends <- function(dt=pb.dat, ip=images.path, cc=cluster.colors[1:8]) {
  # Additional plots to make legends
  # Inputs:
  #   dt - data.table of single cell
  #   ip - path to images directory
  #   cc - vector of cluster hex colors
  # Outputs:
  #   eps and png
  ggplot(dt, aes(tissue, fill=meta)) +
    geom_bar() +
    scale_fill_manual(values=cc)
  ggsave(paste0(ip, "legend_meta.eps"), width=3, height=3)
  
  ggplot(dt, aes(CD45, CD19, color=CD27/max(CD27))) +
    geom_point() +
    scale_color_viridis(option="inferno", name="CD27")
  ggsave(paste0(ip, "legend_inferno.png"), width=3, height=3)
}



##### MAIN #####

if (!dir.exists(images.path)) dir.create(images.path)
dat <- fread(csv.path, stringsAsFactors=T)
dat[meta %in% c("Transitional", "73- Naïve", "73+ Naïve"), meta:="Naive"]
dat[meta %in% c("RB+ 27+ 73- Memory", "RB+ 27+ 73+ Memory"), meta:="RB+ 27+ Memory"]
dat[, `:=`(donor=factor(donor, 1:11),
           isotype=factor(isotype, isotypes),
           tissue=factor(tissue, tissues),
           meta=factor(meta, names(cluster.colors)),
           cluster=factor(cluster, sort(unique(cluster))))]

all.markers <- setdiff(colnames(dat), c(factors, "IgMD", "ND"))

dat[meta!="Naive" & !isotype %in% c("IgD", "IgMD") & CD23>0.32, gate:="MBC2"]
dat[gate=="MBC2", meta:="MBC2"]
mature.dat <- dat[meta!="Naive" & !isotype %in% c("IgD", "IgMD") & tissue=="PB"]
pb.dat <- dat[tissue=="PB" & !meta %in% c("Germinal Center", "39+ Tonsilar")]

### 2A Biaxials
makeCD23Biaxials()

### 2B/D/E UMAP
makeAllergyUmap()

### 2C Heatmap
makeAllergyHeatmap()

### 2D Isotype
makeAllergyIsotype()

### S3B Proportions
makeTissueBars()

### S3C Tissue heatmap
makeTissueHeatmap()

### Additional color legends
makeLegends()
