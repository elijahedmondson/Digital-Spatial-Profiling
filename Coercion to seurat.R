library(GeomxTools)
library(Seurat)
library(SpatialDecon)
library(patchwork)
#https://satijalab.org/seurat/articles/de_vignette.html
library(DESeq2)
library(MAST)


# load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_geoMX.RData")
# assayDataElementNames(target_myData)
# 
# mySeurat <- as.Seurat.NanoStringGeoMxSet(target_myData, normData = "q_norm")
# mySeurat

load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_seurat.RData")

head(mySeurat, 3)
mySeurat@misc[1:8]
head(mySeurat@misc$sequencingMetrics)# sequencing metrics
head(mySeurat@misc$QCMetrics$QCFlags) # QC metrics
head(mySeurat@assays$GeoMx@meta.features) # gene metadata
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)

mySeurat <- as.Seurat.NanoStringGeoMxSet(target_myData, normData = "q_norm", ident = "progression1")
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)




mySeurat <- FindVariableFeatures(mySeurat)
mySeurat <- ScaleData(mySeurat)
mySeurat <- RunPCA(mySeurat, assay = "GeoMx", verbose = FALSE)
mySeurat <- FindNeighbors(mySeurat, reduction = "pca", dims = seq_len(30))
#mySeurat <- FindClusters(mySeurat, verbose = FALSE)
mySeurat <- RunUMAP(mySeurat, reduction = "pca", dims = seq_len(30))

DimPlot(mySeurat, reduction = "umap", pt.size = 5, label = TRUE, group.by = "dx")


levels(mySeurat)

features <- c("Kras", "Trp53", "Cd274", "Cd8a", "Cd68", "Cre",
         "Krt18", #"Notch1", "Notch2", "Notch3", "Notch4","Cldn8",
         "Cdk6","Msh3","Myc","Mastl", "Sox2","Cav1","Fosl1","Gata4",
         #"Cldn18","Capn6",
         "Cpa1","Muc5ac","Epcam",#"Tff1","Smad4","Sox9",
         #"Ptf1a","Pdx1","Neurog3","Bhlha15","Krt19","Dclk1",
         #"Hnf1b","Krt19","Ngn3","Ctrb1", "Hes1", "Smad4",
         #"Onecut1","Onecut2","Onecut3","Cdkn1a","Prss2","Runx1","Gata6",
         "Gata6", "S100a11", "Nr5a2","Agr2", "Foxa2", "Ets2", "Runx3")

# goi.acini <- c("Ctrb1","Cpa1","Gata6","Bhlha15","Nr5a2","Ptf1a")
# goi.duct <- c("Hnf1b","Sox9","Krt19","Gata6","Onecut1")
# features <- c("Cpa1","Gata6","Sox9","Onecut1","Ngn3","Nr5a2","Ptf1a","Pdx1")
# features <- c("Hes1","Dclk1","Sox9","Gata6","Ptf1a","Pdx1")
# goi.PDAC <- c("Dclk1","Pdx1")
# goi.PDACfromDuct <- "Agr2"
# goi.met <- c(c)
# features <- c("Cpa1","Gata6","Sox9","Onecut1","Ngn3","Nr5a2","Ptf1a","Pdx1")


features <- c("Ctrb1","Cpa1","Bhlha15","Nr5a2","Ptf1a",
              "Dclk1","S100a11","Agr2", "Foxa2","Runx3",
              "Cpa1","Muc5ac","Epcam","Cav1","Fosl1",
              "Cd274", "Cd8a", "Cd68","Cdk6","Msh3",
              "Pdzd8", "Mtch2", "Msln", "Prom1", "Vars2")

RidgePlot(mySeurat, features = features, ncol = 5)


levels(mySeurat)

# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
PanINlo.hi.de <- FindMarkers(mySeurat, ident.1 = "4-PanINlo", ident.2 = "5-PanINhi")


PanINlo.hi.de <- FindMarkers(mySeurat, ident.1 = "4-PanINlo", ident.2 = "5-PanINhi", 
                                  test.use = "MAST")

# test.use = 
# "wilcox"
# "bimod"
# "roc"
# "t"
# "poisson"
# "negbinom"
# "LR"
# "MAST" 
# "DESeq2"




results <- PanINlo.hi.de
library(tibble)
results <- tibble::rownames_to_column(results, "Gene")
results <- dplyr::filter(results, p_val < 0.05)

library(ggrepel) 
# Categorize Results based on P-value & FDR for plotting
results$Color <- "NS or FC < 0.5"
results$Color[results$p_val < 0.05] <- "P < 0.05"
results$Color[results$p_val_adj < 0.05] <- "FDR < 0.05"
results$Color[results$p_val_adj < 0.001] <- "FDR < 0.001"
results$Color[abs(results$avg_log2FC) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$p_val)) * sign(results$avg_log2FC)
top_g <- c()
top_g <- c(top_g,
             head(results[, 'Gene'][order(results[,'invert_P'], decreasing= T)[1:30]]),
             head(results[, 'Gene'][order(results[,'invert_P'], decreasing= F)[1:30]]))
  
top_g <- unique(top_g)


# Graph results
ggplot(results,                                                             ###CHANGE
       aes(x = avg_log2FC, y = -log10(p_val),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "___ <- log2(FC) -> ___",                                       ###CHANGE
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% top_g),# & p_val_adj < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") 






