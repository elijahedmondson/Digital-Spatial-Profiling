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

mySeurat <- as.Seurat.NanoStringGeoMxSet(target_myData, normData = "q_norm", ident = "class")
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)




mySeurat <- FindVariableFeatures(mySeurat)
mySeurat <- ScaleData(mySeurat)
mySeurat <- RunPCA(mySeurat, assay = "GeoMx", verbose = FALSE)
mySeurat <- FindNeighbors(mySeurat, reduction = "pca", dims = seq_len(30))
#mySeurat <- FindClusters(mySeurat, verbose = FALSE)
mySeurat <- RunUMAP(mySeurat, reduction = "pca", dims = seq_len(30))

DimPlot(mySeurat, reduction = "umap", pt.size = 5, label = TRUE, group.by = "comps")


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

#Metastases 1
features <- c("Cybrd1","Nr1d1","Bsg","Tmprss4","Tm9sf3",
              "Mmp23","Rhof","Sftpd", "Aqp5","Ccna1",
              "Muc3","Muc5ac","Muc3a","Kif12","Calml4",
              "Dbp", "Mrtfb", "Rplp0","Dnajc10","Rps12",
              "Pdzd8", "Mtch2", "Msln", "Prom1", "Vars2")
#Metastases 2
features <- c("Porcn","Rpl6","Ybx1","Wfdc2","Tpi1",
              "Golim4","Otop3","F3", "Id2","Adamtsl5",
              "Bag1","Rnf186","Glis2","Slc35f5","Tspan12",
              "Slc9a4", "Ephb2", "Tmem45b","Tmprss2","Pdxdc1",
              "Lgals2", "Esrp1", "Tmem54", "Ptprf", "Ccnd2",
              "Ern2","Sult1c2")
#Metastases 3 / CNS
features <- c("Gltp","Spock3","Sgms2","Rasgrf1","St8sia3",
              "Rap1gap","Rbms3","Ccdc92","Ncald","Ppp1r1b",
              "Gabbr2","Nt5c2","Cdkn2a","Atrnl1","Camk2n1",
              "Setbp1","Dennd4c","Hs3st1","Shf")

#Brain/synaptogenesis/neuronal
features <- c("Nfib", "Tuba1b", "Net1", "Ncald","Spock3",
              "Rock2", "Sem1", "Ctnnd1","Adgre5", "Dennd4c",
              "Smad4", "Flna", "Cntn1", "Cntn6","Sgms2",
              "Nrxn1","Nrxn2","Nrxn3","Lamb2","Rasgrf1",
              "Sema3d", "Sema4b","Sema4g","Sema5a","St8sia3",
              "Lama5", "Rtn4", "Picalm","Efnb2", "Rbms3")

list <- "Kras"
list <- c("Ncald", "St8sia3", "Camk2n1", "Atrnl1", "Spock3", "Rasgrf1")

levels(mySeurat)
levels(x = mySeurat) <- c("Acinar", "PanIN","Carcinoma")
fig <- RidgePlot(mySeurat, sort = F, 
          idents = c("Acinar", "PanIN","Carcinoma"), 
          features = list, ncol = 2)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("fig2.tiff", units="in", width=8, height=8, res=300)
fig
dev.off()

# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes

#NormlizeData() before "FindMarkers()
de_markers <- FindMarkers(mySeurat, ident.1 = "4-PanINlo", ident.2 = "6-PDAC", 
                                  test.use = "negbinom")

de_markers <- FindAllMarkers(mySeurat, ident.1 = "5-PanINhi", ident.2 = "4-PanINlo", 
                             test.use = "negbinom")

# test.use = 
# "wilcox"
# "bimod"
# "roc"
# "t"
# "poisson"
# "negbinom" -- appropriate for count
# "LR"
# "MAST" 
# "DESeq2" -- appropriate for count -- need to initialize the mySeurate with count matrix

# 
# I partly figured this out and thought I would update in case someone else 
# comes looking with the same issue.
# 
# It seems that I was getting crazy LogFC's and weird looking volcanoes because 
# of a problem with the normalization of my data. For these analyses I had used 
# the Seurat "SCT integration" pipeline as outlined here: 
# https://satijalab.org/seurat/v3.2/integration.html
# 
# I probably did something wrong, but after triple checking the workflow, I 
# couldn't find it. When I switched back to an integration workflow that 
# includes log normalization, rather than SCT normalization, everything 
# appears to be fixed.
# 
# So in the end, I am unsure of why the SCT integration pipeline failed 
# on me, but switching to log normalization was the solution to my problem.

results <- de_markers
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
             results[, 'Gene'][order(results[,'invert_P'], decreasing= T)[1:50]],
             results[, 'Gene'][order(results[,'invert_P'], decreasing= F)[1:50]])
  
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






