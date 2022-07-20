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

load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_seurat.RData")

head(mySeurat, 3)
mySeurat@misc[1:8]
head(mySeurat@misc$sequencingMetrics)# sequencing metrics
head(mySeurat@misc$QCMetrics$QCFlags) # QC metrics
head(mySeurat@assays$GeoMx@meta.features) # gene metadata
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 0.1)

mySeurat <- as.Seurat.NanoStringGeoMxSet(target_myData, normData = "q_norm", ident = "progression1")
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)




mySeurat <- FindVariableFeatures(mySeurat)
mySeurat <- ScaleData(mySeurat)
mySeurat <- RunPCA(mySeurat, assay = "GeoMx", verbose = FALSE)
mySeurat <- FindNeighbors(mySeurat, reduction = "pca", dims = seq_len(30))
#mySeurat <- FindClusters(mySeurat, verbose = FALSE)
mySeurat <- RunUMAP(mySeurat, reduction = "pca", dims = seq_len(30))

DimPlot(mySeurat, reduction = "umap", pt.size = 5, label = TRUE, group.by = "progression1")


levels(mySeurat)

features <- c("Kras", "Trp53", "Cd274", "Cd8a", "Cd68", "Epcam","Cre",
         "Krt18", #"Notch1", "Notch2", "Notch3", "Notch4","Cldn8",
         "Cdk6","Msh3","Myc","Mastl", "Sox2","Cav1","Fosl1","Gata4",
         #"Cldn18","Capn6","Cpa1","Muc5ac","Tff1","Smad4","Sox9",
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
# goi.met <- c("Pdzd8", "Mtch2", "Spock3", "Serpina3k", "Cybrd1", "Vars2")
# features <- c("Cpa1","Gata6","Sox9","Onecut1","Ngn3","Nr5a2","Ptf1a","Pdx1")


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

de <- PanINlo.hi.de[complete.cases(PanINlo.hi.de), ]

ggplot(data=PanINlo.hi.de, aes(x=avg_log2FC, y=-log10(p_val))) + 
  geom_point() + 
  theme_minimal() 
+
  geom_text()




suppressPackageStartupMessages({
  library(ggplot2)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
})








