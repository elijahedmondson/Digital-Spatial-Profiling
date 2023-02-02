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
#save(mySeurat, target_myData, as.Seurat.NanoStringGeoMxSet, file="C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_seurat_2023.RData")

#load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_seurat.RData")
load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_seurat_2023.RData")

head(mySeurat, 3)
mySeurat@misc[1:8]
head(mySeurat@misc$sequencingMetrics)# sequencing metrics
head(mySeurat@misc$QCMetrics$QCFlags) # QC metrics
head(mySeurat@assays$GeoMx@meta.features) # gene metadata
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)

mySeurat <- as.Seurat.NanoStringGeoMxSet(target_myData, normData = "q_norm", ident = "classes")
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)




mySeurat <- FindVariableFeatures(mySeurat)
mySeurat <- ScaleData(mySeurat)
mySeurat <- RunPCA(mySeurat, assay = "GeoMx", verbose = FALSE)
mySeurat <- FindNeighbors(mySeurat, reduction = "pca", dims = seq_len(30))
#mySeurat <- FindClusters(mySeurat, verbose = FALSE)
mySeurat <- RunUMAP(mySeurat, reduction = "pca", dims = seq_len(30))

DimPlot(mySeurat, reduction = "umap", pt.size = 5, label = TRUE, group.by = "classes")


levels(mySeurat)
levels(x = mySeurat) <- c("Normal", "ADM","PanIN","PDAC", "Lung_met", "Liver_met", "Stroma")


features <- c("Kras","Trp53","Cre","Pdx1",
              "Pdia2","Cel", "Reg1","Pnliprp1","Try4",
              "Hnf1b","Sox9","Krt19","Onecut1")#
         #"Krt18", "Notch1", "Notch2", "Notch3", "Notch4","Cldn8","Cd274", "Cd8a", "Cd68",
         #"Cdk6","Msh3","Mastl", "Sox2","Cav1","Fosl1","Gata4",
         #"Cldn18","Capn6",
         #"Cpa1","Muc5ac","Epcam","Tff1","Sox9",
         #"Ptf1a","Pdx1","Neurog3","Krt19","Dclk1",
         #"Hnf1b","Krt19","Ngn3","Ctrb1", "Hes1", "Smad4",
         #"Onecut1","Onecut2","Onecut3","Cdkn1a","Prss2","Runx1","Gata6",
         #"Gata6", "S100a11", "Nr5a2")

# goi.acini <- c("Ctrb1","Cpa1","Gata6","Bhlha15","Nr5a2","Ptf1a")
# goi.duct <- c("Hnf1b","Sox9","Krt19","Gata6","Onecut1")
# features <- c("Cpa1","Gata6","Sox9","Onecut1","Ngn3","Nr5a2","Ptf1a","Pdx1")
# features <- c("Hes1","Dclk1","Sox9","Gata6","Ptf1a","Pdx1")
# goi.PDAC <- c("Dclk1","Pdx1")
# goi.PDACfromDuct <- "Agr2"
# goi.met <- c(c)
# features <- c("Cpa1","Gata6","Sox9","Onecut1","Ngn3","Nr5a2","Ptf1a","Pdx1")

# #Metastases 1
# features <- c("Cybrd1","Nr1d1","Bsg","Tmprss4","Tm9sf3",
#               "Mmp23","Rhof","Sftpd", "Aqp5","Ccna1",
#               "Muc3","Muc5ac","Muc3a","Kif12","Calml4",
#               "Dbp", "Mrtfb", "Rplp0","Dnajc10","Rps12",
#               "Pdzd8", "Mtch2", "Msln", "Prom1", "Vars2")
# #Metastases 2
# features <- c("Porcn","Rpl6","Ybx1","Wfdc2","Tpi1",
#               "Golim4","Otop3","F3", "Id2","Adamtsl5",
#               "Bag1","Rnf186","Glis2","Slc35f5","Tspan12",
#               "Slc9a4", "Ephb2", "Tmem45b","Tmprss2","Pdxdc1",
#               "Lgals2", "Esrp1", "Tmem54", "Ptprf", "Ccnd2",
#               "Ern2","Sult1c2")
# #Metastases 3 / CNS
# features <- c("Gltp","Spock3","Sgms2","Rasgrf1","St8sia3",
#               "Rap1gap","Rbms3","Ccdc92","Ncald","Ppp1r1b",
#               "Gabbr2","Nt5c2","Cdkn2a","Atrnl1","Camk2n1",
#               "Setbp1","Dennd4c","Hs3st1","Shf")
# 
# #Brain/synaptogenesis/neuronal
# features <- c("Nfib", "Tuba1b", "Net1", "Ncald","Spock3",
#               "Rock2", "Sem1", "Ctnnd1","Adgre5", "Dennd4c",
#               "Smad4", "Flna", "Cntn1", "Cntn6","Sgms2",
#               "Nrxn1","Nrxn2","Nrxn3","Lamb2","Rasgrf1",
#               "Sema3d", "Sema4b","Sema4g","Sema5a","St8sia3",
#               "Lama5", "Rtn4", "Picalm","Efnb2", "Rbms3")

list <- "Kras"
list <- c("Rac1", "St8sia3", "Camk2n1", "Cdc42", "Spock3", "Rasgrf1")



levels(mySeurat)
levels(mySeurat) <- c("Metastasis","Carcinoma", "PanIN","ADM","Bystander","Normal acini",
                      "Normal Islet", "EMT", "Stroma")
levels(mySeurat) <- c("Metastasis","Carcinoma", "PanIN3","PanIN2", "PanIN1", "ADM","Bystander","Normal acini")

levels(x = mySeurat) <- c("1-Normal acini", "2-Bystander","3-ADM","4-PanINlo","5-PanINhi","6-PDAC","7-metastasis")
levels(x = mySeurat) <- c("Metastasis","Carcinoma", "PanIN3","PanIN2", "PanIN1", "ADM","Bystander","Normal acini")


features <- c("Rock2","Ephb2","Efnb2", "Adam10", "Mmp2", "Mmp9", 
              "Nrxn1", "Nrxn2", "Nrxn3", "Nrp2", "Sema3e",
              "Lama5", "Itgb1")

features <- c("Ezr","S100a6", "Gsto1", "Gkn1", 
              "Lypd8l", "Anxa2", "Cdh1", "Prom1", "Myrf", 
              "Flna", "Slc12a2", "Actn1", "Fn1", "Hnf1b",
              "Vasp","Vdac2", "Syncrip", "Rpl5", "Pard3",
              "Dync1i2", "Calm1", "Calm2", "Calm3", "Itgb1")


features <- c("Kras","Trp53","Net1","Nt5c2","Ezr","Clu","S100a6",  
              "Anxa2", "Myrf", "Sema4b","Sema4g","Efnb2",
              "Flna", "Slc12a2", "Actn1", "Actb","Tuba1b",
              "Vasp", "Syncrip", "Pard3","Rock2","Rac1", "Rhoa", "Cdc42",
              "Dync1i2", "Calm1", "Calm2", "Calm3","Lama5", "Itgb1")



features <- c("Sema3f", "Nrp1")

features <- c("Sema7a","Sema3e","Sema4a","Sema4b","Sema4g","Efnb2", "Myo5b")

features <- c("Rock2", "Rhoa","Rhoc","Rac1","Cdc42", "Vcam1", "Ezr")

#LR pairs for Calm
features <- c("Calm1", "Mylk",
              "Calm1", "Ptpra",
              "Calm2", "Mylk",
              "Calm2", "Aqp1",
              "Calm3", "Ar")
#LR pairs for Lama5
features <- c("Lama5", "Itgb1", 
              "Lama5","Itgb4", 
              "Lama5","Itga2",
              "Lama5","Itga6", 
              "Lama5","Sdc1",
              "Lama5","Bcam")
#LR pairs for Efnb2
features <- c("Efnb2","Ephb6",
              "Pecam1",
              "Ephb4",
              "Rhbdl2",
              "Epha3",
              "Ephb1",
              "Ephb3",
              "Ephb2",
              "Epha4")

fig <- RidgePlot(mySeurat, sort = T, #split.by = "dx3.KPC",
         #idents = c("Metastasis","Carcinoma", "PanIN","ADM","Bystander","Normal acini"), 
          idents = c("Normal", "ADM","PanIN","PDAC", "Lung_met", "Liver_met", "Stroma"),
          #idents = c("7-metastasis", "6-PDAC","5-PanINhi","4-PanINlo","3-ADM","2-Bystander","1-Normal acini"),
          features = features, ncol = 2)
fig


setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("fig-LR_pairs_Efnb2.tiff", units="in", width=10, height=17, res=300)
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


































as.Seurat.NanoStringGeoMxSet <- function(x, ident = NULL, normData = NULL, 
                                         coordinates = NULL, 
                                         forceRaw = FALSE, ...){
  
  if (!try(requireNamespace("Seurat", quietly = TRUE))) {
    stop("Please install Seurat from CRAN before converting to a Seurat object")
  }else{
    requireNamespace("Seurat", quietly = TRUE)
  }
  
  
  if(featureType(x) == "Probe"){
    stop("Data must be on Target level before converting to a Seurat Object")
  }
  
  if(is.null(normData)){
    stop("normData must not be NULL, please provide name of normalized counts matrix")
  }
  
  if(!normData %in% assayDataElementNames(x)){
    stop(paste0("The normData name \"", normData, "\" is not a valid assay name. Valid names are: ", 
                paste(assayDataElementNames(x), collapse = ", ")))
  }
  
  normFactor_names <- "normFactors|qFactors|negFactors|hkFactors|hknormFactors"
  
  if(length(grep(pattern = normFactor_names, names(sData(x)))) == 0 & 
     forceRaw == FALSE){
    stop("It is NOT recommended to use Seurat's normalization for GeoMx data. 
             Normalize using GeomxTools::normalize() or set forceRaw to TRUE if you want to continue with Raw data")
  }
  
  
  sequencingMetrics <- c("FileVersion", "SoftwareVersion", "Date", "Plate_ID", 
                         "Well", "SeqSetId", "Raw", "Trimmed", "Stitched", 
                         "Aligned", "umiQ30", "rtsQ30", "DeduplicatedReads", 
                         "NTC_ID", "NTC", "Trimmed (%)", "Stitched (%)", 
                         "Aligned (%)", "Saturated (%)")
  
  QCMetrics <- "QCFlags"
  
  seuratConvert <- suppressWarnings(Seurat::CreateSeuratObject(counts = assayDataElement(x, normData), 
                                                               assay = "GeoMx", 
                                                               project = expinfo(experimentData(x))[["title"]]))
  seuratConvert <- suppressWarnings(Seurat::AddMetaData(object = seuratConvert, 
                                                        metadata = sData(x)[,!colnames(sData(x)) %in% 
                                                                              c(sequencingMetrics,
                                                                                QCMetrics)]))
  seuratConvert@assays$GeoMx <- Seurat::AddMetaData(object = seuratConvert@assays$GeoMx, 
                                                    metadata = fData(x))
  
  if(!is.null(ident)){
    if(!ident %in% colnames(seuratConvert@meta.data)){
      stop(paste0("ident \"", ident, "\" not found in GeoMxSet Object"))
    }
    
    Seurat::Idents(seuratConvert) <- seuratConvert[[ident]]
  }
  
  
  seuratConvert@misc <- otherInfo(experimentData(x)) 
  seuratConvert@misc[["sequencingMetrics"]] <- sData(x)[colnames(sData(x)) %in% 
                                                          sequencingMetrics]
  seuratConvert@misc[["QCMetrics"]] <- sData(x)[colnames(sData(x)) %in% 
                                                  QCMetrics]
  
  if(ncol(seuratConvert@misc[["QCMetrics"]]) == 0){
    seuratConvert@misc[["QCMetrics"]] <- NULL
  }
  
  if(!is.null(coordinates)){
    xcoord <- coordinates[1]
    ycoord <- coordinates[2]
    
    if(xcoord %in% colnames(seuratConvert@meta.data) & 
       ycoord %in% colnames(seuratConvert@meta.data)){
      coord.df <- data.frame(x=seuratConvert@meta.data[[xcoord]], 
                             y=seuratConvert@meta.data[[ycoord]])
      colnames(coord.df) <- coordinates
      seuratConvert@meta.data <- seuratConvert@meta.data[!colnames(seuratConvert@meta.data) %in% 
                                                           coordinates]
    }else{
      if(!xcoord %in% colnames(seuratConvert@meta.data) &
         !ycoord %in% colnames(seuratConvert@meta.data)){
        stop(paste0("xcoord \"", xcoord, "\" and ycoord \"", 
                    ycoord, "\" not found in GeoMxSet Object"))
      }
      
      if(!xcoord %in% colnames(seuratConvert@meta.data)){
        stop(paste0("xcoord \"", xcoord, 
                    "\" not found in GeoMxSet Object"))
      }
      
      if(!ycoord %in% colnames(seuratConvert@meta.data)){
        stop(paste0("ycoord \"", ycoord, 
                    "\" not found in GeoMxSet Object"))
      }
    }
    
    rownames(coord.df) <- rownames(seuratConvert@meta.data)
    
    # need to create DSP specific image class
    seuratConvert@images$image =  new(
      Class = 'SlideSeq',
      assay = "GeoMx",
      key = "image_",
      coordinates = coord.df
    )
  }
  
  return(seuratConvert)
}

#' Convert Object to SpatialExperiment
#' 
#' @param x GeoMxSet object to convert 
#' @param ... arguments to be passed to other methods
#' 
#' @export


