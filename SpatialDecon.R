library(SpatialDecon)
library(GeomxTools)

#https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html
#https://github.com/Nanostring-Biostats/CellProfileLibrary/tree/NewProfileMatrices


# Pancreas - MCA		
# Main Cell Types	Granular
# 1	Acinar	Acinar cell
# 2	B	B cell
# 3	Beta	Beta cell
# 4	Dendritic	Dendritic cell
# 5	Dividing	Dividing cell
# 6	Ductal	Ductal cell
# 7	Endocrine	Endocrine cell
# 8	Endothelial	Endothelial cell_Fabp4 high
# Endothelial cell_Lrg1 high
# Endothelial cell_Tm4sf1 high
# 9	Erythroblast	Erythroblast_Hbb-bt high
# Erythroblast_Igkc high
# 10	Glial	Glial cell
# 11	Granulocyte	Granulocyte
# 12	Macrophage	Macrophage
# 13	Smooth muscle	Smooth muscle cell
# 14	Stromal	Stromal cell_Fn1 high
# Stromal cell_Mfap4 high
# Stromal cell_Smoc2 high
# 15	T	T cell



library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)
# knitr::opts_chunk$set(echo = TRUE)
# output_prefix<-"CPTR474"
# projectname<-"CPTR474"
# datadir<-"C:/Users/edmondsonef/Desktop/DSP GeoMX/data/WTA_04122022/raw_data"
# DCCdir<-"DCC-20220420"
# PKCfilename<-"Mm_R_NGS_WTA_v1.0.pkc"
# WorkSheet<-"final.xlsx"
# final <- read_excel("C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/final.xlsx")
# 
# DCCFiles <- list.files(file.path(datadir , DCCdir), pattern=".dcc$", full.names=TRUE)
# PKCFiles <- file.path(datadir, PKCfilename)
# SampleAnnotationFile <- file.path(datadir, WorkSheet)
# 
# #Shift counts to one to mimic how DSPDA handles zero counts
# myData <- shiftCountsOne(myData, elt="exprs", useDALogic=TRUE) 
# 
# pkcs <- annotation(myData)
# modules <- gsub(".pkc", "", pkcs)
# kable(data.frame(PKCs = pkcs, modules = modules))



# myData <- readRDS(file = "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/my_data.rds")
# target_myData <- readRDS(file = "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/target_myData.rds")
# Pancreas_MCA <- readRDS(file = "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/Pancreas_MCA.RData")

load("C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/Pancreas_MCA.RData")
load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_geoMX_new.RData")

dim(myData)

dim(target_myData)
head(pData(target_myData))
     


target_myData@assayData$exprs[seq_len(5), seq_len(5)]

sampleNames(target_myData) <- paste0(target_myData$dxIPMN2, "-", target_myData$ID)

bg = derive_GeoMx_background(norm = target_myData@assayData$q_norm,  ### COULDN'T FIND exprs_norm?  neg norm negFac
                             probepool = fData(target_myData)$Module,
                             negnames = target_myData@featureData@data$Negative)



signif(profile_matrix[seq_len(3), seq_len(3)], 2)




mousepanc <- download_profile_matrix(species = "Mouse",
                                       age_group = "Adult", 
                                       matrixname = "Pancreas_MCA")
head(mousepanc)
heatmap(sweep(mousepanc, 1, apply(mousepanc, 1, max), "/"),
        labRow = NA, margins = c(10, 5))




res = runspatialdecon(object = target_myData,
                      norm_elt = "q_norm", #neg_norm
                      raw_elt = "exprs",
                      cell_counts = target_myData@phenoData@data$nuclei,
                      X = mousepanc,
                      align_genes = TRUE)

str(pData(res))
names(res@assayData)

heatmap(t(res$beta), 
        #cexCol = 0.5, 
        #cexRow = 0.7, 
        #margins = c(10,7),
        labCol = res$region)

data("cellcols")
cellcols

o = hclust(dist(t(res$cell.counts$cell.counts)))$order
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 3))


TIL_barplot(t(res$cell.counts$cell.counts[, o]), 
            draw_legend = TRUE, 
            cex.names = 0.5)

TIL_barplot(t(res$prop_of_all), 
            draw_legend = TRUE, 
            cex.names = 0.75)





# PCA of the normalized data:
pc = prcomp(t(log2(pmax(res@assayData$q_norm, 1))))$x[, c(1, 2)]

# run florets function:
par(mar = c(5,5,1,1))
layout(mat = (matrix(c(1, 2), 1)), widths = c(6, 2))
florets(x = pc[, 1], y = pc[, 2],
        b = t(res$beta), cex = .5,
        xlab = "PC1", ylab = "PC2")
par(mar = c(0,0,0,0))
frame()
legend("center", 
       #fill = colnames(res$beta), 
       legend = colnames(res$beta), 
       cex = 0.7)








heatmap(sweep(res@experimentData@other$SpatialDeconMatrix, 1, apply(res@experimentData@other$SpatialDeconMatrix, 1, max), "/"),
        labRow = NA, margins = c(10, 5))
colnames(res$beta)
matching = list()
matching$B.cell = c( "B.cell")
matching$DC = "Dendritic.cell"
matching$Granulocyte= "Granulocyte"
matching$Macrophage = "Macrophage"
#matching$erythroblast = c("Erythroblast.Hbb.bt.high", "Erythroblast.Igkc.high")
matching$endothelial = c("Endothelial.cell.Lrg1.high", "Endothelial.cell.Fabp4.high","Endothelial.cell.Tm4sf1.high")
matching$other = c("Erythroblast.Hbb.bt.high", "Erythroblast.Igkc.high",
                   "Glial.cell")
matching$stroma = c("Stromal.cell.Fn1.high", "Stromal.cell.Mfap4.high","Stromal.cell.Smoc2.high")
matching$islet = c("Beta.cell","Endocrine.cell")
matching$acinar = "Acinar.cell"
matching$dividing.cell = "Dividing.cell"
matching$ductal = "Ductal.cell"
matching$T.cell = "T.cell"
matching$smooth.muscle = "Smooth.muscle.cell" 
#matching$glial.cell = "Glial.cell"


collapsed = runCollapseCellTypes(object = res, 
                                 matching = matching)




heatmap(t(collapsed$beta), #cexRow = 0.85, cexCol = 0.75,
        labCol = res$dxIPMN2)


o = hclust(dist(t(collapsed$cell.counts$cell.counts)))$order
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 3))
TIL_barplot(t(collapsed$cell.counts$cell.counts[, o]), draw_legend = TRUE, 
            cex.names = 0.5)

o = hclust(dist(t(collapsed$cell.counts$cell.counts)))$order
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 3))
TIL_barplot(t(collapsed$prop_of_all), 
            draw_legend = TRUE, 
            cex.names = 0.75)



colnames(collapsed$beta)
str(pData(collapsed))

ddd <- pData(collapsed)[ which(collapsed$dxIPMN3=='IPMN'), ]

o = hclust(dist(t(ddd$cell.counts$cell.counts)))$order
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 3))
TIL_barplot(t(ddd$prop_of_all), 
            draw_legend = TRUE, 
            cex.names = 0.75,
            main = "IPMN")



manycols <- c("#8DD3C7", "#FFFFB3", "#BEBADA", 
              "#FB8072", "#80B1D3", "#FDB462", 
              "#B3DE69", "#FCCDE5", "#A6CEE3", 
              "#1F78B4", "#B2DF8A", "#33A02C", 
              "#FB9A99", "#E31A1C", "#FDBF6F", 
              "#FF7F00", "#1B9E77", "#D95F02", 
              "#7570B3", "#E7298A", "#66A61E", 
              "#E6AB02", "#A6761D", "#666666", 
              sample(grDevices::colors(), 99))

# PCA of the normalized data:
pc = prcomp(t(log2(pmax(collapsed@assayData$q_norm, 1))))$x[, c(1, 2)]

# run florets function:
setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("PCA with florets_list.tiff", units="in", width=13, height=8, res=250)
par(mar = c(5,5,1,1))
layout(mat = (matrix(c(1, 2), 1)), widths = c(6, 2))
florets(x = pc[, 1], y = pc[, 2],
        b = t(collapsed$beta), cex = 1.5,
        xlab = "PC1", ylab = "PC2")
par(mar = c(0,0,0,0))
frame()
legend("center", 
       fill = manycols, 
       legend = colnames(collapsed$beta), 
       cex = 0.9)
dev.off()










res = runspatialdecon(object = target_myData,
                      norm_elt = "q_norm", #q_norm
                      raw_elt = "exprs",
                      cell_counts = target_myData@phenoData@data$nuclei,
                      X = mousepanc,
                      align_genes = TRUE)






rdecon = runReverseDecon(object = target_myData,
                         norm_elt = "q_norm",
                         beta = collapsed$beta)
str(fData(rdecon))
#> 'data.frame':    1700 obs. of  12 variables:
#>  $ TargetName            : chr  "ABCF1" "ABL1" "ACVR1B" "ACVR1C" ...
#>  $ HUGOSymbol            : chr  "ABCF1" "ABL1" "ACVR1B" "ACVR1C" ...
#>  $ TargetGroup           : chr  "All Probes;Transport of small molecules" "All Probes;Cell Cycle;Signaling by Rho GTPases;DNA Repair;Factors involved in megakaryocyte development and pla"| __truncated__ "All Probes;Signaling by NODAL;Signaling by TGF-beta family members" "All Probes;Signaling by NODAL;Signaling by TGF-beta family members" ...
#>  $ AnalyteType           : chr  "RNA" "RNA" "RNA" "RNA" ...
#>  $ Codeclass             : chr  "Endogenous" "Endogenous" "Endogenous" "Endogenous" ...
#>  $ Module                : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ CorrelationToNegatives: num  0.597 0.876 0.435 0.928 0.776 ...
#>  $ GlobalOutliers        : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ Negative              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#>  $ coefs                 : num [1:1700, 1:19] 2.13 2.03 1.3 1.13 1.28 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1700] "ABCF1" "ABL1" "ACVR1B" "ACVR1C" ...
#>   .. ..$ : chr [1:19] "(Intercept)" "macrophages" "mast" "B.naive" ...
#>  $ cors                  : num  0.763 0.517 0.804 0.191 0.662 ...
#>  $ resid.sd              : num  0.468 0.315 0.55 0.322 0.383 ...
names(rdecon@assayData)
#> [1] "exprs_norm" "resids"     "exprs"      "yhat"

# look at the residuals:
setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("ReverseDeconvolution.tiff", units="in", width=20, height=20, res=300)
heatmap(pmax(pmin(rdecon@assayData$resids, 2), -2), labCol = res$dxIPMN3)
dev.off()

goi <- c("Kras", "Trp53", "Cd274", "Cd8a", "Cd68", "Epcam","Cre",
                   "Krt18", "Notch1", "Notch2", "Notch3", "Notch4","Cldn8",
                    "Cdk6","Msh3","Myc","Mastl", "Sox2","Cav1","Fosl1","Gata4",
                    "Cldn18","Capn6","Cpa1","Muc5ac","Tff1","Smad4","Sox9",
                    "Ptf1a","Pdx1","Nr5a2","Neurog3","Bhlha15","Krt19","Dclk1",
                    "Elastase","Hnf1b","Krt19","Ngn3","Ctrb1", "Hes1", "Smad4",
                    "Onecut1","Onecut2","Onecut3","Cdkn1a","Prss2","Runx1","Gata6",
                    "Gata6", "S100a11", "Nr5a2","Agr2", "Foxa2", "Fosl1","Ets2", "Runx3")

# look at the two metrics of goodness-of-fit:
plot(fData(rdecon)$cors, fData(rdecon)$resid.sd, col = 0)
showgenes = goi
text(fData(rdecon)$cors[!rownames(fData(rdecon)) %in% showgenes], 
     fData(rdecon)$resid.sd[!rownames(fData(rdecon)) %in% showgenes], 
     setdiff(rownames(fData(rdecon)), showgenes), cex = 0.5)
text(fData(rdecon)$cors[rownames(fData(rdecon)) %in% showgenes], fData(rdecon)$resid.sd[rownames(fData(rdecon)) %in% showgenes], 
     showgenes, cex = 0.75, col = 2)



