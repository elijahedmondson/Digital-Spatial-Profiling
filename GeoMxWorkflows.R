###
###
### TO DO LIST:
### 1. GLMM model to include path grades as continuous variables
###   > fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
###
### 2. Normalization methods comparison:
###   a. Q3 (
###       - Des's methods excluded more genes vs "vignettes/GeoMxWorkflows" 
###   b. Negative normal
####
#### 8/26/2022
####
#### 8/26/2022
####
#### 8/26/2022
####
#### 8/26/2022
####
#### Next: compareCluster()
#### Next: compareCluster()
#### Next: compareCluster()
#### Next: compareCluster()
#### Next: compareCluster()
#### Next: compareCluster()
#### Next: compareCluster()
#### Next: compareCluster()
####
#### 8/26/2022
####
#### 8/26/2022
####
#### 8/26/2022
####
#### 8/26/2022
####
###
###
library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)
library(enrichplot)
library(data.table)
library(fgsea)
library(ggplot2)
library(ggrepel) 
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationHub)
library(GOSemSim)
library(clusterProfiler)
library(GOSemSim)
library(ggwordcloud)
library(ggplot2)
library(cowplot)
library(ReactomePA)
library(DOSE)
library(msigdbr)
library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)

#####


knitr::opts_chunk$set(echo = TRUE)
output_prefix<-"CPTR474"
projectname<-"CPTR474"
datadir<-"C:/Users/edmondsonef/Desktop/DSP GeoMX/data/WTA_04122022/raw_data"
DCCdir<-"DCC-20220420"
PKCfilename<-"Mm_R_NGS_WTA_v1.0.pkc"
WorkSheet<-"final.xlsx"
final <- read_excel("C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/final.xlsx")

DCCFiles <- list.files(file.path(datadir , DCCdir), pattern=".dcc$", full.names=TRUE)
PKCFiles <- file.path(datadir, PKCfilename)
SampleAnnotationFile <- file.path(datadir, WorkSheet)

myData<-readNanoStringGeoMxSet(dccFiles = DCCFiles,
                               pkcFiles = PKCFiles,
                               phenoDataFile = SampleAnnotationFile,
                               phenoDataSheet = "Template",
                               phenoDataDccColName = "Sample_ID",
                               protocolDataColNames = c("aoi", "roi"),
                               experimentDataColNames = c("panel")) 

#Shift counts to one to mimic how DSPDA handles zero counts
myData <- shiftCountsOne(myData, elt="exprs", useDALogic=TRUE) 

pkcs <- annotation(myData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))





#myData <- readRDS(file = "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/my_data.rds")
#target_myData <- readRDS(file = "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/target_myData.rds")




#####
# # select the annotations we want to show, use `` to surround column names with
# # spaces or special symbols
# count_mat <- count(pData(myData), `Position`, Class, Origin, Sex, Age, Strain, Call, dx)
# # simplify the slide names
# count_mat$`core` <- gsub("disease", "d",
#                                gsub("normal", "n", count_mat$`Position`))
# # gather the data and plot in order: class, slide name, region, segment
# test_gr <- gather_set_data(count_mat, 1:7)
# test_gr$x <- factor(test_gr$x,
#                     levels = c("Strain","Sex", "Age", "Position", "Class","Origin", "Call"))
# # plot Sankey
# sampleoverview <- ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
#   geom_parallel_sets(aes(fill = dx), alpha = 0.5, axis.width = 0.1) +
#   geom_parallel_sets_axes(axis.width = 0.2) +
#   geom_parallel_sets_labels(color = "white", size = 4) +
#   theme_classic(base_size = 17) + 
#   theme(legend.position = "bottom",
#         axis.ticks.y = element_blank(),
#         axis.line = element_blank(),
#         axis.text.y = element_blank()) +
#   scale_y_continuous(expand = expansion(0)) + 
#   scale_x_discrete(expand = expansion(0)) +
#   labs(x = "", y = "") +
#   annotate(geom = "segment", x = 7.25, xend = 7.25,
#            y = 0, yend = 20, lwd = 2) +
#   annotate(geom = "text", x = 7.19, y = 7.8, angle = 90, size = 4,
#            hjust = 0.5, label = "20 segments")
# 
# 
# sampleoverview
# 
# setwd("C:/Users/edmondsonef/Desktop/R-plots/")
# tiff("sampleoverview.tiff", units="in", width=19, height=15, res=150)
# sampleoverview
# dev.off()
# 
# 



## ----setqcflagupdated,  eval = TRUE-------------------------------------------
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 4,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
myData <-
  setSegmentQCFlags(myData, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(myData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))


## ----qcflagHistogramsCode, eval = TRUE, warning = FALSE, message = FALSE------
library(ggplot2)

col_by <- "dx"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 200) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 7) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}


## ----plotQCHist, warning = FALSE, message = FALSE-----------------------------
# QC_histogram(sData(myData), "Trimmed (%)", col_by, 80)
# QC_histogram(sData(myData), "Stitched (%)", col_by, 80)
# QC_histogram(sData(myData), "Aligned (%)", col_by, 75)
# QC_histogram(sData(myData), "Saturated (%)", col_by, 50) +
#   labs(title = "Sequencing Saturation (%)",
#        x = "Sequencing Saturation (%)")
# QC_histogram(sData(myData), "area", col_by, 10, scale_trans = "log10")
# QC_histogram(sData(myData), "nuclei", col_by, 10)

#QC_histogram(sData(myData), "Aligned", col_by, 10000)
# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(myData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(myData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(myData)[, negCols] <- sData(myData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(myData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}


# detatch neg_geomean columns ahead of aggregateCounts call
pData(myData) <- pData(myData)[, !colnames(pData(myData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
print("No Template Control (NTC) wells are essential for detecting contamination or non-specific amplification")

kable(table(NTC_Count = sData(myData)$NTC),
      col.names = c("NTC Count", "# of Segments"))



## ----QCSummaryTable, results = "aexprs()## ----QCSummaryTable, results = "asis"-----------------------------------------
kable(QC_Summary, caption = "QC Summary Table for each Segment")

## ----removeQCSampleProbe, eval = TRUE-----------------------------------------
myData <- myData[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(myData)

## ----setbioprobeqcflag,  eval = TRUE------------------------------------------
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
myData <- setBioProbeQCFlags(myData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(myData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

## ----bioprobeQCTable, echo = FALSE, results = "asis"--------------------------
kable(qc_df, caption = "Probes flagged or passed as outliers")


## ----excludeOutlierProbes-----------------------------------------------------
#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(myData, 
         fData(myData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(myData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
myData <- ProbeQCPassed 

## ----aggregateCounts, eval = TRUE---------------------------------------------
# Check how many unique targets the object has
length(unique(featureData(myData)[["TargetName"]]))

# collapse to targets
target_myData <- aggregateCounts(myData)
dim(target_myData)
exprs(target_myData)[100:103, 1:4]


## ----calculateLOQ, eval = TRUE------------------------------------------------
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_myData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_myData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_myData)[, vars[1]] * 
             pData(target_myData)[, vars[2]] ^ cutoff)
  }
}
pData(target_myData)$LOQ <- LOQ

head(pData(target_myData)$LOQ)

## ----LOQMat, eval = TRUE------------------------------------------------------
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_myData)$Module == module
  Mat_i <- t(esApply(target_myData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_myData)$TargetName, ]

## ----segDetectionBarplot------------------------------------------------------
# Save detection rate information to pheno data
pData(target_myData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_myData)$GeneDetectionRate <-
  pData(target_myData)$GenesDetected / nrow(target_myData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_myData)$DetectionThreshold <- 
  cut(pData(target_myData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15,1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_myData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = dx)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

####
####
####
## ----segTable-----------------------------------------------------------------
# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_myData)$DetectionThreshold,
            pData(target_myData)$dx))

## ----filterSegments-----------------------------------------------------------
target_myData <-
  target_myData[, pData(target_myData)$GeneDetectionRate >= .05]                    ########EFE excludes samples with low gene detection
pData(target_myData)[,24:27]

dim(target_myData)

target_myData@phenoData@data$dx

## ----replotSankey, fig.width = 10, fig.height = 8, fig.wide = TRUE, message = FALSE, warning = FALSE----
# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols

# count_mat <- count(pData(myData), `Position`, Class, Origin, Sex, Age, Strain, Call, dx)
# # simplify the slide names
# count_mat$`core` <- gsub("disease", "d",
#                          gsub("normal", "n", count_mat$`Position`))
# # gather the data and plot in order: class, slide name, region, segment
# test_gr <- gather_set_data(count_mat, 1:7)
# test_gr$x <- factor(test_gr$x,
#                     levels = c("Strain","Sex", "Age", "Position", "Class","Origin", "Call"))
# # plot Sankey
# sampleoverview2 <- ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
#   geom_parallel_sets(aes(fill = dx), alpha = 0.5, axis.width = 0.1) +
#   geom_parallel_sets_axes(axis.width = 0.2) +
#   geom_parallel_sets_labels(color = "white", size = 4) +
#   theme_classic(base_size = 17) + 
#   theme(legend.position = "bottom",
#         axis.ticks.y = element_blank(),
#         axis.line = element_blank(),
#         axis.text.y = element_blank()) +
#   scale_y_continuous(expand = expansion(0)) + 
#   scale_x_discrete(expand = expansion(0)) +
#   labs(x = "", y = "") +
#   annotate(geom = "segment", x = 7.25, xend = 7.25,
#            y = 0, yend = 20, lwd = 2) +
#   annotate(geom = "text", x = 7.19, y = 7.8, angle = 90, size = 4,
#            hjust = 0.5, label = "20 segments")
# 
# 
# sampleoverview2

# setwd("C:/Users/edmondsonef/Desktop/R-plots/")
# tiff("sampleoverview2.tiff", units="in", width=19, height=15, res=150)
# sampleoverview2
# dev.off()












## ----goi detection------------------------------------------------------------


library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_myData)]
fData(target_myData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_myData)$DetectionRate <-
  fData(target_myData)$DetectedSegments / nrow(pData(target_myData))

# Gene of interest detection table

goi <- c("Kras", "Trp53", "Cd274", "Cd8a", "Cd68", "Epcam","Cre",
         "Krt18", "Notch1", "Notch2", "Notch3", "Notch4","Cldn8",
         "Cdk6","Msh3","Myc","Mastl", "Sox2","Cav1","Fosl1","Gata4",
         "Cldn18","Capn6","Cpa1","Muc5ac","Tff1","Smad4","Sox9",
         "Ptf1a","Pdx1","Nr5a2","Neurog3","Bhlha15","Krt19","Dclk1",
         "Fap","Hnf1b","Krt19","Ctrb1", "Hes1", "Smad4",
         "Onecut1","Onecut2","Onecut3","Cdkn1a","Prss2","Runx1","Gata6",
         "Gata6", "S100a11", "Nr5a2","Agr2", "Foxa2", "Fosl1","Ets2", "Runx3")
# 
# goi.acini <- c("Ctrb1","Cpa1","Gata6","Bhlha15","Nr5a2","Ptf1a")
# goi.duct <- c("Hnf1b","Sox9","Krt19","Gata6","Onecut1")
# goi.ADM <- c("Cpa1","Gata6","Sox9","Onecut1","Neurog3","Nr5a2","Ptf1a","Pdx1")
# goi.PanIN <- c("Hes1","Dclk1","Sox9","Gata6","Ptf1a","Pdx1")
# goi.PDAC <- c("Dclk1","Pdx1")
# goi.PDACfromDuct <- "Agr2"
# goi.met <- c("Pdzd8", "Mtch2", "Spock3", "Serpina3k", "Cybrd1", "Vars2")
# 


#hnf6 = Onecut1
#Ngn3 = Neurog3
#Mist1 = Bhlha15

goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_myData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_myData)[goi, "DetectionRate"]))

## ----tableGOI, echo = FALSE, results = "asis"---------------------------------
kable(goi_df, caption = "Detection rate for Genes of Interest", align = "c",
      col.names = c("Gene", "Detection, # Segments", "Detection Rate, % of Segments"))

## ----plotDetectionRate, eval = TRUE-------------------------------------------
#Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 3, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_myData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_myData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

## ----subsetGenes, eval = TRUE-------------------------------------------------
# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_myData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_myData <- 
  target_myData[fData(target_myData)$DetectionRate >= 0.035 |                    ########EFE change to include additional genes? 
                    fData(target_myData)$TargetName %in% neg_probes, ]
dim(target_myData)

# retain only detected genes of interest
goi <- goi[goi %in% rownames(target_myData)]

## ----previewNF, fig.width = 8, fig.height = 8, fig.wide = TRUE, eval = TRUE, warning = FALSE, message = FALSE----
library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "dx2"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_myData)),
             Segment = colnames(exprs(target_myData)),
             Annotation = pData(target_myData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_myData), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_myData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) +
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

## ----normalizeObject, eval = TRUE---------------------------------------------
# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_myData <- normalize(target_myData , data_type = "RNA",
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_myData <- normalize(target_myData , data_type = "RNA",
                             norm_method = "neg", 
                             fromElt = "exprs",
                             toElt = "neg_norm")

## ----normplot, fig.small = TRUE-----------------------------------------------
# visualize the first 10 segments with each normalization method
# boxplot(exprs(target_myData)[,1:77],
#         col = "#9EDAE5", main = "Raw Counts",
#         log = "y", names = 1:77, xlab = "Segment",
#         ylab = "Counts, Raw")
# 
# boxplot(assayDataElement(target_myData[,1:77], elt = "q_norm"),
#         col = "#2CA02C", main = "Q3 Norm Counts",
#         log = "y", names = 1:77, xlab = "Segment",
#         ylab = "Counts, Q3 Normalized")
# 
# boxplot(assayDataElement(target_myData[,1:77], elt = "neg_norm"),
#         col = "#FF7F0E", main = "Neg Norm Counts",
#         log = "y", names = 1:77, xlab = "Segment",
#         ylab = "Counts, Neg. Normalized")

## ----dimReduction, eval = TRUE------------------------------------------------
library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_myData , elt = "q_norm"))),  
       config = custom_umap)
pData(target_myData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
ggplot(pData(target_myData),
       aes(x = UMAP1, y = UMAP2, color = comps, shape = Call, label=dsxf)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  theme_bw()+
  theme(legend.position="none")

# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_myData , elt = "q_norm"))),
        perplexity = ncol(target_myData)*.15)
pData(target_myData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
ggplot(pData(target_myData),
       aes(x = tSNE1, y = tSNE2, color = comps, shape = Call, label=dsxf)) +
  geom_point(size = 3) +geom_text(hjust=1.1, vjust=0.2)+
  theme_bw()+
  theme(legend.position="none")


## run PCA
PCAx<-1
PCAy<-2
PCAxy <- c(as.integer( PCAx ),as.integer( PCAy) ) # selected principal components


pca.object <- prcomp(t(log2(assayDataElement(target_myData , elt = "q_norm"))))
pcaData = as.data.frame(pca.object$x[, PCAxy]); 
pData(target_myData)[, c("PC1", "PC2")] <- pcaData[,c(1,2)]
percentVar=round(100*summary(pca.object)$importance[2, PCAxy],0)


ggplot(pData(target_myData),
               aes(x = PC1, y = PC2, color=comps, label=dsxf)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  xlab(paste0("PC", PCAx ,": ", percentVar[1], "% variance")) +
  ylab(paste0("PC", PCAy ,": ", percentVar[2], "% variance")) +

  theme_bw()+
  theme(legend.position="none")





## ----CVheatmap, eval = TRUE, echo = TRUE, fig.width = 8, fig.height = 6.5, fig.wide = TRUE----
library(pheatmap)  # for pheatmap
# create a log2 transform of the data for analysis
assayDataElement(object = target_myData, elt = "log_q") <-
  assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_myData,
                         elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:50]

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.80)]
pheatmap(assayDataElement(target_myData[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
           pData(target_myData)[, c("dx2", "Sex","Strain")])


















## ----Differential Expression----

load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_geoMX_new.RData")
# If comparing structures that co-exist within a given tissue, use an LMM model 
# with a random slope. Diagnosis is our test variable. We control for tissue 
# sub-sampling with slide name using a random slope and intercept; the intercept adjusts for the multiple 
# regions placed per unique tissue, since we have one tissue per slide. If multiple tissues are placed per slide, 
# we would change the intercept variable to the unique tissue name (ex: tissue name, Block ID, etc).

# convert test variables to factors
pData(target_myData)$testRegion <- 
  factor(pData(target_myData)$comps, c("4-PanINlo", "5-PanINhi", "6-PDAC"))                           ###CHANGE
pData(target_myData)[["slide"]] <-                                            ### Control for 
  factor(pData(target_myData)[["MHL Number"]])
assayDataElement(object = target_myData, elt = "log_q") <-
  assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(status in c("Full ROI")) {
  ind <- pData(target_myData)$segment == status
  mixedOutmc <-
    mixedModelDE(target_myData[, ind],
                 elt = "log_q",
                 modelFormula = ~ testRegion + (1 + testRegion | slide),        ### modelFormula =  Reaction ~ Days + (Days || Subject), sleepstudy)
                 #modelFormula = ~ testRegion + (1 | slide),
                 groupVar = "testRegion",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- status
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}

goi <- c("Kras", "Trp53", "Cd274", "Cd8a", "Cd68", "Epcam","Cre",
         "Krt18", "Notch1", "Notch2", "Notch3", "Notch4","Cldn8",
         "Cdk6","Msh3","Myc","Mastl", "Sox2","Cav1","Fosl1","Gata4",
         "Cldn18","Capn6","Cpa1","Muc5ac","Tff1","Smad4","Sox9",
         "Ptf1a","Pdx1","Nr5a2","Neurog3","Bhlha15","Krt19","Dclk1",
         "Elastase","Hnf1b","Krt19","Ngn3","Ctrb1", "Hes1", "Smad4",
         "Onecut1","Onecut2","Onecut3","Cdkn1a","Prss2","Runx1","Gata6",
         "Gata6", "S100a11", "Nr5a2","Agr2", "Foxa2", "Fosl1","Ets2", "Runx3")


Axonogenesis <- c("Cckar","Slit3","Etv1","Sema7a","Cxcr4","Kif5c","Ptprm",
"Ptprs","Bsg","Smad4","Bcl2","Trak1","Ptch1","
Islr2","Taok2","Cdkl3","Actb","Ednra","Pip5k1c",
"Sptbn4","Lama5","Ablim1","Rtn4","Wnt7b","Spg11","
Golga4","Dock7","Ephb2","Cacna1a","Ptpn11","B4gat1",
"Smo","B4galt6","Rab3a","Ntrk3","Neo1","Lrp1",
"Atp5g1","Kif5b","Brsk2","Erbb2","Map1a","Flrt2",
"Tsku","Map1s","Chrnb2","Fstl4","Lrp4","Dag1","Sin3a",
"Mapk8ip3","Dclk1","Adnp","Celsr3","Rpl4","Tubb2b",
"Efna5","Plxnb2","Ngf","Ache","Vim","Flrt3","Fgfr2","
Ephb4","Flot1","Sema4c","Gsk3b","Sema3d","Aatk","Cdh4",
"Tubb3","Agrn","Evl","Brsk1","Notch3","Fzd3",
"Hsp90aa1","Nrn1","Bcl11a","Sema4g","Lama3","Epha8",
"Ntn5","Amigo1","Apbb1","Mgll","Ret","Atp8a2","
Alcam","Unc5a","Grin1","Cntn6","Wnt7a","Pou4f3","Shh")

Synapse_Organization <- c("Abl2", "Ache", "Actb", "Actr3", "Adam10", "Adgrb2", 
                          "Adgre5", "Adgrf1", "Adgrl3", "Adnp", "Arf1", "Arf4", "Arf6", "Arhgap44", 
                          "Baiap2", "Bhlhb9", "C1qa", "C1ql1", "C3", "Cacna1a", "Cacna1s", "Cacnb1", 
                          "Cacnb4", "Camk1", "Camk2b", "Caprin1", "Cel", "Cfl1", "Chchd10", "Chd4", 
                          "Chrnb1", "Clstn1", "Cnksr2", "Cntnap4", "Col4a1", "Col4a5", "Ctnnb1", 
                          "Cttn", "Cttnbp2", "Cyfip1", "Dact1", "Dag1", "Dbn1", "Dctn1", "Dlgap3",
                          "Dock10", "Drd1", "Efna1", "Efnb2", "Eif4g1", "Epha4", "Ephb2", "Erbb4", 
                          "F2r", "Farp1", "Fgfr2", "Flna", "Flrt3", "Fzd5", "Gabrb3", "Ghrl", "Gnpat",
                          "Gphn", "Grm5", "Hnrnpk", "Hspa8", "Igsf9", "Insr", "Itga3", "Itpka", "Klk8",
                          "Lamb2", "Lgi2", "Lrfn2", "Lrfn5", "Lrrc4c", "Lrrk2", "Lrrtm2", "Lzts3",
                          "Magi2", "Marcks", "Mdga1", "Mdga2", "Mef2c", "Mfn2", "Myh10", "Ndrg1", 
                          "Nedd4", "Nfasc", "Nfatc4", "Nfia", "Nlgn1", "Nlgn3", "Nrcam", "Nrg1", 
                          "Nrg2", "Nrp1", "Nrxn1", "Nrxn3", "Ntn1", "Ntrk2", "Numb", "Obsl1", "Ophn1",
                          "Pak3", "Palm", "Pcdh17", "Pcdhgc4", "Pclo", "Pdlim5", "Pdzrn3", "Pfn1", 
                          "Pfn2", "Picalm", "Pik3r1", "Pin1", "Pmp22", "Ppfia2", "Ppfia4", "Prkca", 
                          "Prrt1", "Psen1", "Ptn", "Ptprf", "Ptprt", "Rab17", "Rab29", "Rab39b", 
                          "Rapsn", "Rhoa", "Rims4", "Rock2", "Sdf4", "Sdk1", "Septin7", "Setd5", 
                          "Sez6", "Sez6l", "Shank1", "Shank2", "Shank3", "Sipa1l1", "Six4", "Slitrk6", 
                          "Snta1", "Sorbs1", "Sparc", "Srcin1", "Srgn", "Ssh1", "St8sia2", "Syngap1",
                          "Tanc2", "Taok2", "Tnc", "Tubb5", "Vcp", "Vps35", "Wnt5a", "Wnt7b", "Ywhaz",
                          "Zmynd8", "Lgmn", "Tuba1b")

select_neural <- c("Actb","Tuba1b","Rock2")


# goi.acini <- c("Ctrb1","Cpa1","Gata6","Bhlha15","Nr5a2","Ptf1a")
# goi.duct <- c("Hnf1b","Sox9","Krt19","Gata6","Onecut1")
# goi.ADM <- c("Cpa1","Gata6","Sox9","Onecut1","Ngn3","Nr5a2","Ptf1a","Pdx1")
# goi.PanIN <- c("Hes1","Dclk1","Sox9","Gata6","Ptf1a","Pdx1")
# goi.PDAC <- c("Dclk1","Pdx1")
# goi.PDACfromDuct <- "Agr2"
# goi.met <- c("Pdzd8", "Mtch2", "Spock3", "Serpina3k", "Cybrd1", "Vars2")

#library(biomaRt)

results.sig <- dplyr::filter(results, abs(results$Estimate) > 0.5)
results.sig <- dplyr::filter(results.sig, results.sig$`Pr(>|t|)` < 0.5)
head(results.sig)
#names(results.sig)[6] <- 'Pr(>|t|)'
head(results.sig)

mt_list = split(results.sig, f = results.sig$Contrast)

names(mt_list)

gene <- mt_list[[2]]

gene <- results.sig

kable(subset(results, Gene %in% goi & Subset == "Full ROI"), digits = 3,
      caption = "DE results for Genes of Interest",
      align = "lc", row.names = FALSE)



library(ggrepel) 
# Categorize Results based on P-value & FDR for plotting
gene$Color <- "NS or FC < 0.5"
gene$Color[gene$`Pr(>|t|)` < 0.05] <- "P < 0.05"
gene$Color[gene$FDR < 0.05] <- "FDR < 0.05"
gene$Color[gene$FDR < 0.001] <- "FDR < 0.001"
gene$Color[abs(gene$Estimate) < 0.5] <- "NS or FC < 0.5"
gene$Color <- factor(gene$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
gene$invert_P <- (-log10(gene$`Pr(>|t|)`)) * sign(gene$Estimate)
top_g <- c()
for(cond in c("Full ROI")) {
  ind <- gene$Subset == cond
  top_g <- c(top_g,
             gene[ind, 'Gene'][
               order(gene[ind, 'invert_P'], decreasing = TRUE)[1:30]],
             gene[ind, 'Gene'][
               order(gene[ind, 'invert_P'], decreasing = FALSE)[1:30]])
}
top_g <- unique(top_g)

gene$Contrast

#reverse log fold change to fit with label
gene$Estimate1 <- gene$Estimate*(-1)
# Graph results
ggplot(gene,                                                             ###CHANGE
       aes(x = Estimate1, y = -log10(`Pr(>|t|)`),
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
  geom_text_repel(data = subset(gene, Gene %in% top_g & FDR < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") 








###GO
names(results)[1] <- 'SYMBOL'
eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)
universe <- distinct(results, SYMBOL, .keep_all = T)

head(gene)
names(gene)[1] <- 'SYMBOL'
eg <- bitr(gene$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"),
           OrgDb="org.Mm.eg.db")
gene <- dplyr::left_join(gene, eg, by = "SYMBOL")
rm(eg)


ego <- enrichGO(gene          = gene$ENTREZID,
                keyType       = "ENTREZID",
                universe      = universe$ENTREZID, ##list of all genes?? 
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", #"BP", "MF", and "CC"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

dotplot(ego)
upsetplot(ego)
library(topGO)
plotGOgraph(ego)

head(ego,10)





ggplot(ego[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))

str(ego)

## extract a dataframe with results from object of type enrichResult
egobp.results.df <- ego@result
head(egobp.results.df)

## create a new column for term size from BgRatio
egobp.results.df$term.size <- gsub("/(\\d+)", "", egobp.results.df$BgRatio)

## filter for term size to keep only term.size => 3, gene count >= 5 and subset
egobp.results.df <- egobp.results.df[which(egobp.results.df[,'term.size'] >= 3 & egobp.results.df[,'Count'] >= 5),]
egobp.results.df <- egobp.results.df[c("ID", "Description", "pvalue", "qvalue", "geneID")]

## format gene list column
egobp.results.df$geneID <- gsub("/", ",", egobp.results.df$geneID)

## add column for phenotype
egobp.results.df <- cbind(egobp.results.df, phenotype=1)
egobp.results.df <- egobp.results.df[, c(1, 2, 3, 4, 6, 5)]

## change column headers
colnames(egobp.results.df) <- c("Name","Description", "pvalue","qvalue","phenotype", "genes")

egobp.results.filename <-file.path(getwd(),paste("clusterprofiler_cluster_enr_results.txt",sep="_"))
write.table(egobp.results.df,egobp.results.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

em_command = paste('enrichmentmap build analysisType="generic" ', 
                   'pvalue=',"0.05", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.25",
                   'coeffecients=',"JACCARD",
                   'enrichmentsDataset1=',egobp.results.filename ,
                   sep=" ")

#enrichment map command will return the suid of newly created network.
em_network_suid <- commandsRun(em_command)

renameNetwork("Cluster1_enrichmentmap_cp", network=as.numeric(em_network_suid))









gene1 <- filter(gene, FDR < 0.001)
setwd("C:/Users/edmondsonef/Desktop/R-plots/")
write.csv(gene1, file = '___MHLnumber.csv')



library(patchwork)
setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Volcano.tiff", units="in", width=18, height=15, res=300)
p01 + p02 + p03 + p04 + p05 + p06  
  plot_layout(guides = "collect") 
dev.off()








## ----targetTable, eval = TRUE, as.is = TRUE-----------------------------------

head(gene)
names(gene)[2] <- 'Gene'
head(gene)

kable(subset(gene, Gene %in% c("Pdzd8", "Mtch2", "Spock3", "Serpina3k", "Cybrd1", "Vars2")), 
      row.names = FALSE)
kable(subset(gene, Gene %in% c("Pdzd8")), row.names = FALSE)
## ----targetExprs, eval = TRUE-------------------------------------------------
# show expression for a single target: PDHA1
ggplot(pData(target_myData),
       aes(x = progression1, fill = progression1,
           y = assayDataElement(target_myData["Sem1", ],
                                elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "Sem1 Expression") +
  scale_y_continuous(trans = "log2") +
  #facet_wrap(~class) +
  theme_bw()

## ----targetExprs2, fig.width = 8, fig.wide = TRUE, eval = TRUE----------------
glom <- pData(target_myData)$progression1# == "Metastasis"

# show expression of PDHA1 vs ITGB1
ggplot(pData(target_myData),
       aes(x = assayDataElement(target_myData["Tuba1b", ],
                                elt = "q_norm"),
           y = assayDataElement(target_myData["Mtch2", ],
                                elt = "q_norm"),
           color = progression1)) +
  geom_vline(xintercept =
               max(assayDataElement(target_myData["Dnajc10", glom],
                                    elt = "q_norm")),
             lty = "dashed", col = "darkgray") +
  geom_hline(yintercept =
               max(assayDataElement(target_myData["Mtch2", glom],
                                    elt = "q_norm")),
             lty = "dashed", col = "darkgray") +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  labs(x = "Dnajc10 Expression", y = "Mtch2 Expression") #+
  #facet_wrap(~class)

## ----heatmap, eval = TRUE, fig.width = 8, fig.height = 6.5, fig.wide = TRUE----
# select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(gene, `FDR` < 0.001)$Gene)
pheatmap(log2(assayDataElement(target_myData[GOI, ], elt = "q_norm")),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 3, cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_myData)[, c("progression1", "progression1")])

## ----maPlot, fig.width = 8, fig.height = 12, fig.wide = TRUE, warning = FALSE, message = FALSE----
gene$MeanExp <-
  rowMeans(assayDataElement(target_myData,
                            elt = "q_norm"))

top_g2 <- gene$Gene[gene$Gene %in% top_g &
                         gene$FDR < 0.001 &
                         abs(gene$Estimate) > .5 &
                         gene$MeanExp > quantile(gene$MeanExp, 0.9)]

ggplot(subset(gene, !Gene %in% neg_probes),
       aes(x = MeanExp, y = Estimate,
           size = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
  scale_x_continuous(trans = "log2") +
  geom_point(alpha = 0.5) + 
  labs(y = "Enriched in XXX <- log2(FC) -> Enriched in XXX",
       x = "Mean Expression",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray")) +
  geom_text_repel(data = subset(gene, Gene %in% top_g2),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2) +
  theme_bw(base_size = 16) +
  facet_wrap(~Subset, nrow = 2, ncol = 1)











vennCounts(gene, include="both")



library(ggVennDiagram)

dim(target_myData)
rownames(target_myData)

exprs(target_myData)[1:15, 1:3]


acini_bystander <- dplyr::filter(gene, Contrast == "1 - 2")
head(acini_bystander)

genes <- paste0("gene",1:1000)
set.seed(20210302)
gene_list <- list(A = gene1(genes,100),
                  B = sample(genes,200),
                  C = sample(genes,300),
                  D = sample(genes,200))

ggVennDiagram(acini_bystander,category.names = c("Stage 1","Stage 2","Stage 3", "Stage4"), label = "none")










