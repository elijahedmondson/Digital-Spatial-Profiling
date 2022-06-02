library(NanoStringNCTools)
library(GeomxTools)
knitr::opts_chunk$set(echo = TRUE)
output_prefix<-"CPTR474"
projectname<-"CPTR474"
datadir<-"C:/Users/edmondsonef/Desktop/DSP GeoMX/data/WTA_04122022/raw_data"
DCCdir<-"DCC-20220420"
PKCfilename<-"Mm_R_NGS_WTA_v1.0.pkc"
WorkSheet<-"20220414T0150_efe.xlsx"


DCCFiles <- list.files(file.path(datadir , DCCdir), pattern=".dcc$", full.names=TRUE)
PKCFiles <- file.path(datadir, PKCfilename)
SampleAnnotationFile <- file.path(datadir, WorkSheet)

myData<-readNanoStringGeoMxSet(dccFiles = DCCFiles,
                               pkcFiles = PKCFiles,
                               phenoDataFile = SampleAnnotationFile,
                               phenoDataSheet = "Template",
                               phenoDataDccColName = "Sample_ID",
                               protocolDataColNames = c("aoi", 
                                                        "roi"),
                               experimentDataColNames = c("panel")) 

#Shift counts to one to mimic how DSPDA handles zero counts
myData <- shiftCountsOne(myData, elt="exprs", useDALogic=TRUE) 
pkcs <- annotation(myData)
modules <- gsub(".pkc", "", pkcs)


# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned to known targets (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 10,   # Minimum negative control counts (10)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of cells observed in a segment (100)
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



library(ggplot2)

col_by <- "class"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         xlims = NULL) {
  if(is.null(xlims)) {
    xlims <- range(assay_data[[annotation]])
  }
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    xlim(xlims) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  plt
}

library(knitr)
QC_histogram(sData(myData), "Trimmed (%)", col_by, 80, c(0,101))
QC_histogram(sData(myData), "Stitched (%)", col_by, 80, c(0,101))
QC_histogram(sData(myData), "Aligned (%)", col_by, 75, c(0,101))
QC_histogram(sData(myData), "Saturated (%)", col_by, 50, c(0,101)) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
QC_histogram(sData(myData), "area", col_by, 1000) +
  scale_x_continuous(trans = "log10")
QC_histogram(sData(myData), "nuclei", col_by, 20)

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
  plt <- QC_histogram(pData(myData), ann, col_by, 2) +
    scale_x_continuous(trans = "log10") 
  print(plt)
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(myData) <- pData(myData)[, !colnames(pData(myData)) %in% negCols]


# show all NTC values, Freq = # of Segments with a given NTC count (for WTA samples this is the deduplicated read count from the empty well A01):
kable(table(NTC_Count = sData(myData)$NTC),
      col.names = c("NTC Count", "# of Segments"), caption = "NTC count summary")

kable(QC_Summary, caption = "QC Summary Table for each Segment")

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


kable(qc_df, caption = "Probes flagged or passed as outliers")


#Subset object to exclude all that did not pass Ratio & Global testing
cat("Before excluding probes")
dim(myData)

ProbeQCPassed <- 
  subset(myData, 
         fData(myData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(myData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

cat("After excluding probes")

myData <- ProbeQCPassed 
dim(myData)

# Check how many unique targets the object has
length(unique(featureData(myData)[["TargetName"]]))

# collapse to targets
target_myData <- aggregateCounts(myData)
raw_counts<-data.frame(exprs(target_myData))
raw_counts <- tibble::rownames_to_column(raw_counts, "Gene")
out_raw_counts<-paste0(output_prefix,'_raw_counts.csv')
write.table(raw_counts,out_raw_counts,sep=",",row.names = F,col.names=T,quote=F)

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

# Save detection rate information to pheno data
pData(target_myData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_myData)$GeneDetectionRate <-
  pData(target_myData)$GenesDetected / nrow(target_myData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_myData)$DetectionThreshold <- 
  cut(pData(target_myData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked barplot of different cutpoints (1%, 5%, 10%, 15%)
ggplot(pData(target_myData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = `class`)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment class")

# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_myData)$DetectionThreshold,
            pData(target_myData)$class))


geneDetectionRateThresh<-0.05

target_myData <-
  target_myData[, pData(target_myData)$GeneDetectionRate >= geneDetectionRateThresh]

dim(target_myData)


library(scales) # for precent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_myData)]
fData(target_myData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_myData)$DetectionRate <-
  fData(target_myData)$DetectedSegments / nrow(pData(target_myData))

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
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


# Subset to target genes detected in at least 1 of the samples.
#   Also manually include the negative control probe, for downstream use
targetDetectionRateThresh<-0.10
negativeProbefData <- subset(fData(target_myData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_myData <- 
  target_myData[fData(target_myData)$DetectionRate > targetDetectionRateThresh |
                  fData(target_myData)$TargetName %in% neg_probes, ]

"Remaining targets and ROIs"
dim(target_myData)


library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "class"
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
  facet_wrap(~Annotation, ncol = 3) + 
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

# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_myData <- normalize(target_myData , data_type = "RNA",
                           norm_method = "quant", 
                           desiredQuantile = .75,
                           toElt = "q_norm")

## save q3 data as csv file

norm_data<-data.frame(assayDataElement(target_myData , elt = "q_norm"))
norm_data <- tibble::rownames_to_column(norm_data, "Gene")
out_norm_data<-paste0(output_prefix,'_Q3_normdata.',
                      'geneDetect',
                      geneDetectionRateThresh,
                      '_',
                      'targetDetect',
                      targetDetectionRateThresh,
                      '.csv')

write.table(norm_data,out_norm_data,sep=",",row.names = F,col.names=T,quote=F)
# lets also save the entire targetmyData as an rstructure
save(target_myData, file = paste0(output_prefix,".RData"))

# visualize the first 10 segments with each normalization method
boxplot(exprs(target_myData)[,1:20],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:20, xlab = "Segment",
        ylab = "Counts, Raw",
        ylim=c(1,40000))

boxplot(assayDataElement(target_myData[,1:20], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:20, xlab = "Segment",
        ylab = "Counts, Q3 Normalized",
        ylim=c(1,15000))








library(umap)
library(Rtsne)

myshapes<-c(16,17,18,15,21,22,3,42,4,8,10,19)
myshapes<-c(1,2,16,4,17,6,7,8,18,15,19,12)#3,59,10,11,
# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_myData , elt = "q_norm"))),  
       config = custom_umap)
pData(target_myData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
pUMAP <- ggplot(pData(target_myData),
       aes(x = UMAP1, y = UMAP2, color = dx, label=dxf)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  scale_shape_manual(values=myshapes) +
  theme_bw()+
  theme(legend.position="none")
pUMAP

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("UMAP.tiff", units="in", width=19, height=12, res=150)
pUMAP
dev.off()


# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_myData , elt = "q_norm"))),
        perplexity = ncol(target_myData)*.15)
pData(target_myData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
pTSNE <- ggplot(pData(target_myData),
       aes(x = tSNE1, y = tSNE2, color = dx, label=dxf)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  scale_shape_manual(values=myshapes) +
  theme_bw()+
  theme(legend.position="none")
pTSNE

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("TSNE.tiff", units="in", width=19, height=12, res=150)
pTSNE
dev.off()

## run PCA
PCAx<-1
PCAy<-2
PCAxy <- c(as.integer( PCAx ),as.integer( PCAy) ) # selected principal components


pca.object <- prcomp(t(log2(assayDataElement(target_myData , elt = "q_norm"))))
pcaData = as.data.frame(pca.object$x[, PCAxy]); 
pData(target_myData)[, c("PC1", "PC2")] <- pcaData[,c(1,2)]
percentVar=round(100*summary(pca.object)$importance[2, PCAxy],0)


pPCA <- ggplot(pData(target_myData),
       aes(x = PC1, y = PC2, color=dx, label=dxf)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  xlab(paste0("PC", PCAx ,": ", percentVar[1], "% variance")) +
  ylab(paste0("PC", PCAy ,": ", percentVar[2], "% variance")) +
  scale_shape_manual(values=myshapes) +
  theme_bw()+
  theme(legend.position="none")
pPCA

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("PCA.tiff", units="in", width=19, height=12, res=150)
pPCA
dev.off()

library(pheatmap)  # for pheatmap
# create a log2 transform of the data for analysis
assayDataElement(object = target_myData, elt = "log_q") <-
  assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_myData,
                         elt = "log_q", MARGIN = 1, calc_CV)

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]


pheatmap(assayDataElement(target_myData[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_myData)[, c("class", "dx")])


library(NanoStringNCTools)
library(GeomxTools)
library(EnvStats)
library(ggiraph)
library(ggrepel) 
library(dplyr)
library(ggforce)

# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols
count_mat <- count(pData(demoData), `slide name`, class, region, segment)
# simplify the slide names
count_mat$`slide name` <- gsub("disease", "d",
                               gsub("normal", "n", count_mat$`slide name`))
# gather the data and plot in order: class, slide name, region, segment
test_gr <- gather_set_data(count_mat, 1:4)
test_gr$x <- factor(test_gr$x,
                    levels = c("class", "slide name", "region", "segment"))
# plot Sankey
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = region), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") +
  annotate(geom = "segment", x = 4.25, xend = 4.25,
           y = 20, yend = 120, lwd = 2) +
  annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
           hjust = 0.5, label = "100 segments")

