

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GeomxTools")

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.height = 4,
  dpi=200
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(NanoStringNCTools)
library(GeomxTools)
library(EnvStats)
library(ggiraph)

## ----buildobject--------------------------------------------------------------
datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

demoData <-
  suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                          pkcFiles = PKCFiles,
                                          phenoDataFile = SampleAnnotationFile,
                                          phenoDataSheet = "CW005",
                                          phenoDataDccColName = "Sample_ID",
                                          protocolDataColNames = c("aoi",
                                                                   "cell_line",
                                                                   "roi_rep",
                                                                   "pool_rep",
                                                                   "slide_rep")))

class(demoData)
isS4(demoData)
is(demoData, "ExpressionSet")
demoData

## ----countmatrix--------------------------------------------------------------
# access the count matrix 
assayData(demoData)[["exprs"]][1:3, 1:3]

# access pheno data
pData(demoData)[1:3, ]

# access the protocol data
pData(protocolData(demoData))[1:3, ]

# access the probe information
fData(demoData)[1:3, ]

# check feature type
featureType(demoData)

# access PKC information
annotation(demoData)

## ----accessobject-------------------------------------------------------------
svarLabels(demoData)
head(sData(demoData), 2)

## ----assigndesign-------------------------------------------------------------
design(demoData) <- ~ `segments`
design(demoData)

dimLabels(demoData)
dimLabels(demoData)[2] <- "Sample ID"
dimLabels(demoData)

## ----summaryobject------------------------------------------------------------
head(summary(demoData, MARGIN = 1), 2)
head(summary(demoData, MARGIN = 2), 2)
unique(sData(demoData)$"cell_line")
head(summary(demoData, MARGIN = 2, GROUP = "cell_line")$"HS578T", 2)
head(summary(demoData, MARGIN = 2, GROUP = "cell_line")$"COLO201", 2)
head(summary(demoData, MARGIN = 2, GROUP = "cell_line", log2 = FALSE)$"COLO201", 2)

## ----subsetobject-------------------------------------------------------------
dim(demoData)

## -----------------------------------------------------------------------------
dim(demoData[, demoData$`slide name` == "6panel-old-slide1 (PTL-10891)"])

## -----------------------------------------------------------------------------
dim(subset(demoData, select = phenoData(demoData)[["slide name"]] == "6panel-old-slide1 (PTL-10891)"))

## -----------------------------------------------------------------------------
dim(subset(demoData, TargetName == "ACTA2", `slide name` == "6panel-old-slide1 (PTL-10891)"))
dim(subset(demoData, CodeClass == "Control", `slide name` == "6panel-old-slide1 (PTL-10891)"))

## -----------------------------------------------------------------------------
dim(endogenousSubset(demoData))
dim(negativeControlSubset(demoData))

## -----------------------------------------------------------------------------
dim(endogenousSubset(demoData, 
                     select = phenoData(demoData)[["slide name"]] == "6panel-old-slide1 (PTL-10891)"))

# tally the number of samples according to their protocol or phenodata grouping
with(endogenousSubset(demoData), table(`slide name`))
with(demoData [1:10, 1:10], table(cell_line))    
with(negativeControlSubset(demoData), table(CodeClass))

## ----applyFunctions-----------------------------------------------------------
assayDataElement(demoData, "demoElem") <- 
  assayDataApply(demoData, MARGIN=2, FUN=log, base=10, elt="exprs")
assayDataElement(demoData, "demoElem")[1:3, 1:2]

# loop over the features(1) or samples(2) of the assayData element and get the mean
assayDataApply(demoData, MARGIN=1, FUN=mean, elt="demoElem")[1:5]

# split the data by group column with feature, pheno or protocol data then get the mean
head(esBy(demoData, 
          GROUP = "cell_line", 
          FUN = function(x) { 
            assayDataApply(x, MARGIN = 1, FUN=mean, elt="demoElem") 
          }))

## ----qcobject, eval = TRUE----------------------------------------------------
demoData <- shiftCountsOne(demoData, useDALogic=TRUE)
demoData <- setSegmentQCFlags(demoData)
head(protocolData(demoData)[["QCFlags"]])
demoData <- setBioProbeQCFlags(demoData)
featureData(demoData)[["QCFlags"]][1:5, 1:4]

## ----removeQCSampleProbe,  eval = TRUE----------------------------------------
QCResultsIndex <- which(apply(protocolData(demoData)[["QCFlags"]], 
                              1L , function(x) sum(x) == 0L))
QCPassed <- demoData[, QCResultsIndex]
dim(QCPassed) 

## ---- eval = TRUE-------------------------------------------------------------
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)

## ---- eval = TRUE-------------------------------------------------------------
featureType(target_demoData)
exprs(target_demoData)[1:5, 1:5]

## ----normalizeObject,  eval = TRUE--------------------------------------------
target_demoData <- normalize(target_demoData , data_type="RNA", norm_method="quant", 
                             desiredQuantile = .9, toElt = "q_norm")
target_demoData <- normalize(target_demoData , data_type="RNA", norm_method="neg", fromElt="exprs",  toElt="neg_norm")
target_demoData <- normalize(target_demoData , data_type="RNA", norm_method="hk", fromElt="exprs", toElt="hk_norm")
assayDataElement( target_demoData , elt = "q_norm" )[1:3, 1:2]
assayDataElement( target_demoData , elt = "hk_norm" )[1:3, 1:2]
assayDataElement( target_demoData , elt = "neg_norm" )[1:3, 1:2]

## ----mungeObject--------------------------------------------------------------
neg_set <- negativeControlSubset(demoData)
class(neg_set)
neg_ctrls <- munge(neg_set, ~ exprs)
head(neg_ctrls, 2)
class(neg_ctrls)
head(munge(demoData, ~ exprs), 2)
munge(demoData, mapping = ~`cell_line` + GeneMatrix)

## ----transformObject----------------------------------------------------------
thresh <- assayDataApply(negativeControlSubset(demoData), 2, max)
demoData <-
  transform(demoData,
            negCtrlZeroed = sweep(exprs, 2, thresh),
            log1p_negCtrlZeroed = log1p(pmax(negCtrlZeroed, 0)))
assayDataElementNames(demoData)


## -----------------------------------------------------------------------------
sessionInfo()
