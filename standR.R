

output_prefix<-"Cheng_WTA1"
projectname<-"Cheng_WTA1"

datadir<-"F:/GeoMX KPC/Cheng_WTA1/raw_data"
DCCdir<-"DCC-20230817"
PKCfilename<-"Mm_R_NGS_WTA_v1.0.pkc"
WorkSheet<-"final2_full.xlsx"
final <- read_excel("F:/GeoMX KPC/Cheng_WTA1/raw_data/final2_full.xlsx")
output_dir<-"processed_data/"

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

load("F:/GeoMX KPC/Cheng_WTA1/processed_data/Cheng_WTA1_1_3_2024.RData")




library(standR)
library(SpatialExperiment)
library(limma)
library(edgeR)
library(tidyverse)
library(vissE)
library(GSEABase)
library(msigdb)
library(ggalluvial)




