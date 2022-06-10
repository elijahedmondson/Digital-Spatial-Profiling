library(GeoDiff)
library(dplyr)
library(ggplot2)
library(NanoStringNCTools)
library(GeomxTools)
library(Biobase)
library(reshape2)


library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)


myData <- readRDS(file = "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/my_data.rds")
target_myData <- readRDS(file = "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/target_myData.rds")

head(pData(target_myData))

table(pData(target_myData)$`MHL Number`)
table(pData(target_myData)$progression1)
table(target_myData$`MHL Number`, target_myData$progression1)


#Subset 
kidney <- target_myData[, which(target_myData$`MHL Number` == c("22003074", "22003075"))]

table(kidney$`MHL Number`, kidney$progression1)




object=c(target_myData, final, goi,PKCfilename,PKCFiles,pkcs)


load("C:/Users/edmondsonef/Desktop/myfile.RData")






