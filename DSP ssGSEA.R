# setwd('~/Projects/DSP - Hwang PDAC Round 2/new_DCCs/')
# load('ssGSEA_UpdatedWorkspace_5-23-21.RData')

## Load libraries
library(ggplot2)
library(ggrepel)
library(ggsci)
library(pheatmap)
library(GSVA)
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(clusterProfiler)


load("F:/GeoMX KPC/WTA_04122022/processed_data/KPC_geoMX_exp1.RData")
load("F:/GeoMX KPC/WTA_11232022/processed_data/KPC_geoMX_exp2.RData")

GO_rosetta <- read.csv("C:/Users/edmondsonef/Desktop/GO Terms.csv")


## Create datafram of 
# head(target_myData@assayData$exprs)
# head(target_myData@assayData$log_q)
# head(target_myData@assayData$q_norm)
new_dfs <- target_myData@assayData$log_q
row.names(new_dfs) <- bitr(row.names(new_dfs), fromType="SYMBOL", toType=c("ENTREZID"),OrgDb="org.Mm.eg.db")

eg <- bitr(row.names(new_dfs), fromType="SYMBOL", toType=c("ENTREZID"),OrgDb="org.Mm.eg.db")
row.names(eg) <- eg$SYMBOL
new_dfs <- merge(new_dfs, eg, by = 0)
#make entrezID the rownames
row.names(new_dfs) <- new_dfs$ENTREZID
new_dfs <- as.data.frame(new_dfs)
new_dfs <- subset(new_dfs, select = -c(Row.names, SYMBOL, ENTREZID))
rm(eg)


#Create list object containing a collection of gene sets defined as GO terms with annotated Entrez gene identifiers
goannot <- AnnotationDbi::select(org.Mm.eg.db, keys=keys(org.Mm.eg.db), columns="GO")
columns(org.Mm.eg.db)
goannot <- filter(goannot, ONTOLOGY == "BP")

genesbygo <- split(goannot$ENTREZID, goannot$GO)
length(genesbygo)

colnames(new_dfs)<-target_myData@phenoData@data$dx2
colnames(new_dfs)<-target_myData@phenoData@data$dx
head(new_dfs)

new_dfs <- as.matrix(new_dfs)

ssGSEA <- gsva(new_dfs,
               genesbygo,
               method="ssgsea", 
               min.sz = 5, 
               max.sz=500,
               verbose=FALSE, 
               parallel.sz=1)




colnames(ssGSEA)
rownames(ssGSEA)

#ssGSEA <- select(ssGSEA, -c("KPC(R172H)_Normal_.1"))

head(GO_rosetta)
rownames(GO_rosetta) <- GO_rosetta$GO_ID
ssGSEA <- merge(ssGSEA, GO_rosetta, by = 0)
row.names(ssGSEA) <- ssGSEA$GO_pathway
head(ssGSEA)
ssGSEA <- as.data.frame(ssGSEA)
ssGSEA <- subset(ssGSEA, select = -c(Row.names, GO_ID, GO_pathway))



#zscore the ssgsea output for comparative analysis
mat = (ssGSEA - rowMeans(ssGSEA))/(rowSds(as.matrix(ssGSEA)))[row(ssGSEA)]
mat <- as.data.frame(mat)
colnames(mat)
#mat <- mat[,c(1,7,8,9,15,16,18,19,29,30,32,33,39,41,43,44,45,51,53,55,56,61,62,63,67,69,72,77)]


mat <- as.data.frame(mat)
keywords <- c("axon","DNA", "nerve", "angiogenesis", "dendrite ", "myelin", "synapse", "neuron")
mat1 <- mat %>% filter(grepl(paste(keywords,collapse="|"), rownames(mat)))
rownames(mat1)
mat1 <- as.matrix(mat1)
mat1 <- subset(mat1, (rowMax(mat1) + abs(rowMin(mat1))) > 5.5)
rownames(mat1)

mat <- as.matrix(mat)
mat1 <- subset(mat, (rowMax(mat) + abs(rowMin(mat))) > 5.8)
rownames(mat1)


mat1 <- as.matrix(mat1)
Heatmap(mat1, col = colorRamp2(c(-2,0,2), c("orangered", "white", "purple")))





pheatmap(mat1,
         scale = "row", 
         show_rownames = T, show_colnames = T,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         #breaks = seq(-3, 3, 0.05),
         #color = colorRampPalette(c("orangered", "white", "purple"))(120),
         annotation_col = pData(target_myData)[,c("dx","Strain", "class")])
           #pData(target_myData)[, c("dx", "Sex","Class")])


















