
library(data.table)
library(fgsea)
library(ggplot2)

# load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_geoMX.RData")
# acini_bystander <- dplyr::filter(results, Contrast == "1 - 2")
# head(acini_bystander)




library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationHub)
library(GOSemSim)
library(clusterProfiler)


# The clusterProfiler package provides the bitr() and bitr_kegg() functions for converting ID types. 
# Both bitr() and bitr_kegg() support many species including model and many non-model organisms.

load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_seurat.RData")
head(target_myData@featureData@data$TargetName)


geneuniverse <- target_myData@featureData@data#$TargetName
geneuniverse <- dplyr::filter(geneuniverse, DetectionRate > 0.04)
head(geneuniverse)
geneuniverse <- geneuniverse$TargetName
geneList <- bitr(geneuniverse, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"), 
                 OrgDb="org.Mm.eg.db")
head(geneList)
geneList <- geneList$ENTREZID
head(geneList)

DEgenelist <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/data/progression1_MHLnumber.csv")
head(DEgenelist)
list <- DEgenelist$Gene
head(list)
eg <- bitr(list, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
head(eg)
gene <- eg$ENTREZID
head(gene)

 
# GO analyses (groupGO(), enrichGO() and gseGO()) support organisms that have an 
# OrgDb object available (see also session 2.2).

# If a user has GO annotation data (in a data.frame format with the first column 
# as gene ID and the second column as GO ID), they can use the enricher() and 
# GSEA() functions to perform an over-representation test and gene set enrichment 
# analysis.
# 
# If the genes are annotated by direction annotation, they should also be annotated 
# by their ancestor GO nodes (indirect annotation). If a user only has direct 
# annotation, they can pass their annotation to the buildGOmap function, which 
# will infer indirect annotation and generate a data.frame that is suitable for 
# both enricher() and GSEA().

#"BP" = biological process
#"MF" = molecular function
#"CC" = cellular component

#####groupGO
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Mm.eg.db,
               ont      = "MF", #One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
               level    = 3,
               readable = TRUE)
head(ggo)


#####enrichGO
ego <- enrichGO(gene          = gene,
                keyType = "ENTREZID",
                universe      = geneList, ##list of all genes?? 
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
goplot(ego)
dotplot(ego)

p1 <- dotplot(ego, showCategory = 10, font.size=14)
p2 <- dotplot(ego, showCategory = selected_pathways, font.size=14)


cowplot::plot_grid(p1, p2, labels=LETTERS[1:2])



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GOSemSim")




#####gseGO()
d <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/data/progression1_MHLnumber.csv")
head(d)
d <- dplyr::select(d, Gene, 'Pr...t..')
head(d)

d_list <- bitr(d[,2], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"),
           OrgDb="org.Mm.eg.db")
head(d_list)
names(d)[2] <- 'SYMBOL'

d_new <- dplyr::left_join(d, d_list, by = "SYMBOL")
head(d_new)

geneList = d_new[,6]
names(geneList) = as.character(d_new[,11])
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

ego3 <- gseGO(geneList     = geneList, ##??
              OrgDb        = org.Mm.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

head(ego3)




########
########
######## KEGG 
########
########
library(clusterProfiler)
search_kegg_organism('mmu', by='kegg_code')



load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_seurat.RData")
geneuniverse <- target_myData@featureData@data#$TargetName
geneuniverse <- dplyr::filter(geneuniverse, DetectionRate > 0.04)
head(geneuniverse)
geneuniverse <- geneuniverse$TargetName
geneList <- bitr(geneuniverse, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"), 
                 OrgDb="org.Mm.eg.db")
head(geneList)
geneList <- geneList$UNIPROT
head(geneList)

DEgenelist <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/data/progression1_MHLnumber.csv")
head(DEgenelist)
list <- DEgenelist$Gene
head(list)
eg <- bitr(list, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
head(eg)
gene <- eg$UNIPROT
head(gene)



kk <- enrichKEGG(gene         = gene,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
head(kk)





kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)




kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

########
########





























