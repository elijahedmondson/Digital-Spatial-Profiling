
library(data.table)
library(fgsea)
library(ggplot2)

#####
#####CALCULATE DE
#####CALCULATE DE
load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_geoMX.RData")
# If comparing structures that co-exist within a given tissue, use an LMM model 
# with a random slope. Diagnosis is our test variable. We control for tissue 
# sub-sampling with slide name using a random slope and intercept; the intercept adjusts for the multiple 
# regions placed per unique tissue, since we have one tissue per slide. If multiple tissues are placed per slide, 
# we would change the intercept variable to the unique tissue name (ex: tissue name, Block ID, etc).

# convert test variables to factors
pData(target_myData)$testRegion <- 
  factor(pData(target_myData)$prog4)                         ###CHANGE
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
                 #modelFormula = ~ testRegion + (1 + testRegion | slide),
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


library(ggrepel) 
# Categorize Results based on P-value & FDR for plotting
results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in c("Full ROI")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:30]],
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:30]])
}
top_g <- unique(top_g)
#results <- results[, -1*ncol(results)] # remove invert_P from matrix






#####
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationHub)
library(GOSemSim)
library(clusterProfiler)

library(GOSemSim)

## The clusterProfiler package provides the bitr() and bitr_kegg() functions for converting ID types. 
## Both bitr() and bitr_kegg() support many species including model and many non-model organisms.

# load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_seurat.RData")
# head(target_myData@featureData@data$TargetName)
# 
# 
# geneuniverse <- target_myData@featureData@data#$TargetName
# geneuniverse <- dplyr::filter(geneuniverse, DetectionRate > 0.04)
# head(geneuniverse)
# geneuniverse <- geneuniverse$TargetName
# geneList <- bitr(geneuniverse, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"), 
#                  OrgDb="org.Mm.eg.db")
# head(geneList)
# geneList <- geneList$ENTREZID
# head(geneList)
# 
# DEgenelist <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/data/progression1_MHLnumber.csv")
# head(DEgenelist)
# list <- DEgenelist$Gene
# head(list)
# eg <- bitr(list, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
#            OrgDb="org.Mm.eg.db")
# head(eg)
# gene <- eg$ENTREZID
# head(gene)

gene <- dplyr::filter(results, abs(results$Estimate) > 2)
gene <- dplyr::filter(results, abs(results$FDR) < 0.05) 
gene <- dplyr::filter(results, abs(results$FDR) < 0.001)

head(gene)
names(gene)[1] <- 'SYMBOL'
head(gene)
eg <- bitr(gene$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
head(eg)
gene <- dplyr::left_join(gene, eg, by = "SYMBOL")
head(gene)

ggo <- groupGO(gene     = gene$ENTREZID,
               OrgDb    = org.Mm.eg.db,
               ont      = "CC", #One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
               level    = 3,
               readable = TRUE)

head(ggo)

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

geneList <- bitr(results$Gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
head(geneList)


#####enrichGO
ego <- enrichGO(gene          = gene$ENTREZID,
                keyType = "ENTREZID",
                universe      = geneList$ENTREZID, ##list of all genes?? 
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
goplot(ego)
dotplot(ego)

# p1 <- dotplot(ego, showCategory = 10, font.size=14)
# p2 <- dotplot(ego, showCategory = selected_pathways, font.size=14)
# cowplot::plot_grid(p1, p2, labels=LETTERS[1:2])








#####gseGO()
head(gene)
## assume 1st column is ID
## 2nd column is FC
## feature 1: numeric vector
geneList = gene[,6]
head(geneList)

## feature 2: named vector
names(geneList) = as.character(gene[,10])
head(geneList)
## feature 3: decreasing order
geneList = sort(geneList, decreasing = T)
head(geneList)

#"BP" = biological process
#"MF" = molecular function
#"CC" = cellular component

ego3 <- gseGO(geneList     = geneList, ##??
              OrgDb        = org.Mm.eg.db,
              ont          = "MF",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

head(ego3)
goplot(ego3)
dotplot(ego3)



########
########
######## KEGG 
########
########
library(clusterProfiler)
search_kegg_organism('mmu', by='kegg_code')




kk <- enrichKEGG(gene         = gene$UNIPROT,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
head(kk)





kk2 <- gseKEGG(geneList     = gene$UNIPROT,
               organism     = 'mmu',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)




kk2 <- gseKEGG(geneList     = gene$UNIPROT,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

########
########
########
######## WikiPathways
########
########

get_wp_organisms()

geneList <- bitr(results$Gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
                 OrgDb="org.Mm.eg.db")
head(geneList)
head(gene)

enrichWP(gene$ENTREZID, organism = "Mus musculus") 
gseWP(geneList$ENTREZID, organism = "Mus musculus")











########
########
########
######## Reactome
########
########
library(ReactomePA)

x <- enrichPathway(gene=gene$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
head(x)















