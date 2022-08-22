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

library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)
#####
#####CALCULATE DE
#####CALCULATE DE
#####
load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_geoMX_new.RData")

# If comparing structures that co-exist within a given tissue, use an LMM model 
# with a random slope. Diagnosis is our test variable. We control for tissue 
# sub-sampling with slide name using a random slope and intercept; the intercept 
# adjusts for the multiple regions placed per unique tissue, since we have one 
# tissue per slide. If multiple tissues are placed per slide, we would change the
# intercept variable to the unique tissue name (ex: tissue name, Block ID, etc).

# convert test variables to factors
pData(target_myData)$testRegion <- 
  factor(pData(target_myData)$metastasis3)#, c("Stroma-PanIN","Stroma-nontum"))                           ###CHANGE
pData(target_myData)[["slide"]] <-                                            ### Control for 
  factor(pData(target_myData)[["tissue"]])
assayDataElement(object = target_myData, elt = "log_q") <-
  assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(status in c("Full ROI")) {
  ind <- pData(target_myData)$segment == status
  mixedOutmc <-
    mixedModelDE(target_myData[, ind], elt = "log_q",
                 #modelFormula = ~ testRegion + (1 + testRegion | slide),        ### =Reaction ~ Days + (Days || Subject), sleepstudy)
                 modelFormula = ~ testRegion + (1 | slide),
                 groupVar = "testRegion", nCores = parallel::detectCores(),
                 multiCore = FALSE)
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- status
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
  rm(r_test, tests, mixedOutmc)
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

results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)

###WRITE FILE
head(results)
names(results)[1] <- 'SYMBOL'
head(results)
eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
head(eg)
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)
head(results)

##Change FILENAME
#write.csv(results, "C:/Users/edmondsonef/Desktop/DSP GeoMx/07.08.22_class_MHL_no_int.csv")

#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.06.22_comps_MHL_no.int.csv")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.06.22_comps_MHL_WITH.int.csv")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.08.22_class_MHL_no_int.csv")

results1 <- dplyr::filter(results, abs(results$Estimate) > 0.5)
head(results1)
#names(results1)[6] <- 'Pr(>|t|)'
head(results1)
mt_list = split(results1, f = results1$Contrast)



names(mt_list)
gene <- mt_list[[1]]
head(gene)
top_g <- c()
for(cond in c("Full ROI")) {
  ind <- gene$Subset == cond
  top_g <- c(top_g,
             gene[ind, 'SYMBOL'][
               order(gene[ind, 'invert_P'], decreasing = TRUE)[1:100]],
             gene[ind, 'SYMBOL'][
               order(gene[ind, 'invert_P'], decreasing = FALSE)[1:30]])
}
top_g <- unique(top_g)
top_g

# head(results1)
# head(gene)
# names(gene)[5] <- 'Pr(>|t|)'

gene <- distinct(gene, SYMBOL, .keep_all = T)

pVP <- ggplot(gene,                                                             ###CHANGE
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = SYMBOL)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = " <- log2(FC) -> ",                                       ###CHANGE
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(gene, SYMBOL %in% top_g),# & FDR < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") 
pVP




#####

#gene <- dplyr::filter(results, abs(results$Estimate) > 0.5)
#gene <- dplyr::filter(gene, abs(gene$`Pr(>|t|)`) < 0.01) 
#gene <- dplyr::filter(results, abs(results$FDR) < 0.001)

head(gene)
#names(gene)[1] <- 'SYMBOL'
# head(gene)
# eg <- bitr(gene$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
#            OrgDb="org.Mm.eg.db")
# head(eg)
# gene <- dplyr::left_join(gene, eg, by = "SYMBOL")
# rm(eg)
# head(gene)
# results <- distinct(results, SYMBOL, .keep_all = T)
# universe <- bitr(results$SYMBOL, fromType="SYMBOL", 
#                  toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
#                  OrgDb="org.Mm.eg.db")
# head(universe)

# 
# ggo <- groupGO(gene     = gene$ENTREZID,
#                OrgDb    = org.Mm.eg.db,
#                ont      = "CC", #One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#                level    = 3,
#                readable = TRUE)


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

results <- distinct(results, SYMBOL, .keep_all = T)
universe <- bitr(results$SYMBOL, fromType="SYMBOL", 
                 toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
                 OrgDb="org.Mm.eg.db")
head(universe)
head(gene)

#####enrichGO
ego <- enrichGO(gene          = gene$ENTREZID,
                keyType       = "ENTREZID",
                universe      = universe$ENTREZID, ##list of all genes?? 
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", #"BP", "MF", and "CC"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
#head(ego)

goplot(ego)
pUpS <- upsetplot(ego, 10)
pUpS
pDP <- dotplot(ego)
pDP


wcdf<-read.table(text=ego$GeneRatio, sep = "/")[1]
wcdf$term<-ego[,2]
wcdf <- dplyr::filter(wcdf, V1 > 20)
head(wcdf)
wcdf <- dplyr::top_n(wcdf, 20, V1)
#wcdf = sort(wcdf, decreasing = T)

pWC <- ggplot(wcdf, aes(label = term, size = V1, 
                        color = factor(sample.int(10, nrow(wcdf), replace = TRUE)))) +
  geom_text_wordcloud() +
  #geom_text_wordcloud_area() +
  theme_minimal()

#cowplot::plot_grid(pVP, pDP, pUpS, pWC, ncol=2, labels=LETTERS[1:4])#, rel_widths=c(.8, .8, 1.2))

first_col <- plot_grid(pVP,labels = c('A'))
second_col <- plot_grid(pDP, pWC, ncol=1, labels = c('B','C'))
gg_all = plot_grid(first_col, second_col, labels=c('', ''), ncol=2)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Bystander_Acini.tiff", units="in", width=15, height=10, res=200)
gg_all
dev.off()

#####gseGO()
head(gene)
## assume 1st column is ID
## 2nd column is FC
## feature 1: numeric vector
geneList = gene[,5] #which column? 
head(geneList)

names(geneList) = as.character(gene[,10])
head(geneList)
geneList = sort(geneList, decreasing = T)
head(geneList)

#"BP" = biological process
#"MF" = molecular function
#"CC" = cellular component

ego3 <- gseGO(geneList     = geneList, ##??
              OrgDb        = org.Mm.eg.db,
              ont          = "MF", #"BP", "MF", and "CC"
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

#head(ego3)
goplot(ego3)
dotplot(ego3)
upsetplot(ego3, 10)


########
########
######## KEGG -- network issues
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

#get_wp_organisms()

geneList = gene[,6]
head(geneList)
names(geneList) = as.character(gene[,10])
head(geneList)
geneList = sort(geneList, decreasing = T)
head(geneList)

enrichWP(gene$ENTREZID, organism = "Mus musculus") 
gseWP(geneList$ENTREZID, organism = "Mus musculus")

gmt <- "https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-20220510-gmt-Mus_musculus.gmt"
wp <- read.gmt.wp(gmt)
ewp <- GSEA(geneList, TERM2GENE=wp[,c("wpid",
                                      "gene")], TERM2NAME=wp[,c("wpid", "name")])









########
########
########
######## Reactome
########
########
library(ReactomePA)
head(gene)
de <- gene$ENTREZID
head(de)
x <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)
head(x)




head(gene)
geneList = gene[,6]

names(geneList) = as.character(gene[,10])
head(geneList)
geneList = sort(geneList, decreasing = T)

head(geneList)

y <- gsePathway(geneList, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)
head(y)

viewPathway("E2F mediated regulation of DNA replication", 
            readable = TRUE, 
            foldChange = geneList)







########
########
######## Disease enrichment
########
########

head(gene$ENTREZID)
head(geneList$ENTREZID)

x <- enrichDO(gene          = gene$ENTREZID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = geneList$ENTREZID,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)


########
########
########
######## MeSH
########
########


library(enrichplot)
library(AnnotationHub)
library(MeSHDbi)

# Data source: 
# 1. gendoo
# 2. gene2pubmed
# 3. RBBH
#
# Category:
# C - Diseases
# G - Phenomena and Processes


ah <- AnnotationHub(localHub=TRUE)
mma <- query(ah, c("MeSHDb", "Mus musculus"))
file_mma <- mma[[1]]
db <- MeSHDbi::MeSHDb(file_mma)

de <- names(geneList)[1:100]
x <- enrichMeSH(de, MeSHDb = db, database='gendoo', category = 'C')
head(x)

y <- gseMeSH(geneList, MeSHDb = db, database = 'gene2pubmed', category = "G")
head(y)

########
########
########
######## MSigDb analysis
########
########

# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures

library(msigdbr)
msigdbr_show_species()

m_t2g <- "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/msigdb.v7.5.1.entrez.gmt"

# all gene sets
# m_df <- msigdbr(species = "Mus musculus")
# head(m_df, 2) %>% as.data.frame

m_t2g <- msigdbr(species = "Mus musculus", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

head(gene)
em <- enricher(gene$ENTREZID, TERM2GENE=m_t2g)
head(em)

em2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em2)

dotplot(em, showCategory=20) + ggtitle("dotplot for ORA")
dotplot(em2, showCategory=20) + ggtitle("dotplot for GSEA")


## convert gene ID to Symbol
edox <- setReadable(ego, 'org.Mm.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))


p1 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


########
########
########
######## Biological Theme Comparison
########
########

# Use a named list of gene IDs as the input that passed to the 
# geneCluster parameter.

#gcSample =  list of different samples
results2 <- dplyr::select(results1, Estimate, Contrast, ENTREZID)
results2 <- dplyr::filter(results2, Contrast == c("1-Normal acini - 3-ADM",
                                                 "3-ADM - 4-PanINlo",
                                                 "4-PanINlo - 5-PanINhi",
                                                 "5-PanINhi - 6-PDAC",
                                                 "4-PanINlo - 6-PDAC",
                                                 "5-PanINhi - 7-metastasis"))
                                                 
head(results2)
mt_list = split(results2, f = results2$Contrast)

results2$ENTREZID
str(mt_list)

ids <- unique(results2$Contrast)
mt_list1<-list()
for(i in 1:length(ids)){
  id <- ids[i]
  df <- dplyr::filter(results2, results2$Contrast == id)
  mt_list1[[i]]<-  df$ENTREZID
  }


mt_list1 <- dplyr::select(mt_list, ENTREZID)
str(mt_list1)

ck <- compareCluster(geneCluster = mt_list1, 
                     fun = enrichGO, #"groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
                     OrgDb = org.Mm.eg.db)
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 


dotplot(ck)
# As an alternaitve to using named list, the compareCluster() function also 
# supports passing a formula to describe more complicated experimental designs 
# (e.g., Gene ~ time + tx).



head(results2)
# geneList = results2[,1] #which column? 
# head(geneList)
# 
# names(geneList) = as.character(results2[,3])
# head(geneList)
# geneList = sort(geneList, decreasing = T)
# head(geneList)


# 
# mydf <- data.frame(Entrez=names(geneList), FC=geneList)
# mydf <- mydf[abs(mydf$FC) > 1,]
# mydf$group <- results$Contrast
# mydf$group[mydf$FC < 0] <- "downregulated"
# mydf$othergroup <- "A"
# mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(ENTREZID~Contrast, data=results2, 
                              fun="enrichGO", OrgDb = org.Mm.eg.db, 
                              keyType="ENTREZID")

head(formula_res)
dotplot(formula_res)
formula_res <- setReadable(formula_res, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(formula_res, node_label="category", 
         cex_label_category = 1.2) 
cnetplot(formula_res, node_label="gene", 
         cex_label_category = 1.2) 
cnetplot(formula_res, node_label="all", 
         cex_label_category = 1.2) 



########
########
########
######## Visualization
########
########

library(enrichplot)


universe <- bitr(results$Gene, fromType="SYMBOL", 
                 toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
                 OrgDb="org.Mm.eg.db")
head(universe)

#gene <- dplyr::filter(results, abs(results$Estimate) > 2)
#gene <- dplyr::filter(results, abs(results$FDR) < 0.05) 
gene <- dplyr::filter(results, abs(results$FDR) < 0.001)

head(gene)
names(gene)[1] <- 'SYMBOL'
head(gene)
eg <- bitr(gene$SYMBOL, fromType="SYMBOL", 
           toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
head(eg)
gene <- dplyr::left_join(gene, eg, by = "SYMBOL")
head(gene)

edo <- enrichDGN(gene$ENTREZID, universe = universe$ENTREZID)



## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Mm.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
