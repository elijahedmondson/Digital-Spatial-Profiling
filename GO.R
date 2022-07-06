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

#####
#####CALCULATE DE
#####CALCULATE DE
#####
load("C:/Users/edmondsonef/Desktop/DSP GeoMx/KPC_geoMX_new.RData")

# If comparing structures that co-exist within a given tissue, use an LMM model 
# with a random slope. Diagnosis is our test variable. We control for tissue 
# sub-sampling with slide name using a random slope and intercept; the intercept 
# adjusts for the multiple regions placed per unique tissue, since we have one 
# tissue per slide. If multiple tissues are placed per slide, we would change the
# intercept variable to the unique tissue name (ex: tissue name, Block ID, etc).

# convert test variables to factors
pData(target_myData)$testRegion <- 
  factor(pData(target_myData)$comps)#, c("Stroma-PanIN","Stroma-nontum"))                           ###CHANGE
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
    mixedModelDE(target_myData[, ind], elt = "log_q",
                 modelFormula = ~ testRegion + (1 + testRegion | slide),        ### modelFormula =  Reaction ~ Days + (Days || Subject), sleepstudy)
                 #modelFormula = ~ testRegion + (1 | slide),
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
write.csv(results, "C:/Users/edmondsonef/Desktop/DSP GeoMx/07.06.22_comps_MHL_no.int.csv")


top_g <- c()
for(cond in c("Full ROI")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'SYMBOL'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:30]],
             results[ind, 'SYMBOL'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:30]])
}
top_g <- unique(top_g)
top_g


results1 <- dplyr::filter(results, abs(results$Estimate) > 0.5)
ids <- unique(results1$Contrast)

for(i in 1:length(ids)){
  id <- ids[i]
  mini.df <- data.frame(results1[results1$Contrast == id, ])
  assign(paste("df", id, sep="."), mini.df)
}


######



gene <- `df.4-PanINlo - 5-PanINhi`
head(gene)
top_g <- c()
for(cond in c("Full ROI")) {
  ind <- df$Subset == cond
  top_g <- c(top_g,
             df[ind, 'gene'][
               order(df[ind, 'invert_P'], decreasing = TRUE)[1:30]],
             df[ind, 'gene'][
               order(df[ind, 'invert_P'], decreasing = FALSE)[1:30]])
}
top_g <- unique(top_g)
top_g

head(results1)
head(gene)
names(gene)[5] <- 'Pr(>|t|)'


pVP <- ggplot(gene,                                                             ###CHANGE
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = SYMBOL)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "PanIN-hi <- log2(FC) -> PanIN-lo",                                       ###CHANGE
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(gene, SYMBOL %in% top_g & FDR < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") 





#####

#gene <- dplyr::filter(results, abs(results$Estimate) > 0.5)
#gene <- dplyr::filter(gene, abs(gene$`Pr(>|t|)`) < 0.01) 
#gene <- dplyr::filter(results, abs(results$FDR) < 0.001)

head(gene)
names(gene)[1] <- 'SYMBOL'
head(gene)
eg <- bitr(gene$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
head(eg)
gene <- dplyr::left_join(gene, eg, by = "SYMBOL")
rm(eg)
head(gene)

universe <- bitr(results$Gene, fromType="SYMBOL", 
                 toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
                 OrgDb="org.Mm.eg.db")
head(universe)


ggo <- groupGO(gene     = gene$ENTREZID,
               OrgDb    = org.Mm.eg.db,
               ont      = "CC", #One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
               level    = 3,
               readable = TRUE)


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

universe <- bitr(results$SYMBOL, fromType="SYMBOL", 
                 toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
                 OrgDb="org.Mm.eg.db")
head(universe)


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
pUpS <- upsetplot(ego, 8)

pDP <- dotplot(ego)
pDP


wcdf<-read.table(text=ego$GeneRatio, sep = "/")[1]
wcdf$term<-ego[,2]


pWC <- ggplot(wcdf, aes(label = term, size = V1, 
                        color = factor(sample.int(10, nrow(wcdf), replace = TRUE)))) +
  geom_text_wordcloud_area() +
  theme_minimal()




cowplot::plot_grid(pVP, pDP, pUpS, pWC, ncol=2, labels=LETTERS[1:4])#, rel_widths=c(.8, .8, 1.2))

first_col <- plot_grid(pVP,labels = c('A'))
second_col <- plot_grid(pDP, pWC, ncol=1, labels = c('B','C'))
gg_all = plot_grid(first_col, second_col, labels=c('', ''), ncol=2)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("PanINhi_PanINlow.tiff", units="in", width=12, height=9, res=300)
gg_all
dev.off()

#####gseGO()
head(gene)
## assume 1st column is ID
## 2nd column is FC
## feature 1: numeric vector
geneList = gene[,6]
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
