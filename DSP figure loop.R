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

datadir <-"C:/Users/edmondsonef/Desktop/DSP GeoMX/Results/"
setwd(datadir)

results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.06.22_comps_MHL_no.int.csv")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.06.22_comps_MHL_WITH.int.csv")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.08.22_class_MHL_no_int.csv")

universe <- distinct(results, SYMBOL, .keep_all = T)
m_t2g <- "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/msigdb.v7.5.1.entrez.gmt"

results.sig <- dplyr::filter(results, abs(results$Estimate) > 0.5)
head(results.sig)
names(results.sig)[6] <- 'Pr(>|t|)'
head(results.sig)

mt_list = split(results.sig, f = results.sig$Contrast)
names(mt_list)

#gmt <- "https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-20220510-gmt-Mus_musculus.gmt"
#wp <- read.gmt.wp(gmt)


#options(warn = 0)

names(mt_list)


##FOR LOOP
i = 28
for(i in 28:31){
  suffix <- names(mt_list[i])
  outname <-paste0(suffix, "NEW_comps_MHL_no_int")
  
  gene <- mt_list[[i]]
  gene <- distinct(gene, SYMBOL, .keep_all = T)

  top <- dplyr::filter(gene, gene$FDR < 0.05)
  count = count(top)
  print(paste(suffix, ":",count, "genes with FDR < 0.05."))
  
  write.csv(top, paste0(datadir, outname,"_", count," Genes",".csv"))
  #rm(count, top)

  top_g <- c()
  for(cond in c("Full ROI")) {
    ind <- gene$Subset == cond
    top_g <- c(top_g,
               gene[ind, 'SYMBOL'][
                 order(gene[ind, 'invert_P'], decreasing = TRUE)[1:30]],
               gene[ind, 'SYMBOL'][
                 order(gene[ind, 'invert_P'], decreasing = FALSE)[1:30]])
  }
  top_g <- unique(top_g)
  top_g

  #reverse log fold change to fit with label
  gene$Estimate1 <- gene$Estimate*(-1)
  
  
  pVP <- ggplot(gene,                                                             
                aes(x = Estimate1, y = -log10(`Pr(>|t|)`),
                    color = Color, label = SYMBOL)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x =  suffix,#"<- log2(FC) ->",                                       
         y = "Significance, -log10(P)",
         color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(gene, SYMBOL %in% top_g & `Pr(>|t|)` < 0.05),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,
                    max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") 
  
  volcano <- paste0(datadir, "_", outname, "_volcano.png")
  ggsave(pVP, file=volcano, width = 8, height = 8, units = "in", bg = "white")
  rm(volcano)

  ego <- enrichGO(gene          = gene$ENTREZID,
                  keyType       = "ENTREZID",
                  universe      = universe$ENTREZID, ##list of all genes?? 
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP", #"BP", "MF", and "CC"
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  pDP <- dotplot(ego)
  wcdf<-read.table(text=ego$GeneRatio, sep = "/")[1]
  wcdf$term<-ego[,2]
  wcdf <- dplyr::filter(wcdf, V1 > 20)
  wcdf <- dplyr::top_n(wcdf, 25, V1)
  pWC <- ggplot(wcdf, aes(label = term, size = V1, 
                          color = factor(sample.int(10, nrow(wcdf), replace = TRUE)))) +
    geom_text_wordcloud() +
    theme_minimal()
  
  first_col <- plot_grid(pVP,labels = c('A'))
  second_col <- plot_grid(pDP, pWC, ncol=1, labels = c('B','C'))
  gg_all = plot_grid(first_col, second_col, labels=c('', ''), ncol=2)
  
  multiplot <- paste0(datadir, "_", outname, "_volcano_enrighGO_BP.png")
  # tiff(multiplot, units="in", width=15, height=10, res=175)
  # gg_all
  # dev.off()
  
  ggsave(gg_all, file=multiplot, width = 15, height = 10, units = "in")
  
  
  rm(multiplot, first_col, second_col, gg_all, pDP, wcdf, ego, pVP, pWC)
  
  #####gseGO()
  head(gene)
  gene <- distinct(gene, SYMBOL, .keep_all = T)
  geneList = gene$Estimate1 
  names(geneList) = as.character(gene$ENTREZID)
  geneList = sort(geneList, decreasing = T)
  
  ego <- gseGO(geneList      = geneList, 
                OrgDb        = org.Mm.eg.db,
                ont          = "BP", #"BP", "MF", and "CC"
                minGSSize    = 20,
                maxGSSize    = 1000,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
  p1 <- goplot(ego)
  p2 <- dotplot(ego)
  p3 <- upsetplot(ego, 10)
  gg_all <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
  
  multiplot <- paste0(datadir, "_", outname, "_gseGO_BP.png")
  ggsave(gg_all, file=multiplot, width = 15, height = 15, units = "in", bg = "white")
  rm(ego,gg_all, p1, p2, p3)

  ######## WikiPathways
# 
#   ego <- enrichWP(gene$ENTREZID, organism = "Mus musculus")
#   #ego2 <- gseWP(geneList, organism = "Mus musculus")
# 
# 
#   ewp <- GSEA(geneList, TERM2GENE=wp[,c("wpid","gene")],
#               TERM2NAME=wp[,c("wpid", "name")])
# 
# 
#   p1 <- dotplot(ewp)
#   p2 <- dotplot(ego)
#   gg_all <- cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:3])
# 
#   multiplot <- paste0(datadir, outname, "_wikipathways.png")
#   ggsave(gg_all, file=multiplot, width = 10, height = 10, units = "in", bg = "white")
#   rm(gg_all, ego, ego2, ewp, p1, p2)

  ###MSigDb analysis
  ###MSigDb analysis
  ###MSigDb analysis
  
  
  # H: hallmark gene sets       --good
  # C1: positional gene sets
  # C2: curated gene sets       --good
  # C3: motif gene sets         --good
  # C4: computational gene sets
  # C5: GO gene sets            --good
  # C6: oncogenic signatures    --good
  # C7: immunologic signatures  --good
  msig_list <- c("C2", "C6", "C7", "H", "C5", "C3")
  
  for(j in 1:6){
    Msig <- msig_list[j]
    m_t2g <- msigdbr(species = "Mus musculus", category = Msig) %>% 
      dplyr::select(gs_name, entrez_gene)
    
    edo <- enricher(gene$ENTREZID, TERM2GENE=m_t2g)
    edo2 <- GSEA(geneList, TERM2GENE = m_t2g, eps=0)
    
    p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
    p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
    gg_all <- cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])
    
    multiplot <- paste0(datadir, "_", outname, "_", Msig, "_GSEA_ORA_GSEA.png")
    ggsave(gg_all, file=multiplot, width = 15, height = 10, units = "in", bg = "white")
    rm(multiplot, gg_all, p1, p2)
    
    ## convert gene ID to Symbol
    edox <- setReadable(edo2, 'org.Mm.eg.db', 'ENTREZID')
    p1 <- cnetplot(edox, categorySize="pvalue", node_label="category", foldChange=geneList)
    p2 <- cnetplot(edox, categorySize="pvalue", node_label="gene", foldChange=geneList)
    p3 <- cnetplot(edox, foldChange=geneList, node_label="gene", circular = TRUE, colorEdge = TRUE) 
    gg_all <- cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
      
    multiplot <- paste0(datadir, "_", outname, "_", Msig, "_GSEA_cnetplot.png")
    ggsave(gg_all, file=multiplot, width = 30, height = 10, units = "in", bg = "white")
    rm(p1, p2, p3, gg_all, multiplot)
    
    multiplot <- paste0(datadir, "_", outname, "_", Msig, "_GSEA_heatplot.png")
    heat <- heatplot(edox, foldChange=geneList, showCategory=5)
    ggsave(heat, file=multiplot, width = 20, height = 5, units = "in", bg = "white")
    rm(multiplot, heat)
    
    edox2 <- pairwise_termsim(edox)
    
    multiplot <- paste0(datadir, "_", outname, "_", Msig, "_GSEA_treeplot.png")
  
    tree <- treeplot(edox2, hclust_method = "average")
    ggsave(tree, file=multiplot, width = 20, height = 10, units = "in", bg = "white")
    rm(multiplot, tree)
    
    multiplot <- paste0(datadir, "_", outname, "_", Msig, "_GSEA_gseaplot.png")

    plotty <- gseaplot2(edo2, geneSetID = 1:5)
    ggsave(plotty, file=multiplot, width = 10, height = 5, units = "in", bg = "white")
    rm(multiplot, plotty, edo, edo2, edox, edox2)}
}
  
