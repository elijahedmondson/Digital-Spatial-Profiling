setwd('~/Projects/DSP - Hwang PDAC Round 2/new_DCCs/')
load('ssGSEA_UpdatedWorkspace_5-23-21.RData')

## Load libraries
library(ggplot2)
library(ggrepel)
library(ggsci)
library(pheatmap)

install.packages("BiocManager")
BiocManager::install("GSVA")

library(GSVA)
#https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html#:~:text=Gene%20set%20variation%20analysis%20(GSVA,from%20genes%20to%20gene%20sets.

#### Section 6: ssGSEA analysis & comparison ####
# rerun ssGSEA and see if it holds the same or if it changes dramatically
#

new_dfs <- dfs
new_dfs[[1]] <- dfs[[1]][, c('TargetName', colnames(df_clean))]
new_dfs[[2]] <- cbind(data.frame(Gene = rownames(df_clean)), df_clean)
new_dfs[[3]] <- dfs[[3]][colnames(df_clean), ]

# Convert to entrez id and set as rownames
norm_data <- new_dfs[[2]]
norm_data <- symbol_to_entrez(path_df=norm_data, gene_col="Gene", species = species)
rownames(norm_data) <- norm_data$entrez
norm_data <- as.matrix(norm_data[, -c(1, dim(norm_data)[2])])

# # Convert de results to entrez id
# de <- de_results
# de <- symbol_to_entrez(path_df=de, gene_col="gene", species = species)
gene_ids <- as.character(rownames(norm_data))
ssGSEA_dir <- paste(outdir, "ssGSEA_Detrend", sep = "/")
dir.create(ssGSEA_dir)
ssGSEA_detrend <- gsva(norm_data,
                       pathways,
                       method="ssgsea", min.sz = 5, max.sz=500,
                       verbose=FALSE, parallel.sz=1)

mal_IDs <- colnames(ssGSEA_detrend)[colnames(ssGSEA_detrend) %in% ss_match$Epi_ROIID]
caf_IDs <- colnames(ssGSEA_detrend)[colnames(ssGSEA_detrend) %in% ss_match$CAF_ROIID]

smoothScatter(ssGSEA_results[mal_paths, mal_IDs],
              ssGSEA_detrend[mal_paths, mal_IDs])

as.data.frame(diag(cor(t(ssGSEA_results[mal_paths, mal_IDs]),
                       t(ssGSEA_detrend[mal_paths, mal_IDs]))))

smoothScatter(ssGSEA_results[caf_paths, caf_IDs],
              ssGSEA_detrend[caf_paths, caf_IDs])

diag(cor(t(ssGSEA_results[caf_paths, caf_IDs]),
         t(ssGSEA_detrend[caf_paths, caf_IDs])))
plot(ssGSEA_results['CAF_Immunomodulatory', caf_IDs],
     ssGSEA_detrend['CAF_Immunomodulatory', caf_IDs])

nepi_ROIs <- new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'Epithelial']
ncaf_ROIs <- new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']
use_ROIs <- nepi_ROIs
for(topic in c(mal_paths, caf_paths)) {
  if(topic %in% caf_paths) {
    use_ROIs <- ncaf_ROIs
  }
  mod <- lmer(ssGSEA_results[topic, use_ROIs] ~ treatment +
                (1|Patient), data = dfs[[3]][use_ROIs, ])
  mod_sum <- ls_means(mod, pairwise = TRUE)
  add_p <- FALSE
  if(any(mod_sum[,7] < 0.1)) {
    add_p <- TRUE
    mod_sum$plt <- mod_sum[,7] < 0.1
    mod_sum$sig <- cut(mod_sum[,7], breaks = c(0,0.01,0.05,1),
                       labels = c("**", "*", ""))
    mod_sum$x1 <- c(1,1,1,2,2,3)
    mod_sum$x2 <- c(2,3,4,3,4,4)
    mod_sum <- as.data.frame(mod_sum)
    y_range <- diff(range(scale(ssGSEA_results[topic, use_ROIs])))
    y_max <- max(scale(ssGSEA_results[topic, use_ROIs]))
    mod_sum <- subset(mod_sum, plt)
    for(i in 1:nrow(mod_sum)) {
      mod_sum[i, 'y'] <- y_range*(i * 0.065) + y_max
    }
  }
  
  plt <- ggplot(dfs[[3]][use_ROIs, ],
                aes(x = treatment, fill = treatment)) +
    geom_boxplot(aes(y = scale(ssGSEA_results[topic, use_ROIs]))) +
    theme_bw() +
    scale_fill_jama() +
    labs(x = 'Treatment Group',
         y = paste0(gsub('Malignant_', '', topic), ' Enrichment'),
         title = paste0(gsub('Malignant_', '', topic), ' NMF Score'))
  if(add_p) {
    plt <- plt +
      annotate('segment',
               x = mod_sum$x1, xend = mod_sum$x2,
               y = mod_sum$y, yend = mod_sum$y) +
      annotate('text', x = rowMeans(mod_sum[, c('x1','x2')]),
               y = mod_sum$y + y_range*.01,
               label = paste0("P = ", signif(mod_sum[, 9], 3), mod_sum[, 'sig']),
               hjust = 0.5, vjust = 0, size = 3)
  }
  print(plt)
  ggsave(paste0('ssGSEA_Detrend/', gsub('\\/','-',topic), '_TreatmentAnalysis.png'), plt, 'png',
         width = 5, height = 6)
}

saveRDS(new_dfs, 'dfs_detrendApproach21-6-2.RDS')
saveRDS(ssGSEA_detrend, 'ssGSEA_detrendApproach21-6-2.RDS')


#re-exponentiate
new_dfs_exp <- new_dfs
new_dfs_exp[[2]] <- cbind(new_dfs_exp[[2]][,1], 2^new_dfs_exp[[2]][,-1])

saveRDS(new_dfs_exp, 'dfs_detrendApproach21-6-4_exp.RDS')


test_tsne <- Rtsne(X = t(new_dfs[[2]][,-1]), perplexity = 40, partial_pca = TRUE, PCA = TRUE)
new_dfs[[3]]$Tsne2 <- test_tsne$Y[,2]
new_dfs[[3]]$Tsne1 <- test_tsne$Y[,1]
ggplot(new_dfs[[3]], aes(x = Tsne1, y = Tsne2, color = Segment)) + geom_point() + theme_bw()
ggplot(tsne$samples, aes(x = tsne$X1, y = tsne$X2, color = Segment)) + geom_point() + theme_bw()

test_uamp <- umap::umap(t(new_dfs[[2]][,-1]), preserve.seed = TRUE)
new_dfs[[3]]$Umap1 <- test_uamp$layout[,1]
new_dfs[[3]]$Umap2 <- test_uamp$layout[,2]
ggplot(new_dfs[[3]], aes(x = Umap1, y = Umap2, color = Segment)) + geom_point() + theme_bw() + labs(title = 'Detrended')
ggplot(umap$samples, aes(x = umap$X1, y = umap$X2, color = Segment)) + geom_point() + theme_bw() + labs(title = 'Q3 Norm only')

ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = Umap1, y = Umap2,
           size = ssGSEA_detrend['CAF_Immunomodulatory',
                                 new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']],
           color = ssGSEA_detrend['CAF_Immunomodulatory',
                                  new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']])) +
  geom_point() + theme_bw() + guides(color = FALSE, size = FALSE)

# iCAF vs myCAF: detrend
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_detrend['CAF_Immunomodulatory',
                                    new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']]),
           y = scale(ssGSEA_detrend['CAF_Myofibroblastic',
                                    new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']]),
           color = treatment)) +
  geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
  labs(x = 'iCAF', y = 'myCAF', title = 'Detrended') +
  scale_color_jama() +
  theme(aspect.ratio = 1)

# iCAF vs myCAF: Q3
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_results['CAF_Immunomodulatory',
                                    new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']]),
           y = scale(ssGSEA_results['CAF_Myofibroblastic',
                                    new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']]),
           color = treatment)) +
  geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
  labs(x = 'iCAF', y = 'myCAF', title = 'Q3 Norm') +
  scale_color_jama() +
  theme(aspect.ratio = 1)


table(ss_match$Complete)
# iCAF vs Mesenchymal
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_detrend['CAF_Immunomodulatory',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = scale(ssGSEA_detrend['Malignant_Mesenchymal',
                                    subset(ss_match, Complete)$Epi_ROIID]),
           color = treatment)) +
  geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
  labs(x = 'iCAF, CAF', y = 'Mesenchymal, Epi', title = 'Detrended') +
  scale_color_jama() +
  theme(aspect.ratio = 1)
cor(scale(ssGSEA_detrend['CAF_Immunomodulatory',
                         subset(ss_match, Complete)$CAF_ROIID]),
    scale(ssGSEA_detrend['Malignant_Mesenchymal',
                         subset(ss_match, Complete)$Epi_ROIID]))

ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_results['CAF_Immunomodulatory',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = scale(ssGSEA_results['Malignant_Mesenchymal',
                                    subset(ss_match, Complete)$Epi_ROIID]),
           color = treatment)) +
  geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
  labs(x = 'iCAF, CAF', y = 'Mesenchymal, Epi', title = 'Q3 Norm') +
  scale_color_jama() +
  theme(aspect.ratio = 1)
cor(scale(ssGSEA_results['CAF_Immunomodulatory',
                         subset(ss_match, Complete)$CAF_ROIID]),
    scale(ssGSEA_results['Malignant_Mesenchymal',
                         subset(ss_match, Complete)$Epi_ROIID]))

# myCAF vs Mesenchymal: detrend
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_detrend['CAF_Myofibroblastic',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = scale(ssGSEA_detrend['Malignant_Mesenchymal',
                                    subset(ss_match, Complete)$Epi_ROIID]),
           color = treatment)) +
  geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
  labs(x = 'myCAF, CAF', y = 'Mesenchymal, Epi', title = 'Detrended') +
  scale_color_jama() +
  theme(aspect.ratio = 1)
cor(scale(ssGSEA_detrend['CAF_Myofibroblastic',
                         subset(ss_match, Complete)$CAF_ROIID]),
    scale(ssGSEA_detrend['Malignant_Mesenchymal',
                         subset(ss_match, Complete)$Epi_ROIID]))

ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_results['CAF_Myofibroblastic',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = scale(ssGSEA_results['Malignant_Mesenchymal',
                                    subset(ss_match, Complete)$Epi_ROIID]),
           color = treatment)) +
  geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
  labs(x = 'myCAF, CAF', y = 'Mesenchymal, Epi', title = 'Q3 Norm') +
  scale_color_jama() +
  theme(aspect.ratio = 1)
cor(scale(ssGSEA_results['CAF_Myofibroblastic',
                         subset(ss_match, Complete)$CAF_ROIID]),
    scale(ssGSEA_results['Malignant_Mesenchymal',
                         subset(ss_match, Complete)$Epi_ROIID]))

# CD45 vs iCAF
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_detrend['CAF_Immunomodulatory',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = unlist(new_dfs[[2]]['CD68',
                                   subset(ss_match, Complete)$CAF_ROIID]),
           color = Treatment)) + geom_point() + theme_bw()


ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(fill = scale(ssGSEA_detrend['CAF_Immunomodulatory',
                                       subset(ss_match, Complete)$CAF_ROIID]) > 0,
           y = unlist(new_dfs[[2]]['PTPRC',
                                   subset(ss_match, Complete)$CAF_ROIID]),
           x = Treatment)) +
  geom_boxplot() + theme_bw() +
  labs(x = 'iCaf', y = 'PTPRC (CD45) in CAF ROIs', title = 'Detrended', fill = 'High iCAF') +
  ylim(c(0,7))

ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(fill = scale(ssGSEA_results['CAF_Immunomodulatory',
                                       subset(new_dfs[[3]], Segment == 'CAF')$Sample_ID]) > 0,
           y = unlist(log2(dfs[[2]]['PTPRC',
                                    subset(new_dfs[[3]], Segment == 'CAF')$Sample_ID])),
           x = Treatment)) + geom_boxplot() + theme_bw() +
  labs(x = 'iCaf', y = 'PTPRC (CD45) in CAF ROIs',  title = 'Q3 Norm', fill = 'High iCAF')+
  ylim(c(0,7))