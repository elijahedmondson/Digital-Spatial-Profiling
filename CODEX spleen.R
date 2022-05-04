#https://jef.works/genomic-data-visualization/blog/2022/02/23/nli47/

library(Rtsne)
library(ggplot2)
library(scattermore)
library(gplots)
library(tidyverse)
library(kableExtra)
library(broom)
library(highcharter)
library(magrittr)
library(gridExtra)
library(cowplot)

data<- read.csv("codex_spleen_subset.csv.gz")
pos<-data[,2:3]
#get an overview of what the biopsy looks like
plot(pos,pch=".",cex=1)
area <- data[,4]
pexp <- data[, 5:ncol(data)]
colnames(pexp)
#normalization
numproteins <- rowSums(pexp)
normpexp <- pexp/numproteins*1e6
mat <- log10(normpexp+1)
#do PCA and tSNE based on pcs
pcs <- prcomp(mat)
plot(pcs$sdev[1:30], type="l")
set.seed(0)
emb <- Rtsne(pcs$x[,1:20], dims=2, perplexity=30)$Y

markers<-c("CD4","CD8","CD31","CollagenIV","Podoplanin")
pexp.plot <- lapply(markers, function(g) {
  df1 <- data.frame(x = emb[,1],
                    y = emb[,2],
                    col = mat[,g]) 
  p1 <- ggplot(data = df1,
               mapping = aes(x = x, y = y)) +
    geom_scattermore(mapping = aes(col = col), 
                     pointsize=2) + 
    scale_color_viridis_c(option = "magma" ) 
  plot1 <- p1 + labs(x = "tSNE X" , y = "tSNE Y", title= g ) +
    theme_classic()
  return(plot1)
})


set.seed(0)
com <- kmeans(pcs$x[,1:20], centers = 10)
df1 <- data.frame(x = emb[,1],
                  y = emb[,2],
                  col = as.factor(clus$cluster)) 
ggplot(data = df1,
       mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), 
                   pointsize=2)

# To determine what tissue it is, I want to look at spatial distributions of these proteins.

pos.plot <- lapply(markers, function(g) {
  df2 <- data.frame(x = pos[,1],
                    y = pos[,2],
                    col = mat[,g]) 
  p1 <- ggplot(data = df2,
               mapping = aes(x = x, y = y)) +
    geom_scattermore(mapping = aes(col = col), 
                     pointsize=2) + 
    scale_color_viridis_c(option = "plasma") 
  plot1 <- p1 + labs(x = "tSNE X" , y = "tSNE Y", title= g ) +
    theme_classic()
  return(plot1)
})



row_1<- plot_grid( pexp.plot[[1]],pos.plot[[1]],ncol=2, nrow=1)
row_2<- plot_grid( pexp.plot[[2]],pos.plot[[2]],ncol=2, nrow=1)
row_3<- plot_grid( pexp.plot[[3]],pos.plot[[3]],ncol=2, nrow=1)
row_4<- plot_grid( pexp.plot[[4]],pos.plot[[4]],ncol=2, nrow=1)
row_5<- plot_grid( pexp.plot[[5]],pos.plot[[5]],ncol=2, nrow=1)

plot_grid(row_1,row_2,row_3,row_4,row_5,nrow=5)

ggsave("Stella_Li.png",width=15,height =20 )