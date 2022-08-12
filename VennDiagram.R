

library(VennDiagram)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)
library(tidyverse)
library(plyr)
library(GGally)
library(ggplot2)
library(tidyverse)
library(gapminder)
library(dplyr)

PDAC1 <- PDAC %>% distinct(SYMBOL, .keep_all=T)

#keep only significatn genes with FDR <0.05
lo1 <- lo1 %>% filter(FDR < 0.05)
hi1 <- hi1 %>% filter(FDR < 0.05)
PDAC1 <- PDAC1 %>% filter(FDR < 0.05)

#keep only genes enriched in mets Estimate < -1
lo1 <- lo1 %>% filter(Estimate < -1)
hi1 <- hi1 %>% filter(Estimate < -1)
PDAC1 <- PDAC1 %>% filter(Estimate < -1)

gene_list <- list(PanIN.lo = lo1$SYMBOL, 
                  PanIN.hi = hi1$SYMBOL,
                  PDAC = PDAC1$SYMBOL)
VennDiagram <- venn.diagram(x = gene_list, 
                            fill = c("blue", "red", "green"),
                            cat.col = c("blue", "red", "green"),
                            cex = 2,lty = "blank",
                            cat.cex = 2,
                            filename = NULL)
cowplot::plot_grid(VennDiagram)



list <- get.venn.partitions(gene_list) %>% dplyr::as_tibble()
write.csv(list$..values..$`3`, "C:/Users/edmondsonef/Desktop/MetGenesEnrichedOverlap.csv")
