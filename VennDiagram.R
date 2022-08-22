

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




DF <- ego@result
DF.SO <- DF %>% filter(Description == "axonogenesis")

#Axonogenesis
top_g <- c("Cckar","Slit3","Etv1","Sema7a","Cxcr4","Kif5c","Ptprm",
"Ptprs","Bsg","Smad4","Bcl2","Trak1","Ptch1","
Islr2","Taok2","Cdkl3","Actb","Ednra","Pip5k1c",
"Sptbn4","Lama5","Ablim1","Rtn4","Wnt7b","Spg11","
Golga4","Dock7","Ephb2","Cacna1a","Ptpn11","B4gat1",
"Smo","B4galt6","Rab3a","Ntrk3","Neo1","Lrp1","
Atp5g1","Kif5b","Brsk2","Erbb2","Map1a","Flrt2",
"Tsku","Map1s","Chrnb2","Fstl4","Lrp4","Dag1","Sin3a","
Mapk8ip3","Dclk1","Adnp","Celsr3","Rpl4","Tubb2b",
"Efna5","Plxnb2","Ngf","Ache","Vim","Flrt3","Fgfr2","
Ephb4","Flot1","Sema4c","Gsk3b","Sema3d","Aatk","Cdh4",
"Tubb3","Agrn","Evl","Brsk1","Notch3","Fzd3","
Hsp90aa1","Nrn1","Bcl11a","Sema4g","Lama3","Epha8",
"Ntn5","Amigo1","Apbb1","Mgll","Ret","Atp8a2","
Alcam","Unc5a","Grin1","Cntn6","Wnt7a","Pou4f3","Shh")



#Synapse organization
top_g <- c("Abl2", "Ache", "Actb", "Actr3", "Adam10", "Adgrb2", "Adgre5", 
          "Adgrf1", "Adgrl3", "Adnp", "Arf1","Arf4", "Arf6", "Arhgap44", 
          "Baiap2", "Bhlhb9", "C1qa", "C1ql1", "C3", "Cacna1a", "Cacna1s", 
          "Cacnb1","Cacnb4", "Camk1", "Camk2b", "Caprin1", "Cel", "Cfl1", 
          "Chchd10", "Chd4", "Chrnb1", "Clstn1", "Cnksr2", "Cntnap4", 
          "Col4a1", "Col4a5", 'Ctnnb1', "Cttn", "Cttnbp2", "Cyfip1", 
          "Dact1", "Dag1", "Dbn1", "Dctn1", "Dlgap3", "Dock10", "Drd1", 
          "Efna1", "Efnb2", "Eif4g1", "Epha4", "Ephb2", "Erbb4", "F2r", 
          "Farp1", "Fgfr2", "Flna", "Flrt3", "Fzd5", "Gabrb3", "Ghrl", 
          "Gnpat", "Gphn", "Grm5", "Hnrnpk", "Hspa8", "Igsf9", "Insr", 
          "Itga3", "Itpka", "Klk8", "Lamb2", "Lgi2", "Lrfn2", "Lrfn5", 
          "Lrrc4c", "Lrrk2", "Lrrtm2", "Lzts3", "Magi2", "Marcks", "Mdga1", 
          "Mdga2", "Mef2c", "Mfn2", "Myh10", "Ndrg1", "Nedd4", "Nfasc",
          "Nfatc4", "Nfia", "Nlgn1", "Nlgn3", "Nrcam", "Nrg1", "Nrg2", 
          "Nrp1", "Nrxn1", "Nrxn3","Ntn1", "Ntrk2", "Numb", "Obsl1", "Ophn1", 
          "Pak3", "Palm", "Pcdh17", "Pcdhgc4", "Pclo", "Pdlim5", "Pdzrn3", 
          "Pfn1", "Pfn2", "Picalm", "Pik3r1", "Pin1", "Pmp22", "Ppfia2", 
          "Ppfia4", "Prkca", "Prrt1", "Psen1", "Ptn", "Ptprf", "Ptprt", 
          "Rab17", "Rab29", "Rab39b", "Rapsn", "Rhoa", "Rims4", "Rock2", 
          "Sdf4", "Sdk1", "Septin7", "Setd5", "Sez6", "Sez6l", "Shank1", 
          "Shank2", "Shank3", "Sipa1l1", "Six4", "Slitrk6", "Snta1", 
          "Sorbs1", "Sparc", "Srcin1", "Srgn", "Ssh1", "St8sia2", 
          "Syngap1", "Tanc2", "Taok2", "Tnc", "Tubb5", "Vcp", "Vps35", 
          "Wnt5a", "Wnt7b", "Ywhaz", "Zmynd8", "Lgmn", "Tuba1a")

top_g <- unique(top_g)
top_g

#reverse log fold change to fit with label
#gene$Estimate1 <- gene$Estimate*(-1)


ggplot(gene,                                                             
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
