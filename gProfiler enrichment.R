# Date: Sep 28, 2018
# Author: Shirley Hui, Ruth Isserlin
# Takes supplied list of genes and uses gProfiler to perform enrichment analysis to determine which pathways summarize the genes supplied.
# Input: list of core killing genes (Fig 1G)
# Output: Pathway enrichment output file
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(RCy3)
library(RCurl)

tryCatch(expr = { library("gProfileR")}, 
         error = function(e) { install.packages("gProfileR")}, finally = library("gProfileR"))




# Function to run gprofiler using the gprofiler library
# The function takes the returned gprofiler results and formats it to the generic EM input file
# function returns a data frame in the generic EM file format.
runGprofiler <- function(genes,current_organism = "mmusculus", ###mmusculus for mouse!!
                         significant_only = T, set_size_max = 200, 
                         set_size_min = 3, filter_gs_size_min = 5 , exclude_iea = F){
  
  gprofiler_results <- gprofiler(genes ,
                                 significant=significant_only,ordered_query = F,
                                 exclude_iea=exclude_iea,max_set_size = set_size_max,
                                 min_set_size = set_size_min,
                                 correction_method = "fdr",
                                 organism = current_organism,
                                 src_filter = c("GO:BP","GO:MF","GO:CC","REAC","KEGG","WP","TF","MIRNA","HPA","CORUM"))
  
  # Filter results
  gprofiler_results <- gprofiler_results[which(gprofiler_results[,'term.size'] >= 3
                                               & gprofiler_results[,'overlap.size'] >= filter_gs_size_min ),]
  
  # gProfileR returns corrected p-values only.  Set p-value to corrected p-value
  if(dim(gprofiler_results)[1] > 0){
    
    gprofiler_results_filename <-"C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/Pathway analysis/Output/gprofiler_results/coreCTLgenes_gprofiler.txt"
    write.table(gprofiler_results,gprofiler_results_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)
    
    em_results <- cbind(gprofiler_results[,
                                          c("term.id","term.name","p.value","p.value")], 1,
                        gprofiler_results[,"intersection"])
    colnames(em_results) <- c("Name","Description", "pvalue","qvalue","phenotype","genes")
    
    return(em_results)
  } else {
    return("no gprofiler results for supplied query")
  }
}

genes = read.delim("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/Pathway analysis/Input/coreCTLgenes.txt", header = F)
genes = as.vector(genes[,1])
gprofiler_results = runGprofiler(genes)

# Write out the g:Profiler results
em_results_filename <-"C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/Pathway analysis/Output/em_file/coreCTLgenes_sim0.7.txt"
write.table(gprofiler_results,em_results_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

em_command = paste('enrichmentmap build analysisType="generic" ', 
                   'pvalue=',"0.001", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.7",
                   'coeffecients=',"JACCARD+OVERLAP",
                   'enrichmentsDataset1=',em_results_filename ,
                   sep=" ")

# Enrichment map command will return the suid of newly created network.
em_network_suid <- commandsRun(em_command)
renameNetwork("Cluster1_enrichmentmap", network=as.numeric(em_network_suid))












# Date: Sep 28, 2018
# Author: Shirley Hui
# Takes pathway themes manually identified via a Cytoscape network created using the the core killing gprofiler results (see CoreKillingGProfiler.R).  Pathways were grouped together to form themes if they contained 30% or more similar genes.  Plot the pathway themes into bar plot.
# Input: core killing gprofiler results, core killing enrichment themes
# Output: Bar plot of core killing enriched themes  
gprofilerResults<- read.delim("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/Pathway analysis/Output/gprofiler_results/coreCTLgenes_gprofiler.txt")
emThemes <- read.delim("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/Pathway analysis/Output/enr_file/coreCTLgenes_sim0.7_enrTheme.txt",header=FALSE)
themes <- as.character(unique(emThemes[,1]))
results <- c()
for (ixx in 1:length(themes)) {
  ix <- which(emThemes[,1]==themes[ixx])
  ixs <- c()
  for (i in 1:length(ix)) {
    goid <- as.character(emThemes[ix[i],2])
    ixi <- which(gprofilerResults$term.id==goid)
    ixs <- c(ixs,ixi)
  }
  mean_overlap <- mean(gprofilerResults[ixs,]$overlap.size/gprofilerResults[ixs,]$term.size)*100
  mean_overlap.size <- mean(gprofilerResults[ixs,]$overlap.size)
  mean_term.size <- mean(gprofilerResults[ixs,]$term.size)
  min_pvalue <- -log(min(gprofilerResults[ixs,]$p.value))
  results <- rbind(results,c(mean_overlap,min_pvalue,mean_overlap.size,mean_term.size))
}
rownames(results) <- themes
colnames(results) <- c("overlap","p.value","overlap.size","term.size")

library(ggplot2)
library(RColorBrewer)

#cbPalette <- c("#FED976", "#FD8D3C", "#FC4E2A", "#E31A1C", "#aa0022")
cbPalette <- c("#ededed", "#cccccc", "#969696", "#636363", "#252525") #http://colorbrewer2.org/#type=sequential&scheme=Greys&n=5
cols = cbPalette #<- brewer.pal(6, "YlOrRd")
df = data.frame(results)
df$ratio <- paste(round(df$overlap.size,1), round(df$term.size, 1), sep = "/") #this line adds the overlap/term size ratio, rounds up the term size to xx position after comma
g = ggplot(df, aes(x = reorder (rownames(df),p.value), y = overlap)) +
  ylab("Mean Percentage Overlap") +
  theme(plot.title = element_text(hjust = -0.9)) +
  geom_col(aes(fill = p.value)) +
  geom_text(aes(label = df$ratio, hjust=0))+
  scale_fill_gradientn("-log p", colours = cols, limits=c(min(df$p.value), max(df$p.value))) +
  scale_y_continuous(position = "right") +
  theme(panel.background = element_rect(fill = "white"), axis.line.x = element_line(color="black"), axis.line.y = element_line(color="white"), axis.title.y = element_blank(),axis.ticks.y = element_blank()) +
  coord_flip()

# Adjust aspect_ratio, height and width to output figure to pdf in the desired dimensions
aspect_ratio <- 1.75
ggsave(file="C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/Pathway analysis/Output/coreCTLgenes_sim0.7.pdf",g, height = 4 , width = 5 * aspect_ratio)





