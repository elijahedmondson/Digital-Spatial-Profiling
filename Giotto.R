# library(devtools) # If not installed: install.packages('devtools')
# library(remotes)  #If not installed: install.packages('remotes')
# remotes::install_github("drieslab/Giotto@master")
# 
# remotes::install_github("drieslab/Giotto@cless")
# 
# 
# library(devtools)
# install_github("husson/FactoMineR")

library(scran)
library(smfishHmrf)
library(trendsceek)
#library(SPARK)
library(multinet)
library(RTriangle)
library(FactoMineR)
library(Giotto)

installGiottoEnvironment()
installGiottoEnvironment(force_environment = TRUE)
installGiottoEnvironment(force_miniconda = TRUE)
#removeGiottoEnvironment()

# instrs = createGiottoInstructions(show_plot = FALSE,
#                                   save_plot = TRUE,
#                                   save_dir = 'C:/Users/edmondsonef/desktop/Giotto_Results/',
#                                   python_path = 'C:/Users/edmondsonef/AppData/Local/Microsoft/WindowsApps/python3.exe')
#                                  #python_path = "C:/Users/edmondsonef/AppData/Local/Programs/Python/Python39/")

qp = read.delim("C:/Users/edmondsonef/desktop/measurements.tsv")
head(qp)
qp_expr = qp[,grepl("(Cell|Nucleus|Cytoplasm|Membrane)..",names(qp)) & !grepl("Autofluorescence|DAPI|",names(qp))]
#Area|Circularity|Perimeter|caliper|Eccentricity",names(qp))]
qp_expr = t(qp_expr)
colnames(qp_expr) = rownames(qp)

qp_spatial_loc = qp[,c("Centroid.X.µm","Centroid.Y.µm")]
qp_spatial_loc$Centroid.Y.µm = - qp_spatial_loc$Centroid.Y.µm
qp_spatial_loc$cell_ID = rownames(qp)
qp_spatial_loc = qp_spatial_loc[,c(3,1,2)]

gobj <- createGiottoObject(raw_exprs = qp_expr,
                           spatial_locs = qp_spatial_loc)#,
#instructions = instrs)

# optionally add QuPath metadata such as marker positivity
qp_metadata = qp[,grepl("Class|phenotype",names(qp))]
qp_metadata$cell_ID = rownames(qp)

gobj<-addCellMetadata(gobj, new_metadata = qp_metadata,
                      by_column = T,
                      column_cell_ID = "cell_ID")


## filter
# gobj <- filterGiotto(gobject = gobj,
#                      expression_threshold = 1,
#                      gene_det_in_min_cells = 10,
#                      min_det_genes_per_cell = 2,
#                      expression_values = c('raw'),
#                      verbose = T)

gobj <- normalizeGiotto(gobject = gobj, verbose = T, #scalefactor = 6000, 
                        log_norm = FALSE,library_size_norm = FALSE,
                        scale_genes = FALSE, scale_cells = TRUE)

## add gene & cell statistics
gobj <- addStatistics(gobject = gobj, expression_values = "normalized")

## adjust expression matrix for technical or known variables
gobj <- adjustGiottoMatrix(gobject = gobj, 
                           expression_values = c('normalized'),
                           batch_columns = NULL, 
                           covariate_columns = NULL,
                           return_gobject = TRUE,
                           update_slot = c('custom'))

## visualize
spatPlot(gobject = gobj, point_size = 0.1, 
         coord_fix_ratio = NULL,point_shape = 'no_border',
         save_param = list(save_name = '2_a_spatPlot'))

spatPlot(gobject = gobj, point_size = 0.2,
         coord_fix_ratio = 1, cell_color = 'sample_Xtile_Ytile',
         legend_symbol_size = 3,legend_text = 5,
         save_param = list(save_name = '2_b_spatPlot'))

# PCA
gobj <- runPCA(gobject = gobj, expression_values = 'normalized', scale_unit = T, method = "factominer")
signPCA(gobj, scale_unit = T, scree_ylim = c(0, 3),
        save_param = list(save_name = '3_a_spatPlot'))

plotPCA(gobject = gobj, point_shape = 'no_border', point_size = 0.2,
        save_param = list(save_name = '3_b_PCA'))


# UMAP
gobj <- runUMAP(gobj)#, n_components = 2, n_threads = 12)
plotUMAP(gobject = gobj, point_shape = 'no_border', point_size = 0.2,
         save_param = list(save_name = '3_c_UMAP'))


## sNN network (default)
gobj <- createNearestNetwork(gobject = gobj, k = 20, dimensions_to_use = 11:15)

## 0.1 resolution

gobj <- doLeidenCluster(gobject = gobj, resolution = 1, n_iterations = 100, name = 'leiden')#, python_path = "C:/Users/edmondsonef/Anaconda3/pkgs/python-3.9.12-h6244533_0/")
### CRASH
### CRASH
### CRASH
### CRASH
### CRASH

gobj <- doKmeans(gobject = gobj, dim_reduction_to_use = "umap",dimensions_to_use = 11:15)

qupath_metadata = pDataDT(qupath_test)
leiden_colors = Giotto:::getDistinctColors(length(unique(qupath_metadata$leiden)))
names(leiden_colors) = unique(qupath_metadata$leiden)

plotUMAP(gobject = qupath_test, 
         cell_color = 'leiden', point_shape = 'no_border', point_size = 0.2, cell_color_code = leiden_colors,
         save_param = list(save_name = '4_a_UMAP'))





spatPlot(gobject = gobj, cell_color = 'leiden', point_shape = 'no_border', point_size = 0.2, 
         cell_color_code = leiden_colors, coord_fix_ratio = 1,label_size =2,
         legend_text = 5,legend_symbol_size = 2,
         save_param = list(save_name = '4_b_spatplot'))


spatDimPlot2D(gobject = gobj, cell_color = 'leiden', spat_point_shape = 'no_border', 
              spat_point_size = 0.2, dim_point_shape = 'no_border', dim_point_size = 0.2, 
              cell_color_code = leiden_colors,plot_alignment = c("horizontal"),
              save_param = list(save_name = '5_a_spatdimplot'))


# resolution 0.5
cluster_column = 'leiden'
markers_scran = findMarkers_one_vs_all(gobject=gobj, method="scran",
                                       expression_values="norm", cluster_column=cluster_column, min_genes=3)
markergenes_scran = unique(markers_scran[, head(.SD, 5), by="cluster"][["genes"]])

plotMetaDataHeatmap(gobj, expression_values = "norm", metadata_cols = c(cluster_column), 
                    selected_genes = markergenes_scran,
                    y_text_size = 8, show_values = 'zscores_rescaled',
                    save_param = list(save_name = '6_a_metaheatmap'))



#Once done with the Giotto clustering analysis its results can be saved to a tsv/csv file and fed back to QuPath.


for_qp = qp[,c("Centroid.X.�m","Centroid.Y.�m")]
for_qp$Image = "giotto"
for_qp$leiden = gobj@cell_metadata$cell$rna$leiden
write.table(for_qp,"qp_giotto.tsv",quote = FALSE,sep = "\t",row.names = FALSE)