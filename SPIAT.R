if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SPIAT")

library(SPIAT)

######
###### Cellular neighborhood
######
# 
# The aggregation of cells can result in 'cellular neighbourhoods'. A neighbourhood is defined as a group of cells that cluster together. These can be homotypic, containing cells of a single class (e.g. immune cells), or heterotypic (e.g. a mixture of tumour and immune cells).
# 
# Function identify_neighborhoods() identifies cellular neighbourhoods. Users can select a subset of cell types of interest if desired. SPIAT includes three algorithms for the detection of neighbourhoods.
# 
# Hierarchical Clustering algorithm: Euclidean distances between cells are calculated, and pairs of cells with a distance less than a specified radius are considered to be 'interacting', with the rest being 'non-interacting'. Hierarchical clustering is then used to separate the clusters. Larger radii will result in the merging of individual clusters.
# dbscan
# phenograph
# For Hierarchical Clustering algorithm and dbscan, users need to specify a radius that defines the distance for an interaction. We suggest users to test different radii and select the one that generates intuitive clusters upon visualisation. Cells not assigned to clusters are assigned as Cluster_NA in the output table. The argument min_neighborhood_size specifies the threshold of a neighborhood size to be considered as a neighborhood. Smaller neighbourhoods will be outputted, but will not be assigned a number.
# 
# Rphenograph uses the number of nearest neighbours to detect clusters. This number should be specified by min_neighborhood_size argument. We also encourage users to test different values.
# 
# For this part of the tutorial, we will use the image image_no_markers simulated with the spaSim package. This image contains "Tumour", "Immune", "Immune1" and "Immune2" cells without marker intensities.




######
###### Characterising the distribution of the cells of interest in identified tissue regions
######






