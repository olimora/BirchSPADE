# # this fuction serves for generating clustering analysisi only, does not save any output files
#
# BirchSPADE.clustering.only <- function(input_file_full             # full path to the .fcs file
#                                     ,markers                    # vector of markers to use in clustering analysis
#                                     ,normalization = "minmax"    # c("none", "minmax", "meanstd")
#                                     ,remove_outliers = FALSE
#                                     ,final_cluster_count = 500
#                                     ,kmeans_upsampling_iterations = 1
#                                     ,subcluster_limit = 0
#                                     ,use_density = FALSE
#                                     ,hclust_method = "ward") { #TODO: add BirchTree parameters so they can be set
#
#   # load packages
#   suppressWarnings(library(flowCore))
#   suppressWarnings(library(BirchTree))
#   suppressWarnings(library(Rclusterpp))
#   suppressWarnings(library(igraph))
#
#   whole_analysis_start_time <- Sys.time()
#
#   # set input file names, output directory ...
#   input_file_name = basename(input_file_full)
#
#   # set needed variables
#   comp = TRUE
#   transforms = flowCore::arcsinhTransform(a=0, b=0.2)
#   markers_cout = length(markers)
#
#   ## 1 # read input fcs file, load data, use arcsinh transform
#   message("Loading fcs data, transforming, normalizing ... ")
#   cells_data = BirchSPADE.load_input_data(input_file_full, transforms, markers, normalization)
#
#
#   ## 2 # use BirchTree to reduce data (as downsampling in SPADE)
#   BF_B = 10
#   BF_L = 15
#   if (is.null(subcluster_limit) || subcluster_limit < 1) {
#     subcluster_limit = round(nrow(cells_data)/10)
#   }
#
#   message("BirchTree data reduction ... ")
#   BirchTree_start_time <- Sys.time()
#   birch_out = BirchTree::buildTree(cells_data, BF_B, BF_L, 0, cluster_size_metric = "radius",
#                                    subcluster_limit = subcluster_limit, rebuild_size_factor = 2,
#                                    remove_outliers = remove_outliers)
#   subclusters = as.data.frame(birch_out$subclusters)[,1:markers_cout]
#   colnames(subclusters) = markers
#   if (use_density) {
#     dens_c = ncol(birch_out$subclusters)
#     subclusters$density = birch_out$subclusters[,dens_c]
#   }
#   BirchTree_end_time <- Sys.time()
#   message(paste0("BirchTree data reduction took time (seconds): ",
#                  round(difftime(BirchTree_end_time, BirchTree_start_time, units='secs'), digits = 2)))
#
#   # remove outliers if wanted
#   outliers = NULL
#   rows_without_outliers = NULL
#   if (remove_outliers && nrow(birch_out$outliers) > 0) {
#     message("Outlier removal from fcs data matrix ... ")
#     remove_outliers_start_time <- Sys.time()
#     outliers = as.data.frame(birch_out$outliers)[,1:markers_cout]
#     colnames(outliers) = markers
#     if (use_density) {
#       dens_c = ncol(birch_out$outliers)
#       outliers$density = birch_out$outliers[,dens_c]
#     }
#     removal.result = BirchSPADE.remove_outliers(cells_data, outliers)
#     cells_data = removal.result$data
#     rows_without_outliers = removal.result$rows
#     remove_outliers_end_time <- Sys.time()
#     message(paste0("Outlier removal took time (seconds): ",
#                    round(difftime(remove_outliers_end_time, remove_outliers_start_time, units='secs'), digits = 2)))
#   }
#
#   if (use_density) { # normalize density
#     # normalize loaded data
#     if (normalization == "none") {
#       # print("No normalization.")
#     } else if (normalization == "minmax") {
#       # print("Minmax normalization.")
#       library(caret)
#       pp = preProcess(as.data.frame(subclusters$density), method = "range")
#       density_col = unname(as.vector(predict(pp, as.data.frame(subclusters$density))))
#       subclusters = cbind(subclusters[,1:markers_cout], density = density_col)
#       rm(pp, density_col)
#     } else if (normalization == "meanstd") {
#       # print("Meanstd normalization.")
#       subclusters$density = scale(subclusters$density)
#     }
#   }
#
#   ## 3 # cluster subclusters (reduced data) hierarchicaly
#   message("Hierarchical clustering of subclusters ... ")
#   hclust_start_time <- Sys.time()
#   # methods <- c("ward", "average", "single", "complete")
#   hclust.result <- Rclusterpp.hclust(subclusters, method = hclust_method, distance = "euclidean")
#   subclusters$hier_cluster = cutree(hclust.result, k = final_cluster_count)
#   # get hier_clusters centroids
#   hier_cluster_centroids <- aggregate(subclusters, by=list(subclusters$hier_cluster),FUN=mean)[,-1] # first column is the group
#   hclust_end_time <- Sys.time()
#   message(paste0("Hierarchical clustering took time (seconds): ",
#                  round(difftime(hclust_end_time, hclust_start_time, units='secs'), digits = 2)))
#
#
#   ## 4 # upsampling using K-means
#   # need to find for every point from cells_data closest sub_cluster from birch with k-means, 1 iteration
#   message("Upsampling fcs to clusters using k-means ... ")
#   upsampling_start_time <- Sys.time()
#   # suppress warning that it did not konverge
#   suppressWarnings(kmeans.result <- kmeans(x = cells_data,
#                                            centers = subclusters[,1:markers_cout],
#                                            iter.max = kmeans_upsampling_iterations,
#                                            algorithm = "Lloyd"))
#   #assign clusters to cells by subcluster assigment to hier clusters
#   cluster_col = subclusters$hier_cluster[kmeans.result$cluster]
#   cells_data = cbind(cells_data, "cluster" = cluster_col)
#   # every cell has assigned final cluster
#   upsampling_end_time <- Sys.time()
#   message(paste0("Upsampling took time (seconds): ",
#                  round(difftime(upsampling_end_time, upsampling_start_time, units='secs'), digits = 2)))
#
#
#   whole_analysis_end_time <- Sys.time()
#   message(paste0("Whole analysis took time (seconds): ",
#                  round(difftime(whole_analysis_end_time, whole_analysis_start_time, units='secs'), digits = 2)))
#
#   return (as.data.frame(cells_data))
# }
