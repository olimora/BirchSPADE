## This file hase similar function as driver.R in SPADE.
## Here is declared function that should be called to use this package for SPADE analysis


BirchSPADE.run.analysis <- function(input_file_full             # full path to the .fcs file
                                   ,outputs_dir                # full path where new new folder with output files will be created
                                   ,markers                    # vector of markers to use in clustering analysis
                                   ,normalization = "minmax"    # c("none", "minmax", "meanstd")
                                   ,remove_outliers = FALSE
                                   ,final_cluster_count = 500
                                   ,kmeans_upsampling_iterations = 1
                                   ,plot_trees = TRUE
                                   ,subcluster_limit = 0
                                   ,use_density = FALSE) { #TODO: add BirchTree parameters so they can be set

  # load packages
  suppressWarnings(library(flowCore))
  suppressWarnings(library(BirchTree))
  suppressWarnings(library(Rclusterpp))
  suppressWarnings(library(igraph))

  whole_analysis_start_time <- Sys.time()

  # set input file names, output directory ...
  input_file_name = basename(input_file_full)
  outputs_dir <- paste0(outputs_dir,"/",input_file_name)
  dir.create(outputs_dir, showWarnings = FALSE) #cretae directory if not existing

  # set needed variables
  comp = TRUE
  transforms = flowCore::arcsinhTransform(a=0, b=0.2)
  markers_cout = length(markers)

  ## 1 # read input fcs file, load data, use arcsinh transform
  message("Loading fcs data, transforming, normalizing ... ")
  cells_data = BirchSPADE.load_input_data(input_file_full, transforms, markers, normalization)


  ## 2 # use BirchTree to reduce data (as downsampling in SPADE)
  BF_B = 10
  BF_L = 15
  if (is.null(subcluster_limit) || subcluster_limit < 1) {
    subcluster_limit = round(nrow(cells_data)/10)
  }

  message("BirchTree data reduction ... ")
  BirchTree_start_time <- Sys.time()
  birch_out = BirchTree::buildTree(cells_data, BF_B, BF_L, 0, cluster_size_metric = "radius",
                                   subcluster_limit = subcluster_limit, rebuild_size_factor = 2,
                                   remove_outliers = remove_outliers)
  subclusters = as.data.frame(birch_out$subclusters)[,1:markers_cout]
  colnames(subclusters) = markers
  if (use_density) {
    dens_c = ncol(birch_out$subclusters)
    subclusters$density = birch_out$subclusters[,dens_c]
  }
  BirchTree_end_time <- Sys.time()
  message(paste0("BirchTree data reduction took time (seconds): ",
               round(difftime(BirchTree_end_time, BirchTree_start_time, units='secs'), digits = 2)))

  # remove outliers if wanted
  outliers = NULL
  rows_without_outliers = NULL
  if (remove_outliers && nrow(birch_out$outliers) > 0) {
    message("Outlier removal from fcs data matrix ... ")
    remove_outliers_start_time <- Sys.time()
    outliers = as.data.frame(birch_out$outliers)[,1:markers_cout]
    colnames(outliers) = markers
	if (use_density) {
      dens_c = ncol(birch_out$outliers)
      outliers$density = birch_out$outliers[,dens_c]
    }
    removal.result = BirchSPADE.remove_outliers(cells_data, outliers)
    cells_data = removal.result$data
    rows_without_outliers = removal.result$rows
    remove_outliers_end_time <- Sys.time()
    message(paste0("Outlier removal took time (seconds): ",
                 round(difftime(remove_outliers_end_time, remove_outliers_start_time, units='secs'), digits = 2)))
  }

  if (use_density) { # normalize density
    # normalize loaded data
    if (normalization == "none") {
      # print("No normalization.")
    } else if (normalization == "minmax") {
      # print("Minmax normalization.")
      library(caret)
      pp = preProcess(as.data.frame(subclusters$density), method = "range")
      density_col = unname(as.vector(predict(pp, as.data.frame(subclusters$density))))
      subclusters = cbind(subclusters[,1:markers_cout], density = density_col)
      rm(pp, density_col)
    } else if (normalization == "meanstd") {
      # print("Meanstd normalization.")
      subclusters$density = scale(subclusters$density)
    }
  }

  ## 3 # cluster subclusters (reduced data) hierarchicaly
  message("Hierarchical clustering of subclusters ... ")
  hclust_start_time <- Sys.time()
  # methods <- c("ward", "average", "single", "complete")
  hclust.result <- Rclusterpp.hclust(subclusters, method = "ward", distance = "euclidean")
  subclusters$hier_cluster = cutree(hclust.result, k = final_cluster_count)
  hclust_end_time <- Sys.time()
  message(paste0("Hierarchical clustering took time (seconds): ",
               round(difftime(hclust_end_time, hclust_start_time, units='secs'), digits = 2)))


  ## 4 # upsampling using K-means
  # need to find for every point from cells_data closest sub_cluster from birch/ hier_cluster from hclust with k-means, 1 iteration
  message("Upsampling fcs to clusters using k-means ... ")
  upsampling_start_time <- Sys.time()
  # get hier_clusters centroids
  hier_cluster_centroids <- aggregate(subclusters, by=list(subclusters$hier_cluster),FUN=mean)[,-1] # first column is the group
  # suppress warning that it did not konverge
  suppressWarnings(kmeans.result <- kmeans(x = cells_data,
         centers = hier_cluster_centroids[,1:markers_cout], iter.max = kmeans_upsampling_iterations))
  cells_data <- cbind(cells_data, "cluster" = kmeans.result$cluster)
  # every cell has assigned final cluster
  upsampling_end_time <- Sys.time()
  message(paste0("Upsampling took time (seconds): ",
               round(difftime(upsampling_end_time, upsampling_start_time, units='secs'), digits = 2)))


  ## 5 # add asigned clusters to original not normalized data, bud without outliers if removal was used
  message("Saving orginal fcs data with assigned clusters so medians and others can be computed ... ")
  saving_start_time <- Sys.time()
  saving.returns = BirchSPADE.save_original_data_with_clusters(input_file_full, outputs_dir,
                                                              kmeans.result$cluster, comp, transforms,
                                                              remove_outliers, rows_without_outliers)
  saving_end_time <- Sys.time()
  message(paste0("Saving fcs with clusters took time (seconds): ",
               round(difftime(saving_end_time, saving_start_time, units='secs'), digits = 2)))


  ## 6 # prepare data for visualization: compute graph, layout, medians, etc ...
  message("Computing visualisations = graph, layout, medians and others ... ")
  compute_visualisation_start_time <- Sys.time()
  visualisation = BirchSPADE.compute_visualisation(outputs_dir, hclust.result$merge,
                                  hier_cluster_centroids[1:markers_cout], input_file_full,
                                  saving.returns$file, transforms, markers, comp)
  compute_visualisation_end_time <- Sys.time()
  message(paste0("Computing visualisations took time (seconds): ",
               round(difftime(compute_visualisation_end_time, compute_visualisation_start_time, units='secs'), digits = 2)))


  ## 7 # compute statistic tables
  message("Computing statistics tables ... ")
  statistics_start_time <- Sys.time()
  SPADElike_statistics(outputs_dir, files)
  statistics_end_time <- Sys.time()
  message(paste0("Computing statistics took time (seconds): ",
               round(difftime(statistics_end_time, statistics_start_time, units='secs'), digits = 2)))

  ## 8 # plot trees to pdf
  if (plot_trees) {
    message("Plotting trees to pdf ... ")
    plot_trees_start_time <- Sys.time()
    SPADElike_plotTrees(visualisation$graph, outputs_dir, saving.returns$fcs_params,
                        file_pattern = "*fcs*Rsave",
                        layout = as.matrix(visualisation$layout_table),
                        out_dir = file.path(outputs_dir,"pdf"),
                        size_scale_factor=1.2,
                        cluster_markers = markers)
    plot_trees_end_time <- Sys.time()
    message(paste0("Plotting trees took time (seconds): ",
                 round(difftime(plot_trees_end_time, plot_trees_start_time, units='secs'), digits = 2)))
  }


  whole_analysis_end_time <- Sys.time()
  message(paste0("Whole analysis took time (seconds): ",
               round(difftime(whole_analysis_end_time, whole_analysis_start_time, units='secs'), digits = 2)))

}
