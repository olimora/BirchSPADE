BirchSPADE.compute_visualisation <- function(outputs_dir
                                             ,merge_order
                                             ,hier_cluster_centroids
                                             ,input_file_full
                                             ,data_with_clusters_file
                                             ,transforms
                                             ,cluster_cols
                                             ,comp
                                             ) {
  #from SPADE cluster.R: # Write the DEFAULT merge order
  merge_order_file = paste0(outputs_dir,"/","merge_order.txt")
  write.table(merge_order, file=merge_order_file, sep="\t",quote=F,row.names=F,col.names=F)

  # create / compute graph
  graph_file <- paste0(outputs_dir,"/","graph")
  adjacency  <- dist(hier_cluster_centroids, method="manhattan")
  full_graph <- graph.adjacency(as.matrix(adjacency),mode="undirected",weighted=TRUE)
  graph  <- minimum.spanning.tree(full_graph)
  write.graph(graph, graph_file, format="gml")

  #from driver.R
  # Compute the layout once for the MST, write out for use by other functions
  layout= igraph:::layout.fruchterman.reingold #SPADE.layout.arch #igraph:::layout.fruchterman.reingold #igraph:::layout.kamada.kawai
  layout_table <- layout(graph)
  suppressWarnings(
    if (deparse(substitute(layout)) != "SPADE.layout.arch") { # The igraph internal layouts are much more compact than arch.layout
      layout_table = layout_table * 50
    })
  write.table(layout_table, paste0(outputs_dir, "/layout.table"), row.names = FALSE, col.names = FALSE)

  # compute medians and stuff
  input_file = basename(input_file_full)
  panels <- list(list(panel_files=c(input_file),
                      median_cols=NULL,
                      reference_files=c(input_file),
                      fold_cols=c()))

  # ATTENTION! this code is already prepared for analysis of multiple input files/panles,
  # therefore folowing are vectors and then there is loop through them
  files = c(input_file_full)
  sampled_files = c(data_with_clusters_file)

  # Track all attributes to compute global limits
  attr_values <- list()

  if (is.null(panels)) {  # Initialize panels if NULL
    panels <- list( list(panel_files=basename(files), median_cols=NULL) )
  }

  for (p in panels) {

    reference_medians <- NULL
    if (!is.null(p$reference_files)) {
      # Note assuming the ordering of files and sampled_files are identical...
      reference_files   <- sapply(as.vector(p$reference_files), function(f) { sampled_files[match(f,basename(files))[1]] })
      reference_medians <- SPADElike_markerMedians(reference_files, vcount(graph), cols=p$fold_cols, transforms=transforms, cluster_cols=cluster_cols, comp=comp)
    }

    for (f in as.vector(p$panel_files)) {
      # Note assuming the ordering of files and sampled_files are identical...
      f <- sampled_files[match(f, basename(files))[1]]

      # Compute the median marker intensities in each node, including the overall cell frequency per node
      # message("Computing medians for file: ",f)
      anno <- SPADElike_markerMedians(f, vcount(graph), cols=p$median_cols, transforms=transforms, cluster_cols=cluster_cols, comp=comp)
      if (!is.null(reference_medians)) {	# If a reference file is specified
        # Compute the fold change compared to reference medians
        # message("Computing fold change for file: ", f)
        fold_anno <- SPADElike_markerMedians(f, vcount(graph), cols=p$fold_cols, transforms=transforms, cluster_cols=cluster_cols, comp=comp)
        fold <- fold_anno$medians - reference_medians$medians
        raw_fold <- fold_anno$raw_medians / reference_medians$raw_medians

        ratio <- log10(fold_anno$percenttotal / reference_medians$percenttotal);
        colnames(ratio) <- c("percenttotalratiolog")
        is.na(ratio) <- fold_anno$count == 0 | reference_medians$count == 0

        # Merge the fold-change columns with the count, frequency, and median columns
        anno <- c(anno, list(percenttotalratiolog = ratio, fold = fold, raw_fold=raw_fold))
      }

      SPADElike_writeGraph(SPADElike_annotateGraph(graph, layout=layout_table, anno=anno), paste(f,".medians.gml",sep=""), format="gml")
      # We save an R native version of the annotations to simpify plotting, and other downstream operations
      anno <- SPADElike_flattenAnnotations(anno)
      for (c in colnames(anno)) { attr_values[[c]] <- c(attr_values[[c]], anno[,c]) }
      save(anno, file=paste(f,"anno.Rsave",sep="."))
    }
  }

  # Compute the global limits (cleaning up attribute names to match those in GML files)
  pctile_color=c(0.02,0.98)
  attr_ranges <- t(sapply(attr_values, function(x) { quantile(x, probs=c(0.00, pctile_color, 1.00), na.rm=TRUE) }))
  rownames(attr_ranges) <- sapply(rownames(attr_ranges), function(x) { gsub("[^A-Za-z0-9_]","",x) })
  write.table(attr_ranges, paste0(outputs_dir,"/global_boundaries.table"), col.names=FALSE)

  return(list("graph" = graph, "layout_table" = layout_table))
}
