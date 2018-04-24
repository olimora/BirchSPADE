BirchSPADE.save_original_data_with_clusters <- function(input_file_full
                                                        ,outputs_dir
                                                        ,clusters
                                                        ,comp
                                                        ,transforms
                                                        ,remove_outliers
                                                        ,rows_without_outliers) {

  in_fcs  <- SPADElike_transformFCS(SPADElike_readFCS(input_file_full, comp), transforms);
  cells_data_og = exprs(in_fcs)
  if (remove_outliers) {
    cells_data_og = cells_data_og[rows_without_outliers,]
  }
  cells_data_og = cbind(cells_data_og, "cluster"=clusters)
  params <- parameters(in_fcs)
  pd <- pData(params)

  # Add column named "cluster" to the FCS file
  channel_number <- ncol(in_fcs)+1;
  channel_id     <- paste("$P",channel_number,sep="");
  channel_name   <- "cluster";
  channel_range  <- max(clusters)+1;

  plist <- matrix(c(channel_name,channel_name,channel_range,0,channel_range-1));
  rownames(plist) <- c("name","desc","range","minRange","maxRange");
  colnames(plist) <- c(channel_id);
  pd <- rbind(pd,t(plist))
  pData(params) <- pd

  # and save the file
  out_frame <- flowFrame(cells_data_og, parameters = params, description=description(in_fcs))
  input_file <- basename(input_file_full)
  data_with_clusters_file = paste0(outputs_dir,"/",input_file,".cluster.fcs")
  suppressWarnings(write.FCS(out_frame, data_with_clusters_file))
  return(list("fcs_params" = parameters(in_fcs)@data, "file" = data_with_clusters_file))
}
