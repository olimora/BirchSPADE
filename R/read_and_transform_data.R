BirchSPADE.load_input_data <- function (input_file_full,
                                        transforms,
                                        markers,
                                        normalization) {

  comp = TRUE
  in_fcs  <- SPADElike_transformFCS(SPADElike_readFCS(input_file_full, comp), transforms);
  cells_data <- exprs(in_fcs)
  colnames(cells_data) = parameters(in_fcs)@data$desc # rename colnames from proteins to markers
  cells_data = cells_data[,markers]         # keep only the needed ones
  rm(in_fcs) # remove from memory to save space

  # normalize loaded data
  if (normalization == "none") {
    # print("No normalization.")
  } else if (normalization == "minmax") {
    # print("Minmax normalization.")
    library(caret)
    pp = preProcess(cells_data, method = "range")
    cells_data = predict(pp, cells_data)
    rm(pp)
  } else if (normalization == "meanstd") {
    # print("Meanstd normalization.")
    cells_data = scale(cells_data)
  }

  return(cells_data)
}
