BirchSPADE.remove_outliers <- function(cells_data, outliers) {
  cells_hash = vector(mode = "numeric", length = nrow(cells_data))
  cells_hash = vector(mode = "numeric", length = nrow(outliers))
  for (i in 1:nrow(cells_data)) {
    cells_hash[i] = 0
    for (j in 1:ncol(cells_data)) {
      cells_hash[i] = cells_hash[i] + (cells_data[i,j] * j)
    }
  }
  outliers_hash = vector(mode = "numeric", length = nrow(outliers))
  for (i in 1:nrow(outliers)) {
    outliers_hash[i] = 0
    for (j in 1:ncol(outliers)) {
      outliers_hash[i] = outliers_hash[i] + (outliers[i,j] * j)
    }
  }
  good_rows = which(!(cells_hash %in% outliers_hash))
  cells_data = cells_data[good_rows,]
  return(list("data" = cells_data, "rows" = good_rows))
}
