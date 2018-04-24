SPADElike_statistics <- function(outputs_dir, files) {
  ### Produce statistics tables ###
  # message("Producing tables...")
  dir.create(paste(outputs_dir,'tables',sep='/'),recursive=TRUE,showWarnings=FALSE)
  # Find the files
  files <- dir(outputs_dir,full.names=TRUE,pattern=glob2rx("*.anno.Rsave"))
  # Find all the params
  params <- unique(unlist(c(sapply(files, function(f) { load(f); colnames(anno); })))) #Update to remove redundant file writing

  # Transposition 1: Rows are nodes, cols are files, files are params
  dir.create(paste(outputs_dir,'tables','byAttribute',sep='/'),recursive=TRUE,showWarnings=FALSE)
  for (p in params) {
    pivot <- c()
    names <- c()
    for (f in files){
      load(f)
      f = basename(f)
      if (p %in% colnames(anno)) {
        pivot <- cbind(pivot, anno[,p])
        names <- c(names, f)
      }
    }
    names <- gsub("[[:alnum:][:punct:]]+/output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave", "\\1", names);
    pivot <- cbind(1:nrow(pivot),pivot)
    colnames(pivot) <- c("name", names)
    if (!is.null(pivot) && ncol(pivot) > 0) {
      write.csv(pivot, file=paste(outputs_dir,'/tables/byAttribute/',p,'_table','.csv',sep=''), row.names=FALSE)
    }
  }

  byNodeData = list()
  # Transposition 2: Rows are nodes, cols are params, files are files
  dir.create(paste(outputs_dir,'tables','bySample',sep='/'),recursive=TRUE,showWarnings=FALSE)
  for (f in files) {
    load(f)
    f = basename(f)
    pivot <- anno
    names <- colnames(pivot)
    pivot <- cbind(1:nrow(pivot),pivot)
    colnames(pivot) <- c("ID", names)
    name <- gsub("output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave", "\\1", f)
    write.csv(pivot, file=paste(outputs_dir,'/tables/bySample/',name,'_table','.csv',sep=''), row.names=FALSE)
  }

  # Transposition 3: Rows are params, cols are files, files are nodes
  dir.create(paste(outputs_dir,'tables','byNodeID',sep='/'),recursive=TRUE,showWarnings=FALSE)
  for (node in rownames(pivot)){
    tableData = list()
    for (f in files){
      load(f)
      f = basename(f)
      tableData[[f]]= unlist(anno[node,,drop=T])
    }
    tableData = do.call("cbind",tableData)
    write.csv(tableData, file=paste(outputs_dir,'/tables/byNodeID/',node,'_table','.csv',sep=''), row.names=TRUE,quote=F)
  }
}
