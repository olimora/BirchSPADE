SPADElike_flattenAnnotations <- function(annotations) {
  stopifnot(is.list(annotations))
  flat <- NULL
  for (i in seq_along(annotations)) {	
    df <- data.frame(annotations[[i]])
    
    to_paste <- names(annotations)[i] != colnames(df)
    if (any(to_paste))
      colnames(df) <- paste(names(annotations)[i],colnames(df),sep="")
    
    if (is.null(flat))
      flat <- df
    else
      flat <- cbind(flat, df)
  }
  flat
}

SPADElike_markerMedians <- function(files, num.clusters, cols=NULL, arcsinh_cofactor=NULL, transforms=flowCore::arcsinhTransform(a=0, b=0.2), cluster_cols=NULL, comp=TRUE) {
  
  if (!is.null(arcsinh_cofactor)) {
    warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
    transforms <- flowCore::arcsinhTransform(a=0, b=1/arcsinh_cofactor)
  }
  
  data  <- c()
  
  files <- as.vector(files)
  for (f in files) {
    # Load in FCS file
    in_fcs  <- SPADElike_readFCS(f,comp=comp);
    in_data <- exprs(in_fcs);
    
    params <- parameters(in_fcs);
    pd     <- pData(params);
    
    # Select out the desired columns
    if (is.null(cols)) {
      cols <- as.vector(pd$name) 
    }
    if (!"cluster" %in% cols) {
      cols <- c(cols,"cluster")
    }
    
    idxs <- match(cols,pd$name)
    if (any(is.na(idxs))) { 
      stop("Invalid column specifier") 
    }	
    
    data <- rbind(data, in_data[, idxs,drop=FALSE])
    
  }
  
  
  clst <- data[,"cluster"]
  data <- data[,colnames(data)!="cluster",drop=FALSE]
  data_t <- SPADElike_transformMatrix(data, transforms) 
  
  # TODO: Weird things were being done to the naming, and this breaks that so we can do the transforms cleanly...
  colnames(data) <- sapply(colnames(data),function(x) { 
    if (x %in% cluster_cols)
      x <- paste(x,"clust",sep="_")
    x
  })
  colnames(data_t) = colnames(data)
  
  ids  <- 1:num.clusters
  if (any(is.na(match(unique(clst),ids)))) {
    stop("More clusters in FCS files than indicated")
  }
  
  count   <- matrix(0,  nrow=num.clusters, ncol=1, dimnames=list(ids, "count"))
  medians <- matrix(NA, nrow=num.clusters, ncol=ncol(data), dimnames=list(ids,colnames(data)))
  raw_medians <- matrix(NA, nrow=num.clusters, ncol=ncol(data), dimnames=list(ids,colnames(data)))
  cvs     <- matrix(NA, nrow=num.clusters, ncol=ncol(data), dimnames=list(ids,colnames(data)))
  for (i in ids) {
    data_s  <- subset(data, clst == i)
    data_s_t <- subset(data_t, clst == i)
    
    count[i,1]  <- nrow(data_s_t)
    medians[i,] <- apply(data_s_t, 2, median)
    raw_medians[i,] <- apply(data_s, 2, median)
    cvs[i,]     <- apply(data_s_t, 2, function(d) { 100*sd(d)/abs(mean(d)) })
  } 
  percenttotal <- matrix((count / sum(count)) * 100.0, nrow=num.clusters, ncol=1, dimnames=list(ids, "percenttotal"))
  list(count=count, medians=medians, raw_medians=raw_medians, cvs=cvs, percenttotal=percenttotal)
}

SPADElike_annotateGraph <- function(graph, layout=NULL, anno=NULL) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  if (!is.null(layout) && is.matrix(layout)) {
    if (nrow(layout) != vcount(graph) || ncol(layout) != 2) {
      stop("Ill-formated layout matrix, must 2 columns (x,y) and as many rows as vertices")
    }
    # Over write non-struct graphics attributes if they exist
    v_attr <- list.vertex.attributes(graph)
    for (i in grep("^graphics$",v_attr))
      graph <- remove.vertex.attribute(graph, v_attr[i])
    
    graph <- set.vertex.attribute(graph, "graphics.x", value=layout[,1])
    graph <- set.vertex.attribute(graph, "graphics.y", value=layout[,2])
  }
  
  if (is.null(anno)) {
    return(graph)
  } else if (!is.list(anno)) {
    stop("anno must be a list with named entries");
  }
  
  for (i in seq_len(length(anno))) {
    l <- anno[[i]]
    if (!is.matrix(l)) {
      stop(paste("Argument:",quote(l),"must be a matrix"))
    }
    vt <- V(graph)[match(rownames(l),V(graph)$name)]  # Vertex IDS are 1 indexed
    for (c in colnames(l)) {
      graph <- set.vertex.attribute(graph,ifelse(names(anno)[i] == c,c,paste(names(anno)[i],c,sep="")),index=vt, value=l[,c])
    }
  }
  graph
}

SPADElike_writeGraph <- function(graph, file="", format = c("gml")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  if (file == "") {
    file <- stdout()
  } else if (is.character(file)) {
    file <- file(file, "w")
    on.exit(close(file))
  } else if (!isOpen(file, "w")) {
    open(file, "w")
    on.exit(close(file))
  }
  if (!inherits(file, "connection")) {
    stop("'file' must be a character string or a connection")
  }
  
  write.gml <- function(graph,file) {
    
    write.attr <- function(name, attr) {
      # Strip out non-alphanumeric characters to avoid Cytoscape parsing errors
      name <- gsub("[^A-Za-z0-9_]","",name)
      if (length(grep("^[0-9]",name))) {
        name <- paste("spade",name,sep="")
      }
      if (is.na(attr) || is.nan(attr))
        stop("Unexpected NA or NaN attribute")
      else if (is.character(attr) && nchar(attr) > 0)
        paste(name," \"",attr,"\"",sep="")
      else if (is.integer(attr))
        paste(name,attr)
      else
        paste(name,formatC(attr,format="f"))
    }
    
    writeLines(c("graph [", paste("directed",ifelse(is.directed(graph),1,0))),con=file)
    
    # Identify known "structs"	
    v_attr <- list.vertex.attributes(graph)
    v_attr_g <- v_attr[grep("graphics[.]",v_attr)]  # graphics attributes
    v_attr <- setdiff(v_attr, c(v_attr_g, "id"))  
    
    for (v in V(graph)) {
      writeLines("node [",con=file)
      
      writeLines(paste("id",v),con=file)
      for (a in v_attr) {
        val <- get.vertex.attribute(graph,a,index=v)
        if (!is.na(val) && !is.nan(val))
          writeLines(write.attr(a,val),con=file)
      }
      
      if (length(v_attr_g) > 0) {
        writeLines("graphics [",con=file)
        for (a in v_attr_g) {
          parts <- unlist(strsplit(a,"[.]"))
          val <- get.vertex.attribute(graph,a,index=v)
          if (!is.na(val) && !is.nan(val))
            writeLines(write.attr(parts[2],val),con=file)
        }
        writeLines("]",con=file)
      }
      
      writeLines("]",con=file)
    }
    
    # Identify known "structs"	
    e_attr <- list.edge.attributes(graph)
    if (length(grep("[.]",e_attr)) > 0) {
      stop("Unsupported struct in edge attributes")
    }
    
    for (e in E(graph)) {
      writeLines("edge [",con=file)
      
      pts <- get.edges(graph,e)
      writeLines(c(paste("source",pts[1]), paste("target",pts[2])),con=file)
      for (a in e_attr) {
        val <- get.edge.attribute(graph,a,index=e)
        if (!is.na(val) && !is.nan(val))
          writeLines(write.attr(a,val),con=file)
      }	
      
      writeLines("]",con=file)	
    }
    writeLines("]",con=file)
  }
  
  res <- switch(format,   
                gml = write.gml(graph,file), 
                stop(paste("Unsupported output format:",format)))
  invisible(res)			
}