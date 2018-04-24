SPADElike_readFCS <- function(file, comp=TRUE, verbose=FALSE, ...) {
  if (verbose)
    fcs <- read.FCS(file, ...)
  else 
    fcs <- suppressWarnings(read.FCS(file, ...))
  
  params <- parameters(fcs)
  pd     <- pData(params)
  
  # Replace any null descs with names (for FSC-A, FSC-W, SSC-A)
  bad_col <- grep("^[a-zA-Z0-9]+",pd$desc,invert=TRUE)
  if (length(bad_col) > 0) {
    keyval <- keyword(fcs)
    for (i in bad_col) {
      pd$desc[i] <- pd$name[i]
      keyval[[paste("$P",i,"S",sep="")]] <- pd$name[i]
    }
    pData(params) <- pd;
    fcs <- flowFrame(exprs(fcs),params,description=description(fcs));
    keyword(fcs) <- keyval
  }
  
  # Compensate data if SPILL or SPILLOVER present, stripping compensation matrix 
  # out of the flowFrame, i.e we should only have to do this once
  apply.comp <- function(in_fcs, keyword) {
    comp_fcs <- compensate(in_fcs, description(in_fcs)[[keyword]])
    flowFrame(exprs(comp_fcs), parameters(comp_fcs), description(comp_fcs)[grep("SPILL",names(description(comp_fcs)),invert=TRUE)])
  }
  
  if (comp && !is.null(description(fcs)$SPILL)) {
    fcs <- apply.comp(fcs, "SPILL")	
  } else if (comp && !is.null(description(fcs)$SPILLOVER)) {
    fcs <- apply.comp(fcs, "SPILLOVER")	
  }
  
  fcs
}

SPADElike_transformFCS <- function(ff, tform=NULL) {
  if (is.null(tform)) {
    ff  # No-op
  } else {
    
    if (class(tform) == "transform") {
      new_exprs <- apply(exprs(ff), 2, tform)
      new_range <- apply(range(ff), 2, tform)
    } else {
      new_exprs <- exprs(ff)
      new_range <- range(ff)
      for (name in intersect(colnames(ff), names(tform))) {
        new_exprs[,name] <- tform[[name]](new_exprs[,name])
        new_range[,name] <- tform[[name]](new_range[,name])
      }
      new_range <- as.matrix(new_range)
    }
    
    new_par <- parameters(ff)
    new_par$minRange <- new_range[1,]
    new_par$maxRange <- new_range[2,]
    new_par$range    <- (new_range[2,] - new_range[1,]) + 1
    
    flowFrame(new_exprs, new_par, description(ff))
  } 
}

SPADElike_transformMatrix <- function(mat, tform=NULL) {
  if (is.null(tform)) {
    mat  # No-op
  } else {
    if (class(tform) == "transform") {
      apply(mat, 2, tform)
    } else {
      for (name in intersect(colnames(mat), names(tform))) {
        mat[,name] <- tform[[name]](mat[,name])
      }
      mat
    }
  } 
}