SPADE.layout.arch <-  function(mst_graph) {
  if (!is.igraph(mst_graph)) {
    stop("Input has to be igraph object")
  }
  if (!is.connected(mst_graph)) {
    stop("Cannot handle graph that has disjoint components")
  }

  # Make sure it is a tree, no circles
  if (girth(mst_graph)$girth > 0) {
    stop("Cannot handle graphs with cycles");
  }

  # Find the distance between nodes measured in "hops"
  hops <- c()
  for (v in V(mst_graph)) {
    hops <- rbind(hops, unlist(lapply(get.shortest.paths(mst_graph,v)$vpath, length)))
  }


  # Compute the positions for each vertices
  # --------------------------------------------------------------------------
  v_pos <- array(0,c(vcount(mst_graph),2))

  # The longest path is the backbone arch
  terminals <- which(hops == max(hops), arr.ind=TRUE)[1,]
  back_bone <- unlist(get.shortest.paths(mst_graph, from=terminals["row"], to=terminals["col"]))
  # Layout the backbone arch along an arc
  bb_span <- pi * .55  # span in radians of back bone arch
  bb_unit <- 50.0  # unit distance bewteen nodes
  angles  <- seq(pi/2-bb_span/2, pi/2+bb_span/2, length.out=length(back_bone))
  v_pos[back_bone,] <- bb_unit*length(back_bone)*cbind(cos(angles),-sin(angles))

  # Layout the side chains as trees normal to the backbone
  for (v in back_bone) {
    # Find subset of vertices that compose side chain by deleting links between current
    # backbone vertex and rest of the backbone and then performing subcomponent
    # Note: E(mst_graph,P=c(mapply(c,v,n))) == E(mst_graph)[v %--% n] but is much faster
    n <- intersect(neighbors(mst_graph, v), back_bone)
    side_v <- sort(subcomponent(delete.edges(mst_graph,E(mst_graph,P=c(mapply(c,v,n)))),v))

    # Compute layout for side chains and integrate it into overall layout
    # Note: Relies on side_v being in sorted order
    if (length(side_v) > 1) {
      # Convert side chains to directed graph, deleting edges with decreasing hop distance
      side_h <- hops[v, side_v]  # hops between back_bone node and side chain
      sg <- induced_subgraph(mst_graph, as.vector(side_v))
      side_g <- as.directed(sg, mode="mutual")
      e <- get.edges(side_g,E(side_g))  # edges as a matrix
      side_g <- delete.edges(side_g,subset(E(side_g),side_h[e[,1]] > side_h[e[,2]]))

      # Layout side chain
      # -----------------------------------------------------------------
      root <- which.min(side_h)
      layout <- layout.reingold.tilford(side_g,root=root)

      # rotate tree to be normal to back bone
      polar <- cbind(atan2(layout[,2],layout[,1]), sqrt(rowSums(layout^2)))
      polar[,1] <- polar[,1] + angles[back_bone==v] - pi/2
      layout <- bb_unit*polar[,2]*cbind(cos(polar[,1]),-sin(polar[,1]))
      # translate layout to back_bone
      layout <- layout + matrix(v_pos[v,] - layout[root,], nrow=nrow(layout), ncol=2, byrow=TRUE)
      v_pos[side_v,] <- layout
    }
  }
  v_pos
}

SPADElike_strip_sep <- function(name) {
  ifelse(substr(name,nchar(name),nchar(name))==.Platform$file,substr(name,1,nchar(name)-1),name)
}

SPADElike_normalize_out_dir <- function(out_dir) {
  out_dir <- SPADElike_strip_sep(out_dir)
  out_dir_info <- file.info(out_dir)
  if (is.na(out_dir_info$isdir)) {
    dir.create(out_dir)
  }
  if (!file.info(out_dir)$isdir) {
    stop(paste("out_dir:",out_dir,"is not a directory"))
  }
  out_dir <- paste(SPADElike_strip_sep(out_dir),.Platform$file,sep="")
}

subplot <- function(fun, x, y=NULL, size=c(1,1), vadj=0.5, hadj=0.5,
                    inset=c(0,0), type=c('plt','fig'), pars=NULL){

  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))

  type <- match.arg(type)

  if(missing(x)) x <- locator(2)

  if(is.character(x)) {
    if(length(inset) == 1) inset <- rep(inset,2)
    x.char <- x
    tmp <- par('usr')
    x <- (tmp[1]+tmp[2])/2
    y <- (tmp[3]+tmp[4])/2

    if( length(grep('left',x.char, ignore.case=TRUE))) {
      x <- tmp[1] + inset[1]*(tmp[2]-tmp[1])
      if(missing(hadj)) hadj <- 0
    }
    if( length(grep('right',x.char, ignore.case=TRUE))) {
      x <- tmp[2] - inset[1]*(tmp[2]-tmp[1])
      if(missing(hadj)) hadj <- 1
    }
    if( length(grep('top',x.char, ignore.case=TRUE))) {
      y <- tmp[4] - inset[2]*(tmp[4]-tmp[3])
      if(missing(vadj)) vadj <- 1
    }
    if( length(grep('bottom',x.char, ignore.case=TRUE))) {
      y <- tmp[3] + inset[2]*(tmp[4]-tmp[3])
      if(missing(vadj)) vadj <- 0
    }
  }

  xy <- xy.coords(x,y)

  if(length(xy$x) != 2){
    pin <- par('pin')
    tmp <- cnvrt.coords(xy$x[1],xy$y[1],'usr')$plt

    x <- c( tmp$x - hadj*size[1]/pin[1],
            tmp$x + (1-hadj)*size[1]/pin[1] )
    y <- c( tmp$y - vadj*size[2]/pin[2],
            tmp$y + (1-vadj)*size[2]/pin[2] )

    xy <- cnvrt.coords(x,y,'plt')$fig
  } else {
    xy <- cnvrt.coords(x,y,'usr')$fig
  }

  par(pars)
  if(type=='fig'){
    par(fig=c(xy$x,xy$y), new=TRUE)
  } else {
    par(plt=c(xy$x,xy$y), new=TRUE)
  }
  fun
  tmp.par <- par(no.readonly=TRUE)

  return(invisible(tmp.par))
}

SPADElike_plotTrees <- function(graph, files, params, file_pattern="*anno.Rsave",
                                out_dir=".", layout=SPADE.layout.arch,
                                attr_pattern="percent|medians|fold|cvs", scale=NULL,
                                pctile_color=c(0.02,0.98), normalize="global",
                                size_scale_factor=1, edge.color="grey",
                                bare=FALSE, palette="bluered",
                                cluster_markers = NULL) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.null(scale) && (!is.vector(scale) || length(scale) !=2)) {
    stop("scale must be a two element vector")
  }
  if (!is.vector(pctile_color) || length(pctile_color) != 2) {
    stop("pctile_color must be a two element vector with values in [0,1]")
  }

  if (length(files) == 1 && file.info(SPADElike_strip_sep(files))$isdir) {
    files <- dir(SPADElike_strip_sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))
  }
  out_dir <- SPADElike_normalize_out_dir(out_dir)

  load_attr <- function(save_file) {
    anno <- NULL
    l <- load(save_file)
    stopifnot(l == "anno")
    # Note "anno" populated by load operation
    return(anno)
  }

  boundaries <- NULL
  if (normalize == "global") {
    boundaries <- c()  # Calculate ranges of all encountered attributes with trimmed outliers
    all_attrs  <- c()
    for (f in files) {
      attrs <- load_attr(f)
      for (i in grep(attr_pattern, colnames(attrs))) {
        n <- colnames(attrs)[i]
        all_attrs[[n]] <- c(all_attrs[[n]], attrs[,i])
      }
    }
    for (i in seq_along(all_attrs)) {
      boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]], probs=pctile_color, na.rm=TRUE)
    }
  }

  if (is.function(layout))
    graph_l <- layout(graph)
  else
    graph_l <- layout

  if (palette == "jet")
    palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  else
    if (palette == "bluered")
      palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
  else
    stop("Please use a supported color palette.  Options are 'bluered' or 'jet'")

  colorscale <- palette(100)

  for (f in files) {
    attrs <- load_attr(f)

    vsize <- attrs$percenttotal
    vsize <- vsize/(max(vsize,na.rm=TRUE)^(1/size_scale_factor)) * 3 + 2
    vsize[is.na(vsize) | (attrs$count == 0)] <- 1

    for (i in grep(attr_pattern, colnames(attrs))) {
      name <- colnames(attrs)[i]

      # Compute the color for each vertex using color gradient scaled from min to max
      attr <- attrs[,i]
      if (!is.null(scale))
        boundary <- scale
      else if (normalize == "global") { # Scale to global min/max
        boundary <- boundaries[[name]] # Recall trimmed global boundary for this attribute
      } else # Scale to local min/max
        boundary <- quantile(attr, probs=pctile_color, na.rm=TRUE)  # Trim outliers for this attribtue

      if (length(grep("^medians|percent|cvs", name)))
        boundary <- c(min(boundary), max(boundary))  # Dont make range symmetric for median or percent values
      else
        boundary <- c(-max(abs(boundary)), max(abs(boundary)))  # Make range symmetric for fold-change and ratio values

      boundary <- round(boundary, 2) # Round boundary values to 2 digits of precision so scale looks nice

      # print("boundaries")
      # print(boundary)
      # print(boundary[1])
      # print(boundary[2])

      if (is.na(boundary[1]) && is.na(boundary[2])) {
        boundary = c(0, 0)
      }
      if (boundary[1] == boundary[2]) {  boundary <- c(boundary[1]-1, boundary[2]+1); }  # Prevent "zero" width gradients

      # print("boundaries now")
      # print(boundary)

      grad <- seq(boundary[1], boundary[2], length.out=length(colorscale))

      color <- colorscale[findInterval(attr, grad,all.inside=TRUE)]
      color[is.na(attr) | (attrs$count == 0)] <- "grey"
      if (grepl("^percenttotalratiolog$",name)) {
        # Color nodes with "infinite" ratios with top of color scale
        color[is.na(attr) & attrs$count > 0] <- tail(colorscale,1)
      }

      # Use "empty" circles for nodes with no cells, or otherwise "invalid" information
      fill_color  <- color
      is.na(fill_color) <- is.na(attr)
      frame_color <- color

      # add marker name to the file name too
      protein_name = ""
      marker = ""
      if (grepl("^medians|raw_medians|cvs|raw_fold|fold", name)) {
        protein_name <- gsub(pattern = "^medians|raw_medians|cvs|raw_fold|fold", replacement = "", x = name)
        marker <- params$desc[params$name == protein_name]
        marker <- gsub(pattern = " ", replacement = "", x = marker)
        marker = gsub(pattern = " |/|:|\\?|<|>|\"|\\*|\\||\\\\", replacement = "", x = marker)
        name = gsub(x = name, pattern = protein_name, replacement = paste0(".",marker))
        protein_name = paste0(".", protein_name)
      }
      # Plot the tree, with legend showing the gradient
      pdf(paste(out_dir,basename(f),".",name,protein_name,".pdf",sep=""))
      graph_aspect <- ((max(graph_l[,2])-min(graph_l[,2]))/(max(graph_l[,1])-min(graph_l[,1])))
      plot(graph, layout=graph_l, vertex.shape="circle", vertex.color=fill_color, vertex.frame.color=frame_color, edge.color=edge.color, vertex.size=vsize, vertex.label=NA, edge.arrow.size=.25, edge.arrow.width=1, asp=graph_aspect)

      protein_name = gsub(pattern = "^\\.", replacement = "", x = protein_name)
      name = gsub(pattern = paste0("\\.",marker), replacement = marker, x = name)
      name = paste0(name, " (", protein_name, ")")

      if (!bare) {
        # Substitute pretty attribute names
        if (length(grep("^medians", name)))
          name <- sub("medians", "Median of ", name)
        else if (length(grep("^fold", name)))
          name <- sub("fold", "Arcsinh diff. of ", name)
        else if (grepl("^percenttotal$", name))
          name <- sub("percent", "Percent freq. of ", name)
        else if (grepl("^percenttotalratiolog$", name))
          name <- "Log10 of Ratio of Percent Total of Cells in Each Cluster"
        else if (grepl("^cvs", name))
          name <- sub("cvs", "Coeff. of Variation of ", name)

        # Make parameters used for clustering obvious
        if (!is.null(cluster_markers)) {
          if (marker %in% cluster_markers) {
            name <- paste(name, "\n(Used for tree-building)")
          }
        }

        title(main=paste(strsplit(basename(f),".fcs")[[1]][1], sub=name, sep="\n"))
        subplot(
          image(
            grad, c(1), matrix(1:length(colorscale),ncol=1), col=colorscale,
            xlab=ifelse(is.null(scale),paste("Range:",pctile_color[1],"to",pctile_color[2],"pctile"),""),
            ylab="", yaxt="n", xaxp=c(boundary,1)
          ),
          x="right,bottom",size=c(1,.20)
        )
      }
      dev.off()
    }
  }
}

cnvrt.coords <- function(x,y=NULL,input=c('usr','plt','fig','dev','tdev')) {

  input <- match.arg(input)
  xy <- xy.coords(x,y, recycle=TRUE)

  cusr <- par('usr')
  cplt <- par('plt')
  cfig <- par('fig')
  cdin <- par('din')
  comi <- par('omi')
  cdev <- c(comi[2]/cdin[1],(cdin[1]-comi[4])/cdin[1],
            comi[1]/cdin[2],(cdin[2]-comi[3])/cdin[2])

  if(input=='usr'){
    usr <- xy

    plt <- list()
    plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
    plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])

    fig <- list()
    fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
    fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]

    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

  if(input=='plt') {

    plt <- xy

    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]

    fig <- list()
    fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
    fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]

    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

  if(input=='fig') {

    fig <- xy

    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])

    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]

    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

  if(input=='dev'){
    dev <- xy

    fig <- list()
    fig$x <- (dev$x-cfig[1])/(cfig[2]-cfig[1])
    fig$y <- (dev$y-cfig[3])/(cfig[4]-cfig[3])

    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])

    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

  if(input=='tdev'){
    tdev <- xy

    dev <- list()
    dev$x <- (tdev$x-cdev[1])/(cdev[2]-cdev[1])
    dev$y <- (tdev$y-cdev[3])/(cdev[4]-cdev[3])

    fig <- list()
    fig$x <- (dev$x-cfig[1])/(cfig[2]-cfig[1])
    fig$y <- (dev$y-cfig[3])/(cfig[4]-cfig[3])

    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])

    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

}
