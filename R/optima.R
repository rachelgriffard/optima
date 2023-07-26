#' optima object
#'
#' An optima object contains DNA, protein and CNV for Tapestri platform
#' single cell sequencing data.
#'
#' @slot meta.data User-defined metadata can be kept with the object.
#' @slot cell.ids A vector of cell IDs/barcodes from Tapestri. This
#' vector should contain unique IDs.
#' @slot cell.labels A vector that is used to store the cell type information
#' for each cell.
#' @slot variants A vector of variant IDs.
#' @slot variant.filter A string that keeps track of if optima object is being QC filtered
#' on its variant matrix.
#' @slot vaf.mtx Variant matrix.
#' @slot gt.mtx Genotype matrix, The integers within matrix are reference call GT=0, heterozygous call GT=1, homozygous call GT=2, no calls GT = 3.
#' @slot dp.mtx Sequencing depth matrix.
#' @slot gq.mtx Genotype quality.
#' @slot amps A vector of CNV locus.
#' @slot amp.normalized.method A string that keeps track of if optima object is being normalized
#' on its CNV matrix.
#' @slot amp.mtx CNV matrix.
#' @slot ploidy.mtx Ploidy matrix.
#' @slot proteins A vector of surface protein ID.
#' @slot protein.normalize.method A string that keeps track of if optima object is being normalized
#' on its protein matrix.
#' @slot protein.mtx Protein matrix
#' @return Object containing DNA, protein and CNV single cell sequencing data.
#' @examples setClass()
#' @exportClass optima
optima <- setClass("optima", slots=list(meta.data = "character",
                                        # "id" vectors contains cell barcode
                                        cell.ids = "character",
                                        # cell labels vector contains cell clone info
                                        # could be assigned manually or using clustering method
                                        cell.labels = "character",

                                        # "variants" vector stores variant labels
                                        variants = "character",
                                        variant.filter = "character",
                                        vaf.mtx = "matrix",
                                        gt.mtx = "matrix",
                                        dp.mtx = "matrix",
                                        gq.mtx = "matrix",

                                        # CNVs
                                        amps = "character",
                                        amp.normalize.method = "character",
                                        amp.mtx = "matrix",
                                        ploidy.mtx = "matrix",

                                        # The "protein"  vector stores protein labels
                                        proteins = "character",
                                        protein.normalize.method = "character",
                                        protein.mtx =  "matrix"))

setMethod(
  f = "show",
  signature = "optima",
  function(object) {
    cat("optima object with: \n")
    cat(paste(length(object@cell.ids), "cells\n"))
    cat(paste(length(object@variants), "variants, data", object@variant.filter, "\n"))
    cat(paste(length(object@amps), "CNVs, data", object@amp.normalize.method, "\n"))

    if(object@proteins[1] == "non-protein"){
      cat("no protein data")
    } else {
      cat(paste(length(object@proteins), "proteins, data", object@protein.normalize.method, "\n"))
      cat(paste("Current unique cell labels includes: ", paste(unique(object@cell.labels),  collapse = ", "), "\n"))
    }
  }
)



#' @export
#' @import methods
#' @method names optima
#'
# display all names in the omele object names() function
names.optima <- function(x) {
  slotNames(x)
}

#' @export
# display data in slots using dollar sign
"$.optima" <- function(optima, i) {
  return(slot(object = optima, name = i))
}

#' @export
# set function
"$<-.optima" <- function(optima, i, value) {
  slot(object = optima, name = i) = value
  return(optima)
}


#' Generate elbow plot
#'
#' This function generate an elbow plot. The plot is useful for determining how many PCs will be used for dimension reduction
#' The input matrix can be the DNA matrix in an optima object.
#'
#' @param input.mtx Input matrix
#' @return an elbow plot
#' @examples elbowPlot(example.matrix)
#' @export
elbowPlot <- function(input.mtx){

  pca <- prcomp(input.mtx)

  pca.sd <- pca$sdev

  plot(x = 1:length(pca.sd),
       y = pca.sd,
       xlab = "PC",
       ylab = "Standard Deviation")

}



#' Dimension reduction function
#'
#' This function reduces dimensions for a data matrix, such data matrix
#' can be protein or DNA matrix in an optima object.
#'
#' @param input.mtx Input optima object.
#' @param num.PC The number of PCs being used for preserving variation. Default is 5
#' @import umap
#' @return List containing PCA result and umap result derived from first 5 PCs
#' @examples reduceDim(example.matrix)
#' @export

reduceDim <- function(input.mtx, num.PC = 5){
  pca <- prcomp(input.mtx)
  my_umap <- umap::umap(pca$x[,1:num.PC])
  return(list(pca, my_umap))
}


#' Heatmap function
#'
#' This function creates a heatmap using matrix data. Matrix
#' data can be DNA or protein.
#'
#' @param optima.obj optima object.
#' @param omic.type Type of data for heat map. Potential values "dna" and "protein". If "dna"
#' is specified, then the VAF data will be used for plot. If "protein" is specified, then
#' protein expression values will be used.
#' @import pheatmap
#' @return Heat map visualization.
#' @examples drawHeatMap()
#' @export

drawHeatmap <- function(optima.obj, omic.type){

  # check the number of cluster in the cell label
  num.cluster <- length(unique(optima.obj@cell.labels))

  if(num.cluster == 1){
    stop("only one unique cluster in cell labels")
  } else if(num.cluster > 10){
    stop("More than 10 clusters generated, refine parameters in getClone() function")
  } else {

    # order cells based on clustering result
    cell.order <- order(optima.obj@cell.labels)

    # depending on the omic.type parameter
    if(omic.type == "dna"){

      # get variant matrix
      input.mtx <- optima.obj@vaf.mtx[cell.order,]
      colnames(input.mtx) <- optima.obj@variants

    } else if (omic.type == "protein"){

      # get protein matrix
      input.mtx <- optima.obj@protein.mtx[cell.order,]
      colnames(input.mtx) <- optima.obj@proteins

    } else {
      stop("omic.type parameter can only be dna or protein")
    }

    rownames(input.mtx) <- optima.obj@cell.ids[cell.order]

    # generate dataframe for annotation purpose.
    row.anno <- data.frame(optima.obj@cell.labels[cell.order])
    rownames(row.anno) <- rownames(input.mtx)
    colnames(row.anno) <- "cell type"

    ##########################################################
    # Heatmap
    if(omic.type == "protein"){
      heatmap.title <- paste("Heatmap for",
                             optima.obj$meta.data,
                             optima.obj$protein.normalize.method,
                             omic.type, "expression")
    } else if (omic.type == "dna"){
      heatmap.title <- paste("Heatmap for",
                             optima.obj$meta.data,
                             optima.obj$variant.filter,
                             "DNA VAF")
    }


    pheatmap::pheatmap(input.mtx,
                       treeheight_row = 0,
                       treeheight_col = 0,
                       cluster_rows=FALSE,
                       annotation_row = row.anno,
                       show_rownames=FALSE,
                       main = heatmap.title)

  }
}
