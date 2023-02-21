#' CNV normalization function
#'
#' The function normalizes the CNV matrix to correct for column-wise and row-wise variation and
#' updates the optima object amp.normalize.method from "unnormalized" to "normalized".
#'
#' @param optima.obj optima object.
#' @return optima object with normalized CNV and amp.normalize.method updated to "normalized".
#' @keywords optima.obj
#' @export
#' @examples normalizeCNV(optima.obj)

normalizeCNV <- function(optima.obj){
  cnv.mtx <- optima.obj@amp.mtx

  # conduct normalize
  rowsum.threshold <- sort(rowSums(cnv.mtx), decreasing = TRUE)[10] / 10
  keep.tf <- rowSums(cnv.mtx) > rowsum.threshold
  cnv.mtx <- cnv.mtx / (apply(cnv.mtx, 1, mean) + 1)
  cnv.mtx <- t(t(cnv.mtx) / (apply(cnv.mtx[keep.tf,], 2, median) + 0.05))
  cnv.mtx <- cnv.mtx * 2

  # update object
  optima.obj@amp.mtx <- cnv.mtx
  optima.obj@amp.normalize.method <- "normalized"

  return(optima.obj)
}

#' CNV ploidy calculation function
#'
#' The function uses the normalized CNV matrix to calculate the ploidy for each CNV locus.
#'
#' @param optima.obj optima object.
#' @param diploid.cell The cell type that should be considered as diploid cell. This
#' cell type should be one of cell types in the cell.labels vector.
#' @return optima object with normalized CNV
#' @keywords optima.obj with ploidy.mtx being updated
#' @export
#' @examples calculatePloidy(optima.obj)

calculatePloidy <- function(optima.obj, diploid.cell){

  # check normalization is done
  stopifnot(optima.obj@amp.normalize.method == "normalized")
  stopifnot(length(optima.obj@cell.labels) == nrow(optima.obj@amp.mtx))
  stopifnot(length(unique(optima.obj@cell.labels)) > 1)

  # get cnv matrix
  cnv.mtx <- optima.obj@amp.mtx

  # subset diploid cells
  ref.cnv.mtx <- cnv.mtx[optima.obj@cell.labels == diploid.cell,]
  median.cnv <- apply(ref.cnv.mtx, 2, median)

  median.cnv[median.cnv == 0] <- 1
  ret.ploidy <- t(2*t(cnv.mtx)/median.cnv)

  optima.obj@ploidy.mtx <- ret.ploidy

  return(optima.obj)
}

#' Ploidy scatter plot function
#'
#' For a specified cell type, this function creates a scatter plot indicating
#' ploidy for different CNV loci.
#'
#' @param optima.obj optima object.
#' @param cell.type String that indicates which cell type,
#' @return optima object with normalized CNV.
#' @keywords optima.obj
#' @export
#' @examples
#' plotPloidy()

plotPloidy <- function(optima.obj, cell.type){

  # check data size matches
  stopifnot(dim(optima.obj@amp.mtx) == dim(optima.obj@ploidy.mtx))
  stopifnot(length(unique(optima.obj@cell.labels)) > 1)

  # take median
  num.ploidy <- apply(optima.obj@ploidy.mtx[optima.obj@cell.labels == cell.type,] ,2, median)

  # plot
  par(mar=c(10, 4, 4, 4), xpd=TRUE)
  plot(1:length(num.ploidy), num.ploidy,
       xlab = "",
       xaxt="n",
       ylab = "ploidy",
       main = paste("cell type", cell.type))
  axis(1, at = 1:length(num.ploidy),
       labels = optima.obj@amps,
       las=2,
       cex=0.5)

  title(xlab="CNV loci",
        line=8,
        cex.lab=1.2,
        family="Calibri Light")

}



#' The getter function for CNV matrix
#'
#' This function returns the CNV matrix within the optima object.
#'
#' @param optima.obj optima object.
#' @return A matrix that contains CNV data in the optima object.
#' The row names are cell id, the column names are CNV IDs.
#' @export
#' @examples getCNVmtx(my.obj)

getCNVmtx <- function(optima.obj){
  ret.mtx <- optima.obj@amp.mtx
  colnames(ret.mtx) <- optima.obj@amps
  rownames(ret.mtx) <- optima.obj@cell.ids
  return(ret.mtx)
}


#' The getter function for Ploidy matrix
#'
#' This function returns the Ploidy matrix within the optima object.
#'
#' @param optima.obj optima object.
#' @return A matrix that contains Ploidy data in the optima object.
#' The row names are cell id, the column names are CNV IDs.
#' @export
#' @examples getPloidyMtx(my.obj)

getPloidyMtx <- function(optima.obj){
  ret.mtx <- optima.obj@ploidy.mtx
  colnames(ret.mtx) <- optima.obj@amps
  rownames(ret.mtx) <- optima.obj@cell.ids
  return(ret.mtx)
}
