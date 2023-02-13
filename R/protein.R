#' Protein matrix normalization
#'
#' The function normalizes protein matrix within an optima object using
#' CLR method.
#'
#' @param optima.obj optima object.
#' @return An optima object with protein matrix being normalized and
#'  protein.normalize.method label updated to "normalized"
#' @import compositions
#' @keywords optima.obj
#' @export
#' @examples normalizeProtein(optima.obj)

normalizeProtein <- function(optima.obj) {
  # extract count matrix
  inputMatrix <- optima.obj@protein.mtx
  # apply normalization CLR method
  ret <- (compositions::clr(inputMatrix + 1))

  optima.obj@protein.mtx <- as.matrix(ret)
  optima.obj@protein.normalize.method <- "normalized"
  return(optima.obj)

}

#' Identify signature protein function
#'
#' This function compares protein levels for a input cell type against all other
#' cells using t test. This function returns a data frame ranked by
#' FDR adjusted p-value.
#'
#' @param optima.obj optima object.
#' @param cell.type Input cell type to compare protein level to all other cell types.
#' @return Data frame of all proteins p-values comparing protein levels of input
#' cell type to all other cell types.
#' @keywords optima.obj, cell.type
#' @export
#' @examples findSignature(optima.obj, cell.type)

findSignature <- function(optima.obj, cell.type){
  # check data is normalized
  stopifnot(optima.obj@protein.normalize.method == "normalized")

  # check the cell type is included in dataset is correct
  stopifnot(cell.type %in% unique(optima.obj@cell.labels))

  # extract all proteins
  proteins <- optima.obj@proteins

  p_val_df <- data.frame(matrix(ncol = 2, nrow = length(proteins)))
  colnames(p_val_df) <- c("p_val", "p_val_adj")
  rownames(p_val_df) <- proteins

  for (protein in proteins){

    # subset data
    group.1 <- optima.obj@protein.mtx[optima.obj@cell.labels == cell.type, optima.obj@proteins == protein]
    group.2 <- optima.obj@protein.mtx[optima.obj@cell.labels != cell.type, optima.obj@proteins == protein]

    # conduct analysis
    ret <- t.test(group.1, group.2)$p.value
    p_val_df[protein, "p_val"] <- ret
  }

  p_val_df[,"p_val_adj"] <- p.adjust(p_val_df$p_val, method="fdr")

  order <- order(p_val_df[,"p_val_adj"])
  return(p_val_df[order, ])
}
