#' Protein matrix normalization
#'
#' The function normalizes protein matrix within an optima object using
#' CLR method.
#'
#' @param optima.obj optima object.
#' @import compositions
#' @return An optima object with protein matrix being normalized and
#'  protein.normalize.method label updated to "normalized".
#' @keywords optima.obj
#' @export
#' @examples normalizeProtein(optima.object)

normalizeProtein <- function(optima.object) {
  # extract count matrix
  inputMatrix <- optima.object@protein.mtx
  # apply normalization CLR method
  ret <- (compositions::clr(inputMatrix + 1))

  optima.object@protein.mtx <- as.matrix(ret)
  optima.object@protein.normalize.method <- "normalized"
  return(optima.object)
}

#' Identify signature protein function
#'
#' This function compares protein levels for a input cell type against all other
#' cells using t test. This function returns a data frame ranked by
#' FDR adjusted p-value.
#' If two cell types specified, then this function compares protein expression
#' between the two cell types.
#'
#' @param optima.obj optima object.
#' @param cell.type Input cell type to compare protein level to all other cell types.
#' @param cell.type.2 Second cell type to compare.
#' @return Data frame of all proteins p-values comparing protein levels of input
#' cell type to all other cell types. In addition, it provides the mean difference between
#' Input cell type and other cells. If mean difference is positive, this mean more expression
#' in the cell type.
#' @keywords optima.obj, cell.type
#' @export
#' @examples findSignature(optima.obj, cell.type)

findSignature <- function(optima.obj, cell.type, cell.type.2 = "all"){
  # check data is normalized
  stopifnot(optima.obj@protein.normalize.method == "normalized")

  # check the cell type is included in dataset is correct
  stopifnot(cell.type %in% unique(optima.obj@cell.labels))

  # extract all proteins
  proteins <- optima.obj@proteins

  p_val_df <- data.frame(matrix(ncol = 3, nrow = length(proteins)))
  colnames(p_val_df) <- c("p_val", "p_val_adj", "mean_diff")
  rownames(p_val_df) <- proteins

  for (protein in proteins){

    # subset data
    group.1 <- optima.obj@protein.mtx[optima.obj@cell.labels == cell.type, optima.obj@proteins == protein]

    if(cell.type.2 %in% optima.obj@cell.labels){
      group.2 <- optima.obj@protein.mtx[optima.obj@cell.labels == cell.type.2, optima.obj@proteins == protein]
    } else if (cell.type.2 == "all"){
      group.2 <- optima.obj@protein.mtx[optima.obj@cell.labels != cell.type, optima.obj@proteins == protein]
    }

    mean.diff <- mean(group.1) - mean(group.2)

    # conduct analysis
    ret <- t.test(group.1, group.2)$p.value
    p_val_df[protein, "p_val"] <- ret
    p_val_df[protein, "mean_diff"] <- mean.diff
  }

  p_val_df[,"p_val_adj"] <- p.adjust(p_val_df$p_val, method="fdr")

  order <- order(p_val_df[,"p_val_adj"])
  return(p_val_df[order, ])
}



#' The getter function for Protein matrix
#'
#' This function returns the Protein matrix within the optima object.
#'
#' @param optima.obj optima object.
#' @return A matrix that contains Protein data in the optima object.
#' The row names are cell IDs, the column names are protein IDs.
#' @export
#' @examples getProteinMtx(my.obj)

getProteinMtx <- function(optima.obj){
  ret.mtx <- optima.obj@protein.mtx
  colnames(ret.mtx) <- optima.obj@proteins
  rownames(ret.mtx) <- optima.obj@cell.ids
  return(ret.mtx)
}



#' Plot specific protein levels
#'
#' This function create a plot based on protein level dimension reduction result.
#' Within the plot, Each cell was colored based on the protein level
#'
#' @param optima.obj optima object.
#' @param protein.name the specific protein user interested
#' @param reduceDim.obj dimension reduction result returned by reduceDim() function.
#' @import ggplot2
#' @return A scatter plot based on dimension reduction result. Each cell is colored
#' based on the protein level. The more expression.
#' @export
#' @examples plotProteinFeature(my.obj, "CD11b", protein.reduceDim)


plotProteinFeature <- function(optima.obj,
                               protein.name,
                               reduceDim.obj){

  # get protein expression value
  values = getProteinMtx(optima.obj)[,protein.name]

  # get protein expression value
  values <-  getProteinMtx(optima.obj)[,protein.name]

  my.df <- data.frame(UMAP1 = reduceDim.obj[[2]][[1]][,1],
                      UMAP2 = reduceDim.obj[[2]][[1]][,2],
                      normalized_exp = values)

  ggplot(my.df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = normalized_exp), size = 2)+
    theme_classic() +
    ggtitle(paste(protein.name, "normalized expression"))


}
