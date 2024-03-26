#' DNA Variant filter function
#'
#'
#' This function uses multiple matrices imported from the h5 file to conduct quality filtering.
#' This includes the sequencing depth matrix, genotype matrix, variant allele frequency matrix, genotype quality matrix
#' The function returns an optima object that has been filtered with variant/cells.
#' In addition, the returned optima object's variant.filter label is changed to "filtered".
#' This function is usually applied before protein and CNV analysis.
#'
#' @param optima.obj An optima object with raw data unfiltered.
#' @param min.dp Minimum depth, defaults to 10. Read depth ranges from 0 to positive infinity. When to change: If your sample ended up with more cells being sequenced than planned, then you may have less reads per cell. In such cases, you could try to lower this parameter to 7, 8 or 9. It is not common to increase this parameter above 10.
#' @param min.gq Minimum genotype quality. The default value is 30. The possible range of GQ score is 0-99. This GQ score is derived from GATK.  A higher score means a more confident genotype call. A lower score means genotype call with low confidence. When to change: You may consider decreasing the threshold when the reads sequencing quality is low overall.
#' @param vaf.ref If reference call vaf (GT=0) is larger than vaf.ref, then value in genotype call matrix is converted to GT=3.
#' @param vaf.hom If homozygous call vaf (GT=2) is smaller than vaf.hom, then value in genotype call matrix is converted to GT=3.
#' @param vaf.het If heterozygous call vaf (GT=1) is smaller than vaf.ref, then value in genotype call matrix is converted to GT=3. If heterozygous call vaf (GT=1) is smaller than vaf.ref, then value in genotype call matrix is converted to GT=3. The default value is 35. For heterozygous cell, there should be 50/50 in two alleles, respectively. However, due to sampling error, there are chances when we sequenced more reads in one allele than the other. This parameter allows for tolerating such situations. We don’t normally change this parameter.
#' @param min.cell.pt Minimum threshold for cell percentage that has valid variant call (GT = 0, 1 or 2) after applying the filter. Minimum threshold for cell percentage that has valid variant call (GT = 0, 1 or 2) after applying the filter. The default value is 50%. This means for one variant; we need at least 50% of cells have a valid variant call. When to change: If the variant of interest is in a high GC content region, then PCR amplification is hard. In such cases, you may choose to decrease the percent to 30 or 40 so that your interested variant could come through the filter.
#' @param min.mut.cell.pt Minimum threshold for cell percentage that has mutated genotype (GT = 1 or 2) after applying the filter. The default is 1, corresponds to 1%. This filter is used to remove false positives. When to change: If you know the variant is rare in the data, then you could try lower threshold to try to keep the variant in your dataset.
#' @return An optima object, The DNA data in the object is filtered, the variant.filter label is "filtered".
#' Meanwhile, the protein matrix and CNV matrix is also updated so that only cells withstand DNA variant filter are kept.
#' @keywords filter DNA
#' @export
#' @examples
#' filterVariant(optima.obj)

filterVariant <- function(optima.obj,
                          min.dp=10,
                          min.gq=30,
                          vaf.ref=5,
                          vaf.hom=95,
                          vaf.het=35,
                          min.cell.pt=50,
                          min.mut.cell.pt=1){

  # DNA variant object
  vaf.mtx <- optima.obj@vaf.mtx
  gt.mtx <- optima.obj@gt.mtx
  dp.mtx <- optima.obj@dp.mtx
  gq.mtx <- optima.obj@gq.mtx

  # if True in the matrix, then the value should be dropped
  dp.tf <- dp.mtx < min.dp
  gq.tf <- gq.mtx < min.gq
  vaf.ref.tf <- (vaf.mtx > vaf.ref) & (gt.mtx == 0)
  vaf.hom.tf <- (vaf.mtx < vaf.hom) & (gt.mtx == 2)
  vaf.het.tf <- (vaf.mtx < vaf.het) & (gt.mtx == 1)

  # convert True to False
  keep = !(dp.tf | gq.tf | vaf.ref.tf | vaf.hom.tf | vaf.het.tf)

  # For an value in the matrix,
  # if keep is FALSE, then replace corresponding value in GT matrix with 3
  # if keep is TRUE, then do not replace the value
  gt.mtx[!keep] <- 3

  # assign -1 to filtered genotype.
  vaf.mtx[gt.mtx == 3] <- -1


  num.cells <- nrow(keep)
  num.variants <- ncol(keep)
  # By column, number of cells passed QC should be larger than thresholds
  cell.num.keep.tf <- colSums(apply(gt.mtx, 2, function(x){x %in% 0:2})) > num.cells * min.cell.pt / 100
  mut.cell.num.keep.tf <- colSums(apply(gt.mtx, 2, function(x){x %in% 1:2})) > num.cells * min.mut.cell.pt / 100
  variant.keep.tf <- cell.num.keep.tf & mut.cell.num.keep.tf

  v.names <- optima.obj@variants
  v.names[variant.keep.tf]

  # By row, number of cells passed QC
  cell.variants.keep.tf <- rowSums(gt.mtx != 3) > num.variants * min.cell.pt / 100
  c.names <- optima.obj@cell.ids[cell.variants.keep.tf]

  if(optima.obj@proteins[1] == "non-protein"){
    print("non-protein")
    my.protein.mtx <- optima.obj@protein.mtx
  } else {
    my.protein.mtx <- optima.obj@protein.mtx[cell.variants.keep.tf, ]
  }


  # generate new object
  filtered.obj <- new("optima", meta.data=optima.obj@meta.data,

                      cell.ids = optima.obj@cell.ids[cell.variants.keep.tf],
                      cell.labels = optima.obj@cell.labels[cell.variants.keep.tf],
                      # DNA variant
                      variants = optima.obj@variants[variant.keep.tf],
                      variant.filter = "filtered",
                      vaf.mtx = optima.obj@vaf.mtx[cell.variants.keep.tf, variant.keep.tf],
                      gt.mtx = optima.obj@gt.mtx[cell.variants.keep.tf, variant.keep.tf],
                      dp.mtx = optima.obj@dp.mtx[cell.variants.keep.tf, variant.keep.tf],
                      gq.mtx = optima.obj@gq.mtx[cell.variants.keep.tf, variant.keep.tf],

                      # CNV
                      amps = optima.obj@amps,
                      amp.normalize.method = optima.obj@amp.normalize.method,
                      amp.mtx = optima.obj@amp.mtx[cell.variants.keep.tf, ],

                      # protein
                      proteins = optima.obj@proteins,
                      protein.normalize.method = optima.obj@protein.normalize.method,
                      protein.mtx = my.protein.mtx)

  # print some useful numbers
  cat("Number of cells removed: ")
  cat(length(cell.variants.keep.tf) - sum(cell.variants.keep.tf))
  cat("\nNumber of variants removed: ")
  cat(length(variant.keep.tf) - sum(variant.keep.tf))
  cat("\n")
  return(filtered.obj)
}

#' Single variant ID annotation function
#'
#' Returns annotation from one variant ID. This function is not visable to users.
#'
#' @param variant variant name in a string
#' @import httr
#' @import jsonlite
#' @keywords variant
#' @return Annotation information for one specific variant ID.
#' @examples
#' getInfo(variant_id)
getInfo <- function(variant){

  # helper function for getInfo()
  convertAnnotation <- function(x){
    if(class(x) == "list"){
      return(paste(x[[1]], collapse="-"))
    } else {
      return(x[[1]])
    }
  }

  variant <- gsub(":", "-", variant)
  variant <- gsub("/", "-", variant)
  url <- paste("https://api.missionbio.io/annotations/v1/variants?ids=",
               variant, sep = "")

  # use url to fetch data
  res = httr::GET(url)
  data = jsonlite::fromJSON(rawToChar(res$content))

  # initiate empty vector
  num.annotation <- length(data$annotations)
  ret <- rep(NA, num.annotation)

  # add annotation into the vector
  for (i in 1:num.annotation){
    ret[i] <- convertAnnotation(data$annotations[,i]$value)
  }


  return(ret)
}

#' Variant annotation
#'
#' This function takes variant names as input and
#' returns annotation for annotation table for all variant IDs in a data frame.
#'
#' @param variant.names Input variant IDs, can be a vector.
#' @keywords variant
#' @return A data frame with annotation for all input variant IDs.
#' @export
#' @examples annotateVariant(variants_id)
#'
annotateVariant <- function(variant.names){
  # initiate empty dataframe using NULL
  ret <- NULL
  cat("Retriving data from MissionBio...\n")

  for (i in 1:length(variant.names)){

    # generate result for one variant
    variant.ret <- getInfo(variant.names[i])
    ret <- rbind(ret, variant.ret)
    rownames(ret)[i] <- variant.names[i]
  }

  colnames(ret) <- c("variant_type", "transcript_id", "gene",
                     "protein","cDNA", "protein_coding_impact",
                     "function", "allele_freq","DANN",
                     "dbSNP_id", "clinvar", "cosmic")

  cat("Done!\n")
  return(ret)
}

#' Clustering variant function
#'
#' This function identifies cell clones based on DNA variant data.
#'
#' @param optima.obj optima object.
#' @param eps size/radius of the epsilon neighborhood.
#' This argument will passed to dbscan function.
#' @param minPts number of minimum points required in the eps neighborhood
#' for core points, including the point itself.
#' This argument will passed to dbscan function.
#' @param plot if True, a UMAP plot will be generated based on the dimension reduction
#' result from variant matrix.
#' Default is FALSE.
#' @import dbscan
#' @return A list, the first element in the list is the optima
#' object with labels assigned based on dbscan. the second element
#' in the list is the dimension reduction result based on VAF matrix.
#' @export
#' @examples getClones(my.obj)
getClones <- function(optima.obj,
                      eps = 1,
                      minPts = 100,
                      num.PC = 5,
                      plot = FALSE){
  # dimension reduction using PCA and UMAP
  variant.reduceDim <- reduceDim(optima.obj@vaf.mtx, num.PC = num.PC)

  # Use input from UMAP (two columns) for dbscan clustering
  dbscanRet <- dbscan::dbscan(variant.reduceDim[[2]][[1]],
                              eps = eps,
                              minPts = minPts)

  # convert dbscan cluster labels so that they are correctly ordered
  x <- dbscanRet$cluster
  old <- unique(x)
  new <- 1:length(old)
  x[x %in% old] <- new[match(x, old, nomatch = 0)]

  # assign cell label
  optima.obj@cell.labels <- as.character(x)


  if(plot){

    par(mar=c(5, 4, 4, 8), xpd=TRUE)

    plot(variant.reduceDim[[2]][[1]],
         col=topo.colors(length(unique(optima.obj@cell.labels)))[as.factor(optima.obj@cell.labels)],
         xlab = "UMAP1",
         ylab = "UMAP2",
         main=paste(optima.obj@meta.data, "UMAP plot"))

    legend("topright", inset=c(-0.3, 0),
           title="Cell type",
           unique(optima.obj@cell.labels),
           fill=topo.colors(length(unique(optima.obj@cell.labels))),
           cex=0.8)
  }
  return(list(optima.obj, variant.reduceDim))
}

#' The getter function for VAF matrix
#'
#' This function returns the VAF matrix within the optima object.
#'
#' @param optima.obj optima object.
#' @return A matrix that contains VAF data in the optima object.
#' The row names are cell id, the column names are variant ID.
#' @export
#' @examples getDNAmtx(my.obj)

getDNAmtx <- function(optima.obj){
  ret.mtx <- optima.obj@vaf.mtx
  colnames(ret.mtx) <- optima.obj@variants
  rownames(ret.mtx) <- optima.obj@cell.ids
  return(ret.mtx)
}


#' Plot specific variant VAF levels
#'
#' This function create a plot based on VAF dimension reduction result.
#' Within the plot, Each cell was colored based on the VAF
#'
#' @param optima.obj optima object.
#' @param vaf.name the specific variant name user interested
#' @param reduceDim.obj dimension reduction result returned by reduceDim() function.
#' @import ggplot2
#' @return A scatter plot based on dimension reduction result. Each cell is colored
#' based on the protein level. The more expression.
#' @export
#' @examples plotVariantFeature(my.obj, "chr4:106190862:T/C", vaf.reduceDim)


plotVariantFeature <- function(optima.obj,
                               vaf.name,
                               reduceDim.obj){

  # get variant expression value
  values = getDNAmtx(optima.obj)[,vaf.name]


  my.df <- data.frame(UMAP1 = reduceDim.obj[[2]][[1]][,1],
                      UMAP2 = reduceDim.obj[[2]][[1]][,2],
                      VAF = values)

  ggplot(my.df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = VAF), size = 2)+
    theme_classic() +
    ggtitle(paste(vaf.name, "VAF"))

}
