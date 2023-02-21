#' H5 file to optima object function
#'
#' This function read in a h5 file and return one optima object. The
#' h5 file can be found in the Tapestri pipeline software output. The .h5
#' file contains all necessary data needed for single cell DNA and protein
#' analysis.
#'
#' @param directory Directory for the input h5 file.
#' @param sample.name A sample name that will be used naming visualizations.
#' @param omic.type This parameter indicates if the data set is DNA only or DNA+protein. It has
#' two possible values: "DNA+protein" or "DNA". The default value is "DNA+protein".
#' @import rhdf5
#' @return optima object
#' @keywords directory
#' @export
#' @examples
#' readHdf5("path/to/my/file.h5")

readHdf5 <- function(directory, sample.name, omic.type = "DNA+protein"){

  h5f <- rhdf5::H5Fopen(directory, flags = "H5F_ACC_RDONLY")

  if(omic.type == "DNA+protein"){
    my.proteins = as.character(h5f$assays$protein_read_counts$ca$id)
    my.protein.normalize.method = "unnormalized"
    my.protein.mtx = t(h5f$assays$protein_read_counts$layers$read_counts)
  } else if (omic.type == "DNA") {
    cat("DNA data only, skip reading protein data...\n")
    my.proteins = "non-protein"
    my.protein.normalize.method = "non-protein"
    my.protein.mtx = matrix(NA)
  } else {
    stop("illegale argument for omic.type")
  }


  optima.obj <- new("optima",
                    meta.data=sample.name,

                    cell.ids = as.character(h5f$assays$dna_read_counts$ra$barcode),
                    cell.labels = rep("unassigned", length(h5f$assays$dna_read_counts$ra$barcode)),

                    # DNA variant
                    variants = as.character(h5f$assays$dna_variants$ca$id),
                    variant.filter = "unfiltered",
                    vaf.mtx = t(h5f$assays$dna_variants$layers$AF),
                    gt.mtx = t(h5f$assays$dna_variants$layers$NGT),
                    dp.mtx = t(h5f$assays$dna_variants$layers$DP),
                    gq.mtx = t(h5f$assays$dna_variants$layers$GQ),

                    # CNV
                    amps = as.character(h5f$assays$dna_read_counts$ca$id),
                    amp.normalize.method = "unnormalized",
                    amp.mtx = t(h5f$assays$dna_read_counts$layers$read_counts),
                    ploidy.mtx = matrix(),

                    # protein
                    proteins = my.proteins,
                    protein.normalize.method = my.protein.normalize.method,
                    protein.mtx = my.protein.mtx)

  return(optima.obj)
}
