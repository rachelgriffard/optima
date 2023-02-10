#' H5 file to optima object function
#'
#' This function read in a h5 file and return one optima object. The 
#' h5 file can be found in the Tapestri pipeline software output. The .h5
#' file contains all necessary data needed for single cell DNA and protein 
#' analysis.
#' 
#' @param directory Directory for the input h5 file .
#' @import rhdf5
#' @return optima object
#' @keywords directory
#' @export
#' @examples
#' readHdf5("path/to/my/file.h5")

readHdf5 <- function(directory){
  
  h5f <- rhdf5::H5Fopen(directory, flags = "H5F_ACC_RDONLY")
  
  optima.obj <- new("optima", meta.data="3_cell_mix", 
                   
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
                   proteins = as.character(h5f$assays$protein_read_counts$ca$id),
                   protein.normalize.method = "unnormalized",
                   protein.mtx = t(h5f$assays$protein_read_counts$layers$read_counts))
  
  return(optima.obj)
}