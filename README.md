# optima
optima is an R package for the processing and analysis of data generated from the Tapestri platform providing streamlined functionality for raw data filtering, integration, normalization, transformation, and visualization.

**Manuscript**: https://doi.org/10.1093/bioinformatics/btad611

For questions on installation or usage, please submit an issue via GitHub or contact Dong Pei (dpei@kumc.edu) or Rachel Griffard (rgriffard@kumc.edu).

Please **download** the [vignette file](https://github.com/rachelgriffard/optima/blob/main/vignette.html) file in this repository for a detailed run through of this package functionality.

## Download
To download the optima package, users should use the *install_github* function from the **devtools** package. The full command is as follows:
```
library(devtools)
install_github("rachelgriffard/optima")
```

The latest opitma release is available for download from the [repository](https://github.com/rachelgriffard/optima).

## Basic usage
**Please download the [vignette file](https://github.com/rachelgriffard/optima/blob/main/vignette.html) file in this repository for a detailed run through of this package functionality.**

![optima_flow](https://github.com/rachelgriffard/optima/assets/95938614/62d98e67-621d-413b-8c04-8b15db469894)

### 1. Read in data from h5 format
```
my.obj <- readHdf5(directory = "4-cell-lines-AML-multiomics.dna+protein.h5",
                   sample.name = "Four Cell Mix",![Uploading optima_flow.png…]()

                   omic.type = "DNA+protein")
class(my.obj)
```
The “optima” object is the core object in the optima package. simply run the object to view some summary statistics.

The optima object is defined as a S4 class, a user may also consider define your own optima object within R if they have appropriate data objects.

### 2. Run analyses
   **a. DNA analysis**  
   The first step in DNA analysis is to filter the data set. Cells with too many low quality variants will be discarded.
   Variants with too many low quality data points will be discarded.
   ```
   my.obj.filtered <- filterVariant(my.obj)
   ```
   To get detailed information for each individual variant, The annotateVariant() function can be used. This function will
   generate a dataframe indicating the type of variant, associated genes, protein effect, etc.
   ```
   annotation.table <- annotateVariant(my.obj.filtered$variants)
   head(annotation.table)
   ```
   Users can set clones manually or using the getClones function based on clustering results. The clustering results are stored in the cell.labels vector. Here we print the first 6 cell labels identified from the clustering method.
   ```
   elbowPlot(getDNAmtx(my.obj.filtered))
   set.seed(123)
   clone.ret <- getClones(my.obj.filtered, 
                       minPts = 20, 
                       num.PC = 5,
                       plot = TRUE)
   my.obj.filtered <- clone.ret[[1]]
   vaf.reduceDim <- clone.ret[[2]]
   head(my.obj.filtered$cell.labels)
   ```
   After getting clones labeled, users can use a heatmap to visualize.
   ```
   drawHeatmap(my.obj.filtered, omic.type = "dna")
   ```

   **b. CNV**  
   Users should begin by normalizing their CNV data.
   ```
my.obj.filtered <- normalizeCNV(my.obj.filtered)
my.obj.filtered$amp.normalize.method
```
  Once the cnv data is normalized, we could investigate the CNV for differnt cell populations. This is done by calculating ploidy.
  ```
# diploid cell label
diploid.cell.label = 1

# calculate ploidy
my.obj.filtered <- calculatePloidy(my.obj.filtered, 
                          diploid.cell = diploid.cell.label)
plotPloidy(my.obj.filtered, cell.type = 1)
plotPloidy(my.obj.filtered, cell.type = 2)
```
  **c. Protein**
  The first step in protein analysis is to normalize the read counts using CLR. The input is raw counts, the output is the normalized counts represents protein expression.
  ```
my.obj.filtered <- normalizeProtein(my.obj.filtered)
```
Once the read counts are normalized, we could use dimension reduction to visualize cells in 2D space. The reduceDim() function will conduct PCA then UMAP methods on the protein data. The UMAP result will be used to project cells in a 2D space.
```
# dimension reduction with PCA and UMAP. 
protein.reduceDim <- reduceDim(my.obj.filtered$protein.mtx)
# plot UMAP coordinates after dimension reduction
plot(protein.reduceDim[[2]][[1]][,1],
     protein.reduceDim[[2]][[1]][,2],
     xlab = "UMAP1",
     ylab = "UMAP2",
     col = my.obj.filtered$cell.labels)
```
The drawHeatmap() function will take the protein matrix as input and draw a heatmap based on protein expression data.
```
drawHeatmap(my.obj.filtered, omic.type = "protein")
```
